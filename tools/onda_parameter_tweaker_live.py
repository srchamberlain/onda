#!/usr/bin/env python
#    This file is part of OnDA.
#
#    OnDA is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    OnDA is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with OnDA.  If not, see <http://www.gnu.org/licenses/>.

import sys
import PyQt4.QtCore
import PyQt4.QtGui
import pyqtgraph
import numpy
import random
import signal
import copy
try:
    import configparser
except ImportError:
    import ConfigParser as configparser
import GUI.UI.parameter_tweaker_UI

import collections

from GUI.utils.zmq_gui_utils import ZMQListener
from cfelpyutils.cfeloptarg import parse_parameters
from cfelpyutils.cfelhdf5 import load_nparray_from_hdf5_file
from peakfinder8_extension import peakfinder_8

from cfelpyutils.cfelgeom import (
    pixel_maps_from_geometry_file,
    pixel_maps_for_image_view
)


def check_changed_parameter(param, param_conv_vers, lineedit_element):
        try:
            new_param = param_conv_vers(lineedit_element.text())
            if new_param != param:
                return new_param, True
            else:
                return param, False
        except Exception:
            lineedit_element.setText(str(param))
            return param, False


class MainFrame(PyQt4.QtGui.QMainWindow):
    """
    The main frame of the application
    """
    listening_thread_start_processing = PyQt4.QtCore.pyqtSignal()
    listening_thread_stop_processing = PyQt4.QtCore.pyqtSignal()

    def __init__(self, monitor_params, rec_ip, rec_port):
        super(MainFrame, self).__init__()

        self.monitor_params = monitor_params

        gen_params = monitor_params['General']
        p8pd_params = monitor_params['Peakfinder8PeakDetection']
        
        self.rec_ip, self.rec_port = rec_ip, rec_port
        self.data = collections.deque(maxlen=20)
        self.data_index = 0
        self.image_update_us = 250

        self.zeromq_listener_thread = PyQt4.QtCore.QThread()
        self.zeromq_listener = ZMQListener(self.rec_ip, self.rec_port, u'ondarawdata')
        self.init_listening_thread()

        self.ring_pen = pyqtgraph.mkPen('r', width=2)
        self.circle_pen = pyqtgraph.mkPen('b', width=2)

        pix_maps = pixel_maps_from_geometry_file(monitor_params['General']['geometry_file'])
        self.pixelmap_radius = pix_maps[2]

        self.pixel_maps, self.slab_shape, self.img_shape = pixel_maps_for_image_view(gen_params['geometry_file'])
        self.img_to_draw = numpy.zeros(self.img_shape, dtype=numpy.float32)
        self.mask_to_draw = numpy.zeros(self.img_shape+(3,), dtype=numpy.int16)
        self.max_num_peaks = int(p8pd_params['max_num_peaks'])
        self.asic_nx = int(p8pd_params['asics_nx'])
        self.asic_ny = int(p8pd_params['asics_ny'])
        self.nasics_x = int(p8pd_params['nasics_x'])
        self.nasics_y = int(p8pd_params['nasics_y'])
        self.adc_thresh = float(p8pd_params['adc_threshold'])
        self.minimum_snr = float(p8pd_params['minimum_snr'])
        self.min_pixel_count = int(p8pd_params['min_pixel_count'])
        self.max_pixel_count = int(p8pd_params['max_pixel_count'])
        self.local_bg_radius = int(p8pd_params['local_bg_radius'])
        self.mask_filename = p8pd_params['mask_filename']
        self.mask_hdf5_path = p8pd_params['mask_hdf5_path']
        self.min_res = int(p8pd_params['min_res'])
        self.max_res = int(p8pd_params['max_res'])
        self.loaded_mask = load_nparray_from_hdf5_file(self.mask_filename, self.mask_hdf5_path)
        self.min_num_peaks_for_hit = int(monitor_params['General']['min_num_peaks_for_hit'])
        self.max_num_peaks_for_hit = int(monitor_params['General']['max_num_peaks_for_hit'])

        self.res_mask = numpy.ones(self.slab_shape, dtype=numpy.int8)
        self.res_mask[numpy.where(self.pixelmap_radius < self.min_res)] = 0
        self.res_mask[numpy.where(self.pixelmap_radius > self.max_res)] = 0
        self.mask = self.loaded_mask * self.res_mask

        mask = self.loaded_mask.copy().astype(numpy.float)
        mask = mask * 255./mask.max()
        mask = 255. - mask
        self.mask_to_draw[self.pixel_maps[0], self.pixel_maps[1], 1] = mask.ravel()

        self.mask_image_view = pyqtgraph.ImageItem()
        self.peak_canvas = pyqtgraph.ScatterPlotItem()
        self.circle_canvas = pyqtgraph.ScatterPlotItem()

        self.adc_threshold_label = PyQt4.QtGui.QLabel(self)
        self.adc_threshold_label.setText('adc_threshold')
        self.adc_threshold_lineedit = PyQt4.QtGui.QLineEdit(self)
        self.adc_threshold_lineedit.setText(str(self.monitor_params['Peakfinder8PeakDetection']['adc_threshold']))
        self.adc_threshold_lineedit.editingFinished.connect(self.update_peaks)
        self.hlayout0 = PyQt4.QtGui.QHBoxLayout()
        self.hlayout0.addWidget(self.adc_threshold_label)
        self.hlayout0.addWidget(self.adc_threshold_lineedit)

        self.min_snr_label = PyQt4.QtGui.QLabel(self)
        self.min_snr_label.setText('minmum_snr')
        self.min_snr_lineedit = PyQt4.QtGui.QLineEdit(self)
        self.min_snr_lineedit.setText(str(self.monitor_params['Peakfinder8PeakDetection']['minimum_snr']))
        self.min_snr_lineedit.editingFinished.connect(self.update_peaks)
        self.hlayout1 = PyQt4.QtGui.QHBoxLayout()
        self.hlayout1.addWidget(self.min_snr_label)
        self.hlayout1.addWidget(self.min_snr_lineedit)

        self.min_pixel_count_label = PyQt4.QtGui.QLabel(self)
        self.min_pixel_count_label.setText('min_pixel_count')
        self.min_pixel_count_lineedit = PyQt4.QtGui.QLineEdit(self)
        self.min_pixel_count_lineedit.setText(str(self.monitor_params['Peakfinder8PeakDetection']['min_pixel_count']))
        self.min_pixel_count_lineedit.editingFinished.connect(self.update_peaks)
        self.hlayout2 = PyQt4.QtGui.QHBoxLayout()
        self.hlayout2.addWidget(self.min_pixel_count_label)
        self.hlayout2.addWidget(self.min_pixel_count_lineedit)

        self.max_pixel_count_label = PyQt4.QtGui.QLabel(self)
        self.max_pixel_count_label.setText('max_pixel_count')
        self.max_pixel_count_lineedit = PyQt4.QtGui.QLineEdit(self)
        self.max_pixel_count_lineedit.setText(str(self.monitor_params['Peakfinder8PeakDetection']['max_pixel_count']))
        self.max_pixel_count_lineedit.editingFinished.connect(self.update_peaks)
        self.hlayout3 = PyQt4.QtGui.QHBoxLayout()
        self.hlayout3.addWidget(self.max_pixel_count_label)
        self.hlayout3.addWidget(self.max_pixel_count_lineedit)

        self.local_bg_radius_label = PyQt4.QtGui.QLabel(self)
        self.local_bg_radius_label.setText('local_bg_raidus')
        self.local_bg_radius_lineedit = PyQt4.QtGui.QLineEdit(self)
        self.local_bg_radius_lineedit.setText(str(self.monitor_params['Peakfinder8PeakDetection']['local_bg_radius']))
        self.local_bg_radius_lineedit.editingFinished.connect(self.update_peaks)
        self.hlayout4 = PyQt4.QtGui.QHBoxLayout()
        self.hlayout4.addWidget(self.local_bg_radius_label)
        self.hlayout4.addWidget(self.local_bg_radius_lineedit)

        self.min_res_label = PyQt4.QtGui.QLabel(self)
        self.min_res_label.setText('min_res')
        self.min_res_lineedit = PyQt4.QtGui.QLineEdit(self)
        self.min_res_lineedit.setText(str(self.min_res))
        self.min_res_lineedit.editingFinished.connect(self.update_peaks)
        self.hlayout5 = PyQt4.QtGui.QHBoxLayout()
        self.hlayout5.addWidget(self.min_res_label)
        self.hlayout5.addWidget(self.min_res_lineedit)

        self.max_res_label = PyQt4.QtGui.QLabel(self)
        self.max_res_label.setText('max_res')
        self.max_res_lineedit = PyQt4.QtGui.QLineEdit(self)
        self.max_res_lineedit.setText(str(self.max_res))
        self.max_res_lineedit.editingFinished.connect(self.update_peaks)
        self.hlayout6 = PyQt4.QtGui.QHBoxLayout()
        self.hlayout6.addWidget(self.max_res_label)
        self.hlayout6.addWidget(self.max_res_lineedit)

        self.param_label = PyQt4.QtGui.QLabel(self)
        self.param_label.setText('<b>Peakfinder Parameters:</b>')

        self.ui = GUI.UI.parameter_tweaker_UI.Ui_MainWindow()
        self.ui.setupUi(self)
        self.init_ui()
        self.setWindowTitle('OnDA Live Parameter Tweaker')

        self.proxy = pyqtgraph.SignalProxy(self.ui.imageView.getView().scene().sigMouseClicked,
                                           slot=self.mouse_clicked)
        self.update_peaks()
        self.draw_things()

        self.refresh_timer = PyQt4.QtCore.QTimer()
        self.init_timer()
        self.show()

    def init_ui(self):

        self.ui.imageView.ui.menuBtn.hide()
        self.ui.imageView.ui.roiBtn.hide()

        self.ui.imageView.getView().addItem(self.mask_image_view)

        self.ui.imageView.getView().addItem(self.peak_canvas)
        self.ui.imageView.getView().addItem(self.circle_canvas)
        self.proxy = pyqtgraph.SignalProxy(self.ui.imageView.getView().scene().sigMouseClicked, slot=self.mouse_clicked)

        self.ui.forwardButton.clicked.connect(self.next_event)
        self.ui.backButton.clicked.connect(self.previous_event)
        self.ui.randomButton.clicked.connect(self.play_pause_button_clicked)
        self.ui.randomButton.setText('Play/Pause')

        self.ui.verticalLayout1.insertLayout(0, self.hlayout6)
        self.ui.verticalLayout1.insertLayout(0, self.hlayout5)
        self.ui.verticalLayout1.insertLayout(0, self.hlayout4)
        self.ui.verticalLayout1.insertLayout(0, self.hlayout3)
        self.ui.verticalLayout1.insertLayout(0, self.hlayout2)
        self.ui.verticalLayout1.insertLayout(0, self.hlayout1)
        self.ui.verticalLayout1.insertLayout(0, self.hlayout0)
        self.ui.verticalLayout1.insertWidget(0, self.param_label)
        self.ui.splitter.setStretchFactor(0, 1)
        self.ui.splitter.setStretchFactor(1, 0)

        self.ui.showHidePeaksCheckBox.stateChanged.connect(self.draw_things)
        self.ui.resolutionRingsCheckBox.stateChanged.connect(self.draw_things)

    def init_listening_thread(self):
        self.zeromq_listener.moveToThread(self.zeromq_listener_thread)
        self.zeromq_listener.zmqmessage.connect(self.data_received)
        self.zeromq_listener.start_listening()
        self.listening_thread_start_processing.connect(self.zeromq_listener.start_listening)
        self.listening_thread_stop_processing.connect(self.zeromq_listener.stop_listening)
        self.zeromq_listener_thread.start()
        self.listening_thread_start_processing.emit()

    def init_timer(self):
        self.refresh_timer.timeout.connect(self.draw_things)
        self.refresh_timer.start(self.image_update_us)

    def data_received(self, datdict):
        if self.refresh_timer.isActive():
            self.data.append(copy.deepcopy(datdict))

    def draw_things(self):
        if len(self.data) == 0:
            return None
        
        img = self.data[self.data_index]['raw_data']

        self.img_to_draw[self.pixel_maps[0], self.pixel_maps[1]] = img.ravel()
        self.ui.imageView.setImage(self.img_to_draw.T, autoLevels=False, autoRange=False, autoHistogramRange=False)
        self.mask_image_view.setImage(numpy.transpose(self.mask_to_draw, axes=(1,0,2)), autoLevels=False, autoRange=False, opacity=0.1)

        peak_list = peakfinder_8(
            self.max_num_peaks,
            img.astype(numpy.float32),
            self.mask.astype(numpy.int8),
            self.pixelmap_radius,
            self.asic_nx,
            self.asic_ny,
            self.nasics_x,
            self.nasics_y,
            self.adc_thresh,
            self.minimum_snr,
            self.min_pixel_count,
            self.max_pixel_count,
            self.local_bg_radius)

        if self.ui.showHidePeaksCheckBox.isChecked():

            peak_x = []
            peak_y = []
            for peak_fs, peak_ss in zip(peak_list[0], peak_list[1]):
                peak_in_slab = int(round(peak_ss))*self.slab_shape[1]+int(round(peak_fs))
                try:
                    peak_x.append(self.pixel_maps[0][peak_in_slab])
                    peak_y.append(self.pixel_maps[1][peak_in_slab])
                except Exception:
                    pass
            self.peak_canvas.setData(peak_y, peak_x, symbol='o', size=15, pen=self.ring_pen, brush=(0, 0, 0, 0),
                                     pxMode=False)

            hit = self.min_num_peaks_for_hit < len(peak_list[2]) < self.max_num_peaks_for_hit

            if hit:
                self.ui.hitLabel.setText('Hit [{0}-{1} peaks]: <b>Yes</b> ({2} peaks)'.format(
                    self.min_num_peaks_for_hit, self.max_num_peaks_for_hit, len(peak_list[2])))
            else:
                self.ui.hitLabel.setText('Hit [{0}-{1} peaks]: No ({2} peaks)'.format(
                    self.min_num_peaks_for_hit, self.max_num_peaks_for_hit, len(peak_list[2])))

        else:

            self.ui.hitLabel.setText('Hit [{0}-{1} peaks]: - (- peaks)'.format(self.min_num_peaks_for_hit,
                                                                               self.max_num_peaks_for_hit))
            self.peak_canvas.setData([])

        if self.ui.resolutionRingsCheckBox.isChecked():
            self.circle_canvas.setData([self.img_shape[1]/2, self.img_shape[1]/2],
                                       [self.img_shape[0]/2, self.img_shape[0]/2],
                                       symbol='o', size=[2 * self.min_res, 2 * self.max_res],
                                       pen=self.circle_pen, brush=(0, 0, 0, 0), pxMode=False)

        else:

            self.circle_canvas.setData([])

    def update_peaks(self):

        something_changed = False
        self.adc_thresh, changed = check_changed_parameter(self.adc_thresh, float, self.adc_threshold_lineedit)
        if changed:
            something_changed = True
        self.minimum_snr, changed = check_changed_parameter(self.minimum_snr, float, self.min_snr_lineedit)
        if changed:
            something_changed = True
        self.min_pixel_count, changed = check_changed_parameter(self.min_pixel_count, int,
                                                                self.min_pixel_count_lineedit)
        if changed:
            something_changed = True
        self.max_pixel_count, changed = check_changed_parameter(self.max_pixel_count, int,
                                                                self.max_pixel_count_lineedit)
        if changed:
            something_changed = True
        self.local_bg_radius, changed = check_changed_parameter(self.local_bg_radius, int,
                                                                self.local_bg_radius_lineedit)
        if changed:
            something_changed = True
        self.min_res, changed = check_changed_parameter(self.min_res, int, self.min_res_lineedit)
        if changed:
            something_changed = True
        self.max_res, changed = check_changed_parameter(self.max_res, int, self.max_res_lineedit)
        if changed:
            something_changed = True

        self.res_mask = numpy.ones(self.slab_shape, dtype=numpy.int8)
        self.res_mask[numpy.where(self.pixelmap_radius < self.min_res)] = 0
        self.res_mask[numpy.where(self.pixelmap_radius > self.max_res)] = 0
        self.mask = self.loaded_mask * self.res_mask

        if something_changed:
            self.draw_things()

    def previous_event(self):
        if self.data_index > 0:
            self.data_index -= 1
            self.draw_things()

    def next_event(self):
        if (self.data_index + 1) < len(self.data):
            self.data_index += 1
            self.draw_things()

    def play_pause_button_clicked(self):
        if self.refresh_timer.isActive():
            self.refresh_timer.stop()
            self.data_index = len(self.data) - 1
        else:
            self.refresh_timer.start(self.image_update_us)

    def mouse_clicked(self, event):
        pos = event[0].scenePos()
        if self.ui.imageView.getView().sceneBoundingRect().contains(pos):
            mouse_point = self.ui.imageView.getView().mapSceneToView(pos)
            x_mouse = int(mouse_point.x())
            y_mouse = int(mouse_point.y())
            if 0 < x_mouse < self.img_to_draw.shape[1] and 0 < y_mouse < self.img_to_draw.shape[0]:
                self.ui.lastClickedPositionLabel.setText('Last clicked position: (%g,%g)' % (x_mouse, y_mouse))
                self.ui.lastClickedPixelValueLabel.setText('Pixel Value: %5.1f' % (self.img_to_draw[y_mouse,
                                                                                                    x_mouse]))


def main():
    config = configparser.ConfigParser()

    signal.signal(signal.SIGINT, signal.SIG_DFL)
    app = PyQt4.QtGui.QApplication(sys.argv)
    if len(sys.argv) == 1:
        rec_ip = '127.0.0.1'
        rec_port = 12321
    elif len(sys.argv) == 3:
        rec_ip = sys.argv[1]
        rec_port = int(sys.argv[2])
    else:
        print('Usage: onda_parameter_tweaker_live.py <listening ip> <listening port>')
        sys.exit()

    config.read("monitor.ini")
    monitor_params = parse_parameters(config)

    _ = MainFrame(monitor_params, rec_ip, rec_port)
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()
