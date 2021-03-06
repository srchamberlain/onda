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
import numpy
import signal
import pyqtgraph as pg
import copy
import collections
import cfelpyutils.cfelgeom as cfelgeom
import GUI.UI.onda_hit_viewer_UI

from utils.zmq_gui_utils import ZMQListener

from PyQt4 import (
    QtGui,
    QtCore
)


class MainFrame(QtGui.QMainWindow):
    listening_thread_start_processing = QtCore.pyqtSignal()
    listening_thread_stop_processing = QtCore.pyqtSignal()
    
    def __init__(self, geom_filename, rec_ip, rec_port):
        super(MainFrame, self).__init__()

        self.yx, slab_shape, img_shape = cfelgeom.pixel_maps_for_image_view(geom_filename)

        self.img = numpy.zeros(img_shape, dtype=numpy.float)
        
        self.rec_ip, self.rec_port = rec_ip, rec_port
        self.data = collections.deque(maxlen=20)
        self.data_index = -1
        self.image_update_us = 250
        
        self.zeromq_listener_thread = QtCore.QThread()

        self.zeromq_listener_thread = QtCore.QThread()
        self.zeromq_listener = ZMQListener(self.rec_ip, self.rec_port, u'ondarawdata')
        self.init_listening_thread()

        self.ring_pen = pg.mkPen('r', width=2)
        self.peak_canvas = pg.ScatterPlotItem()
        self.ui = GUI.UI.onda_hit_viewer_UI.Ui_MainWindow()
        self.ui.setupUi(self)
        self.init_ui()

        self.refresh_timer = QtCore.QTimer()
        self.init_timer()
        self.show()

    def init_ui(self):

        self.ui.imageView.ui.menuBtn.hide()
        self.ui.imageView.ui.roiBtn.hide()

        self.ui.imageView.getView().addItem(self.peak_canvas)

        self.ui.backButton.clicked.connect(self.back_button_clicked)
        self.ui.forwardButton.clicked.connect(self.forward_button_clicked)
        self.ui.playPauseButton.clicked.connect(self.play_pause_button_clicked)

    def back_button_clicked(self):
        if self.data_index > 0:
            self.data_index -= 1
            self.update_image_plot()
    
    def forward_button_clicked(self):
        if (self.data_index + 1) < len(self.data):
            self.data_index += 1
            self.update_image_plot()
    
    def play_pause_button_clicked(self):
        if self.refresh_timer.isActive():
            self.refresh_timer.stop()
            self.data_index = len(self.data) - 1
        else:
            self.refresh_timer.start(self.image_update_us)

    def init_listening_thread(self):
        self.zeromq_listener.moveToThread(self.zeromq_listener_thread)
        self.zeromq_listener.zmqmessage.connect(self.data_received)
        self.zeromq_listener.start_listening()
        self.listening_thread_start_processing.connect(self.zeromq_listener.start_listening)
        self.listening_thread_stop_processing.connect(self.zeromq_listener.stop_listening)
        self.zeromq_listener_thread.start()
        self.listening_thread_start_processing.emit()

    def init_timer(self):
        self.refresh_timer.timeout.connect(self.update_image_plot)
        self.refresh_timer.start(self.image_update_us)

    def data_received(self, datdict):
        self.data.append(copy.deepcopy(datdict))

    def update_image_plot(self):
        if len(self.data) > 0:
            data = self.data[self.data_index]

            self.img[self.yx[0], self.yx[1]] = data['raw_data'].ravel().astype(self.img.dtype)

            peak_x = []
            peak_y = []
            for peak_fs, peak_ss in zip(data['peak_list'][0], data['peak_list'][1]):
                peak_in_slab = int(round(peak_ss))*data['raw_data'].shape[1]+int(round(peak_fs))
                peak_x.append(self.yx[1][peak_in_slab])
                peak_y.append(self.yx[0][peak_in_slab])

            self.ui.imageView.setImage(self.img.T, autoLevels=False, autoRange=False, autoHistogramRange=False)
            self.peak_canvas.setData(peak_x, peak_y, symbol='o', size=[5]*len(data['peak_list'][0]),
                                     brush=(255, 255, 255, 0), pen=self.ring_pen,
                                     pxMode=False)
    

def main():
    signal.signal(signal.SIGINT, signal.SIG_DFL)
    app = QtGui.QApplication(sys.argv)
    if len(sys.argv) == 2:
        geom_filename = sys.argv[1]
        rec_ip = '127.0.0.1'
        rec_port = 12321
    elif len(sys.argv) == 4:
        geom_filename = sys.argv[1]
        rec_ip = sys.argv[2]
        rec_port = int(sys.argv[3])
    else:
        print('Usage: onda-hit-viewer-gui.py geometry_filename <listening ip> <listening port>')
        sys.exit()

    _ = MainFrame(geom_filename, rec_ip, rec_port)
    sys.exit(app.exec_())

if __name__ == '__main__':
    main()
