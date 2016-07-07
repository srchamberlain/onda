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
#
#    Added July 7, 2016
#    By Sarah Chamberlain, BioXFEL contribution


import sys
import datetime
import numpy
import scipy.constants
import math
import signal
import pyqtgraph as pg
import copy
import GUI.UI.onda_solutions_UI

from cfelpyutils.cfelgeom import (
    coffset_from_geometry_file,
    res_from_geometry_file,
    pixel_maps_for_image_view
)

from PyQt4 import (
    QtGui,
    QtCore
)

from GUI.utils.zmq_gui_utils import ZMQListener

class QAxis(pg.AxisItem):
    def __init__(self, q_radial_vals, *args, **kwargs):
        pg.AxisItem.__init__(self, *args, **kwargs)
        self.q_radial_vals = numpy.round(q_radial_vals, 2)
    #values = q_radial_vals
    
    def tickStrings(self, values, scale, spacing):
        strings = []
        q_to_resolution = 2*scipy.constants.pi/numpy.round(self.q_radial_vals,2)
        q_to_resolution = numpy.round(q_to_resolution, 2)
        q_and_resolution = numpy.column_stack((self.q_radial_vals, q_to_resolution))
        for v in values[:-1]:
            #strings.append(q_and_resolution[v,:])
            strings.append(self.q_radial_vals[v])
        return strings

class MainFrame(QtGui.QMainWindow):
    """
    The main frame of the application
    """
    listening_thread_start_processing = QtCore.pyqtSignal()
    listening_thread_stop_processing = QtCore.pyqtSignal()

    def __init__(self, geom_filename, rec_ip, rec_port):
        super(MainFrame, self).__init__()

        self.title = 'CSPAD Crystallography Monitor'
        self.data = {}
        self.rec_ip = rec_ip
        self.rec_port = rec_port
        self.geom_filename = geom_filename
        self.local_data = {'time_string': None}
        self.pixel_maps, self.slab_shape, self.img_shape = pixel_maps_for_image_view(self.geom_filename)
        #self.image_center = (self.img_shape[0]/2, self.img_shape[1]/2)
        self.coffset = coffset_from_geometry_file(self.geom_filename)
        self.res = res_from_geometry_file(self.geom_filename)
        #
        self.local_data['radial_average'] = numpy.zeros((500))
        #self.local_data['radial_avg_data'] = numpy.zeros((500))
        self.local_data['qbins'] = numpy.zeros((500))
        self.sum = numpy.zeros((500))
        self.radial = numpy.zeros((500, 2))
        self.q_radial = numpy.zeros((500, 2))
        self.accumulatedRadial = 20
        self.intensity = numpy.zeros((500,400))
        self.Rg_size = 0
        self.Rg = self.Rg_size * [0.0]
        self.intensity_sum_size = 0
        self.intensity_sum = self.intensity_sum_size * [0.0]
        self.count_sum = 0
        self.count = 0
        #

        self.zeromq_listener_thread = QtCore.QThread()
        self.zeromq_listener = ZMQListener(self.rec_ip, self.rec_port, u'ondadata')
        self.init_listening_thread()

        self.vertical_lines = []
        pg.setConfigOption('background', 0.2)
        self.ui = GUI.UI.onda_solutions_UI.Ui_mainWindow()
        self.ui.setupUi(self)
        self.init_ui()

        self.refresh_timer = QtCore.QTimer()
        self.init_timer()
        self.show()

    def init_listening_thread(self):
        self.zeromq_listener.moveToThread(self.zeromq_listener_thread)
        self.zeromq_listener.zmqmessage.connect(self.data_received)
        self.listening_thread_start_processing.connect(self.zeromq_listener.start_listening)
        self.listening_thread_stop_processing.connect(self.zeromq_listener.stop_listening)
        self.zeromq_listener_thread.start()
        self.listening_thread_start_processing.emit()

    def init_timer(self):
        self.refresh_timer.timeout.connect(self.update_image_plot)
        self.refresh_timer.start(500)

    def init_ui(self):

        self.ui.imageView.ui.menuBtn.hide()
        self.ui.imageView.ui.roiBtn.hide()
        
        #
        self.ui.RgPlotWidget.setTitle('Radius of Gyration')
        self.ui.RgPlotWidget.setLabel('left', text = 'Radius of Gyration')
        self.ui.RgPlotWidget.showGrid(True,True)
        self.ui.RgPlotWidget.setYRange(0, 200.0)
        self.Rg_plot = self.ui.RgPlotWidget.plot(self.Rg)
        
        self.ui.sumPlotWidget.setTitle('Cumulative Average Intensity')
        self.ui.sumPlotWidget.setLabel('left', text = 'Cumulative Average Intensity')
        self.ui.sumPlotWidget.showGrid(True,True)
        self.ui.sumPlotWidget.setYRange(0,10.0)
        self.sum_plot = self.ui.sumPlotWidget.plot(self.sum) ########### Plot of Cumulative Radial Intensities

        self.ui.qspacePlotWidget.setTitle('Radial Intensity Plot') ###############
        self.ui.qspacePlotWidget.setLabel('bottom', text = 'pixel bins')
        self.ui.qspacePlotWidget.setLabel('left', text = 'Intenisty')
        self.p1 = self.ui.qspacePlotWidget.plotItem
        self.ui.qspacePlotWidget.showGrid(True,True)
        self.ui.qspacePlotWidget.setYRange(0, 10.0)
        self.ui.qspacePlotWidget.setXRange(0, 500.0)
        self.pixspace_plot = self.p1.plot(self.radial[:,0], self.radial[:,1])
        
        self.ui.intensityPlotWidget.setTitle('unscaled Radial Intensity Sum vs. Events') ############
        self.ui.intensityPlotWidget.setLabel('bottom', text = 'Num of Images')
        self.ui.intensityPlotWidget.setLabel('left', text = 'Radial Intensity Sum')
        self.ui.intensityPlotWidget.showGrid(True,True)
        self.ui.intensityPlotWidget.setYRange(0, 2000.0)
        self.intensity_plot = self.ui.intensityPlotWidget.plot(self.intensity_sum)
        #

        self.ui.savePlotButton.clicked.connect(self.save_plot)
        self.ui.resetPlotsButton.clicked.connect(self.reset_plots)


    def mouse_clicked(self, mouse_evt):
        mouse_pos_in_scene = mouse_evt[0].scenePos()
        if self.ui.intensityPlotWidget.plotItem.sceneBoundingRect().contains(mouse_pos_in_scene):
            if mouse_evt[0].button() == QtCore.Qt.MiddleButton:
                mouse_x_pos_in_data = self.ui.intensityPlotWidget.plotItem.vb.mapSceneToView(mouse_pos_in_scene).x()
                new_vertical_lines = []
                for vert_line in self.vertical_lines:
                    if abs(vert_line.getPos()[0] - mouse_x_pos_in_data) < 5:
                        self.ui.intensityPlotWidget.removeItem(vert_line)
                    else:
                        new_vertical_lines.append(vert_line)
                if len(new_vertical_lines) != len(self.vertical_lines):
                    self.vertical_lines = new_vertical_lines
                    return
                vertical_line = pg.InfiniteLine(mouse_x_pos_in_data, angle=90, movable=False)
                self.vertical_lines.append(vertical_line)
                self.ui.intensityPlotWidget.addItem(vertical_line, ignoreBounds=True)

    def data_received(self, datdict):
        self.data = copy.deepcopy(datdict)


    def reset_plots(self):
        # Now resets sum of intensity plot and starts
        self.intensity_sum = self.intensity_sum_size * [0.0]
        self.intensity_plot.setData(self.intensity_sum)
        self.sum = numpy.zeros((500))
        self.count_sum = 0

    def save_plot(self):
        self.sum = numpy.nan_to_num(self.sum)
        to_save = numpy.column_stack((self.radial[:,0], self.sum))
        for i in range(0,len(self.sum)):
            numpy.savetxt('profile_to_subtract.dat', to_save, delimiter=" ", fmt="%3d %1.4f")


    def update_image_plot(self):

        if len(self.data.keys()) != 0:
            self.local_data = self.data
            self.data = {}
        else:
            return

        QtGui.QApplication.processEvents()

        if math.isnan(self.local_data['intensity_sum']):
            self.intensity_sum.append(0)
        else:
            self.intensity_sum.append(self.local_data['intensity_sum'])

        self.intensity_plot.setData(self.intensity_sum)

        if 'Rg' in self.local_data.keys():
            if math.isnan(self.local_data['Rg']):
                self.Rg.append(0)
            else:
                self.Rg.append(self.local_data['Rg'])
    
        self.Rg_plot.setData(self.Rg)


        QtGui.QApplication.processEvents()

        """if self.local_data['optimized_geometry']:

        else:
        """

        new_vertical_lines = []
        for vline in self.vertical_lines:
            line_pos = vline.getPos()[0]
            line_pos -= 1
            if line_pos > 0.0:
                vline.setPos(line_pos)
                new_vertical_lines.append(vline)
            else:
                self.ui.intensityPlotWidget.removeItem(vline)
        self.ui.vertical_lines = new_vertical_lines

        QtGui.QApplication.processEvents()

        timestamp = self.local_data['timestamp']
        """if timestamp is not None:
            
        else:
        """

        QtGui.QApplication.processEvents()

        #Updates all Radial Intensity plots, including stacked intensities
        if 'radial_average' in self.local_data.keys():
            self.radial[:,0] = numpy.arange(1, 501)
            self.radial[:,1] = self.local_data['radial_average']
            
            self.q_radial[:,0] = self.local_data['qbins']##########
            self.q_radial[:,1] = self.local_data['radial_average']
            
            if self.count == 0:
                self.qaxis = QAxis(self.q_radial[:,0], orientation = 'top', parent = self.p1)
                self.qaxis.setGeometry(self.p1.vb.sceneBoundingRect())
                self.qaxis.linkToView(self.p1.vb)
                self.qaxis.setLabel('qbins')
            else:
                self.qaxis.setGeometry(self.p1.vb.sceneBoundingRect())
                self.qaxis.linkToView(self.p1.vb)
                self.qaxis.setLabel('qbins')
        
            self.pixspace_plot.setData(self.radial[:,0], self.radial[:,1])
            
            ########### Running cumulative average of radial profiles ################
            if self.count_sum == 0:
                self.sum = self.radial[:,1]
            else:
                self.sum = ((self.sum * self.count_sum)+self.radial[:,1])/(self.count_sum+1)
                self.sum_plot.setData(self.sum)
            
            #if 'radial_avg_data' in self.local_data.keys():
            #self.sum = self.local_data['radial_avg_data']
            
            ####### Intensity image Viewer #################
            if self.count == 0:
                self.intensity[:,(400-1)-self.count] = self.radial[:,1]
                self.ui.imageView.setImage(self.intensity, autoHistogramRange = True, autoLevels=False, autoRange=True)
            elif self.count < 400:
                self.intensity[:,(400-1)-self.count] = self.radial[:,1]
                self.ui.imageView.setImage(self.intensity, autoHistogramRange = False, autoLevels=False, autoRange=False)
            else:
                self.intensity = numpy.roll(self.intensity, 1, axis=1)
                self.intensity[:,0]=self.radial[:,1]
                self.ui.imageView.setImage(self.intensity, autoHistogramRange = False, autoLevels=False, autoRange=False)

            self.count += 1
            self.count_sum += 1



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
        print('Usage: onda-solutions-gui.py geometry_filename <listening ip> <listening port>')
        sys.exit()

    _ = MainFrame(geom_filename, rec_ip, rec_port)
    sys.exit(app.exec_())

if __name__ == '__main__':
    main()
