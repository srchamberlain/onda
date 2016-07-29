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
import time
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
        self.coffset = coffset_from_geometry_file(self.geom_filename)
        self.res = res_from_geometry_file(self.geom_filename)
        #
        self.unscaled_radial_profile = numpy.zeros((500))
        self.local_data['radial_average'] = numpy.zeros((500))
        self.local_data['qbins'] = numpy.zeros((500))
        self.intensity_sum_average = 0
        self.std_dev = 0
        self.std_dev_type = 0
        self.N = 1
        self.intensity_threshold = 0
        self.sum = numpy.zeros((500))
        self.radial = numpy.zeros((500, 2))
        self.q_radial = numpy.zeros((500))
        self.profile_to_compare = numpy.zeros((500))
        self.intensity = numpy.zeros((500,400))
        self.Rg_size = 0
        self.Rg = self.Rg_size * [0.0]
        self.I0_size = 0
        self.I0 = self.I0_size * [0.0]
        self.intensity_sum_size = 0
        self.intensity_sum = self.intensity_sum_size * [0.0]
        self.std_dev_1 = numpy.zeros((500))
        self.min_bin = 0
        self.max_bin = 0
        self.count_sum = 0.0
        self.count = 0
        self.count_cumulative = 0.0
        self.percent = 0.0
        self.click = True
        self.click_axis = True
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
        self.refresh_timer.start(10)

    def init_ui(self):

        self.ui.imageView.ui.menuBtn.hide()
        self.ui.imageView.ui.roiBtn.hide()
        
        #
        self.ui.unscaledPlotWidget.setTitle('Unscaled Radial Profile')
        self.ui.unscaledPlotWidget.setLabel('left', text = 'Unscaled Intensity')
        self.ui.unscaledPlotWidget.showGrid(True,True)
        self.ui.unscaledPlotWidget.setXRange(0, 500.0)
        self.ui.unscaledPlotWidget.setXLink(self.ui.scaledPlotWidget)
        self.unscaled_plot = self.ui.unscaledPlotWidget.plot(self.radial[:,0], self.unscaled_radial_profile)
        
        self.ui.scaledPlotWidget.setTitle('Scaled Radial Intensity Profile')
        self.ui.scaledPlotWidget.setLabel('left', text = 'Scaled Intenisty')
        self.ui.scaledPlotWidget.showGrid(True,True)
        self.ui.scaledPlotWidget.setYRange(0, 10.0)
        self.ui.scaledPlotWidget.setXRange(0, 500.0)
        self.scaled_plot = self.ui.scaledPlotWidget.plot(self.radial[:,0], self.radial[:,1])
        
        self.ui.sumPlotWidget.setTitle('Cumulative Average Radial Profile')
        self.ui.sumPlotWidget.setLabel('left', text = 'Cumulative Average Intensity')
        self.ui.sumPlotWidget.showGrid(True,True)
        self.ui.sumPlotWidget.setYRange(0,10.0)
        self.sum_plot = self.ui.sumPlotWidget.plot(self.sum) ########### Plot of Cumulative Radial Intensities
        
        self.ui.intensityPlotWidget.setTitle('unscaled Radial Intensity Sum')
        self.ui.intensityPlotWidget.setLabel('left', text = 'Radial Intensity Sum')
        self.ui.intensityPlotWidget.showGrid(True,True)
        self.ui.intensityPlotWidget.setYRange(0, 2000.0)
        self.intensity_plot = self.ui.intensityPlotWidget.plot(self.intensity_sum)
        
        # Currently not expanded to include specialized SAXS calculations
        #self.ui.RgPlotWidget.setTitle('Radius of Gyration')
        #self.ui.RgPlotWidget.setLabel('left', text = 'Radius of Gyration')
        #self.ui.RgPlotWidget.showGrid(True,True)
        #self.Rg_plot = self.ui.RgPlotWidget.plot(self.Rg)


        self.ui.savePlotButton.clicked.connect(self.save_plot)
        self.ui.resetPlotsButton.clicked.connect(self.reset_plots)
        self.ui.comparePlotsButton.clicked.connect(self.compare_plots)
        self.ui.xAxisTogglePlotsButton.clicked.connect(self.x_toggle_plots)
        #

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
        # Resets Cumulative Average plot and sum of unscaled Radial Intensities plot
        self.intensity_sum_size = 0
        self.intensity_sum = self.intensity_sum_size * [0.0]
        self.intensity_plot.setData(self.intensity_sum)
        self.sum = numpy.zeros((500))
        self.intensity_sum_average = 0
        self.std_dev = 0
        self.count_sum = 0.0
        self.count_cumulative = 0.0


    def save_plot(self):
        ## Saves Cumulative Average radial profile for subtraction or comparison on another run
        self.sum = numpy.nan_to_num(self.sum)
        to_save = numpy.column_stack((self.local_data['qbins'], self.sum))
        numpy.savetxt('profile_to_subtract.dat', to_save, delimiter=" ", fmt="%.5e")
    

    def compare_plots(self):
        ##This function is able to let users compare radial profiles to another profile. Toggles graph on and off. Also checks to see if it should graph with qbins or pixel bins
        ##File fed in through monitor.ini file in process_collect_solutions
        if self.click == True:
            if self.click_axis == True:
                self.profile_to_compare = self.local_data['profile_to_compare']
                xvalues = numpy.arange(0,500)
                self.ui.scaledPlotWidget.plot(xvalues, self.profile_to_compare)
                self.click = False
            else:
                self.profile_to_compare = self.local_data['profile_to_compare']
                self.ui.scaledPlotWidget.plot(self.q_radial, self.profile_to_compare)
                self.click = False
        else:
            self.ui.scaledPlotWidget.clear()
            if self.click_axis == True:
                self.scaled_plot = self.ui.scaledPlotWidget.plot(self.radial[:,0], self.radial[:,1])
                self.click = True
            else:
                self.scaled_plot = self.ui.scaledPlotWidget.plot(self.q_radial, self.radial[:,1])
                self.click = True

    def x_toggle_plots(self):
        ## This Fuction changes variable that indicates whether qbins are on xaxis or pixel bins
        ## Pixel bins are for True
        if self.click_axis == False:
            self.ui.scaledPlotWidget.setLabel('bottom', text = 'pixel bins')
            self.ui.unscaledPlotWidget.setLabel('bottom', text = 'pixel bins')
            self.click_axis = True
        else:
            self.ui.unscaledPlotWidget.setLabel('bottom', text = 'q bins')
            self.ui.scaledPlotWidget.setLabel('bottom', text = 'q bins')
            self.click_axis = False
 
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

        if 'Rg' in self.local_data.keys():
            self.ui.intensityPlotWidget.setTitle('I0') ############
            self.ui.intensityPlotWidget.setLabel('left', text = 'I0')
            if math.isnan(self.local_data['Rg']):
                self.Rg.append(0)
            else:
                self.Rg.append(self.local_data['Rg'])
            if math.isnan(self.local_data['I0']):
                self.I0.append(0)
            else:
                self.I0.append(self.local_data['I0'])
            
            # Not currently expanded to inlude specialized SAXS data
            #self.Rg_plot.setData(self.Rg)
            #self.intensity_plot.setData(self.I0)
            
        else:
            self.intensity_plot.setData(self.intensity_sum)


        QtGui.QApplication.processEvents()

        #if self.local_data['optimized_geometry']:

        #else:

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
        #if timestamp is not None:
            
        #else:

        QtGui.QApplication.processEvents()

        #Updates all Radial Intensity plots, including stacked intensities
        if 'radial_average' in self.local_data.keys():
            # This is defining variables to use
            self.radial[:,0] = numpy.arange(1, 501)
            self.radial[:,1] = self.local_data['radial_average']
            self.unscaled_radial_profile = self.local_data['unscaled_radial_profile']
            self.q_radial = self.local_data['qbins']
            self.N = self.local_data['N']
            self.std_dev_type = self.local_data['std_dev_type']
            self.intensity_threshold = self.local_data['intensity_threshold']
            self.min_bin = self.local_data['min_bin_to_scale']
            self.max_bin = self.local_data['max_bin_to_scale']
            self.num_bins = self.local_data['num_bins']
            self.count_sum += 1
            
        
            ## Incorperate std. dev into cumulative radial profile sum
            # Type 0 corresponds to standard deviation of sum of the intensities
            if self.local_data['intensity_sum'] >= self.intensity_threshold:
                if self.std_dev_type == 0:
                    if self.count_sum == 0:
                        self.sum = self.radial[:,1]
                        self.sum_plot.setData(self.sum)
                        self.count_cumulative += 1
                    else:
                        self.std_dev = numpy.std(self.intensity_sum)
                        self.intensity_sum_average = numpy.mean(self.intensity_sum)
                        if self.local_data['intensity_sum'] >= (self.intensity_sum_average-(self.N*self.std_dev)) and (self.local_data['intensity_sum'] <= self.intensity_sum_average+(self.N*self.std_dev)):
                            self.sum = ((self.sum * self.count_sum)+self.radial[:,1])/(self.count_sum+1)
                            self.sum_plot.setData(self.sum)
                            self.count_cumulative += 1
                            self.percent = numpy.round(self.count_cumulative/(self.count_sum)*100.0, 1)
                            self.ui.NumPlotLabel.setText("Number Averaged: {}".format(self.count_cumulative))
                            self.ui.NumTotalLabel.setText("Number Processed: {}".format(self.count_sum))
                            self.ui.percentPlotLabel.setText("{}% Plotted".format(self.percent))
                #Type 1 indicates standard deviation on each value of scaled region in radial profile
                elif self.std_dev_type == 1:
                    if self.count_sum == 0:
                        self.sum = self.radial[:,1]
                        self.sum_plot.setData(self.sum)
                        self.count_cumulative += 1
                    elif self.count_sum == 1:
                        self.start_std = numpy.column_stack((self.sum, self.radial[:,1]))
                        self.sum = ((self.sum * self.count_sum)+self.radial[:,1])/(self.count_sum+1)
                        self.sum_plot.setData(self.sum)
                        self.std_dev_1 = numpy.std(self.start_std, axis=1)
                        self.count_cumulative += 1
                    else:
                        self.std_dev_1 = numpy.sqrt(((self.count_sum-1)*numpy.square(self.std_dev_1)+(self.radial[:,1]-((self.count_sum*self.sum+self.radial[:,1])/(self.count_sum+1)))*(self.radial[:,1]-self.sum))/self.count_sum)
                        if numpy.all(self.radial[self.min_bin:self.max_bin,1] >= (self.sum[self.min_bin:self.max_bin]-(self.N*self.std_dev_1[self.min_bin:self.max_bin]))) and numpy.all(self.radial[self.min_bin:self.max_bin,1] <= (self.sum[self.min_bin:self.max_bin]+(self.N*self.std_dev_1[self.min_bin:self.max_bin]))):
                            self.sum = ((self.sum * self.count_sum)+self.radial[:,1])/(self.count_sum+1)
                            self.sum_plot.setData(self.sum)
                            self.count_cumulative += 1
                            self.percent = numpy.round(self.count_cumulative/(self.count_sum)*100.0, 1)
                            self.ui.NumPlotLabel.setText("Number Averaged: {}".format(self.count_cumulative))
                            self.ui.NumTotalLabel.setText("Number Processed: {}".format(self.count_sum))
                            self.ui.percentPlotLabel.setText("{}% Plotted".format(self.percent))
                    # None or 2 specifies no std. dev filter
                else:
                    if self.count_sum == 0:
                        self.sum = self.radial[:,1]
                        self.sum_plot.setData(self.sum)
                        self.count_cumulative += 1
                    else:
                        self.sum = ((self.sum * self.count_sum)+self.radial[:,1])/(self.count_sum+1)
                        self.sum_plot.setData(self.sum)
                        self.count_cumulative += 1
                        self.percent = numpy.round(self.count_cumulative/(self.count_sum)*100.0, 1)
                        self.ui.NumPlotLabel.setText("Number Averaged: {}".format(self.count_cumulative))
                        self.ui.NumTotalLabel.setText("Number Processed: {}".format(self.count_sum))
                        self.ui.percentPlotLabel.setText("{}% Plotted".format(self.percent))
            
            
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


            if self.click_axis == True:
                self.unscaled_plot.setData(self.radial[:,0], self.unscaled_radial_profile)
                self.scaled_plot.setData(self.radial[:,0], self.radial[:,1])
            else:
                self.unscaled_plot.setData(self.q_radial, self.unscaled_radial_profile)
                self.scaled_plot.setData(self.q_radial, self.radial[:,1])

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
