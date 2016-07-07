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


import time
import sys
import zmq
import numpy as np
import scipy.ndimage

from cfelpyutils.cfelgeom import pixel_maps_from_geometry_file, coffset_from_geometry_file

from parallelization_layer.utils import (
    global_params as gp,
    zmq_monitor_utils as zmq_mon,
    dynamic_import as dyn_imp
)
from processing_layer.algorithms.generic_algorithms import (
    DarkCalCorrection
)
from processing_layer.algorithms.solution_algorithms import(
ScaleRadial, QSpace, RpixbinsCalc, RadiusG
)

par_layer = dyn_imp.import_layer_module('parallelization_layer', gp.monitor_params)

MasterWorker = getattr(par_layer, 'MasterWorker')



def radial_intensity(img, rpixbins):
    """
        Averages radial intensities
        """
    
    radial_average = scipy.ndimage.mean(img, labels=rpixbins, index=np.arange(0, rpixbins.max()))
    
    
    return radial_average


def profile_to_subtract(fnam):
    """Extracts profile to be subtracted from the radial average profile from a file, interpolates to rpixbins.
        
        Args
        
        fnam (str): file with profile to subtract. The profile should be same experimental set up!
        
        rpixbins (array-like): bins used to interpolate profile.
        
        Returns
        
        sub_profile (array-like): profile to be subtracted from radial average profiles.
        """
    
    profile = np.loadtxt(fnam, usecols=(0,1))
    nbins = np.arange(1,501)
    
    #x = np.linspace(1,501,501)
    #y = np.sin(0.01*x)**4
    #sub_profile = np.interp(nbins, x, y)
    sub_profile = np.interp(nbins, profile[:,0], profile[:,1])
    #print subtract.shape
    
    return sub_profile



class Onda(MasterWorker):

    def __init__(self, source, monitor_params):

        super(Onda, self).__init__(map_func=self.process_data,
                                   reduce_func=self.collect_data,
                                   source=source, monitor_params=monitor_params)

        gen_params = monitor_params['General']
        dkc_params = monitor_params['DarkCalCorrection']
        rad_params = monitor_params['RadialAveraging']

        self.max_saturated_peaks = gen_params['max_saturated_peaks']
        self.min_num_peaks_for_hit = gen_params['min_num_peaks_for_hit']
        self.max_num_peaks_for_hit = gen_params['max_num_peaks_for_hit']
        self.saturation_value = gen_params['saturation_value']
        self.optimized_geometry = gen_params['geometry_is_optimized']
        self.hit_sending_interval = gen_params['hit_sending_interval']

        self.dkc_filename = dkc_params['filename']
        self.dkc_hdf5_group = dkc_params['hdf5_group']
        
        self.rad_scale = rad_params['scale']
        self.rad_num_rpixbins = rad_params['num_bins']
        self.rad_subtract_profile_filename = rad_params['subtract_profile_filename']
        self.rad_Rg_Calc = rad_params['rg_calc']
        self.rad_ignore_bins = rad_params['ignore_bins']
        self.rad_min_rbin = rad_params['min_bin_to_scale']
        self.rad_max_rbin = rad_params['max_bin_to_scale']

        self.coffset = coffset_from_geometry_file(gen_params['geometry_file'])
        pix_maps = pixel_maps_from_geometry_file(gen_params['geometry_file'])

        self.pixelmap_radius = pix_maps[2]

        self.dkc_apply_mask = 'mask' in dkc_params.keys() and dkc_params['mask'] is True

        self.dkc_mask_filename = None
        self.dkc_mask_hdf5_group = None
        if self.dkc_apply_mask is True:
            if 'mask_filename' in dkc_params.keys():
                self.dkc_mask_filename = dkc_params['mask_filename']

        if self.dkc_apply_mask is True:
            if 'mask_hdf5_group' in dkc_params.keys():
                self.dkc_mask_hdf5_group = dkc_params['mask_hdf5_group']

        self.dkc_gain_map_correction = 'gain_map' in dkc_params.keys() and dkc_params['gain_map'] is True

        self.dkc_gain_map_filename = None
        self.dkc_gain_map_hdf5_group = None
        if self.dkc_gain_map_correction is True:
            if 'gain_map_filename' in dkc_params.keys():
                self.dkc_gain_map_filename = dkc_params['gain_map_filename']

        if self.dkc_gain_map_correction is True:
            if 'gain_map_hdf5_group' in dkc_params.keys():
                self.dkc_gain_map_hdf5_group = dkc_params['gain_map_hdf5_group']


        self.dark_cal_correction = DarkCalCorrection(self.role,
                                                     self.dkc_filename,
                                                     self.dkc_hdf5_group,
                                                     self.dkc_apply_mask,
                                                     self.dkc_mask_filename,
                                                     self.dkc_mask_hdf5_group,
                                                     self.dkc_gain_map_correction,
                                                     self.dkc_gain_map_filename,
                                                     self.dkc_gain_map_hdf5_group)


        
        ######This is for Solution ONLY##########
        self.rpixbins_calculator = RpixbinsCalc(self.role, self.pixelmap_radius)
        self.rpixbins, self.dr = self.rpixbins_calculator.rpixbins_calc(self.rad_num_rpixbins)
        
        
        # Radial Scaling
        if self.rad_scale == True:
            self.radial_data_scaling = ScaleRadial(self.role, self.rad_min_rbin, self.rad_max_rbin)
        
        # Defines Profile to be subtracted from Radial Intensity profile
        if self.rad_subtract_profile_filename == None:
            self.sub_profile = np.zeros((500))
        else:
            self.sub_profile = profile_to_subtract(self.rad_subtract_profile_filename)
            if self.rad_scale == True:
                self.sub_profile = self.radial_data_scaling.scale_profile(self.sub_profile)
        
        # Qspace Calculator
        self.q_space_convert = QSpace(self.role, self.rad_num_rpixbins, self.coffset, self.dr)
        
        # Radial Gyration calcuator #
        if self.rad_Rg_Calc == True:
            self.radius_gyration = RadiusG(self.role)
        #######################################

        if self.role == 'master':

            self.collected_data = {}
            self.publish_ip = gen_params['publish_ip']
            self.publish_port = gen_params['publish_port']
            self.speed_rep_int = gen_params['speed_report_interval']

            running_win_size = gen_params['running_average_size']

            self.hit_rate_running_w = [0.0] * running_win_size
            self.saturation_rate_running_w = [0.0] * running_win_size

            print('Starting the monitor...')
            sys.stdout.flush()

            zmq_mon.init_zmq_to_gui(self, self.publish_ip, self.publish_port)

            self.num_events = 0
            self.old_time = time.time()

            self.time = None
            self.hit_rate = 0
            self.sat_rate = 0

        if self.role == 'worker':

            self.results_dict = {}

            self.hit_sending_counter = 0

            print('Starting worker: {0}.'.format(self.mpi_rank))
            sys.stdout.flush()

        return

    def process_data(self):

        self.results_dict = {}
        
        #corr_raw_data = self.raw_data
        corr_raw_data = self.dark_cal_correction.apply_darkcal_correction(self.raw_data)

        self.results_dict['timestamp'] = self.event_timestamp
        self.results_dict['detector_distance'] = self.detector_distance
        self.results_dict['beam_energy'] = self.beam_energy
        
        
        ###### Creates radial profiles ######
        self.radial_average = radial_intensity(corr_raw_data, self.rpixbins)
        self.results_dict['radial_average']=self.radial_average
        
        ###### converting to q-space  (rpixbins to qbins)######
        self.qbins = self.q_space_convert.convert_to_q(self.role, self.detector_distance, self.beam_energy)
        self.results_dict['qbins'] = self.qbins
        
        ###### Calculates Radius of Gyration ######
        if self.rad_Rg_Calc == True:
            self.Rg = self.radius_gyration.RgCalculator(self.role, self.qbins, self.radial_average, self.ignore_bins)
            self.results_dict['Rg'] = self.Rg
        

        return self.results_dict, self.mpi_rank

    def collect_data(self, new):

        self.collected_data = {}
        self.collected_rawdata = {}

        self.results_dict, _ = new
        self.num_events += 1

        ####### For Solution data analysis #########
        if self.results_dict['radial_average'] is not None:
            self.collected_data['num_cum_avg'] = 20
            self.collected_data['qbins'] = self.results_dict['qbins']
            
            if self.rad_Rg_Calc == True:
                self.collected_data['Rg'] = self.results_dict['Rg']
        
            # Scaling radial intensities, if scaling set to false, then just stores radial average intensity data unscaled
            if self.rad_scale == True:
                scale_this = self.results_dict['radial_average']
                self.radial_average = self.radial_data_scaling.scale_profile(scale_this)
                self.radial_average = self.radial_average-self.sub_profile
                self.collected_data['radial_average'] = self.radial_average
                #### This sums the intesity of the radial profiles after scaling ######
                intensity_sum = np.nansum(scale_this)
                self.collected_data['intensity_sum'] = intensity_sum
            else:
                self.collected_data['radial_average'] = self.results_dict['radial_average']
                self.radial_average = self.radial_average-self.sub_profile
                ### If no Scaling, sums unscaled radial average intesity data #####
                intensity_sum = np.nansum(self.radial_average)
                self.collected_data['intensity_sum'] = intensity_sum
                
            self.collected_data['detector_distance'] = self.results_dict['detector_distance']
            self.collected_data['beam_energy'] = self.results_dict['beam_energy']
            self.collected_data['optimized_geometry'] = self.optimized_geometry
            self.collected_data['timestamp'] = self.results_dict['timestamp']

        self.zmq_publish.send(b'ondadata', zmq.SNDMORE)
        self.zmq_publish.send_pyobj(self.collected_data)

        if 'raw_data' in self.results_dict.keys():
            self.collected_rawdata['raw_data'] = self.results_dict['raw_data']

            self.zmq_publish.send(b'ondarawdata', zmq.SNDMORE)
            self.zmq_publish.send_pyobj(self.collected_rawdata)

        if self.num_events % self.speed_rep_int == 0:
            self.time = time.time()
            print('Processed: {0} in {1:.2f} seconds ({2:.2f} Hz)'.format(
                self.num_events,
                self.time - self.old_time,
                float(self.speed_rep_int)/float(self.time-self.old_time)))
            sys.stdout.flush()
            self.old_time = self.time

        return
