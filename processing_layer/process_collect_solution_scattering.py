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
from processing_layer.algorithms.solution_scattering_algorithms import (
  PixelSpaceQSpaceConversion, RadiusOfGyration
)

par_layer = dyn_imp.import_layer_module('parallelization_layer', gp.monitor_params)

MasterWorker = getattr(par_layer, 'MasterWorker')


def calculate_rpixbins (pixelmap_radius, nbins):
    """Calculates radius bin pixel map

    Calculates a pixelmap containing radius bin labels for each pixel

        Args:

            pixelmap_radius (numpy.ndarray): pixel map containing absolute radius value for each pixel

            nbins (int): number if bins required by the user

        Returns:

            rpixbins (numpy.ndarray int): array of labels cooresponding to each value in pixelmap_radius
            for radial intensity function

            dr (int): size in pixels of each bin
     """

    rpixbins = np.zeros(pixelmap_radius.shape, dtype=int)

    dr = float(np.max(pixelmap_radius)) / nbins

    for i in range(0, nbins - 1):

        rpixbins[(pixelmap_radius >= i * dr) & (pixelmap_radius < (i + 1) * dr)] = i
        rpixbins[pixelmap_radius >= (nbins) * dr] = nbins

    return rpixbins, dr


def scale_radial_profile(radial_int, min_rpixbin, max_rpixbin):
    """Scales a radial profile

    Scales a radial profile based on the average intensity value in a region defined by user

    Args:

       radial_int (numpy.ndarray): radial profile of intensities to scale

       min_rpixbin (int): Minimum bin number for scaling region, readi in from monitor params

       max_rpixbin (int): Maximum bin number for scaling region, readi in from monitor params

    Returns:

        self.radial_int_new (numpy.ndarray): array of intensity values scaled using the region specified by the user
    """

    region_int = radial_int[min_rpixbin:max_rpixbin]
    average = np.average(region_int)
    radial_int_new = radial_int / average
    return radial_int_new


def calculate_average_radial_intensity(img, rpixbins):
    """Calculates average radial intensities.

    Calculates radial average of input data. Input data is subdivided into radius bins according to a radius bin pixel
    map. An average intensity is then computed for each radius bin.

    Args:

        img (ndarray): raw image

        rpixbins (ndarray): pixel map describing the radius bin each pixels falls into

    Returns:

        radial_average (ndarray): average intensity values for each radius bin
    """

    radial_average = scipy.ndimage.mean(img, labels=rpixbins, index=np.arange(0, rpixbins.max()))
    
    
    return radial_average


def load_and_interpolate_radial_profile(fnam, num_of_bins):
    """Loads a radial profile and interpolates to a different number of radius bins.

    Loads a radial profile from a text file and interpolates intensity values if the required number of radius bins
    differs from the one in the file.

    Args:

        fnam (str): name of the file which contains the radial profile.

        num_of_bins (int): required number of radius bins.

    Returns:

        interpolated_profile (array-like): profile to be subtracted from radial average profiles.
    """
    
    profile = np.loadtxt(fnam, usecols=(0,1))
    nbins = np.arange(0,num_of_bins)
    
    #x = np.linspace(1,501,501)
    #y = np.sin(0.01*x)**4
    #sub_profile = np.interp(nbins, x, y)
    sub_profile = np.interp(nbins, profile[:,0], profile[:,1], 0, 0)

    
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
        
        # Reads in Radial Section of monitor.ini file
        self.rad_pump_laser_on = rad_params['pump_laser_on']
        self.rad_pump_laser_exr_code = rad_params['pump_laser_evr_code']
        self.rad_scale = rad_params['scale']
        self.rad_num_rpixbins = rad_params['num_bins']
        self.rad_subtract_profile_filename = rad_params['subtract_profile_filename']
        self.rad_rg_calc = rad_params['rg_calc']
        self.rad_ignore_bins = rad_params['ignore_bins']
        self.rad_min_rbin = rad_params['min_bin_to_scale']
        self.rad_max_rbin = rad_params['max_bin_to_scale']
        self.rad_profile_compare_filename = rad_params['profile_compare_filename']
        self.rad_threshold = rad_params['intensity_threshold']
        self.rad_std_dev_type = rad_params['std_dev_type']
        self.rad_n = rad_params['n_sigma']

        self.coffset = coffset_from_geometry_file(gen_params['geometry_file'])
        _, _, self.pixelmap_radius = pixel_maps_from_geometry_file(gen_params['geometry_file'])

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

        # Calculates radius bin pixelmap
        self.rpixbins, self.dr = calculate_rpixbins(self.pixelmap_radius, self.rad_num_rpixbins)

        # Initialize Q bins calculator
        self.q_space_convert = PixelSpaceQSpaceConversion(self.role, self.rad_num_rpixbins, self.coffset, self.dr)

        # Radial Gyration calculator
        # Not yet ready for use, set to false in Monitor.ini file
        if self.rad_rg_calc == True:
            self.radius_gyration = RadiusOfGyration(self.role)

        if self.role == 'master':

            self.collected_data = {}
            self.publish_ip = gen_params['publish_ip']
            self.publish_port = gen_params['publish_port']
            self.speed_rep_int = gen_params['speed_report_interval']

            print('Starting the monitor...')
            sys.stdout.flush()

            zmq_mon.init_zmq_to_gui(self, self.publish_ip, self.publish_port)

            self.num_events = 0
            self.old_time = time.time()

            self.time = None

            # Defines profile to be subtracted from Radial Intensity profile
            if self.rad_subtract_profile_filename is None:
                self.sub_profile = np.zeros((500))
                # TODO: Why hardcode 500? Shouldn't it be np.zeros((self.rad_num_rpixbins))?
            else:
                self.sub_profile = load_and_interpolate_radial_profile(self.rad_subtract_profile_filename,
                                                                       self.rad_num_rpixbins)
                if self.rad_scale == True:
                    self.sub_profile = scale_radial_profile(self.sub_profile, self.rad_min_rbin, self.rad_max_rbin)

            if self.rad_scale == True:
                self.scale_profile = scale_radial_profile
            else:
                self.scale_profile = lambda radial_int, min_rpixbin, max_rpixbin: radial_int

        if self.role == 'worker':

            self.results_dict = {}

            self.hit_sending_counter = 0

            # Profile to compare with Radial Intensity Profile
            if self.rad_profile_compare_filename is not None:
                self.profile_compare = np.loadtxt(self.rad_profile_compare_filename, usecols=(0, 1))
                self.profile_compare_x = self.profile_compare[:, 0]
                self.profile_compare_y = self.profile_compare[:, 1]
            else:
                self.profile_to_compare = np.zeros((self.rad_num_rpixbins))

            print('Starting worker: {0}.'.format(self.mpi_rank))
            sys.stdout.flush()

        return

    def process_data(self):

        self.results_dict = {}
        
        corr_raw_data = self.dark_cal_correction.apply_darkcal_correction(self.raw_data)

        # Defines variables to be sent to the master
        self.results_dict['timestamp'] = self.event_timestamp
        self.results_dict['detector_distance'] = self.detector_distance
        self.results_dict['beam_energy'] = self.beam_energy
        self.results_dict['std_dev_type'] = self.rad_std_dev_type
        self.results_dict['N'] = self.rad_n
        self.results_dict['min_bin_to_scale'] = self.rad_min_rbin
        self.results_dict['max_bin_to_scale'] = self.rad_max_rbin
        self.results_dict['num_bins'] = self.rad_num_rpixbins
        self.results_dict['intensity_threshold'] = self.rad_threshold
        self.results_dict['pump_laser_state'] = self.pump_laser_state
        self.results_dict['pump_laser_on'] = self.rad_pump_laser_on
        
        # Calculate radial pofile
        self.results_dict['radial_average'] = calculate_average_radial_intensity(corr_raw_data, self.rpixbins)
        
        # Convertion to q-space
        self.results_dict['qbins'] = self.q_space_convert.convert_to_q(self.role, self.detector_distance,
                                                                       self.beam_energy)

        # Compare radial average to provided profile
        if self.rad_profile_compare_filename is not None:
            self.profile_to_compare = np.interp(self.qbins, self.profile_compare_x, self.profile_compare_y, 0, 0)
            ## Scales profile to same region as radial profile, only uses nonzero values of profile
            data_to_scale = self.profile_to_compare[self.rad_min_rbin:self.rad_max_rbin]
            self.profile_to_compare /= np.average(data_to_scale) #, weights=data_to_scale.astype(bool))
        self.results_dict['profile_to_compare'] = self.profile_to_compare
        
        # Calculates Radius of Gyration
        # Not yet ready for use, set to false in Monitor.ini file
        if self.rad_rg_calc == True:
            self.rg, self.izero = self.radius_gyration.calculate_rog(self.role, self.qbins, self.radial_average,
                                                                 self.rad_ignore_bins)
            self.results_dict['rg'] = self.rg
            self.results_dict['izero'] = self.izero
    
        return self.results_dict, self.mpi_rank

    def collect_data(self, new):

        self.collected_data = {}
        self.collected_rawdata = {}

        self.results_dict, _ = new
        self.num_events += 1

        self.collected_data['qbins'] = self.results_dict['qbins']

        # Not yet ready for use, set to false in Monitor.ini file
        if self.rad_rg_calc == True:
            self.collected_data['rg'] = self.results_dict['rg']
            self.collected_data['izero'] = self.results_dict['izero']

        # Scaling radial intensities and subtracts profile
        self.collected_data['unscaled_radial_profile'] = self.results_dict['radial_average']
        self.collected_data['radial_average'] = (self.scale_profile(self.results_dict['radial_average']) -
                                                 self.sub_profile)

        # This sums the intesity of the radial profiles before scaling
        self.collected_data['intensity_sum'] = np.nansum(self.results_dict['radial_average'])

        self.collected_data['min_bin_to_scale'] = self.results_dict['min_bin_to_scale']
        self.collected_data['max_bin_to_scale'] = self.results_dict['max_bin_to_scale']
        self.collected_data['num_bins'] = self.results_dict['num_bins']
        self.collected_data['std_dev_type'] = self.results_dict['std_dev_type']
        self.collected_data['N'] = self.results_dict['N']
        self.collected_data['intensity_threshold'] = self.results_dict['intensity_threshold']
        self.collected_data['profile_to_compare'] = self.results_dict['profile_to_compare']
        self.collected_data['detector_distance'] = self.results_dict['detector_distance']
        self.collected_data['beam_energy'] = self.results_dict['beam_energy']
        self.collected_data['optimized_geometry'] = self.optimized_geometry
        self.collected_data['timestamp'] = self.results_dict['timestamp']
        self.collected_data['pump_laser_state'] = self.results_dict['pump_laser_state']
        self.collected_data['pump_laser_on'] = self.results_dict['pump_laser_on']

        self.zmq_publish.send(b'ondadata', zmq.SNDMORE)
        self.zmq_publish.send_pyobj(self.collected_data)

        if self.num_events % self.speed_rep_int == 0:
            self.time = time.time()
            print('Processed: {0} in {1:.2f} seconds ({2:.2f} Hz)'.format(
                self.num_events,
                self.time - self.old_time,
                float(self.speed_rep_int)/float(self.time-self.old_time)))
            sys.stdout.flush()
            self.old_time = self.time

        return
