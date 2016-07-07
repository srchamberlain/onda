"""
    This file is part of OnDA.

    OnDA is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OnDA is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with OnDA.  If not, see <http://www.gnu.org/licenses/>.
    
    Added July 7, 2016
    By Sarah Chamberlain, BioXFEL contribution
"""
import numpy
import scipy.ndimage
import scipy.constants
import math


#################
# RPIXBINS CALC #
#################
class RpixbinsCalc:
    """
        Calculates rpixbins for doing radial intensity averaging in pixel space
        
    """

    def __init__(self, role, pixelmap_radius):
        """Initializes rpixbins algorithm
            
        Args: 
            pixelmap_radius (numpy.ndarray): array of radius values calculated using the raw data image
            
        """
    
        self.rpixbins = numpy.zeros(pixelmap_radius.shape, dtype=int)
        self.pixelmap_radius = pixelmap_radius
    
    def rpixbins_calc(self, nbins):
        """
        rpixbins for radial intensity averaging
        
        Args:
            nbins (int): number if bins wanted by the user
            
        Returns:
            self.rpixbins (numpy.ndarray): array of labels cooresponding to each value in pixelmap_radius for radial intensity function
            
            dr (int): length in pixels of one rpixbin, needed for qSpace algorithm later
        """
        #rpixbins = np.zeros(pixelmap_radius.shape, dtype=int)
        #nbins = 500
        dr = float(numpy.max(self.pixelmap_radius))/nbins
    
        for i in range(0, nbins-1):
            self.rpixbins[(self.pixelmap_radius>=i*dr) & (self.pixelmap_radius<(i+1)*dr)] = i
        
            self.rpixbins[self.pixelmap_radius>=(nbins)*dr] = nbins
    
    
        return self.rpixbins, dr



##################
# RADIAL CUM AVG #
##################
class RadialCumAvg:
    """
        
        Master algorithm, without end function
        
        Accumulates raw data images from the detector data in slab format, and when the required number of shots have
        been collected, puts the average in the the collected_data dict.
        
        Required monitor ini parameters:
        
        accumulated_shots:        the number of accumulated shots
        
        """
    
    def __init__(self, role, monitor_params):
        
        self.accumulated_shots = monitor_params['RadialAveraging']['accumulated_shots']
        
        # Initialized on master
        if role == 'master':
            
            self.apply = self.master
            
            self.num_radial_int = 0
            self.sum_radial_int = numpy.zeros((500))

    def master(self, rawavg, rpixbins):
        """
        Args: 
            radial_int: radial intensity to be added to the accumulator
        
        Returns: 
            cum_radial_average (numpy.ndarray): array of averaged intesities calculated using Averaged Raw Data as input
        """

        """
        Averages radial intensities
        """
        
        ###This uses avg raw data for cumulative radial average
        cum_radial_average = scipy.ndimage.mean(rawavg, labels=rpixbins, index=numpy.arange(0, rpixbins.max()))
    
    
        return cum_radial_average
        
        
        ###This part below does the cumulative radial average one at a time
        """if self.num_radial_int == self.accumulated_shots:
            self.num_radial_int = 0
            self.sum_radial_int.fill(0)
                
        self.sum_radial_int += (rawavg / self.accumulated_shots)
        self.num_radial_int += 1
        
        if self.num_radial_int == self.accumulated_shots:
            return self.sum_radial_int
        return None"""

###########
# SCALING #
###########
class ScaleRadial:
    """
       Scaling function, scales profile to region defined by user
       
       Inputs:  self.min_rpixbin: Minimum bin number for scaling region, readi in from monitor params
       
                self.max_rpixbin: Maximum bin number for scaling region, readi in from monitor params
    """

    def __init__(self, role, min_rpixbin, max_rpixbin):
    
        self.min_rpixbin = min_rpixbin
        self.max_rpixbin = max_rpixbin
        self.region = self.max_rpixbin - self.min_rpixbin
        self.region_int = numpy.zeros((self.region))

    def scale_profile(self, radial_int):
        """
        Args: 
            radial_int (numpy.ndarray): radial profile of intensities to scale
           
        Returns:
            self.radial_int_new (numpy.ndarray): array of intensity values scaled using the region specified by the user
        
        """

        self.region_int = radial_int[self.min_rpixbin:self.max_rpixbin:1]
        self.average = numpy.average(self.region_int)
        self.radial_int_new = radial_int/(self.average)

        return self.radial_int_new


###########
# Q-SPACE #
###########

class QSpace:
    """
       Converts Pixel space bins to q-space
       

    """
    
    def __init__(self, role, num_of_bins, coffset, dr):
        """ Initializes QSpace bin calculator algorithm
            
            Args:
            num_of_bins (int): number of pixel bins
            
            coffset (int) : coffset read in process_collect_soultions.py from the geometry file
                
            dr (int): defined in rbins calculator, number of pixels long a rpixbin is
        """
    
    
        self.num_of_bins = num_of_bins
        self.nbins = numpy.arange(1, num_of_bins+1)
        self.theta = numpy.zeros([self.num_of_bins])
        self.qbins = numpy.zeros([self.num_of_bins])
        self.coffset = coffset
        self.dr = dr

    def convert_to_q(self, role, detector_distance, beam_energy):
        """qbin calculator
            
        Calculates Qbins associated with rpixbins
           
        Args:
            detector_distance: distance from solution to detector, read in process_collect_soultions.py, from the geometry file
           
            beam_energy: beam_energy, read in process_collect_soultions.py from the geometry file
            
        Returns:
            self.qbins (numpy.ndarray): array of q values associated with each averaged intensity in the radial profile
        """

        self.detector_distance = detector_distance
        self.beam_energy = beam_energy
        self.lambd = scipy.constants.h * scipy.constants.c /(scipy.constants.e * beam_energy)
        
        for i in self.nbins:
            self.theta[i-1] = 0.5*numpy.arctan((i*self.dr*0.00011)/(self.detector_distance*10e-4+self.coffset))
        #print self.theta
        ##theta is a list of angles
                                               
        for i in self.nbins:
            self.qbins[i-1] = 4*scipy.constants.pi*numpy.sin(self.theta[i-1])/(self.lambd/10e-11)
    
        return self.qbins
        #return qbins
        ### qbins should be a 500 item long column array of q values ###


####################
# RADIUSOFGYRATION #
####################

class RadiusG:
    """
       Calculate Radius of Gyration
       
       for each radial intensity profile, calculates Rg
       
       ***Currently not called in process_collect_solutions because not finished/optimized for time*****
    """

    def __init__(self, role):
        """ Initializes the Rg algorithm

        
        """
        pass
    

    def RgCalculator(self, role, qbins, radial_int, ignore_bins):
        """ Radius of Gyration Algorithm
        
        Args:
            qbins(numpy.ndarray): values of q associated with each bin of intensities
        
            radial_int (numpy.ndarray): values of averaged intensities for defined bins
        
        
        """
        q = qbins[20:]
        self.rmax = scipy.constants.pi/(qbins[1]-qbins[0])

        r = numpy.linspace(0,self.rmax,(500-ignore_bins))
        
        qr = q[:,None]*r[None,:]
        P = (1/(2*scipy.constants.pi**2))*numpy.trapz(qr*numpy.sin(qr), q)
        
        topRg = numpy.trapz(P*r**2, r)
        botRg = 2*numpy.trapz(P, r)
        
        self.Rg = numpy.sqrt(topRg/botRg)

        
        return self.Rg











