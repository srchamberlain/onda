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


###################################
# PIXEL SPACE / QSPACE CONVERSION #
###################################

class PixelSpaceQSpaceConversion:
    """Pixel space to q-space conversion.

    Implements conversion between pixel space and q-space.
    """
    
    def __init__(self, role, num_of_bins, coffset, dr):
        """ Initializes the QSpace bin calculator algorithm.
            
         Args:

            num_of_bins (int): number of pixel bins
            
            coffset (int) : detector distance's offset (coffset)

            dr (int): size of each bin in pixel
        """
    
    
        self.num_of_bins = num_of_bins
        self.nbins = numpy.arange(1, num_of_bins+1)
        self.theta = numpy.zeros([self.num_of_bins])
        self.qbins = numpy.zeros([self.num_of_bins])
        self.coffset = coffset
        self.dr = dr

    def convert_to_q(self, role, detector_distance, beam_energy):
        """QSpace bin conversion
            
        Calculates bins in QSpace associated with bins in pixel space

        Args:

            detector_distance (float): detector distance

            beam_energy (float): beam_energy

        Returns:

            self.qbins (numpy.ndarray float): array of q values associated with each radius bin
        """

        self.detector_distance = detector_distance
        self.beam_energy = beam_energy
        self.lambd = scipy.constants.h * scipy.constants.c /(scipy.constants.e * beam_energy)
        
        for i in self.nbins:
            self.theta[i-1] = 0.5*numpy.arctan((i*self.dr*0.00011)/(self.detector_distance*10e-4+self.coffset))

                                               
        for i in self.nbins:
            self.qbins[i-1] = 4*scipy.constants.pi*numpy.sin(self.theta[i-1])/(self.lambd/10e-11)
    
        return self.qbins



######################
# RADIUS OF GYRATION #
######################

class RadiusOfGyration:
    """Radius of Gyration calculation.
       
    Calculates radius of gyration for a radial intensity profile.
    """

    def __init__(self, role):
        """ Initializes the algorithm

        """

        pass
    

    def calculate_rog(self, role, qbins, radial_int, ignore_bins):
        """Calculates radius of gyration.

        Calculates radius of gyration for a radial profile.
        
        Args:

            qbins(numpy.ndarray): values of q associated with each bin in the radial profile
        
            radial_int (numpy.ndarray): radial profile of averaged intensities

        Returns:

            rg (float): radius of gyration

            i0: # TODO
        """
        
        q = qbins[ignore_bins:-1]
        rmax = scipy.constants.pi/(q[1]-q[0])

        r = numpy.linspace(0,self.rmax,(500-(ignore_bins+1)))
        
        qr = q[:,None]*r[None,:]
        P = (1/(2*scipy.constants.pi**2))*numpy.trapz(radial_int[ignore_bins:-1, None]*qr*numpy.sin(qr), q)

        topRg = numpy.trapz(P*r**2, r)
        botRg = 2*numpy.trapz(P, r)
        
        izero = 2*numpy.pi*botRg
        
        rg = numpy.sqrt(topRg/botRg)

        
        return rg, izero











