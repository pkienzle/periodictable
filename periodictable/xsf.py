# This program is public domain
r"""
This module has one class and nine fuctions.

   :class:`Xray`
       X-ray scattering properties for the elements.

The following attributes are added to each element:

   :func:`Xray.sftable`
       Three column table of energy vs. scattering factors f1, f2.

   :func:`Xray.scattering_factors`
       Returns f1, f2, the X-ray scattering factors for the given wavelengths
       interpolated from sftable.

   :func:`Xray.f0`
       Returns f0 for the given vector Q, with q_i in [0,24\pi] inv Ang.

   :func:`Xray.sld`
       Returns scattering length density (*real*, *imaginary*) for the
       given wavelengths or energies.

The following functions are available for X-ray scatting information processing:

     :func:`xray_wavelength`
         Finds X-ray wavelength in angstroms given energy in keV.

     :func:`xray_energy`
         Finds X-ray energy in keV given wavelength in angstroms.
         
     :func:`init`
         Initializes a periodic table with the Lawrence Berkeley Laboratory
         Center for X-Ray Optics xray scattering factors.

     :func:`init_spectral_lines`
         Sets the K_alpha and K_beta1 wavelengths for select elements.

     :func:`sld_table`
         Prints the xray SLD table for the given wavelength.

     :func:`xray_sld`
         Computes xray scattering length densities for molecules.

     :func:`xray_sld_from_atoms`
         The underlying scattering length density calculator. This works with
         a dictionary of atoms and quanties directly.

     :func:`emission_table`
         Prints a table of emission lines.

K_alpha, K_beta1 (Angstrom):
    X-ray emission lines for various elements, including Ag, Pd, Rh, Mo,
    Zn, Cu, Ni, Co, Fe, Mn, Cr and Ti. K_alpha is the average of
    K_alpha1 and K_alpha2 lines.

X-ray scattering factors:
    Low-Energy X-ray Interaction Coefficients: Photoabsorption, scattering
    and reflection for E in 30 to 30,000 eV, and Z in 1 to 92.

.. Note::

    For :ref:`custom tables <custom-table>`, use :func:`init` and 
    :func:`init_spectral_lines` to set the data.


X-ray f1 and f2 tables
======================
The data for the tables is stored in the ``periodictable/xsf``. 
directory.  The following information is from ``periodictable/xsf/read.me``, 
with minor formatting changes.
These ``[*.nff]`` files were used to generate the tables published in 
reference [#Henke1993]_. The files contain three columns of data::

    Energy(eV), f_1, f_2,

where *f_1* and *f_2* are the atomic (forward) scattering factors.
There are 500+ points on a uniform logarithmic mesh with points
added 0.1 eV above and below "sharp" absorption edges. The
tabulated values of *f_1* contain a relativistic, energy 
independent, correction given by::

    Z* = Z - (Z/82.5)^(2.37)

.. Note::
    Below 29 eV *f_1* is set equal to -9999.

The atomic photoabsorption cross section, mu_a, may be readily obtained
from the values of *f_2* using the relation::

    mu_a = 2*r_0*lambda*f_2

where *r_0* is the classical electron radius, and lambda is the wavelength.
The index of refraction for a material with *N* atoms per unit volume
is calculated by::

    n = 1 - N*r_0*(lambda)^2*(f_1+if_2)/(2*pi).

These (semi-empirical) atomic scattering factors are based upon
photoabsorption measurements of elements in their elemental state.
The basic assumption is that condensed matter may be modeled as a
collection of non-interacting atoms.  This assumption is in general
a good one for energies sufficiently far from absorption thresholds.
In the threshold regions, the specific chemical state is important
and direct experimental measurements must be made.

These tables are based on a compilation of the available experimental
measurements and theoretical calculations.  For many elements there is
little or no published data and in such cases it was necessary to
rely on theoretical calculations and interpolations across Z.
In order to improve the accuracy in the future considerably more
experimental measurements are needed.

Please send any comments about the tables to EMGullikson@lbl.gov.

.. table:: Note that the following elements have been updated since the 
           publication of Ref. [#Henke1993]_ in July 1993.

           ========  ==========   ==============   
           Element   Updated      Energy Range       
           ========  ==========   ==============   
           Mg        1/15/94      30-50 eV
           Al        1/15/94      30-73 eV
           Si        1/15/94      30-100 eV
           Au        11/7/94      2000-6500 eV
           Li        11/15/94     2000-30000 eV
           Si        6/95         30-500 eV
           Fe        10/95        600-800 eV
           Mo        11/97        10-930 eV
           Be        8/04         40-250 eV
           Mo        8/04         25-60 eV
           W         8/04         35-250 eV
           Ru        8/04         40-1300 eV
           Ti        8/04         20-150 eV
           Sc        4/06         50-1300 eV
           Gd        6/07         12-450 eV
           La        6/07         14-440 eV
           ========  ==========   ==============  
          

Data available at:
  #. http://henke.lbl.gov/optical_constants/asf.html
  #. http://henke.lbl.gov/optical_constants/update.html

.. [#Henke1993] B. L. Henke, E. M. Gullikson, and J. C. Davis.  "X-ray interactions:
       photoabsorption, scattering, transmission, and reflection at E=50-30000 eV,
       Z=1-92", Atomic Data and Nuclear Data Tables 54 no.2, 181-342 (July 1993).

"""
from __future__ import with_statement
__all__ = ['Xray', 'init', 'init_spectral_lines',
           'xray_energy','xray_wavelength',
           'xray_sld','xray_sld_from_atoms',
           'emission_table','sld_table','plot_xsf',
           ]
import os.path
import glob
import numpy
from numpy import nan

from . import core
from .core import Element, Ion, default_table
from .constants import (avogadro_number, plancks_constant, speed_of_light,
                        electron_radius)

import logging
from distutils.filelist import findall
def xray_wavelength(energy):
    """
    Convert X-ray energy to wavelength.

    :Parameters:
        *energy* : float or vector | keV

    :Returns:

    :Algorithm:

        Use the formula:

            lambda = h c / E

        where:

            h = planck's constant in eV s
            c = speed of light in m/s
    """
    return plancks_constant*speed_of_light/numpy.asarray(energy)*1e7

def xray_energy(wavelength):
    """
    Convert X-ray wavelength to energy.

    :Parameters:
        *wavelength* : float or vector | A

    :Returns:
        *energy* : float or vector | keV

    """
    return plancks_constant*speed_of_light/numpy.asarray(wavelength)*1e7

class Xray(object):
    """
    X-ray scattering properties for the elements. Refer help(periodictable.xsf)
    from command prompt for details.
    """
    _nff_path = core.get_data_path('xsf')
    sftable_units = ["eV","",""]
    scattering_factors_units = ["",""]
    sld_units = ["10^-6 inv A^2","10^-6 invA^2"]
    _table = None
    def __init__(self, element):
        self.element = element

    def _gettable(self):
        if self._table is None:
            # Load table when necessary; note there is no table for
            # neutrons (n), and lowercase nitrogen=> n.nff, so it must
            # be checked for explicitly.
            filename = os.path.join(self._nff_path,
                                    self.element.symbol.lower()+".nff")
            if self.element.symbol != 'n' and os.path.exists(filename):
                xsf = numpy.loadtxt(filename,skiprows=1).T
                xsf[1,xsf[1]==-9999.] = numpy.NaN
                xsf[0] *= 0.001  # Use keV in table rather than eV
                self._table = xsf
        return self._table
    sftable = property(_gettable,doc="X-ray scattering factor table (E,f1,f2)")

    def scattering_factors(self, energy):
        """
        X-ray scattering factors f',f''.

        :Parameters: 
            *energy* : float or vector | keV
                X-ray energy.
               
        :Returns:
            *scattering_factors* : (float, float)
                Values outside the range return NaN.

        :Algorithm:

            Linear interpolation within the Henke Xray scattering factors 
            database at the Lawrence Berkeley Laboratory Center for X-ray 
            Optics.
        """
        xsf = self.sftable
        if xsf is None:
            return None,None

        scalar = numpy.isscalar(energy)
        if scalar:
            energy = numpy.array([energy])
        f1 = numpy.interp(energy,xsf[0],xsf[1],left=nan,right=nan)
        f2 = numpy.interp(energy,xsf[0],xsf[2],left=nan,right=nan)
        if scalar:
            f1,f2 = f1[0],f2[0]
        return f1,f2

    def f0(self, Q):
        r"""
        Isotropic X-ray scattering factors f0 for the input Q.

        :Parameters: 
            *Q* : float or vector in [0, 24*pi] | inv A
                X-ray scattering properties for the elements.
               
        :Returns:
            *f0* : float
                Values outside the valid range return NaN.
      

        .. Note:: f0 is often given as a function of sin(theta)/lambda
                  whereas we are using  Q = 4*pi*sin(theta)/lambda, or
                  in terms of energy Q = 4*pi*sin(theta)*E/(h c).
        

        Reference:
             D. Wassmaier, A. Kerfel, Acta Crystallogr. A51 (1995) 416.
             http://dx.doi.org/10.1107/S0108767394013292
        """
        from . import cromermann
        f = cromermann.fxrayatq(Q=Q, 
                                symbol=self.element.symbol,
                                charge=self.element.charge)
        return f

    def sld(self, wavelength=None, energy=None):
        """
        X ray scattering length density.

        :Parameters: 

            *wavelength* : float or vector | A
                Wavelength of the X-ray.

            *energy* : float or vector | keV
                Energy of the X-ray (if *wavelength* not specified).

            .. note:
                Only one of *wavelength* and *energy* is needed.
               
        :Returns:
            *sld* : (float, float) | inv A^2
                (*real*, *imaginary*) X-ray scattering length density.

        :Raises:
            *TypeError* : neither *wavelength* nor *energy* was specified.

        :Algorithm:
            The element SLD is r_eN(f1+1jf2), where *r_e* is the electron
            radius and *N* is number density = density/mass * Avogadro's Number.

            The constants are available directly::
     
                *r_e* = periodictable.xsf.electron_radius
                *N_A* = periodictable.constants.avogadro_number

        Data comes from the Henke Xray scattering factors database at the
        Lawrence Berkeley Laboratory Center for X-ray Optics.
        """
        if wavelength is not None:
            energy = xray_energy(wavelength)
        if energy is None:
            raise TypeError('X-ray SLD needs wavelength or energy')
        f1,f2 = self.scattering_factors(energy)
        if f1 is None or self.element.number_density is None:
            return None,None
        rho = f1*electron_radius*self.element.number_density*1e-8
        irho = f2*electron_radius*self.element.number_density*1e-8
        return rho,irho

# Note: docs and function prototype are reproduced in __init__
def xray_sld(compound, density=None, 
             wavelength=None, energy=None, natural_density=None):
    """
    Compute xray scattering length densities for molecules. 

    :Parameters: 
        *compound* : Formula initializer
            Chemical formula initializer.
        *density* : float | g/cm^3
            Density of the compound. 
        *wavelength* : float | A
            Wavelength of the X-ray.
        *energy* : float | keV
            Energy of the X-ray, if *wavelength* is not specified.

    :Returns:
        *sld* : (float, float) | 10^-6 inv A^2
            (*real*, *imaginary*) scattering length density. 

    :Raises:
        *AssertionError* :  *density* or *wavelength*/*energy* is missing.
    """
    import formulas
    compound = formulas.Formula(compound)
    if density is None:
        if natural_density is not None: 
            density = natural_density/compound.natural_mass_ratio()
        else:
            density = compound.density # defaults to molecule density
    assert density is not None, "scattering calculation needs density"

    if wavelength is not None: energy = xray_energy(wavelength)
    assert energy is not None, "scattering calculation needs energy or wavelength"

    mass,sum_f1,sum_f2 = 0,0,0
    for element,quantity in compound.atoms.iteritems():
        mass += element.mass*quantity
        f1,f2 = element.xray.scattering_factors(energy)
        #print element,f1,f2,wavelength
        sum_f1 += f1*quantity
        sum_f2 += f2*quantity
    if mass == 0: # because the formula is empty
        rho,irho = 0,0
    else:
        N = (density/mass*avogadro_number*1e-8)
        rho = N*sum_f1*electron_radius
        irho = N*sum_f2*electron_radius
    return rho,irho

def xray_sld_from_atoms(*args, **kw):
    """
    .. deprecated:: 0.91
        :func:`xray_sld` now accepts dictionaries of {atom: count} directly.
    """
    return xray_sld(*args, **kw)


def init_spectral_lines(table):
    """
    Sets the K_alpha and K_beta1 wavelengths for select elements
    """
    Element.K_alpha_units = "angstrom"
    Element.K_beta1_units = "angstrom"
    table.Ag.K_alpha = 0.5608
    table.Ag.K_beta1 = 0.4970
    table.Pd.K_alpha = 0.5869
    table.Pd.K_beta1 = 0.5205
    table.Rh.K_alpha = 0.6147
    table.Rh.K_beta1 = 0.5456
    table.Mo.K_alpha = 0.7107
    table.Mo.K_beta1 = 0.6323
    table.Zn.K_alpha = 1.4364
    table.Zn.K_beta1 = 1.2952
    table.Cu.K_alpha = 1.5418
    table.Cu.K_beta1 = 1.3922
    table.Ni.K_alpha = 1.6591
    table.Ni.K_beta1 = 1.5001
    table.Co.K_alpha = 1.7905
    table.Co.K_beta1 = 1.6208
    table.Fe.K_alpha = 1.9373
    table.Fe.K_beta1 = 1.7565
    table.Mn.K_alpha = 2.1031
    table.Mn.K_beta1 = 1.9102
    table.Cr.K_alpha = 2.2909
    table.Cr.K_beta1 = 2.0848
    table.Ti.K_alpha = 2.7496
    table.Ti.K_beta1 = 2.5138

def init(table, reload=False):
    
    if 'xray' in table.properties and not reload: return
    table.properties.append('xray')
    
    # Create an xray object for the particular element/ion.  Note that
    # we must not use normal attribute tests such as "hasattr(el,'attr')" 
    # or "try: el.attr; except:" since the delegation methods on Ion will
    # just return the attribute from the base element.  Instead we check
    # for an instance specific xray object for the particular ion prior
    # to delegating.
    # TODO: is there a better way to set up delegation on a field by
    # field basis?
    def _cache_xray(el): 
        if '_xray' not in el.__dict__ and isinstance(el, (Element,Ion)):
            el._xray = Xray(el)
        return el._xray
    Element.xray = property(_cache_xray)
    Ion.xray = property(_cache_xray)
    
    ## Note: the simple approach below fails for e.g., Ni[58].ion[3].xray
    #for el in table:
    #    for charge in el.ions:
    #        el.ion[charge].xray = Xray(el.ion[charge])
    #    el.xray = Xray(el)

def plot_xsf(el):
    """
    Plots the xray scattering factors for the given element.

    :Parameters: 
        *el* : Element
               
    :Returns: None
    """
    import pylab
    xsf = el.xray.sftable
    pylab.title('X-ray scattering factors for '+el.name)
    pylab.plot(xsf[0],xsf[1])
    pylab.plot(xsf[0],xsf[2])
    pylab.xlabel('Energy (keV)')
    pylab.ylabel('Scattering factor')
    pylab.legend(['f1','f2'])
    pylab.show()

def sld_table(wavelength, table=None):
    """
    Prints the xray SLD table for the given wavelength.

    :Parameters: 
        *wavelength* : float | A
            X-ray wavelength.
        *table* : PeriodicTable
            The default periodictable unless a specific table has been requested.
               
    :Returns: None
    """
    table = default_table(table)
    ## Cu K-alpha Zeff and dZ table
    #for el in table:
    #    f1,f2 = el.xsf(table.Cu.K_alpha)
    #    if f1 is not None:
    #        print el.Z,el.symbol,"%.1f"%f1,"%.4f"%(f1-el.Z)

    # NBCU spreadsheet format
    print "X-ray scattering length density for",wavelength,"A"
    print "%3s %6s %6s"%('El','rho','irho')
    for el in table:
        rho,irho = el.xray.sld(table.Cu.K_alpha)
        if rho is not None:
            print "%3s %6.2f %6.2f"%(el.symbol,rho,irho)

def emission_table(table=None):
    """
    Prints a table of emission lines.

    :Parameters: 
        *table* : PeriodicTable.
            The default periodictable unless a specific table has been requested.
               
    :Returns: None
    """
    table = default_table(table)
    print "%3s %7s %7s"%('El','Kalpha','Kbeta1')
    for el in table:
        if hasattr(el,'K_alpha'):
            print "%3s %7.4f %7.4f"%(el.symbol,el.K_alpha,el.K_beta1)


