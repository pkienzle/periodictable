# This program is public domain
"""
X-ray scatting information.

    xray.sftable (eV,?,?)
        Three column table of energy vs. scattering factors f1,f2.
    xray.scattering_factors(wavelength)
        Returns f1,f2, the X-ray scattering factors for the given wavelengths
        interploted from sftable.
    xray.sld(wavelength)
        Returns scattering length density and absorption for the
        given wavelengths.

    K_alpha, K_beta1 (Angstrom)

        X-ray emission lines for various elements:
           Ag, Pd, Rh, Mo, Zn, Cu, Ni, Co, Fe, Mn, Cr and Ti
        K_alpha is the average of K_alpha1 and K_alpha2 lines.
   

K-alpha and K-beta1 lines:

   D.C. Creagh and J.H. Hubbell (1999), "4.2.4 X-ray absorption 
   (or attenuation) coefficients" in A.J.C. Wilson, E.Prince Eds., 
   International Tables for Crystallography Vol C p220-236.  
   Kluwer Academic Publishers, Dordrecht, The Netherlands.

X-ray scattering factors:

                Low-Energy X-ray Interaction Coefficients:
                Photoabsorption, Scattering, and Reflection
                        E = 30-30,000 eV, Z = 1-92

                B. L. Henke, E. M. Gullikson, and J. C. Davis
                        Center for X-Ray Optics, 2-400
                        Lawrence Berkeley Laboratory
                        Berkeley, California 94720

These files were used to generate the tables published in reference [1].
The files contain three columns of data: Energy(eV), f_1, f_2,
where f_1 and f_2 are the atomic (forward) scattering factors.
There are 500+ points on a uniform logarithmic mesh with points
added 0.1 eV above and below "sharp" absorption edges.
(Note: below 29 eV f_1 is set equal to -9999.) 
The tabulated values of f_1 contain a relativistic, energy independent,
correction given by, Z* = Z - (Z/82.5)^(2.37).


The atomic photoabsorption cross section, mu_a, may be readily obtained
from the values of f_2 using the relation,

                        mu_a = 2*r_0*lambda*f_2 

where r_0 is the classical electron radius, and lambda is the wavelength.
The index of refraction for a material with N atoms per unit volume
is calculated by,

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

Please send any comments to EMGullikson@lbl.gov.

Note that the following elements have been updated since the publication
of Ref. [1] in July 1993.

Element         Updated         Energy Range
Mg              1/15/94         30-50 eV
Al              1/15/94         30-73 eV
Si              1/15/94         30-100 eV
Au              11/7/94         2000-6500 eV
Li              11/15/94        2000-30000 eV
Si                6/95          30-500 eV
Fe               10/95          600-800 eV
Mo               11/97          10-930 eV
Be                8/04          40-250 eV
Mo                8/04          25-60 eV
W                 8/04          35-250 eV
Ru                8/04          40-1300 eV
Ti                8/04          20-150 eV
Sc                4/06          50-1300 eV
Gd                6/07          12-450 eV
La                6/07          14-440 eV

Data available at: 
  http://henke.lbl.gov/optical_constants/asf.html
  http://henke.lbl.gov/optical_constants/update.html

[1] B.L. Henke, E.M. Gullikson, and J.C. Davis.  "X-ray interactions: 
photoabsorption, scattering, transmission, and reflection at E=50-30000 eV, 
Z=1-92", Atomic Data and Nuclear Data Tables 54 no.2, 181-342 (July 1993).

[2] J.H. Hubbell, W.J. Veigele, E.A. Briggs, R.T. Brown, D.T. Cromer, 
and R.J. Howerton. "Atomic Form Factors, Incoherent Scattering Functions, 
and Photon Scattering Cross Sections",  J. Phys. Chem. Ref. Data 4, 
471-538 (1975); erratum in 6, 615-616 (1977).

"""
from __future__ import with_statement
import os.path
import numpy

import elements
from elements import periodic_table
from density import avogadro_number
import mass,density

# Need X-ray energy/wavelength conversions to lookup wavelengths in scattering
# function tables.
plancks_constant = 4.13566733e-15 #(10) eV s
speed_of_light = 299792458 # m/s (exact)
electron_radius = 2.8179402894e-15 #(58) m  

def xray_wavelength(energy):
    """
    Find X-ray wavelength in angstroms given energy in keV.

    lambda = h c / E

    h = planck's constant in eV s
    c = speed of light in m/s
    """
    return plancks_constant*speed_of_light/numpy.asarray(energy)*1e7

def xray_energy(wavelength):
    """
    Find X-ray energy in keV given wavelength in angstroms.
    """
    return plancks_constant*speed_of_light/numpy.asarray(wavelength)*1e7


class Xray(object):
    """
    X-ray scattering properties for the elements.

    Center for X-Ray Optics.
    L. Henke, E. M. Gullikson, and J. C. Davis

    K_alpha and K_beta1 emission lines for selected elements.

    D. C. Creagh and J. H. Hubbell
    International Tables for Crystallography Volume C.  
    
    See help(elements.xsf) for details.
    """
    sftable_units = ["eV","",""]
    scattering_factors_units = ["",""]
    sld_units = ["Nb","Nb"]
    _table = None
    def __init__(self, element):
       self._symbol = element.symbol
       self._number_density = element.number_density

    def _gettable(self):
        if self._table is None:
            # Load table when necessary; note there is no table for
            # neutrons (n), and lowercase nitrogen=> n.nff, so it must
            # be checked for explicitly.
            filename = os.path.join(os.path.dirname(__file__),'xsf',
                                    self._symbol.lower()+".nff")
            if self._symbol != 'n' and os.path.exists(filename):
                xsf = numpy.loadtxt(filename,skiprows=1).T
                xsf[1,xsf[1]==-9999.] = numpy.NaN
                xsf[0] *= 0.001  # Use keV in table rather than eV
                self._table = xsf
        return self._table
    sftable = property(_gettable,doc="X-ray scattering factor table (E,f1,f2)")

    def scattering_factors(self, energy):
        """
        Return the X-ray scattering factors f1,f2 for the given energy (keV).
        Energy can be a scalar or a vector.

        Data comes from the Henke Xray scattering factors database at the
        Lawrence Berkeley Laboratory Center for X-ray Optics.  
        """
        xsf = self.sftable
        if xsf is None:
            return None,None

        scalar = numpy.isscalar(energy)
        if scalar:
            energy = numpy.array([energy])
        f1 = numpy.interp(energy,xsf[0],xsf[1])
        f2 = numpy.interp(energy,xsf[0],xsf[2])
        if scalar:
            f1,f2 = f1[0],f2[0]
        return f1,f2

    def sld(self, wavelength=None, energy=None):
        """
        Return the X ray scattering length density and absorption in
        units of 10**-6 angstrom**-2.

        If wavelength is given it is converted to energy, ignoring whatever
        energy was specified.  Wavelength is in angstroms. Energy is in keV.
        Energy/wavelength can be a scalar or vector.

        The element SLD is r_e * N * (f1 + 1j*f2), where r_e is the electron
        radius and N is number density = density/mass * Avogadro's Number.

        The constants are available directly::
            elements.xsf.electron_radius
            elements.density.avogadro_number

        Data comes from the Henke Xray scattering factors database at the
        Lawrence Berkeley Laboratory Center for X-ray Optics.
        
        Raises TypeError if neither wavelength or energy is specified.
        """
        if wavelength is not None: 
            energy = xray_energy(wavelength)
        if energy is None: 
            raise TypeError('X-ray SLD needs wavelength or energy')
        f1,f2 = self.scattering_factors(energy)
        if f1 is None or self._number_density is None:
            return None,None
        scattering = f1*electron_radius*self._number_density*1e-8
        absorption = f2*electron_radius*self._number_density*1e-8
        return scattering,absorption

# Note: docs and function prototype are reproduced in __init__
def xray_sld(input,density=None,wavelength=None,energy=None):
    """
    Compute xray scattering length densities for molecules.
    Returns the scattering length density and absorption
    in units of 10**-6 angstrom**-2.
    
    Wavelength is in angstroms.  Energy is in keV.
    
    Raises AssertionError if density or wavelength/energy is not provided
    """
    import molecules
    return xray_sld_from_atoms(molecules.Molecule(input).atoms,
                               density=density,wavelength=wavelength,
                               energy=energy)

def xray_sld_from_atoms(atoms,density=None,wavelength=None,energy=None):
    """
    The underlying scattering length density calculator.  This
    works with a dictionary of atoms and quanties directly, such
    as returned by molecule.atoms.
    
    Wavelength is in angstroms.  Energy is in keV.
    
    Raises AssertionError if density or wavelength/energy is not provided
    """

    if wavelength is not None: energy = xray_energy(wavelength)
    assert density is not None, "xray_sld needs density"
    assert energy is not None, "xray_sld needs energy or wavelength"

    mass,scattering,absorption = 0,0,0
    for element,quantity in atoms.iteritems():
        mass += element.mass*quantity
        f1,f2 = element.xray.scattering_factors(energy)
        #print element,f1,f2,wavelength
        scattering += f1*quantity
        absorption += f2*quantity
    N = (density/mass*avogadro_number*1e-8)
    rho = N*scattering*electron_radius
    mu = N*absorption*electron_radius
    #print "ro",electron_radius
    #print "scattering",scattering,"absorption",absorption,"1/N",1/N
    return rho,mu



def _init():
    if 'xray' in periodic_table.properties: return
    periodic_table.properties.append('xray')
    for el in periodic_table:
        el.xray = Xray(el)
    elements.Element.K_alpha_units = "angstrom"
    elements.Element.K_beta1_units = "angstrom"
        
    # Sets the K_alpha and K_beta1 wavelengths for select elements
    periodic_table.Ag.K_alpha = 0.5608
    periodic_table.Ag.K_beta1 = 0.4970
    periodic_table.Pd.K_alpha = 0.5869
    periodic_table.Pd.K_beta1 = 0.5205
    periodic_table.Rh.K_alpha = 0.6147
    periodic_table.Rh.K_beta1 = 0.5456
    periodic_table.Mo.K_alpha = 0.7107
    periodic_table.Mo.K_beta1 = 0.6323
    periodic_table.Zn.K_alpha = 1.4364
    periodic_table.Zn.K_beta1 = 1.2952
    periodic_table.Cu.K_alpha = 1.5418
    periodic_table.Cu.K_beta1 = 1.3922
    periodic_table.Ni.K_alpha = 1.6591
    periodic_table.Ni.K_beta1 = 1.5001
    periodic_table.Co.K_alpha = 1.7905
    periodic_table.Co.K_beta1 = 1.6208
    periodic_table.Fe.K_alpha = 1.9373
    periodic_table.Fe.K_beta1 = 1.7565
    periodic_table.Mn.K_alpha = 2.1031
    periodic_table.Mn.K_beta1 = 1.9102
    periodic_table.Cr.K_alpha = 2.2909
    periodic_table.Cr.K_beta1 = 2.0848
    periodic_table.Ti.K_alpha = 2.7496
    periodic_table.Ti.K_beta1 = 2.5138

def plot_xsf(el):
    import pylab
    xsf = el.xray.sftable
    pylab.title('X-ray scattering factors for '+el.name)
    pylab.plot(xsf[0],xsf[1])
    pylab.plot(xsf[0],xsf[2])
    pylab.xlabel('Energy (keV)')
    pylab.ylabel('Scattering factor')
    pylab.legend(['f1','f2'])
    pylab.show()

def sld_table(wavelength):
    ## Cu K-alpha Zeff and dZ table
    #for el in periodic_table:
    #    f1,f2 = el.xsf(periodic_table.Cu.K_alpha)
    #    if f1 is not None:
    #        print el.Z,el.symbol,"%.1f"%f1,"%.4f"%(f1-el.Z)

    # NBCU spreadsheet format
    print "X-ray scattering length density and absorption for",wavelength,"A"
    print "%3s %6s %6s"%('El','rho','mu')
    for el in periodic_table:
        rho,mu = el.xray.sld(periodic_table.Cu.K_alpha)
        if rho is not None:
            print "%3s %6.2f %6.2f"%(el.symbol,rho,mu)

def emission_table():
    """
    Print a table of emission lines.
    """
    print "%3s %7s %7s"%('El','Kalpha','Kbeta1')
    for el in periodic_table:
        if hasattr(el,'K_alpha'):
            print "%3s %7.4f %7.4f"%(el.symbol,el.K_alpha,el.K_beta1)

_init()

def test():
    from molecules import Molecule
    from elements import Cu,Mo,Ni,Fe,Si

    # Check some K_alpha and K_beta1 values
    assert Cu.K_alpha == 1.5418
    assert Cu.K_beta1 == 1.3922
    
    # Check scalar scattering factor lookup
    f1,f2 = Ni.xray.scattering_factors(xray_energy(Cu.K_alpha))
    assert abs(f1-25.0229)<0.0001
    assert abs(f2-0.5249)<0.0001

    # Check array scattering factor lookup
    f1,f2 = Ni.xray.scattering_factors(xray_energy(Cu.K_alpha))
    m1,m2 = Ni.xray.scattering_factors(xray_energy(Mo.K_alpha))
    B1,B2 = Ni.xray.scattering_factors(xray_energy([Cu.K_alpha,Mo.K_alpha]))
    assert (B1==[f1,m1]).all() and (B2==[f2,m2]).all()

    # Check that we can lookup sld by wavelength and energy
    Fe_rho,Fe_mu = Fe.xray.sld(wavelength=Cu.K_alpha)
    assert abs(Fe_rho-59.45) < 0.01
    Si_rho,Si_mu = Si.xray.sld(energy=8.050)
    assert abs(Si_rho-20.0701) < 0.0001
    assert abs(Si_mu-0.4572) < 0.0001

    # Check that wavelength is the default
    Fe_rho_default,Fe_mu_default = Fe.xray.sld(Cu.K_alpha)
    assert Fe_rho == Fe_rho_default and Fe_mu == Fe_mu_default

    # Check array form of sld lookup
    f1,f2 = Si.xray.sld(wavelength=Cu.K_alpha)
    m1,m2 = Si.xray.sld(wavelength=Mo.K_alpha)
    B1,B2 = Si.xray.sld(wavelength=[Cu.K_alpha,Mo.K_alpha])
    assert (B1==[f1,m1]).all() and (B2==[f2,m2]).all()

    #print Cu.xray.sftable
    #plot_xsf(Cu)
    #emission_table()
    #sld_table(table,Cu.K_alpha)

    """
    # Table of scattering length densities for various molecules
    for molecule,density in [('SiO2',2.2),('B4C',2.52)]:
        atoms = Molecule(molecule).atoms
        rho,mu = xray_sld(atoms,density,wavelength=Cu.K_alpha)
        print "sld for %s(%g g/cm**3)  rho=%.4g mu=%.4g"\
            %(molecule,density,rho,mu)
    """

    # Cross check against mo
    rho,mu = xray_sld_from_atoms({Si:1},density=Si.density,wavelength=1.54)
    rhoSi,muSi = Si.xray.sld(wavelength=1.54)
    assert abs(rho - rhoSi) < 1e-14
    assert abs(mu - muSi) < 1e-14

    # Check that neutron_sld works as expected
    atoms = Molecule('SiO2').atoms
    rho,mu = xray_sld_from_atoms(atoms,2.2,energy=xray_energy(Cu.K_alpha))
    assert abs(rho-18.87)<0.1
    atoms = Molecule('B4C').atoms
    rho,mu = xray_sld_from_atoms(atoms,2.52,energy=xray_energy(Cu.K_alpha))
    assert abs(rho-20.17)<0.1
    

if __name__ == "__main__": test()
