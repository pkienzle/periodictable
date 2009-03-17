# This program is public domain
"""
X-ray scatting information.

    xray.sftable (eV,?,?)
        Three column table of energy vs. scattering factors f1,f2.
    xray.scattering_factors(wavelength)
        Returns f1,f2, the X-ray scattering factors for the given wavelengths
        interploted from sftable.
    xray.sld(wavelength=A or energy=keV)
        Returns scattering length density and absorption for the
        given wavelengths or energies.

    K_alpha, K_beta1 (Angstrom)

        X-ray emission lines for various elements:
           Ag, Pd, Rh, Mo, Zn, Cu, Ni, Co, Fe, Mn, Cr and Ti
        K_alpha is the average of K_alpha1 and K_alpha2 lines.

For private tables, use init(table) and spectral_lines(table) to set
the data.

K-alpha and K-beta1 lines::

    D.C. Creagh and J.H. Hubbell (1999), "4.2.4 X-ray absorption
    (or attenuation) coefficients" in A.J.C. Wilson, E.Prince Eds.,
    International Tables for Crystallography Vol C p220-236.
    Kluwer Academic Publishers, Dordrecht, The Netherlands.

X-ray scattering factors::

    Low-Energy X-ray Interaction Coefficients: Photoabsorption,
    Scattering, and Reflection E = 30-30,000 eV, Z = 1-92

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

from .core import elements, Element, Isotope
#from . import mass,density
from .constants import (avogadro_number, plancks_constant, speed_of_light,
                        electron_radius)


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

    See help(periodictable.xsf) for details.
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
            periodictable.xsf.electron_radius
            periodictable.constants.avogadro_number

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

    Element.K_alpha_units = "angstrom"
    Element.K_beta1_units = "angstrom"

def init_spectral_lines(table):
    """
    Sets the K_alpha and K_beta1 wavelengths for select elements
    """
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
    for el in table:
        el.xray = Xray(el)

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

def sld_table(wavelength, table=elements):
    ## Cu K-alpha Zeff and dZ table
    #for el in table:
    #    f1,f2 = el.xsf(table.Cu.K_alpha)
    #    if f1 is not None:
    #        print el.Z,el.symbol,"%.1f"%f1,"%.4f"%(f1-el.Z)

    # NBCU spreadsheet format
    print "X-ray scattering length density and absorption for",wavelength,"A"
    print "%3s %6s %6s"%('El','rho','mu')
    for el in table:
        rho,mu = el.xray.sld(table.Cu.K_alpha)
        if rho is not None:
            print "%3s %6.2f %6.2f"%(el.symbol,rho,mu)

def emission_table(table=elements):
    """
    Print a table of emission lines.
    """
    print "%3s %7s %7s"%('El','Kalpha','Kbeta1')
    for el in table:
        if hasattr(el,'K_alpha'):
            print "%3s %7.4f %7.4f"%(el.symbol,el.K_alpha,el.K_beta1)
