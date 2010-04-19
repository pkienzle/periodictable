# This program is public domain

"""
Neutron scattering factors for the elements and isotopes, from the
ILL Neutron Data Booklet[1].

This module adds the neutron property to the periodic table.  For
details of neutron scattering factor values, see `periodictable.nsf.Neutron`.
The property is set to None if there is no neutron scattering information
for the element.  Individual isotopes may have their own scattering
information

Example
=======

The following plots the symbol versus scattering length density for
all elements::

    import periodictable,pylab
    SLDs = [(el.number,el.neutron.sld()[0],el.symbol)
            for el in periodictable.elements
            if el.neutron.has_sld()]
    for Z,sld,sym in SLDs:
        if sld is not None: pylab.text(Z,sld,sym)
    pylab.axis([0,100,-4,10])

Similarly, you can print a table of scattering length densities for isotopes
of a particular element::

    for iso in periodictable.Ni:
        if iso.neutron.has_sld():
            print iso,iso.neutron.sld()[0]

Details
=======

There are a number of functions available in periodictable.nsf

    neutron_sld(molecule)
        computes sld for a molecule (available as periodictable.neutron_sld()).
    neutron_sld_from_atoms({isotope:quantity}, density)
        computes sld from a set of isotopes, a density and a wavelength
    energy_dependent_table()
        lists energy dependent isotopes
    sld_table()
        lists all elements in natural abundance
    absorption_comparison_table()
        compares el.neutron.b_c_i and el.neutron.absorption
    coherent_comparison_table()
        compares el.neutron.b_c and el.neutron.coherent

For private tables use init(table) to set the data.

The neutron scattering information table is reproduced from the Atomic
Institute for Austrian Universities (2007 version)::

    http://www.ati.ac.at/~neutropt/scattering/table.html

The above site has references to the published values for every entry in
the table.  We have included these in the documentation directory
associated with the periodictable package.

[1] H. Rauch and W. Waschkowski (2003). Neutron Scattering Lengths
in ILL Neutron Data Booklet (second edition), A.-J. Dianoux, G. Lander, Eds.
Old City Publishing, Philidelphia, PA. pp 1.1-1 to 1.1-17.

[2] H. Rauch, W. Waschkowski (2000). Landolt-Boernstein, New Series I/16A.
H.Schopper Ed. Chap. 6. Springer: Berlin.

[3] L. Koester, H. Rauch, E. Seymann. Atomic Data Nuclear
Data Tables 49 (1991) 65
"""
from numpy import sqrt, pi
from .core import Element, Isotope, default_table
from .constants import avogadro_number

__all__ = ['init', 'Neutron',
           'absorption_comparison_table', 'coherent_comparison_table',
           'energy_dependent_table',
           'neutron_sld', 'neutron_sld_from_atoms', 'sld_plot', 'sld_table'
           ]

class Neutron(object):
    """
    Neutron scattering factors.

    Neutron scattering factors are attached to each element in the periodic
    table for which values are available.  If no information is available,
    then the neutron field of the element will be None.  Even when neutron
    information is available, it may not be complete, so individual fields
    may be None.

    The following fields are used

    b_c (fm)
        Bound coherent scattering length.
    b_c_i (fm)
        Imaginary part of bound coherent scattering length.  This is
        related to absorption cross section by 2 pi/k where k = 2 pi / lambda
        with a factor of 100 for converting between barns and fm.  b_c_i
        is not available for all isotopes for which absorption cross sections
        have been measured.
    bp,bm (fm)
        Spin-dependent scattering for I+1/2 and I-1/2 (not always available).
        Incoherent scattering arises from the spin-dependent scattering b+
        and b-.  The Neutron Data Booklet[1] gives formulas for calculating
        coherent and incoherent scattering from b+ and b- alone.
    bp_i,bm_i (fm)
        Imaginary portion.  See the Neutron Data Booklet[1] for details.
    is_energy_dependent (boolean)
        Do not use this data if scattering is energy dependent.
    coherent (barns)
        Coherent scattering cross section.  In theory coherent scattering
        is related to bound coherent scattering by :math:`4 \pi b_c^2 / 100`.
        In practice, these values are different, with the following examples
        the worst

        +--------+--------+-------+--------+-------+-------+-------+---------+
        +========+========+=======+========+=======+=======+=======+=========+
        | Sc  3% | Ti  4% | V 34% | Mn  1% | Cd 4% | Te 4% | Xe 9% | Sm 100% |
        +--------+--------+-------+--------+-------+-------+-------+---------+
        | Eu 46% | Gd 61% | Tb 1% | Ho 11% | W  4% | Au 7% | Hg 2% |         |
        +--------+--------+-------+--------+-------+-------+-------+---------+
    incoherent (barns)
        Incoherent scattering cross section.
    total (barns)
        Total scattering cross section.  This is just coherent+incoherent.
    absorption (barns)
        Absorption cross section at 1.798 angstroms.  Scale to your energy
        by dividing by 1.798 and multiplying by your wavelength.

    For elements, the scattering cross-sections are based on the natural
    abundance of the individual isotopes.  Individual isotopes may have
    additional information as follows

    abundance (%)
        Abundance used in elemental measurements.
    nuclear_spin (string)
        Spin on the nucleus: '0', '1/2', '3/2', etc.

    Each field above has a corresponding \*_units field with the string
    value for the units.

    For scattering calculations, the scattering length density is the
    value of interest.  This is computed from the number_density of the
    individual elements, as derived from the element density and atomic
    mass.

    Note: 1 barn = 100 fm^2

    [1] H. Rauch and W. Waschkowski (2003). Neutron Scattering Lengths
    in ILL Neutron Data Booklet (second edition), A.-J. Dianoux, G. Lander, Eds.
    Old City Publishing, Philidelphia, PA. pp 1.1-1 to 1.1-17.
    """
    b_c = None
    b_c_i = None
    b_c_units = "fm"
    bp = None
    bp_i = None
    bp_units = "fm"
    bm = None
    bm_i = None
    bm_units = "fm"
    coherent = None
    coherent_units = "barns"
    incoherent = None
    incoherent_units = "barns"
    total = None
    total_units = "barns"
    absorption = None
    absorption_units = "barns"
    abundance = 0.
    abundance_units = "%"
    is_energy_dependent = False
    def __init__(self):
        self._number_density = None
    def __str__(self):
        return "b_c=%.3g coh=%.3g inc=%.3g abs=%.3g"\
            %(self.b_c,self.coherent,self.incoherent,self.absorption)
    def has_sld(self):
        """Return True if sld is defined for this element/isotope"""
        return None not in [self.b_c, self._number_density]
    def sld(self,wavelength=1.798):
        """
        Returns scattering length density (real, imaginary, incoherent)
        for the element at natural density.

        The coherent scattering returned is N * b_c*10, where N is
        the number density computed from the bulk density of the element
        and b_c is the bound coherent scattering length for the isotope.
        For most elements, the coherent scattering cross section is
        independent of energy.  Those elements and isotopes for which this 
        does not hold have neutron.is_energy_dependent set to True.

        The incoherent scattering is returned as N * incoherent.

        The absorption cross sections are tabulated at wavelength 1.798 A.
        In the thermal neutron energy range the absorption cross section
        is assumed to scale linearly with wavelength.  See [1] for
        details.

        The absorption cross section is related to the imaginary portion
        of the bound coherent scattering cross section by the formula[2]
        absorption = 4*pi/k b'', where k is 2*pi/lambda, and there is an 
        additional factor of 100 converting from barns to fm^2.
        The value we return for the SLD is thus::

            Im(b_c) = 0.01 * absorption*N/(2 * 1.798) * wavelength

        Note: There is a factor of 10 unaccounted for, but required in order
        to match the b_c_i values given in the underlying tables.  Run::

            periodictable.nsf.absorption_comparison_table()

        to show how well b_c_i corresponds to absorption.

        [1] Lynn, JE and Seeger, PA (1990). Resonance effects in neutron
        scattering lengths of rare-earth nuclides. Atomic Data and
        Nuclear Data Tables 44, 191-207

        [2] Sears, VF (1999) 4.4.4 Scattering lengths for neutrons.
        In Wilson & Prince eds. Intl. Tables for Crystallography C
        Kluwer Academic Publishers. pp 448-449
        """

        # Compute number and absorption density assuming isotope has
        # same structure as the bulk element
        if not self.has_sld(): return None,None,None
        N = self._number_density*1e-24
        bp = self.b_c*10*N
        bpp = self.absorption/(2*1.798)*N*0.01
        binc = self.incoherent*N*0.01
        return bp,bpp,binc

def init(table, reload=False):
    """
    Load the Rauch table from the neutron data book.
    """
    if 'neutron' in table.properties and not reload: return
    table.properties.append('neutron')
    assert ('density' in table.properties and
        'mass' in table.properties), \
        "Neutron table requires mass and density properties"

    # Defaults for missing neutron information
    missing = Neutron()
    Isotope.neutron = missing
    Element.neutron = missing

    for line in nsftable.split('\n'):
        columns = line.split(',')

        nsf = Neutron()
        p = columns[1]
        spin = columns[2]
        nsf.b_c,nsf.bp, nsf.bm = [fix_number(a) for a in columns[3:6]]
        nsf.is_energy_dependent = (columns[6] == 'E')
        nsf.coherent,nsf.incoherent,nsf.total,nsf.absorption \
            = [fix_number(a) for a in columns[7:]]

        parts = columns[0].split('-')
        Z = int(parts[0])
        symbol = parts[1]
        isotope_number = int(parts[2]) if len(parts)==3 else 0

        # Fetch element from the table and check that the symbol matches
        element = table[Z]
        assert element.symbol == symbol, \
            "Symbol %s does not match %s"%(symbol,element.symbol)

        # Plug the default number density for the element into the nsf so
        # it can calculate sld.
        nsf._number_density = element.number_density


        # For new elements, clear out 'neutron' attribute for isotopes
        # This protects against isotope using the element data when
        # they don't have any specific neutron data.
        #if isotope_number == 0 or not hasattr(element,'neutron'):
        #    for iso in element: iso.neutron = None

        if isotope_number == 0:
            # Bulk values using laboratory abundances of isotopes
            element.neutron = nsf
        else:
            # Values for the individual isotope
            isotope = element.add_isotope(isotope_number)
            isotope.neutron = nsf
            isotope.nuclear_spin = spin
            # p column contains either abundance(uncertainty) or "half-life Y"
            isotope.neutron.abundance = fix_number(p) if ' ' not in p else 0

            # If the element is not yet initialized, copy info into the atom.
            # This serves to set the element info for elements with only
            # one isotope.
            if element.neutron is missing:
                element.neutron = nsf

    for line in nsftableI.split('\n'):
        columns = line.split(',')

        # Fetch the nsf record
        parts = columns[0].split('-')
        Z = int(parts[0])
        symbol = parts[1]
        isotope_number = int(parts[2]) if len(parts)==3 else 0
        element = table[Z]
        if isotope_number == 0:
            nsf = element.neutron
        else:
            nsf = element[isotope_number].neutron

        # Read imaginary values
        nsf.b_c_i,nsf.bp_i,nsf.bm_i = [fix_number(a) for a in columns[1:]]

def fix_number(str):
    """
    Convert strings of the form e.g., 35.24(2)* into numbers without
    uncertainty.  Also accepts a limited range, e.g., <1e-6, which is
    converted as 1e-6.  Missing values are set to 0.
    """
    if str == '': return None
    idx = str.find('(')
    if idx >= 0: str = str[0:idx]
    if str[0] == '<': str = str[1:]
    return float(str)

def sld_table(wavelength=1, table=None, isotopes=True):
    """
    Scattering length density table for wavelength 4.75 A
    """
    table = default_table(table)
    # Table for comparison with scattering length density calculators
    # b_c for Sc, Te, Xe, Sm, Eu, Gd, W, Au, Hg are different from Neutron News
    # The Rauch data have cited references to back up the numbers
    # (see doc directory), though it is not clear what criteria are
    # used to select amongst the available measurements.
    print " Neutron scattering length density table"
    print "%-7s %7s %7s %7s %7s %7s"%('atom','mass','density',
                                         'sld','imag','incoh')
    for el in table:
        if el.neutron.has_sld():
            coh,jcoh,inc = el.neutron.sld(wavelength)
            print "%-7s %7.3f %7.3f %7.3f %7.3f %7.3f %s"\
                %(el,el.mass,el.density,coh,jcoh,inc,
                  '*' if el.neutron.is_energy_dependent else '')
            if isotopes:
                isos = [iso for iso in el if iso.neutron != None and iso.neutron.has_sld()]
            else:
                isos = []
            for iso in isos:
                    coh,jcoh,inc = iso.neutron.sld(wavelength)
                    print "%-7s %7.3f %7.3f %7.3f %7.3f %7.3f %s"\
                        %(iso,iso.mass,iso.density,coh,jcoh,inc,
                          '*' if iso.neutron.is_energy_dependent else '')
    print "* Energy dependent cross sections"

def energy_dependent_table(table=None):
    """
    Print a table of energy dependent isotopes.
    """
    table = default_table(table)
    # List of energy dependent elements and isotopes
    print "Elements and isotopes with energy dependent absorption:"
    for el in table:
        if not hasattr(el,'neutron'): continue
        dep = []
        if el.neutron.is_energy_dependent:
            dep += [str(el)]
        dep += [str(el)+'-'+str(iso.isotope)
                for iso in el
                if iso.neutron != None and iso.neutron.is_energy_dependent]
        if len(dep) > 0: print "   "," ".join(dep)


# Note: docs and function prototype are reproduced in __init__
def neutron_sld(input, density=None, wavelength=1):
    """
    Compute neutron scattering length densities for molecules.

    Returns the scattering length density (real, imaginary incoherent).

    Raises AssertionError if density is missing.
    """
    import formulas
    mol = formulas.Formula(input)
    if density is None: density = mol.density # defaults to molecule density
    return neutron_sld_from_atoms(mol.atoms,
                                  density=density,
                                  wavelength=wavelength)

def neutron_sld_from_atoms(atoms, density=None, wavelength=1):
    """
    The underlying scattering length density calculator.  This
    works with a dictionary of atoms and quanties directly, such
    as returned by molecule.atoms.

    Returns the scattering length density (real, imaginary, incoherent).

    Raises AssertionError if density is missing.
    """
    assert density is not None, "neutron_sld needs density"

    b_c = 0
    absorption = 0
    incoherent = 0
    mass = 0
    is_energy_dependent = False
    for element,quantity in atoms.iteritems():
        #print element,quantity,element.neutron.b_c,element.neutron.absorption,element.neutron.incoherent
        mass += element.mass*quantity
        b_c += element.neutron.b_c*quantity
        absorption += element.neutron.absorption*quantity*wavelength/1.798
        incoherent += element.neutron.incoherent*quantity
        is_energy_dependent |= element.neutron.is_energy_dependent

    if mass == 0:  # for empty formula
        return 0,0,0
    else:
        cell_volume = mass/(density*avogadro_number*1e-24) # (10^8 A/cm)^3
        bp = 10*b_c/cell_volume
        bpp = 0.01*absorption/(2*wavelength)/cell_volume  # b'' = sigma_a/(2*lambda)
        binc = 0.01*incoherent/cell_volume
        return bp,bpp,binc

# We are including the complete original table here in case somebody in
# future wants to extract uncertainties or other information.
#
# Z-Symbol-A
#   This is the atomic number, the symbol and the isotope.
#   If Z-Symbol only, the line represents an element with scattering determined
#   by the natural abundance of the isotopes in laboratory samples.  If there
#   is only one isotope, then there is no corresponding element definition.
# concentration/half-life
#   This is the natural abundance of the isotope expressed as a percentage, or
#   it is the half-life in years (number Y) or seconds (number S).
# spin I
#   For isotopes, the nuclear spin.
# b_c, bp, bm
#   Bound coherent scattering length in fm
#   b+/b- if present are spin dependent scattering for I+1/2 and I-1/2
#   respectively
# c
#   'E' if there is a strong energy dependency.
#   '+/-' if separate b+/b- values are available [doesn't seem true -PAK]
# coherent, incoherent, total
#   The coherent and incoherent scattering cross-sections in barns.
# absorption
#   The thermal absorption cross section in barns at 1.798 Angstroms/25.30 meV.
#
# Numbers in parenthesis represents uncertainty.
# Numbers followed by '*' are estimated.
# Numbers may be given as limit, e.g., <1.0e-6
#
# Formatting corrections by Paul Kienzle

nsftable="""\
0-n-1,618 S,1/2,-37.0(6),0,-37.0(6),,43.01(2),,43.01(2),0
1-H,,,-3.7409(11),,,,1.7568(10),80.26(6),82.02(6),0.3326(7)
1-H-1,99.985,1/2,-3.7423(12),10.817(5),-47.420(14),+/-,1.7583(10),80.27(6),82.03(6),0.3326(7)
1-H-2,0.0149,1,6.674(6),9.53(3),0.975(60),,5.592(7),2.05(3),7.64(3),0.000519(7)
1-H-3,12.26 Y,1/2,4.792(27),4.18(15),6.56(37),,2.89(3),0.14(4),3.03(5),<6.0E-6
2-He,,,3.26(3),,,,1.34(2),0,1.34(2),0.00747(1)
2-He-3,0.013,1/2,5.74(7),4.7(5),8.8(1.4),E,4.42(10),1.6(4),6.0(4),5333.0(7.0)
2-He-4,99.987,0,3.26(3),,,,1.34(2),0,1.34(2),0
3-Li,,,-1.90(3),,,,0.454(10),0.92(3),1.37(3),70.5(3)
3-Li-6,7.5,1,2.0(1),0.67(14),4.67(17),+/-,0.51(5),0.46(5),0.97(7),940.0(4.0)
3-Li-7,92.5,3/2,-2.22(2),-4.15(6),1.00(8),+/-,0.619(11),0.78(3),1.40(3),0.0454(3)
4-Be-9,100,3/2,7.79(1),,,,7.63(2),0.0018(9),7.63(2),0.0076(8)
5-B,,,5.30(4),,,,3.54(5),1.70(12),5.24(11),767.0(8.0)
5-B-10,19.4,3,-0.2(4),-4.2(4),5.2(4),,0.144(6),3.0(4),3.1(4),3835.0(9.0)
5-B-11,80.2,3/2,6.65(4),5.6(3),8.3(3),,5.56(7),0.21(7),5.77(10),0.0055(33)
6-C,,,6.6484(13),,,,5.551(2),0.001(4),5.551(3),0.00350(7)
6-C-12,98.89,0,6.6535(14),,,,5.559(3),0,5.559(3),0.00353(7)
6-C-13,1.11,1/2,6.19(9),5.6(5),6.2(5),+/-,4.81(14),0.034(11),4.84(14),0.00137(4)
7-N,,,9.36(2),,,,11.01(5),0.50(12),11.51(11),1.90(3)
7-N-14,99.635,1,9.37(2),10.7(2),6.2(3),,11.03(5),0.50(12),11.53(11),1.91(3)
7-N-15,0.365,1/2,6.44(3),6.77(10),6.21(10),,5.21(5),0.00005(10),5.21(5),0.000024(8)
8-O,,,5.805(4),,,,4.232(6),0.000(8),4.232(6),0.00019(2)
8-O-16,99.75,0,5.805(5),,,,4.232(6),0,4.232(6),0.00010(2)
8-O-17,0.039,5/2,5.6(5),5.52(20),5.17(20),,4.20(22),0.004(3),4.20(22),0.236(10)
8-O-18,0.208,0,5.84(7),,,,4.29(10),0,4.29(10),0.00016(1)
9-F-19,100,1/2,5.654(12),5.632(10),5.767(10),+/-,4.017(14),0.0008(2),4.018(14),0.0096(5)
10-Ne,,,4.566(6),,,,2.620(7),0.008(9),2.628(6),0.039(4)
10-Ne-20,90.5,0,4.631(6),,,,2.695(7),0,2.695(7),0.036(4)
10-Ne-21,0.27,3/2,6.66(19),,,,5.6(3),0.05(2),5.7(3),0.67(11)
10-Ne-22,9.2,0,3.87(1),,,,1.88(1),0,1.88(1),0.046(6)
11-Na-23,100,3/2,3.63(2),6.42(4),-1.00(6),+/-,1.66(2),1.62(3),3.28(4),0.530(5)
12-Mg,,,5.375(4),,,,3.631(5),0.08(6),3.71(4),0.063(3)
12-Mg-24,78.99,0,5.49(18),,,,4.03(4),0,4.03(4),0.050(5)
12-Mg-25,10,5/2,3.62(14),4.73(30),1.76(20),+/-,1.65(13),0.28(4),1.93(14),0.19(3)
12-Mg-26,11,0,4.89(15),,,,3.00(18),0,3.00(18),0.0382(8)
13-Al-27,100,5/2,3.449(5),3.67(2),3.15(2),,1.495(4),0.0082(6),1.503(4),0.231(3)
14-Si,,,4.15071(22),,,,2.1633(10),0.004(8),2.167(8),0.171(3)
14-Si-28,92.2,0,4.106(6),,,,2.120(6),0,2.120(6),0.177(3)
14-Si-29,4.7,1/2,4.7(1),4.50(15),4.7(4),+/-,2.78(12),0.001(2),2.78(12),0.101(14)
14-Si-30,3.1,0,4.58(8),,,,2.64(9),0,2.64(9),0.107(2)
15-P-31,100,1/2,5.13(1),,,+/-,3.307(13),0.005(10),3.312(16),0.172(6)
16-S,,,2.847(1),,,,1.0186(7),0.007(5),1.026(5),0.53(1)
16-S-32,95,0,2.804(2),,,,0.9880(14),0,0.9880(14),0.54(4)
16-S-33,0.74,3/2,4.74(19),,,+/-,2.8(2),0.3(6),3.1(6),0.54(4)
16-S-34,4.2,0,3.48(3),,,,1.52(3),0,1.52(3),0.227(5)
16-S-36,0.02,0,3.0(1.0)*,,,,1.1(8),0,1.1(8),0.15(3)
17-Cl,,,9.5792(8),,,,11.528(2),5.3(5),16.8(5),33.5(3)
17-Cl-35,75.77,3/2,11.70(9),16.3(2),4.0(3),+/-,17.06(6),4.7(6),21.8(6),44.1(4)
17-Cl-37,24.23,3/2,3.08(6),3.10(7),3.05(7),+/-,1.19(5),0.001(3),1.19(5),0.433(6)
18-Ar,,,1.909(6),,,,0.458(3),0.225(5),0.683(4),0.675(9)
18-Ar-36,0.34,0,24.9(7),,,,77.9(4),0,77.9(4),5.2(5)
18-Ar-38,0.07,0,3.5(3.5),,,,1.5(3.1),0,1.5(3.1),0.8(5)
18-Ar-40,99.59,0,1.7,,,,0.421(3),0,0.421(3),0.660(9)
19-K,,,3.67(2),,,,1.69(2),0.27(11),1.96(11),2.1(1)
19-K-39,93.3,3/2,3.79(2),5.15,1.51,+/-,1.76(2),0.25(11),2.01(11),2.1(1)
19-K-40,0.012,4,3.1(1.0)*,,,,1.1(6),0.5(5)*,1.6(9),35.0(8.0)
19-K-41,6.7,3/2,2.69(8),,,,0.91(5),0.3(6),1.2(6),1.46(3)
20-Ca,,,4.70(2),,,,2.78(2),0.05(3),2.83(2),0.43(2)
20-Ca-40,96.94,0,4.78(5),,,,2.90(2),0,2.90(2),0.41(2)
20-Ca-42,0.64,0,3.36(10),,,,1.42(8),0,1.42(8),0.68(7)
20-Ca-43,0.13,7/2,-1.56(9),,,,0.31(4),0.5(5),0.8(5),6.2(6)
20-Ca-44,2.13,0,1.42(6),,,,0.25(2),0,0.25(2),0.88(5)
20-Ca-46,0.003,0,3.55(21),,,,1.6(2),0,1.6(2),0.74(7)
20-Ca-48,0.18,0,0.39(9),,,,0.019(9),0,0.019(9),1.09(14)
21-Sc-45,100,7/2,12.1(1),6.91(22),18.99(28),+/-,19.0(3),4.5(3),23.5(6),27.5(2)
22-Ti,,,-3.370(13),,,,1.485(2),2.87(3),4.35(3),6.09(13)
22-Ti-46,8,0,4.72(5),,,,3.05(7),0,3.05(7),0.59(18)
22-Ti-47,7.5,5/2,3.53(7),0.46(23),7.64(13),,1.66(11),1.5(2),3.2(2),1.7(2)
22-Ti-48,73.7,0,-5.86(2),,,,4.65(3),0,4.65(3),7.84(25)
22-Ti-49,5.5,7/2,0.98(5),2.6(3),-1.2(4),,0.14(1),3.3(3),3.4(3),2.2(3)
22-Ti-50,5.3,0,5.88(10),,,,4.80(12),0,4.80(12),0.179(3)
23-V,,,-0.443(14),,,,0.01838(12),5.08(6),5.10(6),5.08(4)
23-V-50,0.25,6,7.6(6)*,,,,7.3(1.1),0.5(5)*,7.8(1.0),60.0(40.0)
23-V-51,99.75,7/2,-0.402(2),4.93(25),-7.58(28),+/-,0.0203(2),5.07(6),5.09(6),4.9(1)
24-Cr,,,3.635(7),,,,1.660(6),1.83(2),3.49(2),3.05(6)
24-Cr-50,4.35,0,-4.50(5),,,,2.54(6),0,2.54(6),15.8(2)
24-Cr-52,83.8,0,4.914(15),,,,3.042(12),0,3.042(12),0.76(6)
24-Cr-53,9.59,3/2,-4.20(3),1.16(10),-13.0(2),,2.22(3),5.93(17),8.15(17),18.1(1.5)
24-Cr-54,2.36,0,4.55(10),,,,2.60(11),0,2.60(11),0.36(4)
25-Mn-55,100,5/2,-3.750(18),-4.93(46),-1.46(33),,1.75(2),0.40(2),2.15(3),13.3(2)
26-Fe,,,9.45(2),,,,11.22(5),0.40(11),11.62(10),2.56(3)
26-Fe-54,5.8,0,4.2(1),,,,2.2(1),0,2.2(1),2.25(18)
26-Fe-56,91.7,0,10.1(2),,,,12.42(7),0,12.42(7),2.59(14)
26-Fe-57,2.19,1/2,2.3(1),,,,0.66(6),0.3(3)*,1.0(3),2.48(30)
26-Fe-58,0.28,0,15(7),,,,28.0(26.0),0,28.0(26.0),1.28(5)
27-Co-59,100,7/2,2.49(2),-9.21(10),3.58(10),+/-,0.779(13),4.8(3),5.6(3),37.18(6)
28-Ni,,,10.3(1),,,,13.3(3),5.2(4),18.5(3),4.49(16)
28-Ni-58,67.88,0,14.4(1),,,,26.1(4),0,26.1(4),4.6(3)
28-Ni-60,26.23,0,2.8(1),,,,0.99(7),0,0.99(7),2.9(2)
28-Ni-61,1.19,3/2,7.60(6),,,,7.26(11),1.9(3),9.2(3),2.5(8)
28-Ni-62,3.66,0,-8.7(2),,,,9.5(4),0,9.5(4),14.5(3)
28-Ni-64,1.08,0,-0.37(7),,,,0.017(7),0,0.017(7),1.52(3)
29-Cu,,,7.718(4),,,,7.485(8),0.55(3),8.03(3),3.78(2)
29-Cu-63,69.1,3/2,6.477(13),,,+/-,5.2(2),0.006(1),5.2(2),4.50(2)
29-Cu-65,30.9,3/2,10.204(20),,,+/-,14.1(5),0.40(4),14.5(5),2.17(3)
30-Zn,,,5.680(5),,,,4.054(7),0.077(7),4.131(10),1.11(2)
30-Zn-64,48.9,0,5.23(4),,,,3.42(5),0,3.42(5),0.93(9)
30-Zn-66,27.8,0,5.98(5),,,,4.48(8),0,4.48(8),0.62(6)
30-Zn-67,4.1,5/2,7.58(8),5.8(5),10.1(7),+/-,7.18(15),0.28(3),7.46(15),6.8(8)
30-Zn-68,18.6,0,6.04(3),,,,4.57(5),0,4.57(5),1.1(1)
30-Zn-70,0.62,0,6.9(1.0)*,,,,4.5(1.5),0,4.5(1.5),0.092(5)
31-Ga,,,7.288(2),,,,6.675(4),0.16(3),6.83(3),2.75(3)
31-Ga-69,60,3/2,8.043(16),6.3(2),10.5(4),+/-,7.80(4),0.091(11),7.89(4),2.18(5)
31-Ga-71,40,3/2,6.170(11),5.5(6),7.8(1),+/-,5.15(5),0.084(8),5.23(5),3.61(10)
32-Ge,,,8.185(20),,,,8.42(4),0.18(7),8.60(6),2.20(4)
32-Ge-70,20.7,0,10.0(1),,,,12.6(3),0,12.6(3),3.0(2)
32-Ge-72,27.5,0,8.51(10),,,,9.1(2),0,9.1(2),0.8(2)
32-Ge-73,7.7,9/2,5.02(4),8.1(4),1.2(4),,3.17(5),1.5(3),4.7(3),15.1(4)
32-Ge-74,36.4,0,7.58(10),,,,7.2(2),0,7.2(2),0.4(2)
32-Ge-76,7.7,0,8.2(1.5),,,,8.0(3.0),0,8.0(3.0),0.16(2)
33-As-75,100,3/2,6.58(1),6.04(5),7.47(8),+/-,5.44(2),0.060(10),5.50(2),4.5(1)
34-Se,,,7.970(9),,,,7.98(2),0.32(6),8.30(6),11.7(2)
34-Se-74,0.9,0,0.8(3.0),,,,0.1(6),0,0.1(6),51.8(1.2)
34-Se-76,9,0,12.2(1),,,,18.7(3),0,18.7(3),85.0(7.0)
34-Se-77,7.5,0,8.25(8),,,,8.6(2),0.05(25),8.65(16),42.0(4.0)
34-Se-78,23.5,0,8.24(9),,,,8.5(2),0,8.5(2),0.43(2)
34-Se-80,50,0,7.48(3),,,,7.03(6),0,7.03(6),0.61(5)
34-Se-82,8.84,0,6.34(8),,,,5.05(13),0,5.05(13),0.044(3)
35-Br,,,6.79(2),,,,5.80(3),0.10(9),5.90(9),6.9(2)
35-Br-79,50.49,3/2,6.79(7),,,+/-,5.81(2),0.15(6),5.96(13),11.0(7)
35-Br-81,49.31,3/2,6.78(7),,,+/-,5.79(12),0.05(2),5.84(12),2.7(2)
36-Kr,,,7.81(2),,,,7.67(4),0.01(14),7.68(13),25.0(1.0)
36-Kr-78,0.35,0,,,,,,0,,6.4(9)
36-Kr-80,2.5,0,,,,,,0,,11.8(5)
36-Kr-82,11.6,0,,,,,,0,,29.0(20.0)
36-Kr-83,11.5,9/2,,,,,,,,185.0(30.0)
36-Kr-84,57,0,,,,,,0,6.6,0.113(15)
36-Kr-86,17.3,0,8.07(26),,,,8.2(4),0,8.2(4),0.003(2)
37-Rb,,,7.08(2),,,,6.32(4),0.5(4),6.8(4),0.38(1)
37-Rb-85,72.17,5/2,7.07(10),,,,6.2(2),0.5(5)*,6.7(5),0.48(1)
37-Rb-87,27.83,3/2,7.27(12),,,,6.6(2),0.5(5)*,7.1(5),0.12(3)
38-Sr,,,7.02(2),,,,6.19(4),0.06(11),6.25(10),1.28(6)
38-Sr-84,0.56,0,5.0(2.0),,,,6.0(2.0),0,6.0(2.0),0.87(7)
38-Sr-86,9.9,0,5.68(5),,,,4.04(7),0,4.04(7),1.04(7)
38-Sr-87,7,9/2,7.41(7),,,,6.88(13),0.5(5)*,7.4(5),16.0(3.0)
38-Sr-88,82.6,0,7.16(6),,,,6.42(11),0,6.42(11),0.058(4)
39-Y-89,100,1/2,7.75(2),8.4(2),5.8(5),+/-,7.55(4),0.15(8),7.70(9),1.28(2)
40-Zr,,,7.16(3),,,,6.44(5),0.02(15),6.46(14),0.185(3)
40-Zr-90,51.48,0,6.5(1),,,,5.1(2),0,5.1(2),0.011(59
40-Zr-91,11.23,5/2,8.8(1),7.9(2),10.1(2),+/-,9.5(2),0.15(4),9.7(2),1.17(10)
40-Zr-92,17.11,0,7.5(2),,,,6.9(4),0,6.9(4),0.22(6)
40-Zr-94,17.4,0,8.3(2),,,,8.4(4),0,8.4(4),0.0499(24)
40-Zr-96,2.8,0,5.5(1),,,,3.8(1),0,3.8(1),0.0229(10)
41-Nb-93,100,9/2,7.054(3),7.06(4),7.35(4),+/-,6.253(5),0.0024(3),6.255(5),1.15(6)
42-Mo,,,6.715(20),,,,5.67(3),0.04(5),5.71(4),2.48(4)
42-Mo-92,15.48,0,6.93(8),,,,6.00(14),0,6.00(14),0.019(2)
42-Mo-94,9.1,0,6.82(7),,,,5.81(12),0,5.81(12),0.015(2)
42-Mo-95,15.72,5/2,6.93(7),,,,6.00(10),0.5(5)*,6.5(5),13.1(3)
42-Mo-96,16.53,0,6.22(6),,,,4.83(9),0,4.83(9),0.5(2)
42-Mo-97,9.5,5/2,7.26(8),,,,6.59(15),0.5(5)*,7.1(5),2.5(2)
42-Mo-98,23.78,0,6.60(7),,,,5.44(12),0,5.44(12),0.127(6)
42-Mo-100,9.6,0,6.75(7),,,,5.69(12),0,5.69(12),0.4(2)
43-Tc-99,210000 Y,9/2,6.8(3),,,,5.8(5),0.5(5)*,6.3(7),20.0(1.0)
44-Ru,,,7.02(2),,,,6.21(5),0.4(1),6.6(1),2.56(13)
44-Ru-96,5.8,0,,,,,,0,,0.28(2)
44-Ru-98,1.9,0,,,,,,0,,<8.0
44-Ru-99,12.7,5/2,,,,,,,,6.9(1.0)
44-Ru-100,12.6,0,,,,,,0,,4.8(6)
44-Ru-101,17.07,5/2,,,,,,,,3.3(9)
44-Ru-102,31.61,0,,,,,,0,,1.17(7)
44-Ru-104,18.58,0,,,,,,0,,0.31(2)
45-Rh-103,100,1/2,5.90(4),8.15(6),6.74(6),,4.34(6),0.3(3)*,4.6(3),144.8(7)
46-Pd,,,5.91(6),,,,4.39(9),0.093(9),4.48(9),6.9(4)
46-Pd-102,1,0,7.7(7)*,,,,7.5(1.4),0,7.5(1.4),3.4(3)
46-Pd-104,11,0,7.7(7)*,,,,7.5(1.4),0,7.5(1.4),0.6(3)
46-Pd-105,22.33,5/2,5.5(3),,,+/-,3.8(4),0.8(1.0),4.6(1.1),20.0(3.0)
46-Pd-106,27.33,0,6.4(4),,,,5.1(6),0,5.1(6),0.304(29)
46-Pd-108,26.71,0,4.1(3),,,,2.1(3),0,2.1(3),8.5(5)
46-Pd-110,11.8,0,7.7(7)*,,,,7.5(1.4),0,7.5(1.4),0.226(31)
47-Ag,,,5.922(7),,,,4.407(10),0.58(3),4.99(3),63.3(4)
47-Ag-107,51.8,1/2,7.555(11),8.14(9),5.8(3),+/-,7.17(2),0.13(3),7.30(4),37.6(1.2)
47-Ag-109,48.2,1/2,4.165(11),3.24(8),6.9(2),+/-,2.18(1),0.32(5),2.50(5),91.0(1.0)
48-Cd,,,4.83(5),,,E,3.04(6),3.46(13),6.50(12),2520.0(50.0)
48-Cd-106,1.2,0,5.0(2.0)*,,,,3.1(2.5),0,3.1(2.5),1.0(2.0)
48-Cd-108,0.9,0,5.31(24),,,,3.7(1),0,3.7(1),1.1(3)
48-Cd-110,12.39,0,5.78(8),,,,4.4(1),0,4.4(1),11.0(1.0)
48-Cd-111,12.75,1/2,6.47(8),,,,5.3(2),0.3(3)*,5.6(4),24.0(5.0)
48-Cd-112,24.07,0,6.34(6),,,,5.1(2),0,5.1(2),2.2(5)
48-Cd-113,12.36,1/2,-8.0(1),,,E,12.1(4),0.3(3)*,12.4(5),20600.0(400.0)
48-Cd-114,28.86,0,7.48(5),,,,7.1(2),0,7.1(2),0.34(2)
48-Cd-116,7.58,0,6.26(9),,,,5.0(2),0,5.0(2),0.075(13)
49-In,,,4.065(20),,,,2.08(2),0.54(11),2.62(11),193.8(1.5)
49-In-113,4.28,9/2,5.39(6),,,,3.65(8),0.000037(5),3.65(8),12.0(1.1)
49-In-115,95.72,9/2,4.00(3),2.1(1),6.4(4),,2.02(2),0.55(11),2.57(11),202.0(2.0)
50-Sn,,,6.225(2),,,,4.871(3),0.022(5),4.892(6),0.626(9)
50-Sn-112,1,0,6.0(1.0)*,,,,4.5(1.5),0,4.5(1.5),1.00(11)
50-Sn-114,0.66,0,6.0(3),,,,4.8(5),0,4.8(5),0.114(30)
50-Sn-115,0.35,1/2,6.0(1.0)*,,,,4.5(1.5),0.3(3)*,4.8(1.5),30.0(7.0)
50-Sn-116,14.3,0,6.10(1),,,,4.42(7),0,4.42(7),0.14(3)
50-Sn-117,7.61,1/2,6.59(8),0.22(10),-0.23(10),,5.28(8),0.3(3)*,5.6(3),2.3(5)
50-Sn-118,24.03,0,6.23(4),,,,4.63(8),0,4.63(8),0.22(5)
50-Sn-119,8.58,1/2,6.28(3),0.14(10),0.0(1),,4.71(8),0.3(3)*,5.0(3),2.2(5)
50-Sn-120,32.86,0,6.67(4),,,,5.29(8),0,5.29(8),0.14(3)
50-Sn-122,4.72,0,5.93(3),,,,4.14(7),0,4.14(7),0.18(2)
50-Sn-124,5.94,0,6.15(3),,,,4.48(8),0,4.48(8),0.133(5)
51-Sb,,,5.57(3),,,,3.90(4),0.00(7),3.90(6),4.91(5)
51-Sb-121,57.25,5/2,5.71(6),5.7(2),5.8(2),,4.10(9),0.0003(19),4.10(19),5.75(12)
51-Sb-123,42.75,7/2,5.38(7),5.2(2),5.4(2),,3.64(9),0.001(4),3.64(9),3.8(2)
52-Te,,,5.68(2),,,,4.23(4),0.09(6),4.32(5),4.7(1)
52-Te-120,0.09,0,5.3(5),,,,3.5(7),0,3.5(7),2.3(3)
52-Te-122,2.4,0,3.8(2),,,,1.8(2),0,1.8(2),3.4(5)
52-Te-123,0.87,1/2,-0.05(25),-1.2(2),3.5(2),,0.002(3),0.52(5),0.52(5),418.0(30.0)
52-Te-124,4.61,0,7.95(10),,,,8.0(2),0,8.0(2,6.8(1.3)
52-Te-125,6.99,1/2,5.01(8),4.9(2),5.5(2),,3.17(10),0.008(8),3.18(10),1.55(16)
52-Te-126,18.71,0,5.55(7),,,,3.88(10),0,3.88(10),1.04(15)
52-Te-128,31.79,0,5.88(8),,,,4.36(10),0,4.36(10),0.215(8)
52-Te-130,34.48,0,6.01(7),,,,4.55(11),0,4.55(11),0.29(6)
53-I-127,100,5/2,5.28(2),6.6(2),3.4(2),,3.50(3),0.31(6),3.81(7),6.15(6)
54-Xe,,,4.69(4),,,,3.04(4),0,,23.9(1.2)
54-Xe-124,0.1,0,,,,,,0,,165.0(20.0)
54-Xe-126,0.09,0,,,,,,0,,3.5(8)
54-Xe-128,1.9,0,,,,,,0,,<8.0
54-Xe-129,26.14,1/2,,,,,,,,21.0(5.0)
54-Xe-130,3.3,0,,,,,,0,,<26.0
54-Xe-131,21.18,3/2,,,,,,,,85.0(10.0)
54-Xe-132,26.89,0,,,,,,0,,0.45(6)
54-Xe-134,10.4,0,,,,,,0,,0.265(20)
54-Xe-136,8.9,0,,,,,,0,,0.26(2)
55-Cs-133,100,7/2,5.42(2),,,+/-,3.69(15),0.21(5),3.90(6),29.0(1.5)
56-Ba,,,5.07(3),,,,3.23(4),0.15(11),3.38(10),1.1(1)
56-Ba-130,0.1,0,-3.6(6),,,,1.6(5),0,1.6(5),30.0(5.0)
56-Ba-132,0.09,0,7.8(3),,,,7.6(6),0,7.6(6),7.0(8)
56-Ba-134,2.4,0,5.7(1),,,,4.08(14),0,4.08(14),2.0(1.6)
56-Ba-135,6.59,3/2,4.66(10),,,,2.74(12),0.5(5)*,3.2(5),5.8(9)
56-Ba-136,7.81,0,4.90(8),,,,3.03(10),0,3.03(10),0.68(17)
56-Ba-137,11.32,3/2,6.82(10),,,,5.86(17),0.5(5)*,6.4(5),3.6(2)
56-Ba-138,71.66,0,4.83(8),,,,2.94(10),0,2.94(19),0.27(14)
57-La,,,8.24(4),,,,8.53(8),1.13(19),9.66(17),8.97(2)
57-La-138,0.09,5,8.0(2.0)*,,,,8.0(4.0),0.5(5)*,8.5(4.0),57.0(6.0)
57-La-139,99.91,7/2,8.24(4),11.4(3),4.5(4),+/-,8.53(8),1.13(15),9.66(17),8.93(4)
58-Ce,,,4.84(2),,,,2.94(2),0.00(10),2.94(10),0.63(4)
58-Ce-136,0.19,0,5.76(9),,,,4.23(13),0,4.23(13),7.3(1.5)
58-Ce-138,0.26,0,6.65(9),,,,5.64(15),0,5.64(15),1.1(3)
58-Ce-140,88.48,0,4.81(9),,,,2.94(11),0,2.94(11),0.57(4)
58-Ce-142,11.07,0,4.72(9),,,,2.84(11),0,2.84(11),0.95(5)
59-Pr-141,100,5/2,4.58(5),,,+/-,2.64(6),0.015(3),2.66(6),11.5(3)
60-Nd,,,7.69(5),,,,7.43(19),9.2(8),16.6(8),50.5(1.2)
60-Nd-142,27.11,0,7.7(3),,,,7.5(6),0,7.5(6),18.7(7)
60-Nd-143,12.17,7/2,14.0(2.0)*,,,,25.0(7.0),55.0(7.0),80.0(2.0),337.0(10.0)
60-Nd-144,23.85,0,2.8(3),,,,1.0(2),0,1.0(2),3.6(3)
60-Nd-145,8.5,7/2,14.0(2.0)*,,,,25.0(7.0),5.0(5.0)*,30.0(9.0),42.0(2.0)
60-Nd-146,17.22,0,8.7(2),,,,9.5(4),0,9.5(4),1.4(1)
60-Nd-148,5.7,0,5.7(3),,,,4.1(4),0,4.1(4),2.5(2)
60-Nd-150,5.6,0,5.28(20),,,,3.5(3),0,3.5(3),1.2(2)
61-Pm-147,2.62 Y,7/2,12.6(4),,,,20.0(1.3),1.3(2.0),21.3(1.5),168.4(3.5)
62-Sm,,,0.00(5),,,E,0.422(9),39.0(3.0),39.4(3.0),5922.0(56.0)
62-Sm-144,3.1,0,-3.0(4.0)*,,,,1.0(3.0),0,1.0(3.0),0.7(3)
62-Sm-147,15,7/2,14.0(3.0),,,,25.0(11.0),14.0(19.0.),39.0(16.0),57.0(3.0)
62-Sm-148,11.2,0,-3.0(4.0)*,,,,1.0(3.0),0,1.0(3.0),2.4(6)
62-Sm-149,13.8,7/2,18.7(28),,,E,63.5(6),137.0(5.0),200.0(5.0),42080.0(400.0)
62-Sm-150,7.4,0,14.0(3.0),,,,25.0(11.0),0,25.0(11.0),104.0(4.0)
62-Sm-152,26.7,0,-5.0(6),,,,3.1(8),0,3.1(8),206.0(6.0)
62-Sm-154,22.8,0,8.0(1.0),,,,11.0(2.0),0,11.0(2.0),8.4(5)
63-Eu,,,5.3(3),,,E,6.57(4),2.5(4),9.2(4),4530.0(40.0)
63-Eu-151,47.8,5/2,,,,E,5.5(2),3.1(4),8.6(4),9100.0(100.0)
63-Eu-153,52.8,5/2,8.22(12),,,,8.5(2),1.3(7),9.8(7),312.0(7.0)
64-Gd,,,9.5(2),,,E,29.3(8),151.0(2.0),180.0(2.0),49700.0(125.0)
64-Gd-152,0.2,0,10.0(3.0)*,,,,13.0(8.0),0,13.0(8.0),735.0(20.0)
64-Gd-154,2.2,0,10.0(3.0)*,,,,13.0(8.0),0,13.0(8.0),85.0(12.0)
64-Gd-155,14.9,3/2,13.8(3),,,E,40.8(4),25.0(6.0),66.0(6.0),61100.0(400.0)
64-Gd-156,20.6,0,6.3(4),,,,5.0(6),0,5.0(6),1.5(1.2)
64-Gd-157,15.7,3/2,4.0(2.0),,,E,650.0(4.0),394.0(7.0),1044.0(8.0),259000.0(700.0)
64-Gd-158,24.7,0,9.0(2.0),,,,10.0(5.0),0,10.0(5.0),2.2(2)
64-Gd-160,21.7,0,9.15(5),,,,10.52(11),0,10.52(11),0.77(2)
65-Tb-159,100,3/2,7.34(2),6.8(2),8.1(2),+/-,6.84(6),0.004(3),6.84(6),23.4(4)
66-Dy,,,16.9(3),,,,35.9(8),54.4(1.2),90.3(9),994.0(13.0)
66-Dy-156,0.06,0,6.1(5),,,,4.7(8),0,4.7(8),33.0(3.0)
66-Dy-158,0.1,0,6.0(4.0)*,,,,5.0(6.0),0,5.(6.),43.0(6.0)
66-Dy-160,2.3,0,6.7(4),,,,5.6(7),0,5.6(7),56.0(5.0)
66-Dy-161,18.9,5/2,10.3(4),,,,13.3(1.0),3.0(1.0),16.0(1.0),600.0(25.0)
66-Dy-162,25.5,0,-1.4(5),,,,0.25(18),0,0.25(18),194.0(10.0)
66-Dy-163,24.9,5/2,5.0(4),6.1(5),3.5(5),,3.1(5),0.21(19),3.3(5),124.0(7.0)
66-Dy-164,28.2,0,49.4(5),,,,307.0(3.0),0,307.0(3.0),2840.0(40.0)
67-Ho-165,100,7/2,8.44(3),6.9(2),10.3(2),+/-,8.06(8),0.36(3),8.42(16),64.7(1.2)
68-Er,,,7.79(2),,,,7.63(4),1.1(3),8.7(3),159.0(4.0)
68-Er-162,0.14,0,9.01(11),,,,9.7(4),0,9.7(4),19.0(2.0)
68-Er-164,1.6,0,7.95(14),,,,8.4(4),0,8.4(4),13.0(2.0)
68-Er-166,33.4,0,10.51(19),,,,14.1(5),0,14.1(5),19.6(1.5)
68-Er-167,22.9,7/2,3.06(5),5.3(3),0.0(3),,1.1(2),0.13(6),1.2(2),659.0(16.0)
68-Er-168,27,0,7.43(8),,,,6.9(7),0,6.9(7),2.74(8)
68-Er-170,15,0,9.61(6),,,,11.6(1.2),0,11.6(1.2),5.8(3)
69-Tm-169,100,1/2,7.07(3),,,+/-,6.28(5),0.10(7),6.38(9),100.0(2.0)
70-Yb,,,12.41(3),,,,19.42(9),4.0(2),23.4(2),34.8(8)
70-Yb-168,0.14,0,-4.07(2),,,E,2.13(2),0,2.13(2),2230.0(40.0)
70-Yb-170,3,0,6.8(1),,,,5.8(2),0,5.8(2),11.4(1.0)
70-Yb-171,14.3,1/2,9.7(1),6.5(2),19.4(4),,11.7(2),3.9(2),15.6(3),48.6(2.5)
70-Yb-172,21.9,0,9.5(1),,,,11.2(2),0,11.2(2),0.8(4)
70-Yb-173,16.3,5/2,9.56(10),2.5(2),13.3(3),,11.5(2),3.5,15,17.1(1.3)
70-Yb-174,31.8,0,19.2(1),,,,46.8(5),0,46.8(5),69.4(5.0)
70-Yb-176,12.7,0,8.7(1),,,,9.6(2),0,9.6(2),2.85(5)
71-Lu,,,7.21(3),,,,6.53(5),0.7(4),7.2(4),74.0(2.0)
71-Lu-175,97.4,7/2,7.28(9),,,,6.59(5),0.6(4),7.2(4),21.0(3.0)
71-Lu-176,2.6,7,6.1(2),,,,4.7(2),1.2(3),5.9,2065.(35.)
72-Hf,,,7.77(14),,,,7.6(3),2.6(5),10.2(4),104.1(5)
72-Hf-174,0.184,0,10.9(1.1),,,,15.0(3.0),0,15.0(3.0),561.0(35.0)
72-Hf-176,5.2,0,6.61(18),,,,5.5(3),0,5.5(3),23.5(3.1)
72-Hf-177,18.5,0,0.8(1.0)*,,,,0.1(2),0.1(3),0.2(2),373.0(10.0)
72-Hf-178,27.2,0,5.9(2),,,,4.4(3),0,4.4(3),84.0(4.0)
72-Hf-179,13.8,9/2,7.46(16),,,,7.0(3),0.14(2),7.1(3),41.0(3.0)
72-Hf-180,35.1,0,13.2(3),,,,21.9(1.0),0,21.9(1.0),13.04(7)
73-Ta,,,6.91(7),,,,6.00(12),0.01(17),6.01(12),20.6(5)
73-Ta-180,0.012,9,7.0(2.0)*,,,,6.2(3.5),0.5(5)*,7.0(4.0),563.0(60.0)
73-Ta-181,99.98,7/2,6.91(7),,,+/-,6.00(12),0.011(2),6.01(12),20.5(5)
74-W,,,4.755(18),,,,2.97(2),1.63(6),4.60(6),18.3(2)
74-W-180,0.13,0,5.0(3.0)*,,,,3.0(4.0),0,3.0(4.0),30.0(20.0)
74-W-182,26.3,1/2,7.04(4),,,,6.10(7),0,6.10(7),20.7(5)
74-W-183,14.3,1/2,6.59(4),6.3(4),7.0(4),,5.36(7),0.3(3)*,5.7(3),10.1(3)
74-W-184,30.7,0,7.55(6),,,,7.03(11),0,7.03(11),1.7(1)
74-W-186,28.6,0,-0.73(4),,,,0.065(7),0,0.065(7),37.9(6)
75-Re,,,9.2(2),,,,10.6(5),0.9(6),11.5(3),89.7(1.0)
75-Re-185,37.5,5/2,9.0(3),,,,10.2(7),0.5(9),10.7(6),112.0(2.0)
75-Re-187,62.5,5/2,9.3(3),,,,10.9(7),1.0(6),11.9(4),76.4(1.0)
76-Os,,,10.7(2),,,,14.4(5),0.3(8),14.7(6),16.0(4.0)
76-Os-184,0.02,0,10.0(2.0)*,,,,13.0(5.0),0,13.0(5.0),3000.0(150.0)
76-Os-186,1.6,0,12.0(1.7),,,,17.0(5.0),0,17.0(5.0),80.0(13.0)
76-Os-187,1.6,1/2,10.0(2.0)*,,,,13.0(5.0),0.3(3)*,13.0(5.0),320.0(10.0)
76-Os-188,13.3,0,7.8(3),,,,7.3(6),0,7.3(6),4.7(5)
76-Os-189,16.1,3/2,11.0(3),,,,14.4(8),0.5(5)*,14.9(9),25.0(4.0)
76-Os-190,26.4,0,11.4(3),,,,15.2(8),0,15.2(8),13.1(3)
76-Os-192,41,0,11.9(4),,,,16.6(1.2),0,16.6(1.2),2.0(1)
77-Ir,,,10.6(3),,,,14.1(8),0.0(3.0),14.0(3.0),425.0(2.0)
77-Ir-191,37.4,3/2,,,,,,,,954.0(10.0)
77-Ir-193,62.6,3/2,,,,,,,,111.0(5.0)
78-Pt,,,9.60(1),,,,11.58(2),0.13(11),11.71(11),10.3(3)
78-Pt-190,0.01,0,9.0(1.0),,,,10.0(2.0),0,10.0(2.0),152.0(4.0)
78-Pt-192,1.78,0,9.9(5),,,,12.3(1.2),0,12.3(1.2),10.0(2.5)
78-Pt-194,32.9,0,10.55(8),,,,14.0(2),0,14.0(2),1.44(19)
78-Pt-195,33.8,1/2,8.91(9),9.5(3),7.2(3),+/-,9.8(2),0.13(4),9.9(2),27.5(1.2)
78-Pt-196,25.3,0,9.89(8),,,,12.3(2),0,12.3(2),0.72(4)
78-Pt-198,7.2,0,7.8(1),,,,7.6(2),0,7.6(2),3.66(19)
79-Au-197,100,3/2,7.90(7),6.26(10),9.90(14),+/-,7.32(12),0.43(5),7.75(13),98.65(9)
80-Hg,,,12.595(45),,,,20.24(5),6.6(1),26.8(1),372.3(4.0)
80-Hg-196,0.15,0,30.3(1.0),,,E,115.0(8.0),0,115.0(8.0),3080.0(180.0)
80-Hg-198,10.1,0,,,,,,0,,2.0(3)
80-Hg-199,16.9,0,16.9(4),,,E,36.0(2.0),30.0(3.0),66.0(2.0),2150.0(48.0)
80-Hg-200,23.1,0,,,,,,0,,<60.0
80-Hg-201,13.2,3/2,,,,,,,,7.8(2.0)
80-Hg-202,29.7,0,11.002(43),,,,15.2108(2),0,15.2108(2),4.89(5)
80-Hg-204,6.8,0,,,,,,0,,0.43(10)
81-Tl,,,8.776(5),,,,9.678(11),0.21(15),9.89(15),3.43(6)
81-Tl-203,29.5,1/2,8.51(8),9.08(10),6.62(10),,6.14(28),0.14(4),6.28(28),11.4(2)
81-Tl-205,70.5,1/2,8.87(7),5.15(10),9.43(10),+/-,11.39(17),0.007(1),11.40(17),0.104(17)
82-Pb,,,9.401(2),,,,11.115(7),0.0030(7),11.118(7),0.171(2)
82-Pb-204,1.4,0,10.893(78),,,,12.3(2),0,12.3(2),0.65(7)
82-Pb-206,24.1,0,9.221(78),,,,10.68(12),0,10.68(12),0.0300(8)
82-Pb-207,22.1,1/2,9.286(16),,,+/-,10.82(9),0.002(2),10.82(9),0.699(10)
82-Pb-208,52.4,0,9.494(30),,,,11.34(5),0,11.34(5),0.00048(3)
83-Bi-209,100,9/2,8.532(2),8.26(1),8.74(1),,9.148(4),0.0084(19),9.156(4),0.0338(7)
88-Ra-226,1620 Y,0,10.0(1.0),,,,13.0(3.0),0,13.0(3.0),12.8(1.5)
90-Th-232,100,0,10.31(3),,,,13.36(8),0,13.36(8),7.37(6)
91-Pa-231,32500 Y,3/2,9.1(3),,,,10.4(7),0.1(3.3),10.5(3.2),200.6(2.3)
92-U,,,8.417(5),,,,8.903(11),0.005(16),8.908(11),7.57(2)
92-U-233,159000 Y,5/2,10.1(2),,,,12.8(5),0.1(6),12.9(3),574.7(1.0)
92-U-234,0.005,0,12.4(3),,,,19.3(9),0,19.3(9),100.1(1.3)
92-U-235,0.72,7/2,10.50(3),,,,13.78(11),0.2(2),14.0(2),680.9(1.1)
92-U-238,99.27,0,8.407(7),,,,8.871(11),0,8.871(11),2.68(2)
93-Np-237,2140000 Y,5/2,10.55(10),,,,14.0(3),0.5(5)*,14.5(6),175.9(2.9)
94-Pu-239,24400 Y,1/2,7.7(1),,,,7.5(2),0.2(6),7.7(6),1017.3(2.1)
94-Pu-240,6540 Y,0,3.5(1),,,,1.54(9),0,1.54(9),289.6(1.4)
94-Pu-242,376000 Y,0,8.1(1),,,,8.2(2),0,8.2(2),18.5(5)
95-Am-243,7370 Y,5/2,8.3(2),,,,8.7(4),0.3(2.6),9.0(2.6),75.3(1.8)
96-Cm-244,17.9 Y,0,9.5(3),,,,11.3(7),0,11.3(7),16.2(1.2)
96-Cm-246,4700 Y,0,9.3(2),,,,10.9(5),0,10.9(5),1.36(17)
96-Cm-248,340000 Y,0,7.7(2),,,,7.5(4),0,7.5(4),3.00(26)\
"""

# Imaginary values for select isotopes
# isotope, b_c_i, bp_i, bm_i
nsftableI="""\
2-He-3,-1.48,,-5.925
3-Li-6,-0.26,-0.08(1),-0.62(2)
5-B,-0.21,,
47-Ag-107,-0.01,,
47-Ag-109,-0.025,,
48-Cd,-1.2,,
48-Cd-113,-12,,
49-In,-0.054,,
49-In-115,-0.056,,
52-Te-123,-0.1,,
62-Sm,-1.5,,
62-Sm-149,-11,,
64-Gd,-13.6,,
64-Gd-155,-10.3,,
71-Lu-176,-0.57(2),,
80-Hg-196,-0.8,,\
"""
# Excluding the following because the measurements for the real parts
# were not used in nsftable table.
# 63-Eu-151,-2.46,,
# 64-Gd-157,-47,-75,
def _diff(s,a,b):
    print "%10s %8.2f %8.2f %3.0f%%"%(s, a, b, 100*abs((a-b)/b))

def absorption_comparison_table(table=None):
    """
    Print a table of b_c_i and -absorption/(2*1.798*1000) for each isotope.
    
    The factor of 100 is required to go from barn to fm^2.
    
    This author does not know where the remaining factor of 10 comes from.

    This is useful for checking the integrity of the data and formula.
    """
    table = default_table(table)
    print "Comparison of b_c_i and absorption where b_c_i exists"
    for el in table:
        if el.neutron.b_c_i is not None:
            _diff(el, el.neutron.b_c_i, -el.neutron.absorption/2/1.798/1000)
        for iso in el:
            if iso.neutron.b_c_i is not None:
                _diff(iso, iso.neutron.b_c_i, -iso.neutron.absorption/2/1.798/1000)

def coherent_comparison_table(table=None):
    """
    Print a table of 4*pi*b_c**2/100 and coherent for each isotope.

    This is useful for checking the integrity of the data and formula.
    """
    table = default_table(table)
    import numpy
    print "Comparison of b_c and coherent where b_c exists"
    for el in table:
        if el.neutron.b_c is not None:
            _diff(el, el.neutron.b_c**2*4*numpy.pi/100, el.neutron.coherent)
        for iso in el:
            if iso.neutron.b_c is not None:
                _diff(iso, iso.neutron.b_c**2*4*numpy.pi/100, iso.neutron.coherent)


def sld_plot(table=None):
    """
    Plot SLD as a function of element number.
    """
    table = default_table(table)
    import pylab

    SLDs = [(el.number,el.neutron.sld()[0],el.symbol)
            for el in table
            if el.neutron.has_sld()]
    for Z,sld,sym in SLDs:
        if sld is not None: pylab.text(Z,sld,sym)
    pylab.axis([0,100,-4,10])
    pylab.xlabel('Element number')
    pylab.ylabel('Scattering length density (10**-6 Nb)')
    pylab.title('Neutron SLD for elements in natural abundance')
