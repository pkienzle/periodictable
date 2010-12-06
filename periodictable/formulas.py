# This program is public domain
# Author: Paul Kienzle
"""
Chemical formula parser.
"""
from __future__ import division
from copy import copy
from math import pi, sqrt

# Requires that the pyparsing module is installed.

from pyparsing import (Literal, Optional, White, Regex,
                       ZeroOrMore, OneOrMore, Forward, StringEnd)

from .core import default_table, isatom, isisotope, change_table
from .constants import avogadro_number
from .util import require_keywords, cell_volume

PACKING_FACTORS = dict(cubic=pi/6, bcc=pi*sqrt(3)/8, hcp=pi/sqrt(18),
                       fcc=pi/sqrt(18), diamond=pi*sqrt(3)/16)


def mix_by_weight(*args, **kw):
    """
    Generate a mixture which apportions each formula by weight.

    :Parameters:
        *formula1* : Formula OR string
            Material

        *quantity1* : float
            Relative quantity of that material

        *formula2* : Formula OR string
            Material

        *quantity2* : float
            Relative quantity of that material

        ...

        *density* : float
            Density of the mixture, if known

        *natural_density* : float
            Density of the mixture with natural abundances, if known.

        *name* : string
            Name of the mixture

        *table* : PeriodicTable
            Private table to use when parsing string formulas.

    :Returns:

        *formula* : Formula

    If density is not given, then it will be computed from the density
    of the components, assuming the components take up no more nor less
    space because they are in the mixture.  If component densities are
    not available, then the resulting density will not be computed.
    """
    table = default_table(kw.pop('table',None))
    density = kw.pop('density',None)
    natural_density = kw.pop('natural_density',None)
    name = kw.pop('name',None)
    if kw: raise TypeError("Unexpected keyword "+kw.keys()[0])

    if len(args)%2 != 0:
        raise ValueError("Need a quantity for each formula")
    pairs = [(formula(args[i],table=table),args[i+1])
             for i in range(0, len(args), 2)]

    # Drop pairs with zero quantity
    pairs = [(f,q) for f,q in pairs if q > 0]

    result = Formula()
    if len(pairs) > 0:
        # cell mass = mass
        # target mass = q
        # cell mass * n = target mass
        #   => n = target mass / cell mass
        #        = q / mass
        # scale this so that n = 1 for the smallest quantity
        scale = min(q/f.mass for f,q in pairs)
        for f,q in pairs:
            result += ((q/f.mass)/scale) * f
    else:
        scale = 0

    if natural_density: result.natural_density = natural_density
    if density: result.density = density
    if scale and not result.density and all(f.density for f,_ in pairs):
        volume = sum(q/f.density for f,q in pairs)/scale
        result.density = result.mass/volume

    if name: result.name = name
    return result

def mix_by_volume(*args, **kw):
    """
    Generate a mixture which apportions each formula by volume.

    :Parameters:
        *formula1* : Formula OR string
            Material

        *quantity1* : float
            Relative quantity of that material

        *formula2* : Formula OR string
            Material

        *quantity2* : float
            Relative quantity of that material

        ...

        *density* : float
            Density of the mixture, if known

        *natural_density* : float
            Density of the mixture with natural abundances, if known.

        *name* : string
            Name of the mixture

        *table* : PeriodicTable
            Private table to use when parsing string formulas.

    :Returns:

        *formula* : Formula

    If density is not given, then it will be computed from the density
    of the components, assuming the components take up no more nor less
    space because they are in the mixture.  If component densities are
    not available, then the resulting density will not be computed.
    """
    table = default_table(kw.pop('table',None))
    density = kw.pop('density',None)
    natural_density = kw.pop('natural_density',None)
    name = kw.pop('name',None)
    if kw: raise TypeError("Unexpected keyword "+kw.keys()[0])

    if len(args)%2 != 0:
        raise ValueError("Need a quantity for each formula")
    pairs = [(formula(args[i],table=table),args[i+1])
             for i in range(0, len(args), 2)]

    if not all(f.density for f,_ in pairs):
        raise ValueError("Need a density for each formula")

    # Drop pairs with zero quantity
    pairs = [(f,q) for f,q in pairs if q > 0]

    result = Formula()
    if len(pairs) > 0:
        # cell volume = mass/density
        # target volume = q
        # cell volume * n = target volume
        #   => n = target volume / cell volume
        #        = q / (mass/density)
        #        = q * density / mass
        # scale this so that n = 1 for the smallest quantity
        scale = min(q*f.density/f.mass for f,q in pairs)
        for f,q in pairs:
            result += ((q*f.density/f.mass)/scale) * f
    else:
        scale = 0

    if natural_density: result.natural_density = natural_density
    if density: result.density = density
    if scale and not result.density:
        volume = sum(q for _,q in pairs)/scale
        result.density = result.mass/volume

    if name: result.name = name
    return result

def formula(value=None, density=None, natural_density=None,
            name=None, table=None):
    r"""
    Construct a chemical formula representation from a string, a
    dictionary of atoms or another formula.

    :Parameters:
        *formula* : see below
            Chemical formula.

        *density* : float | |g/cm^3|
            Material density.  Not needed for single element formulas.

        *natural_density* : float | |g/cm^3|
            Material density assuming naturally occurring isotopes and no
            change in cell volume.

        *name* : string
            Common name for the molecule.

        *table* : PeriodicTable
            Private table to use when parsing string formulas.

    :Exceptions:
        *ValueError* : invalid formula initializer

    After creating a formula, a rough estimate of the density can be
    computed using:

         formula.density = formula.molecular_mass/formula.volume(packing_factor=...)

    The volume() calculation uses the covalent radii of the components and
    the known packing factor or crystal structure name.  If the lattice
    constants for the crystal are known, then they can be used instead:

        formula.density = formula.molecular_mass/formula.volume(a,b,c,alpha,beta,gamma)

    Formulas are designed for calculating quantities such as molar mass and
    scattering length density, not for representing bonds or atom positions.
    The representations are simple, but preserve some of the structure for
    display purposes.
    """
    if value == None or value == '':
        structure = tuple()
    elif isinstance(value,Formula):
        structure = value.structure
        if not density and not natural_density: density = value.density
        if not name: name = value.name
    elif isatom(value):
        structure = ((1,value),)
    elif isinstance(value,dict):
        structure = _convert_to_hill_notation(value)
    elif _is_string_like(value):
        try:
            structure = _immutable(parse_formula(value, table=table))
        except ValueError,exception:
            raise ValueError(str(exception))
            #print "parsed",value,"as",self
    else:
        try:
            structure = _immutable(value)
        except:
            raise ValueError("not a valid chemical formula: "+str(value))
    return Formula(structure=structure, name=name, density=density,
                   natural_density=natural_density)

class Formula(object):
    """
    Simple chemical formula representation.

    """
    def __init__(self, structure=tuple(), density=None, natural_density=None,
                 name=None):
        self.structure = structure
        self.name = name

        # If natural_density or density are specified, use them.
        # If only one element in the formula, use its density.
        # Otherwise, leave density unspecified, and let the user
        # assign it separately if they need it.
        if natural_density:
            self.natural_density = natural_density
        elif density:
            self.density = density
        elif len(self.atoms) == 1:
            self.density = self.atoms.keys()[0].density
        else:
            self.density = None

    def _atoms(self):
        """
        { *atom*: *count*, ... }

        Composition of the molecule.  Referencing this attribute computes
        the *count* as the total number of each element or isotope in the
        chemical formula, summed across all subgroups.
        """
        return _count_atoms(self.structure)
    atoms = property(_atoms,doc=_atoms.__doc__)


    def _hill(self):
        """
        Formula

        Convert the formula to a formula in Hill notation.  Carbon appears
        first followed by hydrogen then the remaining elements in alphabetical
        order.
        """
        return formula(self.atoms)
    hill = property(_hill, doc=_hill.__doc__)

    def natural_mass_ratio(self):
        """
        Natural mass to isotope mass ratio.

        :Returns:
            *ratio* : float

        The ratio is computed from the sum of the masses of the individual
        elements using natural abundance divided by the sum of the masses
        of the isotopes used in the formula.  If the cell volume is
        preserved with isotope substitution, then the ratio of the masses
        will be the ratio of the densities.
        """
        total_natural_mass = total_isotope_mass = 0
        for el,count in self.atoms.items():
            try:
                natural_mass = el.element.mass
            except AttributeError:
                natural_mass = el.mass
            total_natural_mass += count * natural_mass
            total_isotope_mass += count * el.mass
        return total_natural_mass/total_isotope_mass
    def _get_natural_density(self):
        """
        |g/cm^3|

        Density of the formula with specific isotopes of each element
        replaced by the naturally occurring abundance of the element
        without changing the cell volume.
        """
        return self.density/self.natural_mass_ratio()
    def _set_natural_density(self, natural_density):
        self.density = natural_density / self.natural_mass_ratio()
    natural_density = property(_get_natural_density, _set_natural_density,
                               doc=_get_natural_density.__doc__)
    def _mass(self):
        """
        atomic mass units u (C[12] = 12 u)

        Molar mass of the molecule.  Use molecular_mass to get the mass in
        grams.
        """
        mass = 0
        for el,count in self.atoms.iteritems():
            mass += el.mass*count
        return mass
    mass = property(_mass,doc=_mass.__doc__)

    @property
    def molecular_mass(self):
        """
        g

        Mass of the molecule in grams.
        """
        return self.mass/avogadro_number

    def _pf(self):
        """
        packing factor  | unitless

        packing factor estimated from density.
        """
        return self.density

    def volume(self, packing_factor='hcp', *args, **kw):
        r"""
        Estimate unit cell volume.

        The crystal volume can be estimated from the element covalent radius
        and the atomic packing factor using::

            packing_factor = N_atoms V_atom / V_crystal

        Packing factors for a number of crystal lattice structures are defined.

        .. table:: Crystal lattice names and packing factors

         ======== ======================= ====================== ==============
         Code     Description             Formula                Packing factor
         ======== ======================= ====================== ==============
         cubic    simple cubic            $\pi/6$                0.52360
         bcc      body-centered cubic     $\pi\sqrt{3/8}$        0.68017
         hcp      hexagonal close-packed  $\pi/\sqrt{18}$        0.74048
         fcc      face-centered cubic     $\pi/\sqrt{18}$        0.74048
         diamond  diamond cubic           $\pi\sqrt{3/16}$       0.34009
         ======== ======================= ====================== ==============

        :Parameters:
            *packing_factor*  = 'hcp' : float or string
                Atomic packing factor.  If *packing_factor* is the name of
                a crystal lattice, use the *lattice* packing factor.
            *a*, *b*, *c* : float | |Ang|
                Lattice spacings. *b* and *c* default to *a*.
            *alpha*, *beta*, *gamma* : float | |deg|
                Lattice angles.  These default to 90\ |deg|

        :Returns:

            *volume* : float | |cm^3|
                Molecular volume.

        :Raises:

            *KeyError* : unknown lattice type

            *TypeError* : missing or bad lattice parameters

        Using the cell volume, mass density can be set with:

            formula.density = formula.molecular_mass/formula.volume()

        """
        # Let cell_volume sort out its own parameters.
        if args or kw:
            return cell_volume(*args, **kw)*1e-24

        # Compute atomic volume
        V = 0
        for el,count in self.atoms.items():
            V += el.covalent_radius**3*count
        V *= 4.*pi/3

        # Translate packing factor from string
        try:
            _ = packing_factor + ""
        except:
            pass
        else:
            packing_factor = PACKING_FACTORS[packing_factor.lower()]
        return V/packing_factor*1e-24

    @require_keywords
    def neutron_sld(self, wavelength=None, energy=None):
        """
        Neutron scattering information for the molecule.

        :Parameters:
            *wavelength* : float | |Ang|
                Wavelength of the neutron beam.

        :Returns:

            *sld* : (float, float, float) | |1e-6/Ang^2|
                Neutron scattering length density is returned as the tuple
                (*real*, *imaginary*, *incoherent*), or as (None, None, None)
                if the mass density is not known.

        .. deprecated:: 0.95
            Use periodictable.neutron_sld(formula) instead.
        """
        from .nsf import neutron_sld
        if self.density is None: return None,None,None
        return neutron_sld(self.atoms,density=self.density,
                           wavelength=wavelength, energy=energy)

    @require_keywords
    def xray_sld(self, energy=None, wavelength=None):
        """
        X-ray scattering length density for the molecule.

        :Parameters:
            *energy* : float | keV
                Energy of atom.

            *wavelength* : float | |Ang|
                Wavelength of atom.

            .. Note: One of *wavelength* or *energy* is required.

        :Returns:

            *sld* : (float, float) | |1e-6/Ang^2|
                X-ray scattering length density is returned as the tuple
                    (*real*, *imaginary*), or as (None, None) if the mass
                    density is not known.

        .. deprecated:: 0.95
            Use periodictable.xray_sld(formula) instead.
        """
        from .xsf import xray_sld
        if self.density is None: return None,None
        return xray_sld(self.atoms,density=self.density,
                        wavelength=wavelength,energy=energy)

    def change_table(self, table):
        """
        Replace the table used for the components of the formula.
        """
        self.structure = _change_table(self.structure, table)
        return self

    def __eq__(self,other):
        """
        Return True if two formulas represent the same structure. Note
        that they may still have different names and densities.
        Note: use hill representation for an order independent comparison.
        """
        if not isinstance(other,Formula): return False
        return self.structure==other.structure

    def __add__(self,other):
        """
        Join two formulas.
        """
        #print "adding",self,other
        if not isinstance(other,Formula):
            raise TypeError("expected formula+formula")
        ret = Formula()
        ret.structure = tuple(list(self.structure) + list(other.structure))
        return ret

    def __iadd__(self, other):
        """
        Extend a formula with another.
        """
        self.structure = tuple(list(self.structure) + list(other.structure))
        return self

    def __rmul__(self,other):
        """
        Provide a multiplier for formula.
        """
        #print "multiplying",self,other
        try:
            other += 0
        except TypeError:
            raise TypeError("n*formula expects numeric n")
        ret = copy(self)
        if other != 1 and self.structure:
            if len(self.structure) == 1:
                q,f = self.structure[0]
                ret.structure = ((other*q,f),)
            else:
                ret.structure = ((other,ret.structure),)
        return ret

    def __str__(self):
        return self.name if self.name else _str_atoms(self.structure)

    def __repr__(self):
        return "formula('%s')"%(str(self))


def formula_grammar(table):
    """
    Construct a parser for molecular formulas.

    :Parameters:

        *table* = None : PeriodicTable
             If table is specified, then elements and their associated fields
             will be chosen from that periodic table rather than the default.

    :Returns:
        *parser* : pyparsing.ParserElement.
            The ``parser.parseString()`` method returns a list of
            pairs (*count,fragment*), where fragment is an *isotope*,
            an *element* or a list of pairs (*count,fragment*).

    """
    # Recursive
    formula = Forward()

    # Lookup the element in the element table
    symbol = Regex("[A-Z][a-z]*")
    symbol = symbol.setParseAction(lambda s,l,t: table.symbol(t[0]))

    # Translate isotope
    openiso = Literal('[').suppress()
    closeiso = Literal(']').suppress()
    isotope = Optional(~White()+openiso+Regex("[1-9][0-9]*")+closeiso,
                       default='0')
    isotope = isotope.setParseAction(lambda s,l,t: int(t[0]) if t[0] else 0)

    # Translate counts
    fract = Regex("(0|[1-9][0-9]*|)([.][0-9]*)")
    fract = fract.setParseAction(lambda s,l,t: float(t[0]) if t[0] else 1)
    whole = Regex("[1-9][0-9]*")
    whole = whole.setParseAction(lambda s,l,t: int(t[0]) if t[0] else 1)
    count = Optional(~White()+(fract|whole),default=1)

    # Convert symbol,isotope,count to (count,isotope)
    element = symbol+isotope+count
    def convert_element(string,location,tokens):
        #print "convert_element received",tokens
        symbol,isotope,count = tokens[0:3]
        if isotope != 0: symbol = symbol[isotope]
        return (count,symbol)
    element = element.setParseAction(convert_element)

    # Convert "count elements" to a pair
    implicit_group = count+OneOrMore(element)
    def convert_implicit(string,location,tokens):
        #print "convert_implicit received",tokens
        count = tokens[0]
        fragment = tokens[1:]
        return fragment if count==1 else (count,fragment)
    implicit_group = implicit_group.setParseAction(convert_implicit)

    # Convert "(formula) count" to a pair
    opengrp = Literal('(').suppress()
    closegrp = Literal(')').suppress()
    explicit_group = opengrp + formula + closegrp + count
    def convert_explicit(string,location,tokens):
        #print "convert_group received",tokens
        count = tokens[-1]
        fragment = tokens[:-1]
        return fragment if count == 1 else (count,fragment)
    explicit_group = explicit_group.setParseAction(convert_explicit)

    group = implicit_group | explicit_group
    separator = Optional(Literal('+').suppress()) + Optional(White().suppress())
    formula << group + ZeroOrMore(Optional(White().suppress())+separator+group)
    grammar = Optional(formula) + StringEnd()

    grammar.setName('Chemical Formula')
    return grammar

_PARSER_CACHE = {}
def parse_formula(str, table = None):
    """
    Parse a chemical formula, returning a structure with elements from the
    given periodic table.
    """
    table = default_table(table)
    if table not in _PARSER_CACHE:
        _PARSER_CACHE[table] = formula_grammar(table)
    return _PARSER_CACHE[table].parseString(str)

def _count_atoms(seq):
    """
    Traverse formula structure, counting the total number of atoms.
    """
    total = {}
    for count,fragment in seq:
        if isinstance(fragment,(list,tuple)):
            partial = _count_atoms(fragment)
        else:
            partial = {fragment: 1}
        for el,elcount in partial.iteritems():
            if el not in total: total[el] = 0
            total[el] += elcount*count
    return total

def _immutable(seq):
    """
    Traverse formula structure, checking that the counts are numeric and
    units are atoms.  Returns an immutable copy of the structure, with all
    lists replaced by tuples.
    """
    if isatom(seq):
        return seq
    else:
        return tuple((count+0,_immutable(fragment)) for count,fragment in seq)

def _change_table(seq, table):
    """Converts lists to tuples so that structure is immutable."""
    if isatom(seq):
        return change_table(seq, table)
    else:
        return tuple((count,_change_table(fragment, table))
                     for count,fragment in seq)

def _hill_compare(a,b):
    """
    Compare elements in standard order.
    """
    if a.symbol == b.symbol:
        a = a.isotope if isisotope(a) else 0
        b = b.isotope if isisotope(b) else 0
        return cmp(a,b)
    elif a.symbol in ("C", "H"):
        if b.symbol in ("C", "H"):
            return cmp(a.symbol, b.symbol)
        else:
            return -1
    else:
        return cmp(a.symbol, b.symbol)

def _convert_to_hill_notation(atoms):
    """
    Return elements listed in standard order.
    """
    return [(atoms[el], el) for el in sorted(atoms.keys(), cmp=_hill_compare)]


def _str_atoms(seq):
    """
    Convert formula structure to string.
    """
    #print "str",seq
    ret = ""
    for count,fragment in seq:
        if isatom(fragment):
            # Isotopes are Sym[iso] except for D and T
            if isisotope(fragment) and 'symbol' not in fragment.__dict__:
                ret += "%s[%d]"%(fragment.symbol,fragment.isotope)
            else:
                ret += fragment.symbol
            if count!=1:
                ret += "%g"%count
        else:
            if count == 1:
                piece = _str_atoms(fragment)
            else:
                piece = "(%s)%g"%(_str_atoms(fragment),count)
            #ret = ret+" "+piece if ret else piece
            ret += piece

    return ret

def _is_string_like(val):
    """Returns True if val acts like a string"""
    try: val+''
    except: return False
    return True
