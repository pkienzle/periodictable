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
    not available, then the resulting density will not be computed.  The
    density calculation assumes the cell volume remains constant for the
    original materials, which is not in general the case.
    """
    table = default_table(kw.pop('table',None))
    density = kw.pop('density',None)
    natural_density = kw.pop('natural_density',None)
    name = kw.pop('name',None)
    if kw: raise TypeError("Unexpected arguments "+", ".join(kw.keys()))

    if len(args)%2 != 0:
        raise ValueError("Need a quantity for each formula")
    pairs = [(formula(args[i],table=table),args[i+1])
             for i in range(0, len(args), 2)]
    result = _mix_by_weight_pairs(pairs)
    if natural_density: result.natural_density = natural_density
    if density: result.density = density
    if name: result.name = name
    return result

def _mix_by_weight_pairs(pairs):

    # Drop pairs with zero quantity
    # Note: must be first statement in order to accept iterators
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
        if all(f.density for f,_ in pairs):
            volume = sum(q/f.density for f,q in pairs)/scale
            result.density = result.mass/volume
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
    not available, then a ValueError is raised. The  density calculation
    assumes the cell volume remains constant for the original materials,
    which is not in general the case.
    """
    table = default_table(kw.pop('table',None))
    density = kw.pop('density',None)
    natural_density = kw.pop('natural_density',None)
    name = kw.pop('name',None)
    if kw: raise TypeError("Unexpected arguments "+", ".join(kw.keys()))

    if len(args)%2 != 0:
        raise ValueError("Need a quantity for each formula")
    pairs = [(formula(args[i],table=table),args[i+1])
             for i in range(0, len(args), 2)]
    result = _mix_by_volume_pairs(pairs)
    if natural_density: result.natural_density = natural_density
    if density: result.density = density
    if name: result.name = name
    return result

def _mix_by_volume_pairs(pairs):

    # Drop pairs with zero quantity
    # Note: must be first statement in order to accept iterators
    pairs = [(f,q) for f,q in pairs if q > 0]

    if not all(f.density for f,_ in pairs):
        raise ValueError("Need a density for each formula")

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

        volume = sum(q for _,q in pairs)/scale
        result.density = result.mass/volume

    return result

def _mix_by_layer_triples(pairs, substrate, sub_vol):

    tunits2mm = { 'nm': 1e-6,
              'um': 1e-3,
              'mm': 1}
    # Drop pairs with zero quantity
    # Note: must be first statement in order to accept iterators
    pairs = [(f,t,u) for f,t,u in pairs if t > 0]

    result = Formula()
    if len(pairs) > 0:
        # layer volume: substrae area*thickness
        # layer mass: layer volume * density
        # mass fraction: layer mass/ total mass
        # scale this so that n = 1 for the substrate
        if not all(f.density for f,_,_ in pairs):
            raise ValueError("Need a density for each formula")
        if not substrate.density:
            raise ValueError("Need a density for each formula")

        v_sub = sub_vol[0]*sub_vol[1]*sub_vol[2]
        a_sub = sub_vol[0]*sub_vol[1]
        m_sub = v_sub * substrate.density
        m_total = m_sub
        p2 = []
        for f,t,u in pairs:
            v = a_sub*t*tunits2mm[u]
            m = v * f.density
            m_total += m
            p2.append((m,f))

        scale = 1/((m_sub)/m_total)
        for (m,f) in p2:
            result += scale*(m/m_total)*f
        result += substrate
    return result



def formula(compound=None, density=None, natural_density=None,
            name=None, table=None):
    r"""
    Construct a chemical formula representation from a string, a
    dictionary of atoms or another formula.

    :Parameters:
        *compound* : Formula initializer
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
    if compound == None or compound == '':
        structure = tuple()
    elif isinstance(compound,Formula):
        structure = compound.structure
        if density is None and natural_density is None: 
            density = compound.density
        if not name: name = compound.name
    elif isatom(compound):
        structure = ((1,compound),)
    elif isinstance(compound,dict):
        structure = _convert_to_hill_notation(compound)
    elif _is_string_like(compound):
        try:
            formula = parse_formula(compound, table=table)
            if name: formula.name = name
            if density is not None: formula.density = density
            elif natural_density is not None: formula.natural_density = natural_density
            return formula
        except ValueError as exception:
            raise ValueError(str(exception))
            #print "parsed",compound,"as",self
    else:
        try:
            structure = _immutable(compound)
        except:
            raise ValueError("not a valid chemical formula: "+str(compound))
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
        if natural_density is not None:
            self.natural_density = natural_density
        elif density is not None:
            self.density = density
        elif len(self.atoms) == 1:
            self.density = list(self.atoms.keys())[0].density
        else:
            self.density = None

    @property
    def atoms(self):
        """
        { *atom*: *count*, ... }

        Composition of the molecule.  Referencing this attribute computes
        the *count* as the total number of each element or isotope in the
        chemical formula, summed across all subgroups.
        """
        return _count_atoms(self.structure)

    @property
    def hill(self):
        """
        Formula

        Convert the formula to a formula in Hill notation.  Carbon appears
        first followed by hydrogen then the remaining elements in alphabetical
        order.
        """
        return formula(self.atoms)

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

    @property
    def natural_density(self):
        """
        |g/cm^3|

        Density of the formula with specific isotopes of each element
        replaced by the naturally occurring abundance of the element
        without changing the cell volume.
        """
        return self.density*self.natural_mass_ratio()

    @natural_density.setter
    def natural_density(self, natural_density):
        self.density = natural_density / self.natural_mass_ratio()

    @property
    def mass(self):
        """
        atomic mass units u (C[12] = 12 u)

        Molar mass of the molecule.  Use molecular_mass to get the mass in
        grams.
        """
        mass = 0
        for el,count in self.atoms.items():
            mass += el.mass*count
        return mass

    @property
    def molecular_mass(self):
        """
        g

        Mass of the molecule in grams.
        """
        return self.mass/avogadro_number

    @property
    def charge(self):
        """
        Net charge of the molecule.
        """
        return sum([m*a.charge for a,m in self.atoms.items()])

    @property
    def mass_fraction(self):
        """
        Fractional mass representation of each element/isotope/ion.
        """
        total_mass = self.mass
        return dict((a,m*a.mass/total_mass) for a,m in self.atoms.items())

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

        Using the cell volume, mass density can be set with::

            formula.density = n*formula.molecular_mass/formula.volume()

        where n is the number of molecules per unit cell.
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
    composite = Forward()
    mixture = Forward()
    thinlayer = Forward()
    layers = Forward()

    # whitespace and separators
    space = Optional(White().suppress())
    separator = space+Literal('+').suppress()+space

    # Lookup the element in the element table
    symbol = Regex("[A-Z][a-z]*")
    symbol = symbol.setParseAction(lambda s,l,t: table.symbol(t[0]))

    # Translate isotope
    openiso = Literal('[').suppress()
    closeiso = Literal(']').suppress()
    isotope = Optional(~White()+openiso+Regex("[1-9][0-9]*")+closeiso,
                       default='0')
    isotope = isotope.setParseAction(lambda s,l,t: int(t[0]) if t[0] else 0)

    # Translate ion
    openion = Literal('{').suppress()
    closeion = Literal('}').suppress()
    ion = Optional(~White()
                     +openion
                     +Regex("([1-9][0-9]*)?[+-]")
                     +closeion,
                   default='0+')
    ion = ion.setParseAction(lambda s,l,t: int(t[0][-1]+(t[0][:-1] if len(t[0])>1 else '1')))

    # Translate counts
    fract = Regex("(0|[1-9][0-9]*|)([.][0-9]*)")
    fract = fract.setParseAction(lambda s,l,t: float(t[0]) if t[0] else 1)
    whole = Regex("[1-9][0-9]*")
    whole = whole.setParseAction(lambda s,l,t: int(t[0]) if t[0] else 1)
    count = Optional(~White()+(fract|whole),default=1)

    # Convert symbol,isotope,ion,count to (count,isotope)
    element = symbol+isotope+ion+count
    def convert_element(string,location,tokens):
        #print "convert_element received",tokens
        symbol,isotope,ion,count = tokens[0:4]
        if isotope != 0: symbol = symbol[isotope]
        if ion != 0: symbol = symbol.ion[ion]
        return (count,symbol)
    element = element.setParseAction(convert_element)

    # Convert "count elements" to a pair
    implicit_group = count+OneOrMore(element)
    def convert_implicit(string,location,tokens):
        #print "implicit",tokens
        count = tokens[0]
        fragment = tokens[1:]
        return fragment if count==1 else (count,fragment)
    implicit_group = implicit_group.setParseAction(convert_implicit)

    # Convert "(composite) count" to a pair
    opengrp = space + Literal('(').suppress() + space
    closegrp = space + Literal(')').suppress() + space
    explicit_group = opengrp + composite + closegrp + count
    def convert_explicit(string,location,tokens):
        #print "explicit",tokens
        count = tokens[-1]
        fragment = tokens[:-1]
        return fragment if count == 1 else (count,fragment)
    explicit_group = explicit_group.setParseAction(convert_explicit)

    # Build composite from a set of groups
    group = implicit_group | explicit_group
    implicit_separator = separator | space
    composite << group + ZeroOrMore(implicit_separator + group)

    density = Literal('@').suppress() + count + Optional(Regex("[ni]"),default='i')
    compound = composite + Optional(density,default=None)
    def convert_compound(string,location,tokens):
        #print "compound",tokens
        if tokens[-1] is None:
            return Formula(structure=_immutable(tokens[:-1]))
        elif tokens[-1] == 'n': 
            return Formula(structure=_immutable(tokens[:-2]), natural_density=tokens[-2])
        else:
            return Formula(structure=_immutable(tokens[:-2]), density=tokens[-2])
    compound = compound.setParseAction(convert_compound)

    partsep = space + Literal('//').suppress() + space
    percent = Literal('%').suppress()

    weight_percent = Regex("%(w((eigh)?t)?|m(ass)?)").suppress() + space
    by_weight = count + weight_percent + mixture + ZeroOrMore(partsep+count+(weight_percent|percent)+mixture) + partsep + mixture
    def convert_by_weight(string,location,tokens):
        #print "by weight",tokens
        piece = tokens[1:-1:2] + [tokens[-1]]
        fract = [float(v) for v in tokens[:-1:2]]
        fract.append(100-sum(fract))
        #print piece, fract
        if len(piece) != len(fract): raise ValueError("Missing base component of mixture")
        if fract[-1] < 0: raise ValueError("Formula percentages must sum to less than 100%")
        return _mix_by_weight_pairs(zip(piece,fract))
    mixture_by_weight = by_weight.setParseAction(convert_by_weight)

    volume_percent = Regex("%v(ol(ume)?)?").suppress() + space
    by_volume = count + volume_percent + mixture + ZeroOrMore(partsep+count+(volume_percent|percent)+mixture) + partsep + mixture 
    def convert_by_volume(string,location,tokens):
        #print "by volume",tokens
        piece = tokens[1:-1:2] + [tokens[-1]]
        fract = [float(v) for v in tokens[:-1:2]]
        fract.append(100-sum(fract))
        #print piece, fract
        if len(piece) != len(fract): raise ValueError("Missing base component of mixture "+string)
        if fract[-1] < 0: raise ValueError("Formula percentages must sum to less than 100%")
        return _mix_by_volume_pairs(zip(piece,fract))
    mixture_by_volume = by_volume.setParseAction(convert_by_volume)

    layer_mul = Literal("*").suppress() + count + Optional(space)
    layer_base = Literal('#').suppress() + count + layer_mul + layer_mul + Optional(space)
    layer_thick = count +Regex("(nm|um|mm)") + space
    layer_start = Literal('[').suppress() + Optional(space)
    layer_end = layer_base + Literal(']').suppress()

    layer_spec = layer_thick + mixture + ZeroOrMore(partsep + layer_thick + mixture)
    layer_group = opengrp + layer_spec + closegrp + count
    def convert_layercount(string,location,tokens):
        # print "layercount",tokens
        count = tokens[-1]
        fragment = tokens[:-1]
        if count == 1:
            return fragment
        else:
            return fragment * count

    layer_group = layer_group.setParseAction(convert_layercount)

    layers << (layer_spec|layer_group) + ZeroOrMore(partsep + (layer_spec|layer_group))
    by_layer = layer_start + layers + partsep + mixture + layer_end
    def convert_by_layer(string,location,tokens):
        try:
            # print "by layers",tokens
            piece = tokens[2:-4:3]
            fract = [float(v) for v in tokens[:-4:3]]
            tunit = tokens[1:-4:3]
            substrat =  tokens[-4]
            sub_vol =  tokens[-3:]
            if len(piece) != len(fract) : raise ValueError("Missing substrate in layer definition"+string)
            return _mix_by_layer_triples(zip(piece,fract,tunit), substrat, sub_vol)
        except Exception as e:
            print e
            raise
    mixture_by_layer = by_layer.setParseAction(convert_by_layer)

    mixture << (compound | (opengrp + (mixture_by_weight | mixture_by_volume | mixture_by_layer) + closegrp))
    formula = compound | mixture_by_weight | mixture_by_volume | mixture_by_layer
    grammar = Optional(formula, default=Formula()) + StringEnd()

    grammar.setName('Chemical Formula')
    return grammar

_PARSER_CACHE = {}
def parse_formula(formula_str, table = None):
    """
    Parse a chemical formula, returning a structure with elements from the
    given periodic table.
    """
    table = default_table(table)
    if table not in _PARSER_CACHE:
        _PARSER_CACHE[table] = formula_grammar(table)
    return _PARSER_CACHE[table].parseString(formula_str)[0]

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
        for el,elcount in partial.items():
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
        if b.symbol in ("C","H"):
            return 1
        else:
            return cmp(a.symbol, b.symbol)

def _hill_key(a):
    return "".join((("0" if a.symbol in ("C","H") else "1"),
                    a.symbol,
                    "%4d"%(a.isotope if isisotope(a) else 0)))

def _convert_to_hill_notation(atoms):
    """
    Return elements listed in standard order.
    """
    #return [(atoms[el], el) for el in sorted(atoms.keys(), cmp=_hill_compare)]
    return [(atoms[el], el) for el in sorted(atoms.keys(), key=_hill_key)]


def _str_atoms(seq):
    """
    Convert formula structure to string.
    """
    #print "str",seq
    ret = ""
    for count,fragment in seq:
        if isatom(fragment):
            # Normal isotope string form is #-Yy, but we want Yy[#]
            if isisotope(fragment) and 'symbol' not in fragment.__dict__:
                ret += "%s[%d]"%(fragment.symbol,fragment.isotope)
            else:
                ret += fragment.symbol
            if fragment.charge != 0:
                sign = '+' if fragment.charge > 0 else '-'
                value = str(abs(fragment.charge)) if abs(fragment.charge)>1 else ''
                ret += '{'+value+sign+'}'
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
