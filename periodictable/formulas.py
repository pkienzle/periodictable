# This program is public domain
"""

As of this writing, this module includes formula for parsing, computing molar mass and
computing scattering length density.
"""

from copy import copy
from math import pi, sqrt

# Requires that the pyparsing module is installed.

from pyparsing import (Literal, Optional, White, Regex,
                       ZeroOrMore, OneOrMore, Forward, StringEnd)

from .core import Element, Isotope, default_table

PACKING_FACTORS = dict(cubic=pi/6, bcc=pi*sqrt(3)/8, hcp=pi/sqrt(18),
                       fcc=pi/sqrt(18), diamond=pi*sqrt(3)/16)


class Formula(object):
    r"""
    Simple chemical formula representation.
    This is designed for calculating molar mass and scattering
    length density, not for representing bonds or atom positions.
    We preserve the structure of the formula so that it can
    be used as a basis for a rich text representation such as
    `matplotlib TeX markup <http://matplotlib.sourceforge.net/users/mathtext.html>`_.

    :Parameters:
        *formula* : see below
            Chemical formula.

        *density* : float | g/cm**3
            Material density.

        *name* : string
            Common name for the molecule.
    
    :Exceptions:
        *ValueError* : invalid formula initializer

    Formula initializers can have a variety of forms:

    * string:
        m = Formula( "CaCO3+6H2O" )

        For full details see :ref:`Formula grammar <formula>`
    
    * structure:
        m = Formula( [(1,Ca),(2,C),(3, O),(6,[(2,H),(1,O)]] )

    * formula math:
        m = Formula( "CaCO3" ) + 6*Formula( "H2O" )

    * another formula (makes a copy):
        m = Formula( Formula("CaCO3+6H2O") )

    * an atom:
        m = Formula(Ca)

    * nothing:
        m = Formula()
    """
    def __init__(self,value=None,density=None,name=None):
        self.density,self.name = None,None
        if value == None or value == '':
            self.structure = []
        elif isinstance(value,Formula):
            self.__dict__ = copy(value.__dict__)
        elif _isatom(value):
            self.structure = ((1,value),)
        elif isinstance(value,dict):
            self.structure = _convert_to_hill_notation(value)
        elif _is_string_like(value):
            try:
                self._parse_string(value)
            except ValueError,msg:
                raise ValueError,msg
            #print "parsed",value,"as",self
        else:
            try:
                self._parse_structure(value)
            except:
                raise ValueError("not a valid chemical formula")
            
        if density: self.density = density
        if name: self.name = name

        if self.density is None:
            # Density defaults to element density if only one element.
            # Users can estimate the density as:
            #    self.density = self.mass/self.volume()
            # if they feel so inclined, but it is too misleading to
            # make that assumption by default.
            if len(self.atoms) == 1:
                self.density = self.atoms.keys()[0].density

    def _parse_structure(self,input):
        """
        Set the formula to the given structure, checking that it is valid.
        """
        _check_atoms(input)
        self.structure = _immutable(input)

    def _parse_string(self,input):
        """
        Parse the formula string.
        """
        self.structure = _immutable(parser.parseString(input))

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
        return Formula(self.atoms)
    hill = property(_hill, doc=_hill.__doc__)

    def _mass(self):
        """
        atomic mass units u (C[12] = 12 u)

        Atomic mass of the molecule.
        
        Referencing this attribute computes the mass of the chemical formula.
        """
        mass = 0
        for el,count in self.atoms.iteritems():
            mass += el.mass*count
        return mass
    mass = property(_mass,doc=_mass.__doc__)

    def _mass(self):
        """
        atomic mass units u (C[12] = 12 u)

        Atomic mass of the molecule.
        
        Referencing this attribute computes the mass of the chemical formula.
        """
        mass = 0
        for el,count in self.atoms.iteritems():
            mass += el.mass*count
        return mass
    mass = property(_mass,doc=_mass.__doc__)


    def volume(self, packing_factor='hcp'):
        """
        Estimate molecular volume. 

        The crystal volume can be estimated from the element covalent radius 
        and the atomic packing factor using::
        
            packing_factor = N_atoms V_atom / V_crystal

        Packing factors for a number of crystal lattice structures are defined.

        .. table:: Crystal lattice names and packing factors
            
            ========   =======================   =============    ==============
            Code       Description               Formula          Packing factor
            ========   =======================   =============    ==============
            cubic      simple cubic              pi 1/6           0.52360
            bcc        body-centered cubic       pi sqrt(3)/8     0.68017
            hcp        hexagonal close-packed    pi 1/sqrt(18)    0.74048
            fcc        face-centered cubic       pi 1/sqrt(18)    0.74048
            diamond    diamond cubic             pi sqrt(3)/16    0.34009
            ========   =======================   =============    ==============
         
        :Parameters:
            *packing_factor*  = 'hcp' : float or string
                Atomic packing factor.  If *packing_factor* is the name of
                a crystal lattice, use the *lattice* packing factor.
        
        :Returns: *volume* : float | A^3 
             Molecular volume. 

        :Raises:

             *ValueError* : *lattice* is not defined
        """
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
        return V/packing_factor

    # TODO: move neutron/xray sld to extension
    def neutron_sld(self, wavelength=1):
        """
        Neutron scattering information for the molecule.

        :Parameters:
            *wavelength* : float | A
                Wavelength of the neutron beam.
        
        :Returns: *sld* : (float, float, float) | 10^-6 inv A^2
 
             Neutron scattering length density is returned as the tuple
             (*real*, *imaginary*, *incoherent*), or as (None, None, None) 
             if the mass density is not known.

        """
        if self.density is None: return None,None,None
        from nsf import neutron_sld_from_atoms
        return neutron_sld_from_atoms(self.atoms,density=self.density,
                                      wavelength=wavelength)

    def xray_sld(self, energy=None, wavelength=None):
        """
        X-ray scattering length density for the molecule.

        :Parameters:

            *energy* : float | keV
                Energy of atom.

            *wavelength* : float | A
                Wavelength of atom.

            .. Note: One of *wavelength* or *energy* is required.

        :Returns: *sld* : (float, float) | inv A^2
                X-ray scattering length density is return as the tuple
                (*real*, *imaginary*), or as (None, None) if the mass 
                density is not known.

        """
        if self.density is None: return None,None
        from xsf import xray_sld_from_atoms
        return xray_sld_from_atoms(self.atoms,density=self.density,
                                   wavelength=wavelength,energy=energy)

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

    def __rmul__(self,other):
        """
        Provide a multiplier for formula.
        """
        #print "multiplying",self,other
        try:
            other += 0
        except TypeError:
            raise TypeError("n*formula expects numeric n")
        ret = Formula()
        if self.structure != []:
            ret.structure = ((other,self.structure),)
        return ret

    def __str__(self):
        return self.name if self.name else _str_atoms(self.structure)

    def __repr__(self):
        return "formula('%s')"%(str(self))

    def __getstate__(self):
        """
        Pickle formula structure using a string representation since
        elements can't be pickled.
        """
        output = copy(self.__dict__)
        output['structure'] = _str_atoms(self.structure)
        return output

    def __setstate__(self,input):
        """
        Unpickle formula structure from a string representation since
        elements can't be pickled.
        """
        self.__dict__ = input
        self._parse_string(input['structure'])

def formula_grammar(table=None):
    """
    Construct a parser for molecular formulas. 

    :Parameters:
        *table* = None : PeriodicTable
             If table is specified, then elements and their associated fields
             will be chosen from that periodic table rather than the default.
    :Returns: *parser* : pyparsing.ParserElement
         The ``parser.parseString()`` method returns a list of 
         pairs (*count,fragment*), where fragment is an *isotope*, 
         an *element* or a list of pairs (*count,fragment*).    
                
    """
    table = default_table(table)

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
        if isotope != 0:
            try:
                symbol = symbol[isotope]
            except KeyError:
                raise ValueError("%d is not an isotope of %s"%(isotope,symbol))
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
    separator = Literal('+').suppress()
    formula << group + ZeroOrMore(Optional(separator)+group)
    grammar = Optional(formula) + StringEnd()

    return grammar

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

def _check_atoms(seq):
    """
    Traverse formula structure, checking that the counts are numeric and
    units are elements or isotopes.  Raises an error if this is not the
    case.
    """
    for count,fragment in seq:
        count += 0 # Fails if not numeric
        if not _isatom(fragment):
            _check_atoms(fragment) # Fails if not sequence of pairs

def _immutable(seq):
    """Converts lists to tuples so that structure is immutable."""
    if _isatom(seq): return seq
    return tuple( (count,_immutable(fragment)) for count,fragment in seq )

def _hill_compare(a,b):
    """
    Compare elements in standard order.
    """
    if a.symbol == b.symbol:
        a = a.isotope if isinstance(a, Isotope) else 0
        b = b.isotope if isinstance(b, Isotope) else 0
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
        if _isatom(fragment):
            # Isotopes are Sym[iso] except for D and T
            if _isisotope(fragment) and 'symbol' not in fragment.__dict__:
                ret += "%s[%d]"%(fragment.symbol,fragment.isotope)
            else:
                ret += fragment.symbol
            if count!=1:
                ret += "%g"%count
        else:
            if count == 1:
                ret += _str_atoms(fragment)
            else:
                ret += "(%s)%g"%(_str_atoms(fragment),count)
    return ret

def _isatom(val):
    """Return true if value is an element or isotope"""
    return isinstance(val,(Element,Isotope))

def _isisotope(val):
    """Return true if value is an isotope"""
    return isinstance(val,Isotope)

def _is_string_like(val):
    """Returns True if val acts like a string"""
    try: val+''
    except: return False
    return True

parser = formula_grammar()
parser.setName('Chemical Formula')
