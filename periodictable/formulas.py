# This program is public domain
"""
Operations on chemical formulas.

As of this writing, this includes parsing, computing molar mass and
computing scattering length density.

Requires that the pyparsing module is installed.
"""
from copy import copy

from pyparsing import (Literal, Optional, White, Regex,
                       ZeroOrMore, OneOrMore, Forward, StringEnd)

from .core import elements, Element, Isotope

class Formula(object):
    """
    Simple chemical formula representation.

    This is designed for calculating molar mass and scattering
    length density, not for representing bonds or atom positions.
    We preserve the structure of the formula so that it can
    be used as a basis for a rich text representation such as
    matplotlib TeX markup.

    Example initializers::

       string:
          m = Formula( "CaCO3+6H2O" )
       structure:
          m = Formula( [(1,Ca),(2,C),(3,O),(6,[(2,H),(1,O)]] )
       formula math:
          m = Formula( "CaCO3" ) + 6*Formula( "H2O" )
       another formula (makes a copy):
          m = Formula( Formula("CaCO3+6H2O") )
       an atom:
          m = Formula( Ca )
       nothing:
          m = Formula()

    Additional information can be provided::

       density (g / cm**3)   material density
       name (string) common name for the molecule

    Operations::

        m.atoms returns a dictionary of {isotope: count} for each isotope
          in the formula
        m.mass (u) returns the molar mass of the molecule

    Formula strings consist of counts and atoms such as "CaCO3+6H2O".

    Groups can be separated by '+' or space, so "CaCO3 6H2O" works as well.

    Groups and be defined using parentheses, such as "CaCO3(H2O)6".

    Parentheses can nest: "(CaCO3(H2O)6)1"

    Isotopes are represented by index, e.g., "CaCO[18]3+6H2O".

    Counts can be integer or decimal, e.g. "CaCO3+(3HO0.5)2".

    For full details see help(periodictable.formulas.formula_grammar)

    """
    def __init__(self,value=None,density=None,name=None):
        """
        Initialize formula with a string, a structure or another formula.
        """
        self.density,self.name = None,None
        if value==None:
            self.structure = []
        elif isinstance(value,Formula):
            self.__dict__ = copy(value.__dict__)
        elif _isatom(value):
            self.structure = ((1,value),)
            self.density = value.density
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
                raise TypeError("not a valid chemical formula")
        if density: self.density = density
        if name: self.name = name

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
        Dictionary of {isotope: count} for the formula.
        """
        return _count_atoms(self.structure)
    atoms = property(_atoms,doc=_atoms.__doc__)

    def _mass(self):
        """
        Molar mass of the molecule.
        """
        mass = 0
        for el,count in self.atoms.iteritems():
            mass += el.mass*count
        return mass
    mass = property(_mass,doc=_mass.__doc__)

    def volume(self,packing_fraction=0.68):
        """
        Estimate molecular volume.

        The volume is estimated from the element covalent
        radius and a given packing fraction.  Density is
        approximately chem.mass/chem.volume().

        The default packing fraction is 0.68.

        Common fractions are::

           simple cubic: 0.52
           body-centered cubic: 0.68
           hexabonal close-packed: 0.74
           face-centered cubic: 0.74
           diamond cubic: 0.34

        Volume can be more accurately estimated from the lattice
        parameters for the crystalline form, but that is beyond
        the scope of this package.
        """
        V = 0
        for el,count in self.atoms.items():
            V += el.covalent_radius**3*count
        V *= 4./3*math.pi
        return V*packing_fraction

    def neutron_sld(self, wavelength=1):
        """
        Neutron scattering information for the molecule.
        Returns (sld,absorption,incoherent scattering).
        Returns None if the density is not known.
        """
        if self.density is None: return None,None,None
        from nsf import neutron_sld_from_atoms
        return neutron_sld_from_atoms(self.atoms,density=self.density,
                                      wavelength=wavelength)

    def xray_sld(self, energy=None, wavelength=None):
        """
        X-ray scattering information for the molecule.
        Returns (sld,absorption,incoherent scattering).
        Returns None if the density is not known.
        """
        if self.density is None: return None,None
        from xsf import xray_sld_from_atoms
        return xray_sld_from_atoms(self.atoms,density=self.density,
                                   wavelength=wavelength,energy=energy)

    def __eq__(self,other):
        """
        Return True if two formulas represent the same structure.  Note
        that they may still have different names and densities.
        Note: doesn't check order.
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

def formula_grammar():
    """
    Return a parser for molecular formulas.

    The parser.parseString() method returns a list of
    pairs (count,fragment),  where fragment is an isotope, an
    element or a list of pairs (count,fragment).

    Formulas are parsed from strings using the following grammar::

        number    :: [1-9][0-9]*
        fraction  :: ( | '0' | number) '.' [0-9]*
        count     :: number | fraction | ''
        symbol    :: [A-Z][a-z]*
        isotope   :: '[' number ']' | ''
        element   :: symbol isotope count
        separator :: '+' | ' '
        group     :: count element+ | '(' formula ')' count
        grammar   :: group separator formula | group
    """
    # Recursive
    formula = Forward()

    # Lookup the element in the element table
    symbol = Regex("[A-Z][a-z]*")
    symbol = symbol.setParseAction(lambda s,l,t: elements.symbol(t[0]))

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
    grammar = formula + StringEnd()

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
