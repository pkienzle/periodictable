# This program is public domain
"""
Operations on molecules.

As of this writing, this includes parsing, computing molar mass and 
computing scattering length density.

Requires that the pyparsing module is installed.
"""
from copy import copy
from pyparsing import Literal, Optional, White, Regex, \
    ZeroOrMore, OneOrMore, Forward, StringEnd
import elements
from elements import periodic_table

class Molecule(object):
    """
    Simple molecule representation.  

    This is designed for calculating molar mass and scattering 
    length density, not for representing bonds or atom positions.  
    We preserve the structure of the formula so that it can 
    be used as a basis for a rich text representation such as 
    matplotlib TeX markup.

    Example initializers::

       string: 
          m = Molecule( "CaCO3+6H2O" )
       structure:
          m = Molecule( [(1,Ca),(2,C),(3,O),(6,[(2,H),(1,O)]] )
       molecular math:
          m = Molecule( "CaCO3" ) + 6*Molecule( "H2O" )
       another molecule (makes a copy):
          m = Molecule( Molecule("CaCO3+6H2O") )
       an atom:
          m = Molecule( Ca )
       nothing:
          m = Molecule()

    Additional information can be provided::

       density (g / cm**3)   material density
       name (string) common name for the molecule

    Operations::
    
        m.atoms returns a dictionary of {isotope: count} for each isotope
          in the molecule
        m.mass (u) returns the molar mass of the molecule

    Molecule strings consist of counts and atoms such as "CaCO3+6H2O".  

    Groups can be separated by '+' or space, so "CaCO3 6H2O" works as well. 
    
    Groups and be defined using parentheses, such as "CaCO3(H2O)6".
    
    Parentheses can nest: "(CaCO3(H2O)6)1"
    
    Isotopes are represented by index, e.g., "CaCO[18]3+6H2O". 
    
    Counts can be integer or decimal, e.g. "CaCO3+(3HO0.5)2".

    For full details see help(elements.molecule.molecule_grammar)

    """
    def __init__(self,value=None,density=None,name=None):
        """
        Initialize molecule with a string, a structure or another molecule.
        """
        self.density,self.name = None,None
        if value==None:
            self.structure = []
        elif isinstance(value,Molecule):
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
                raise TypeError("not a valid molecule structure")
        if density: self.density = density
        if name: self.name = name

    def _parse_structure(self,input):
        """
        Set the molecule to the given structure, checking that it is valid.
        """
        _check_atoms(input)
        self.structure = _immutable(input)

    def _parse_string(self,input):
        """
        Parse the molecule string.
        """
        self.structure = _immutable(parser.parseString(input))

    def _atoms(self):
        """
        Dictionary of {isotope: count} for the molecule.
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

    def neutron_sld(self, wavelength=1):
        """
        Neutron scattering information for the molecule.  
        Returns (sld,absorption,incoherent scattering).
        Returns None if the density is not known.
        """
        if self.density is None: return None,None,None
        from nsf import neutron_sld_from_atoms
        return neutron_sld_from_atoms(self.atoms,self.density,wavelength)

    def xray_sld(self, wavelength):
        """
        X-ray scattering information for the molecule.
        Returns (sld,absorption,incoherent scattering).
        Returns None if the density is not known.
        """
        if self.density is None: return None,None
        from xsf import xray_sld_from_atoms
        return xray_sld_from_atoms(self.atoms,self.density,wavelength)

    def __eq__(self,other):
        """
        Return True if two molecules represent the same structure.  Note
        that they may still have different names and densities.
        Note: doesn't check order.
        """
        if not isinstance(other,Molecule): return False
        return self.structure==other.structure

    def __add__(self,other):
        """
        Join two molecules.
        """
        #print "adding",self,other
        if not isinstance(other,Molecule):
           raise TypeError("expected molecule+molecule")
        ret = Molecule()
        ret.structure = tuple(list(self.structure) + list(other.structure))
        return ret

    def __rmul__(self,other):
        """
        Provide a multiplier for molecule.
        """
        #print "multiplying",self,other
        try:
            other += 0
        except TypeError:
            raise TypeError("n*molecule expects numeric n")
        ret = Molecule()
        if self.structure != []:
            ret.structure = ((other,self.structure),)
        return ret

    def __str__(self):
        return self.name if self.name else _str_atoms(self.structure)
    
    def __repr__(self):
        return "molecule('%s')"%(str(self))
    
    def __getstate__(self):
        """
        Pickle molecule structure using a string representation since
        elements can't be pickled.
        """
        output = copy(self.__dict__)
        output['structure'] = _str_atoms(self.structure)
        return output

    def __setstate__(self,input):
        """
        Unpickle molecule structure from a string representation since
        elements can't be pickled.
        """
        self.__dict__ = input
        self._parse_string(input['structure'])
    

def molecule_grammar():
    """
    Return a parser for molecular formulae.
      
    The parser.parseString() method returns a list of 
    pairs (count,fragment),  where fragment is an isotope, an 
    element or a list of pairs (count,fragment).

    Molecules are parsed from strings using the following grammar::
    
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
    symbol = symbol.setParseAction(lambda s,l,t: periodic_table.symbol(t[0]))

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
    Traverse molecule structure, counting the total number of atoms.
    """
    total = {}
    for count,fragment in seq:
         if isinstance(fragment,list):
             partial = _count_atoms(fragment)
         else:
             partial = {fragment: 1}
         for el,elcount in partial.iteritems():
             if el not in total: total[el] = 0
             total[el] += elcount*count 
    return total

def _check_atoms(seq):
    """
    Traverse molecule structure, checking that the counts are numeric and
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
    Convert molecule structure to string.
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
    return isinstance(val,(elements.Element,elements.Isotope))

def _isisotope(val):
    """Return true if value is an isotope"""
    return isinstance(val,elements.Isotope)

def _is_string_like(val):
    """Returns True if val acts like a string"""
    try: val+''
    except: return False
    return True

parser = molecule_grammar()
parser.setName('Chemical Formula')

def test():
    from copy import deepcopy
    from elements import Ca,C,O,H,Fe,Ni,Si
    import mass, density
    ikaite=Molecule()
    # Note: this should be a tuple of tuples
    ikaite.structure = ((1,Ca),(1,C),(3,O), (6,((2,H),(1,O))))

    # Test print
    assert str(ikaite)=="CaCO3(H2O)6"

    # Test constructors
    assert ikaite==Molecule( [(1,Ca),(1,C),(3,O),(6,[(2,H),(1,O)])] )
    assert ikaite==Molecule( ikaite )
    assert ikaite is not Molecule(ikaite)
    assert ikaite.structure is Molecule(ikaite).structure

    # Test parsers
    assert Molecule("Ca") == Molecule([(1,Ca)])
    assert Molecule("Ca") == Molecule(Ca)
    assert Molecule("CaCO3") == Molecule([(1,Ca),(1,C),(3,O)])
    assert ikaite==Molecule("CaCO3+6H2O")
    assert ikaite==Molecule("(CaCO3+6H2O)1")
    assert ikaite==Molecule("CaCO3 6H2O")
    assert ikaite==Molecule("CaCO3(H2O)6")
    assert ikaite==Molecule("(CaCO3(H2O)6)1")
    assert Molecule([(0.75,Fe),(0.25,Ni)])==Molecule("Fe0.75Ni0.25")

    # Test composition
    assert ikaite==Molecule( "CaCO3" ) + 6*Molecule( "H2O" )
    
    # Check the mass calculator
    assert Molecule('H2O').mass == 2*H.mass+O.mass

    # Test isotopes; make sure this is last since it changes ikaite!
    assert ikaite!=Molecule("CaCO[18]3+6H2O")
    assert Molecule("O[18]").mass == O[18].mass
    
    # Check x-ray and neutron sld
    import density,xsf,nsf
    rho,mu,inc = Molecule('Si',Si.density).neutron_sld()
    rhoSi,muSi,incSi = Si.neutron.sld()
    assert abs(rho - rhoSi) < 1e-14
    assert abs(mu - muSi) < 1e-14
    assert abs(inc - incSi) < 1e-14
    
    # Check that names work
    permalloy = Molecule('Ni8Fe2',8.692,name='permalloy')
    assert str(permalloy)=='permalloy'

    # Check that get/restore state works
    assert deepcopy(permalloy).__dict__ == permalloy.__dict__

    # Check that copy constructor works
    assert Molecule(permalloy).__dict__ == permalloy.__dict__
    assert Molecule('Si',name='Silicon').__dict__ != Molecule('Si').__dict__

if __name__ == "__main__": test()
