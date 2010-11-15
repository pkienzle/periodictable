from periodictable import elements, formula, mix_by_weight, mix_by_volume
from periodictable.core import PeriodicTable
from periodictable import (mass, density, xsf, nsf, formulas,
                           crystal_structure, covalent_radius)

def test():
    # Create private table
    private = PeriodicTable("mine")
    mass.init(private)
    density.init(private)
    xsf.init(private)
    nsf.init(private)
    crystal_structure.init(private)
    covalent_radius.init(private)

    # Test that it matches public table
    assert elements.Cm.mass == private.Cm.mass
    assert elements.Cm.density == private.Cm.density
    assert elements.Cm.crystal_structure == private.Cm.crystal_structure
    assert elements.Cm.covalent_radius == private.Cm.covalent_radius
    slde,abse = elements.Cm.xray.sld(energy=8)
    sldp,absp = private.Cm.xray.sld(energy=8)
    assert slde == sldp
    assert abse == absp

    # Check that modifications to private table don't change public table.
    private.Cm._mass = elements.Cm.mass+1
    assert elements.Cm.mass != private.Cm.mass
    private.Cm._density = elements.Cm.density+1
    assert elements.Cm.density != private.Cm.density
    private.Cm.covalent_radius = elements.Cm.covalent_radius+1
    assert elements.Cm.covalent_radius != private.Cm.covalent_radius
    private.Cm.xray.newfield = 5
    assert not hasattr(elements.Cm.xray,'newfield')

    # Check that the formula parser can use a private table
    privateH2O = formula('H2O', table=private)
    publicH2O = formula('H2O', table=elements)
    for el in privateH2O.atoms.keys():
        assert id(el) == id(private[el.number])
    assert privateH2O != publicH2O
    assert str(privateH2O) == str(publicH2O)
    assert privateH2O == publicH2O.change_table(private)
    private.H._mass = 1
    assert formula('H2',table=private).mass == 2
