from periodictable import elements
from periodictable.core import PeriodicTable
from periodictable import (mass, density, xsf, nsf,
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
    assert elements.Cm.crystal_structure == private.cm.crystal_structure
    assert elements.Cm.covalent_radius == private.cm.covalent_radius
    slde,abse = elements.Cm.xray.sld(energy=8)
    sldp,absp = private.Cm.xray.sld(energy=8)
    assert slde == sldp
    assert abse == absp

    # Check that mods don't change public table.
    private.Cm.mass = elements.Cm.mass+1
    assert elements.Cm.mass != private.Cm.mass
    private.Cm.density = elements.Cm.density+1
    assert elements.Cm.density != private.Cm.density
    private.Cm.covalent_radius = elements.Cm.covalent_radius+1
    assert elements.Cm.covalent_radius != private.Cm.covalent_radius
    private.Cm.xray.newfield = 5
    assert not hasattr(elements.Cm.xray,'newfield')

    # Check that the formula parser can use a private table
    formula = periodictable.formulas.formula_grammar(private)
    H2O = formula('H2O')
    for el in H2O.atoms.keys():
        assert id(el) == private[el.number]

    # A more user friendly approach may be to use a formula parsed against
    # the default table, but somehow evaluate it against the private
    # table within a context.  This will likely require a restructured
    # periodic table class for this to be practical.
