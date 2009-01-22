from copy import deepcopy

from elements import Ca,C,O,H,Fe,Ni,Si
from elements import molecule

def test():
    ikaite=molecule()
    # Note: this should be a tuple of tuples
    ikaite.structure = ((1,Ca),(1,C),(3,O), (6,((2,H),(1,O))))

    # Test print
    assert str(ikaite)=="CaCO3(H2O)6"

    # Test constructors
    assert ikaite==molecule( [(1,Ca),(1,C),(3,O),(6,[(2,H),(1,O)])] )
    assert ikaite==molecule( ikaite )
    assert ikaite is not molecule(ikaite)
    assert ikaite.structure is molecule(ikaite).structure

    # Test parsers
    assert molecule("Ca") == molecule([(1,Ca)])
    assert molecule("Ca") == molecule(Ca)
    assert molecule("CaCO3") == molecule([(1,Ca),(1,C),(3,O)])
    assert ikaite==molecule("CaCO3+6H2O")
    assert ikaite==molecule("(CaCO3+6H2O)1")
    assert ikaite==molecule("CaCO3 6H2O")
    assert ikaite==molecule("CaCO3(H2O)6")
    assert ikaite==molecule("(CaCO3(H2O)6)1")
    assert molecule([(0.75,Fe),(0.25,Ni)])==molecule("Fe0.75Ni0.25")

    # Test composition
    assert ikaite==molecule( "CaCO3" ) + 6*molecule( "H2O" )

    # Check atom count
    assert molecule("Fe2O4+3H2O").atoms == {Fe:2,O:7,H:6}

    # Check the mass calculator
    assert molecule('H2O').mass == 2*H.mass+O.mass
    assert molecule("Fe2O4+3H2O").mass == 2*Fe.mass+7*O.mass+6*H.mass
    assert molecule("Fe2O[18]4+3H2O").mass == 2*Fe.mass+4*O[18].mass+3*O.mass+6*H.mass

    # Test isotopes; make sure this is last since it changes ikaite!
    assert ikaite!=molecule("CaCO[18]3+6H2O")
    assert molecule("O[18]").mass == O[18].mass
    
    # Check x-ray and neutron sld
    import density,xsf,nsf
    rho,mu,inc = molecule('Si',Si.density).neutron_sld(wavelength=4.5)
    rhoSi,muSi,incSi = Si.neutron.sld(wavelength=4.5)
    assert abs(rho - rhoSi) < 1e-14
    assert abs(mu - muSi) < 1e-14
    assert abs(inc - incSi) < 1e-14

    rho,mu = molecule('Si',Si.density).xray_sld(wavelength=1.54)
    rhoSi,muSi = Si.xray.sld(wavelength=1.54)
    assert abs(rho - rhoSi) < 1e-14
    assert abs(mu - muSi) < 1e-14
    
    # Check that names work
    permalloy = molecule('Ni8Fe2',8.692,name='permalloy')
    assert str(permalloy)=='permalloy'

    # Check that get/restore state works
    assert deepcopy(permalloy).__dict__ == permalloy.__dict__

    # Check that copy constructor works
    assert molecule(permalloy).__dict__ == permalloy.__dict__
    assert molecule('Si',name='Silicon').__dict__ != molecule('Si').__dict__

if __name__ == "__main__": test()
