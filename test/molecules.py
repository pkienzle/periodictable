from copy import deepcopy

from periodictable import Ca,C,O,H,Fe,Ni,Si
from periodictable import formula

def test():
    ikaite=formula()
    # Note: this should be a tuple of tuples
    ikaite.structure = ((1,Ca),(1,C),(3,O), (6,((2,H),(1,O))))

    # Test print
    assert str(ikaite)=="CaCO3(H2O)6"

    # Test constructors
    assert ikaite==formula( [(1,Ca),(1,C),(3,O),(6,[(2,H),(1,O)])] )
    assert ikaite==formula( ikaite )
    assert ikaite is not formula(ikaite)
    assert ikaite.structure is formula(ikaite).structure

    # Test parsers
    assert formula("Ca") == formula([(1,Ca)])
    assert formula("Ca") == formula(Ca)
    assert formula("CaCO3") == formula([(1,Ca),(1,C),(3,O)])
    assert ikaite==formula("CaCO3+6H2O")
    assert ikaite==formula("(CaCO3+6H2O)1")
    assert ikaite==formula("CaCO3 6H2O")
    assert ikaite==formula("CaCO3(H2O)6")
    assert ikaite==formula("(CaCO3(H2O)6)1")
    assert formula([(0.75,Fe),(0.25,Ni)])==formula("Fe0.75Ni0.25")

    # Test composition
    assert ikaite==formula( "CaCO3" ) + 6*formula( "H2O" )

    # Check atom count
    assert formula("Fe2O4+3H2O").atoms == {Fe:2,O:7,H:6}

    # Check the mass calculator
    assert formula('H2O').mass == 2*H.mass+O.mass
    assert formula("Fe2O4+3H2O").mass == 2*Fe.mass+7*O.mass+6*H.mass
    assert (formula("Fe2O[18]4+3H2O").mass
            == 2*Fe.mass+4*O[18].mass+3*O.mass+6*H.mass)

    # Test isotopes; make sure this is last since it changes ikaite!
    assert ikaite!=formula("CaCO[18]3+6H2O")
    assert formula("O[18]").mass == O[18].mass

    # Check x-ray and neutron sld
    import density,xsf,nsf
    rho,mu,inc = formula('Si',Si.density).neutron_sld(wavelength=4.5)
    rhoSi,muSi,incSi = Si.neutron.sld(wavelength=4.5)
    assert abs(rho - rhoSi) < 1e-14
    assert abs(mu - muSi) < 1e-14
    assert abs(inc - incSi) < 1e-14

    rho,mu = formula('Si',Si.density).xray_sld(wavelength=1.54)
    rhoSi,muSi = Si.xray.sld(wavelength=1.54)
    assert abs(rho - rhoSi) < 1e-14
    assert abs(mu - muSi) < 1e-14

    # Check that names work
    permalloy = formula('Ni8Fe2',8.692,name='permalloy')
    assert str(permalloy)=='permalloy'

    # Check that get/restore state works
    assert deepcopy(permalloy).__dict__ == permalloy.__dict__

    # Check that copy constructor works
    assert formula(permalloy).__dict__ == permalloy.__dict__
    assert formula('Si',name='Silicon').__dict__ != formula('Si').__dict__

if __name__ == "__main__": test()
