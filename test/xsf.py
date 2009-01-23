from elements import molecule
from elements import Cu,Mo,Ni,Fe,Si
from elements.xsf import xray_energy, xray_sld_from_atoms

def test():

    # Check some K_alpha and K_beta1 values
    assert Cu.K_alpha == 1.5418
    assert Cu.K_beta1 == 1.3922

    # Check scalar scattering factor lookup
    f1,f2 = Ni.xray.scattering_factors(xray_energy(Cu.K_alpha))
    assert abs(f1-25.0229)<0.0001
    assert abs(f2-0.5249)<0.0001

    # Check array scattering factor lookup
    f1,f2 = Ni.xray.scattering_factors(xray_energy(Cu.K_alpha))
    m1,m2 = Ni.xray.scattering_factors(xray_energy(Mo.K_alpha))
    B1,B2 = Ni.xray.scattering_factors(xray_energy([Cu.K_alpha,Mo.K_alpha]))
    assert (B1==[f1,m1]).all() and (B2==[f2,m2]).all()

    # Check that we can lookup sld by wavelength and energy
    Fe_rho,Fe_mu = Fe.xray.sld(wavelength=Cu.K_alpha)
    assert abs(Fe_rho-59.45) < 0.01
    Si_rho,Si_mu = Si.xray.sld(energy=8.050)
    assert abs(Si_rho-20.0701) < 0.0001
    assert abs(Si_mu-0.4572) < 0.0001

    # Check that wavelength is the default
    Fe_rho_default,Fe_mu_default = Fe.xray.sld(Cu.K_alpha)
    assert Fe_rho == Fe_rho_default and Fe_mu == Fe_mu_default

    # Check array form of sld lookup
    f1,f2 = Si.xray.sld(wavelength=Cu.K_alpha)
    m1,m2 = Si.xray.sld(wavelength=Mo.K_alpha)
    B1,B2 = Si.xray.sld(wavelength=[Cu.K_alpha,Mo.K_alpha])
    assert (B1==[f1,m1]).all() and (B2==[f2,m2]).all()

    #print Cu.xray.sftable
    #plot_xsf(Cu)
    #emission_table()
    #sld_table(table,Cu.K_alpha)

    """
    # Table of scattering length densities for various molecules
    for molecule,density in [('SiO2',2.2),('B4C',2.52)]:
        atoms = molecule(molecule).atoms
        rho,mu = xray_sld(atoms,density,wavelength=Cu.K_alpha)
        print "sld for %s(%g g/cm**3)  rho=%.4g mu=%.4g"\
            %(molecule,density,rho,mu)
    """

    # Cross check against mo
    rho,mu = xray_sld_from_atoms({Si:1},density=Si.density,wavelength=1.54)
    rhoSi,muSi = Si.xray.sld(wavelength=1.54)
    assert abs(rho - rhoSi) < 1e-14
    assert abs(mu - muSi) < 1e-14

    # Check that neutron_sld works as expected
    atoms = molecule('SiO2').atoms
    rho,mu = xray_sld_from_atoms(atoms,2.2,energy=xray_energy(Cu.K_alpha))
    assert abs(rho-18.87)<0.1
    atoms = molecule('B4C').atoms
    rho,mu = xray_sld_from_atoms(atoms,2.52,energy=xray_energy(Cu.K_alpha))
    assert abs(rho-20.17)<0.1


if __name__ == "__main__": test()
