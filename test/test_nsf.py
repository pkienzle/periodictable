import numpy
import periodictable
from periodictable import elements, formula, nsf
from periodictable.nsf import neutron_scattering, neutron_sld
from math import sqrt, pi

def test():
    H,He,D,O = elements.H,elements.He,elements.D,elements.O
    assert H.neutron.absorption == 0.3326
    assert H.neutron.total == 82.02
    assert H.neutron.incoherent == 80.26
    assert H.neutron.coherent == 1.7568
    assert elements.Ru[101].neutron.bp == None
    assert H[1].nuclear_spin == '1/2'
    assert H[2].nuclear_spin == '1'
    assert not H[6].neutron.has_sld()

    assert He[3].neutron.b_c_i == -1.48
    assert He[3].neutron.bm_i == -5.925

    Nb = elements.Nb
    assert Nb.neutron.absorption == Nb[93].neutron.absorption

    # Check that b_c values match abundance-weighted values
    # Note: Currently they do not for match within 5% for Ar,V,Sm or Gd
    for el in elements:
        if not hasattr(el,'neutron'): continue
        b_c = 0
        complete = True
        for iso in el:
            if iso.neutron != None:
                if iso.neutron.b_c == None:
                    complete = False
                else:
                    b_c += iso.neutron.b_c*iso.neutron.abundance/100.
        if complete and b_c != 0 and abs((b_c-el.neutron.b_c)/b_c) > 0.05:
            err = abs((b_c-el.neutron.b_c)/b_c)
            ## Printing suppressed for the release version
            #print("%2s %.3f % 7.3f % 7.3f"%(el.symbol,err,b_c,el.neutron.b_c))

    # Check neutron_sld and neutron_xs against NIST calculator
    # Note that we are using different tables, so a general comparison with
    # NIST numbers is not possible, but ^30Si and ^18O are the same in both.
    M = formula('Si[30]O[18]2',density=2.2)
    sld,xs,depth = neutron_scattering(M,wavelength=4.75)
    sld2 = neutron_sld(M,wavelength=4.75)
    assert all(abs(v-w)<1e-10 for v,w in zip(sld,sld2))
    #_summarize(M)
    #_summarize(formula('O2',density=1.14))
    # Alan's numbers:
    assert abs(sld[0] - 3.27) < 0.01
    assert abs(sld[1] - 0) < 0.01
    #assert abs(xs[2] - 0.00292) < 0.00001   # TODO fix test
    assert abs(xs[1] - 0.00569) < 0.00001
    #assert abs(depth - 4.329) < 0.001       # TODO fix test

    # Cu/Mo K-alpha = 1.89e-5 + 2.45e-7i / 1.87e-5 + 5.16e-8i

    Ni,Si = elements.Ni, elements.Si

    # Make sure molecular calculation corresponds to direct calculation
    sld = neutron_sld('Si',density=Si.density,wavelength=4.75)
    sld2 = Si.neutron.sld(wavelength=4.75)
    assert all(abs(v-w)<1e-10 for v,w in zip(sld,sld2))

    sld,_,_ = Si.neutron.scattering(wavelength=4.75)
    sld2 = Si.neutron.sld(wavelength=4.75)
    assert all(abs(v-w)<1e-10 for v,w in zip(sld,sld2))

    sld,xs,depth = neutron_scattering('Si',density=Si.density,wavelength=4.75)
    sld2,xs2,depth2 = Si.neutron.scattering(wavelength=4.75)
    assert all(abs(v-w)<1e-10 for v,w in zip(sld,sld2))
    assert all(abs(v-w)<1e-10 for v,w in zip(xs,xs2))
    assert abs(depth-depth2) < 1e-14

    # incoherent cross sections for Ni[62] used to be negative
    sld,xs,depth = neutron_scattering('Ni[62]',density=Ni[62].density,
                                      wavelength=4.75)
    assert sld[2] == 0 and xs[2] == 0
    sld,xs,depth = Ni[62].neutron.scattering(wavelength=4.75)
    assert sld[2] == 0 and xs[2] == 0
    assert Ni[62].neutron.sld()[2] == 0


    # Test call from periodictable
    sld,xs,depth = periodictable.neutron_scattering('H2O',density=1,wavelength=4.75)
    sld2,xs2,depth2 = neutron_scattering('H2O',density=1,wavelength=4.75)
    assert all(abs(v-w)<1e-10 for v,w in zip(sld,sld2))
    assert all(abs(v-w)<1e-10 for v,w in zip(xs,xs2))
    assert depth==depth2
    sld = periodictable.neutron_sld('H2O',density=1,wavelength=4.75)
    assert all(abs(v-w)<1e-10 for v,w in zip(sld,sld2))

    # Check empty formula
    sld,xs,depth = neutron_scattering('',density=0,wavelength=4.75)
    assert all(v == 0 for v in sld)
    assert all(v == 0 for v in xs)
    assert numpy.isinf(depth)

    # Check density == 0 works
    sld,xs,depth = neutron_scattering('Si',density=0,wavelength=4.75)
    assert all(v == 0 for v in sld)
    assert all(v == 0 for v in xs)
    assert numpy.isinf(depth)

    # Test natural density
    D2O_density = (2*D.mass + O.mass)/(2*H.mass + O.mass)
    sld,xs,depth = neutron_scattering('D2O',natural_density=1,wavelength=4.75)
    sld2,xs2,depth2 = neutron_scattering('D2O',density=D2O_density,wavelength=4.75)
    assert all(abs(v-w)<1e-14 for v,w in zip(sld,sld2))
    assert all(abs(v-w)<1e-14 for v,w in zip(xs,xs2))
    assert abs(depth-depth2)<1e-14

    # Test that sld depends on density not on the size of the unit cell
    D2O_density = (2*D.mass + O.mass)/(2*H.mass + O.mass)
    sld,xs,depth = neutron_scattering('D2O',natural_density=1,wavelength=4.75)
    sld2,xs2,depth2 = neutron_scattering('2D2O',natural_density=1,wavelength=4.75)
    assert all(abs(v-w)<1e-14 for v,w in zip(sld,sld2))
    assert all(abs(v-w)<1e-14 for v,w in zip(xs,xs2))
    assert abs(depth-depth2)<1e-14

    # Test energy <=> velocity <=> wavelength
    assert abs(nsf.neutron_wavelength_from_velocity(2200) - 1.7981972618436388) < 1e-14
    assert abs(nsf.neutron_wavelength(25) - 1.8) < 0.1
    assert abs(nsf.neutron_energy(nsf.neutron_wavelength(25)) - 25) < 1e-14

    # Confirm scattering functions accept energy and wavelength
    sld,xs,depth = neutron_scattering('H2O',density=1,wavelength=4.75)
    sld2,xs2,depth2 = neutron_scattering('H2O',density=1,energy=nsf.neutron_energy(4.75))
    assert all(abs(v-w)<1e-14 for v,w in zip(sld,sld2))
    assert all(abs(v-w)<1e-14 for v,w in zip(xs,xs2))
    assert abs(depth-depth2)<1e-14

def test_formula():
    M = formula('B4C', density=2.52)
    sld,xs,depth = neutron_scattering(M,wavelength=4.75)
    # Compare to Alan Munter's numbers:
    #   SLD=7.65e-6 - 2.34e-7i /A^2
    #   inc,abs XS = 0.193, 222.4 / cm
    #   1/e = 0.004483 cm
    #   Cu/Mo K-alpha = 2.02e-5 + 1.93e-8i / 2.01e-5 + 4.64e-9i
    # Using lambda=1.798 rather than 1.8
    #   abs XS => 222.6
    #   1/e => 0.004478
    assert abs(sld[0]-7.649)<0.001
    assert abs(sld[1]-0.234)<0.001
    assert abs(xs[1]-222.6)<0.1
    #assert abs(xs[2]-0.193)<0.001   # TODO: fix test
    #assert abs(depth-0.004478)<0.000001 # TODO: fix test
    # Check that sld_inc and coh_xs are consistent
    #   cell_volume = (molar_mass/density) / N_A * 1e24
    #   number_density = num_atoms/cell_volume
    #   sigma_i = inc_xs/number_density
    #   sld_inc = number_density * sqrt ( 100/(4*pi) * sigma_i ) * 10
    #   sld_re = number_density * b_c * 10
    #   sigma_c = 4*pi/100*b_c**2
    #   coh_xs = sigma_c * number_density
    Nb = 0.13732585020640778
    sld_inc = Nb*sqrt(100/(4*pi)*xs[2]/Nb)*10
    coh_xs = Nb*4*pi/100*(sld[0]/(10*Nb))**2
    assert abs(sld[2] - sld_inc) < 1e-14
    assert abs(xs[0] - coh_xs) < 1e-14

def test_composite():
    from periodictable.nsf import neutron_composite_sld
    calc = neutron_composite_sld([formula(s) for s in ('HSO4','H2O','CCl4')],
                                 wavelength=4.75)
    sld = calc(numpy.array([3,1,2]),density=1.2)
    sld2 = neutron_sld('3HSO4+1H2O+2CCl4',density=1.2,wavelength=4.75)
    #print(sld)
    #print(sld2)
    assert all(abs(v-w)<1e-14 for v,w in zip(sld,sld2))

def time_composite():
    from periodictable.nsf import neutron_composite_sld
    import time
    calc = neutron_composite_sld([formula(s) for s in ('HSO4','H2O','CCl4')],
                                 wavelength=4.75)
    q = numpy.array([3,1,2])
    N = 1000
    bits = [formula(s) for s in ('HSO4','H2O','CCl4')]
    tic = time.time()
    for i in range(N):
        sld = calc(q,density=1.2)
    toc = time.time()-tic
    print("composite %.1f us"%(toc/N*1e6))
    tic = time.time()
    for i in range(N):
        sld = neutron_sld(q[0]*bits[0]+q[1]*bits[1]+q[2]*bits[2],
                          density=1.2,wavelength=4.75)
    toc = time.time()-tic
    print("direct %.1f us"%(toc/N*1e6))


def test_abundance():
    # Check abundance totals to 0% or 100%
    for el in elements:
        if not hasattr(el,'neutron'): continue
        abundance=0
        for iso in el:
            if iso.neutron == None: continue
            if not hasattr(iso.neutron,'abundance'):
                print("abundance missing for %s"%iso)
            if iso.neutron.abundance == None:
                print("%s abundance=None"%iso)
            else:
                abundance += iso.neutron.abundance
        # TODO: abundance tables are not very good
        assert abs(abundance-100) < 1.1 or abundance==0,\
            "total abundance for %s is %.15g%%"%(el.symbol,abundance)



def _summarize(M):
    from periodictable.nsf import neutron_sld, neutron_xs
    sld = neutron_sld(M,wavelength=4.75)
    xs = neutron_xs(M,wavelength=4.75)
    print("%s sld %s"%(M,sld))
    print("%s xs %s 1/e %s"%(M,xs,1/sum(xs)))
    #return
    for el in list(M.atoms.keys()):
        print("%s density %s"%(el,el.density))
        print("%s sld %s"%(el,el.neutron.sld(wavelength=4.75)))
        print("%s xs"%el + " %.15g %.15g %.15g"%el.neutron.xs(wavelength=4.75))
        print("%s 1/e %s"%(el,1./sum(el.neutron.xs(wavelength=4.75))))

def molecule_table():
    # Table of scattering length densities for various molecules
    print("SLDS for some molecules")
    for molecule,density in [('SiO2',2.2),('B4C',2.52)]:
        atoms = formula(molecule).atoms
        rho,mu,inc = neutron_sld(atoms,density,wavelength=4.75)
        print("%s(%g g/cm**3)  rho=%.4g mu=%.4g inc=%.4g"
              %(molecule,density,rho,mu,inc))

def show_tables():
    molecule_table()
    periodictable.nsf.sld_table(4.75)
    periodictable.nsf.energy_dependent_table()
    periodictable.nsf.total_comparison_table()
    periodictable.nsf.coherent_comparison_table()
    periodictable.nsf.incoherent_comparison_table()

    print("""\
Specific elements with b_c values different from Neutron News 1992.
This is not a complete list.""")
    for sym in ['Sc','Te','Xe','Sm','Eu','Gd','W','Au','Hg']:
        el = getattr(elements,sym)
        print("%s %s %s %s %s"%(el.symbol,el.neutron.b_c,el.neutron.coherent,
              el.neutron.incoherent,el.neutron.absorption))



if __name__ == "__main__":
    time_composite()
