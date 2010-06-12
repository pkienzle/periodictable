import periodictable
from periodictable import elements, formula

def test():
    H,He = elements.H,elements.He
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

    # Check abundance totals to 0% or 100%
    for el in elements:
        if not hasattr(el,'neutron'): continue
        abundance=0
        for iso in el:
            if iso.neutron == None: continue
            if not hasattr(iso.neutron,'abundance'):
                print "abundance missing for",iso
            if iso.neutron.abundance == None:
                print iso,"abundance=None"
            else:
                abundance += iso.neutron.abundance
        # TODO: abundance tables are not very good
        assert abs(abundance-100) < 1.1 or abundance==0,\
            "total abundance for %s is %.15g%%"%(el.symbol,abundance)

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
            #print "%2s %.3f % 7.3f % 7.3f"%(el.symbol,err,b_c,el.neutron.b_c)

    neutron_sld = periodictable.nsf.neutron_sld
    #periodictable.nsf.sld_table(4.75)
    #periodictable.nsf.energy_dependent_table()

    """
    # Specific elements with b_c values different from Neutron News 1992.
    # This is not a complete list.
    for sym in ['Sc','Te','Xe','Sm','Eu','Gd','W','Au','Hg']:
        el = getattr(elements,sym)
        print el.symbol,el.neutron.b_c,el.neutron.coherent,\
            el.neutron.incoherent,el.neutron.absorption
    """

    # Check that neutron_sld works as expected
    atoms = formula('SiO2').atoms
    coh,absorp,inc = neutron_sld(atoms,2.2,wavelength=4.75)
    assert abs(coh-3.475)<0.001
    #assert abs(absorp-0.0001)<0.00001
    atoms = formula('B4C').atoms
    coh,absorp,inc = neutron_sld(atoms,2.52,wavelength=4.75)
    assert abs(coh-7.649)<0.001
    #assert abs(absorp-2.226)<0.001
    Si = elements.Si

    # Make sure molecular calculation corresponds to direct calculation
    atoms = formula('Si').atoms
    coh,absorp,inc = neutron_sld(atoms,Si.density,wavelength=4.75)
    coh2,absorp2,inc2 = Si.neutron.sld(wavelength=4.75)
    assert abs(coh-coh2)<0.001
    assert abs(absorp-absorp2)<0.001

    """
    # Table of scattering length densities for various molecules
    for molecule,density in [('SiO2',2.2),('B4C',2.52)]:
        atoms = formula(molecule).atoms
        rho,mu,inc = neutron_sld(atoms,density,wavelength=4.75)
        print "sld for %s(%g g/cm**3)  rho=%.4g mu=%.4g inc=%.4g"\
            %(molecule,density,rho,mu,inc)
    """

    coh,absorp,inc = periodictable.nsf.neutron_sld('', density=0)
    assert coh == absorp == inc == 0

if __name__ == "__main__": test()
