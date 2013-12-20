import periodictable

def test():

    assert periodictable.Be[12].mass == 12.026921
    assert periodictable.Be.mass == 9.012182
    assert abs(periodictable.Be[12].ion[2].mass - (12.026921 - 2*periodictable.constants.electron_mass))<1e-12
    assert abs(periodictable.Be.ion[2].mass - (9.012182 - 2*periodictable.constants.electron_mass))<1e-12
    assert periodictable.Pb[206].abundance == 24.1
    assert periodictable.Pb[209].abundance == 0
    assert periodictable.Pb.mass == 207.2
    assert periodictable.n.mass == 1.00866491597

    # Check abundance totals to 0% or 100%
    for el in periodictable.elements:
        abundance=0
        for iso in el:
            if iso.abundance == None:
                print("%s abundance=None"%iso)
            else:
                abundance += iso.abundance
        assert abs(abundance-100) < 1e-4 or abundance==0,\
            "total abundance for %s is %.15g%%"%(el.symbol,abundance)


    # Check average mass corresponds to abundance information
    # Note: should check that this is true within uncertainty, but
    # uncertainties are not being loaded.
    for el in periodictable.elements:
        abundance=0
        mass=0
        for iso in el:
            if iso.abundance == None:
                print("%s abundance=None"%iso)
            else:
                abundance += iso.abundance
                mass += iso.mass*iso.abundance/100.
        assert abundance==0 or abs(mass - el.mass)/el.mass  < 1e-3,\
            "avg mass for %s is %g != %g"%(el.symbol,el.mass,mass)


if __name__ == "__main__": test()
