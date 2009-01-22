import elements

def test():

    assert elements.table.Be[12].mass == 12.026921
    assert elements.table.Be.mass == 9.012182
    assert elements.table.Pb[206].abundance == 24.1
    assert elements.table.Pb[209].abundance == 0
    assert elements.table.Pb.mass == 207.2
    assert elements.table.n.mass == 1.00866491597

    # Check abundance totals to 0% or 100%
    for el in elements.table:
        abundance=0
        for iso in el:
            if iso.abundance == None:
                print iso,"abundance=None"
            else:
                abundance += iso.abundance
        assert abs(abundance-100) < 1e-4 or abundance==0,\
            "total abundance for %s is %.15g%%"%(el.symbol,abundance)


    # Check average mass corresponds to abundance information
    # Note: should check that this is true within uncertainty, but
    # uncertainties are not being loaded.
    for el in elements.table:
        abundance=0
        mass=0
        for iso in el:
            if iso.abundance == None:
                print iso,"abundance=None"
            else:
                abundance += iso.abundance
                mass += iso.mass*iso.abundance/100.
        assert abundance==0 or abs(mass - el.mass)/el.mass  < 1e-3,\
            "avg mass for %s is %g != %g"%(el.symbol,el.mass,mass)

    
if __name__ == "__main__": test()
