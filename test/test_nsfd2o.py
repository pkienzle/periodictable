#!/usr/bin/env python

import periodictable

def test():
    # Regression test: make sure neutron sld is an isotope property
    # Note: it is important that this test be run separately
    D2O = periodictable.formula('D2O',density=20./18)
    D2O_coh,D2O_abs,D2O_inc = periodictable.neutron_sld(D2O)
    D2O_coh2,D2O_abs2,D2O_inc2 = periodictable.neutron_sld(D2O)
    H2O = periodictable.formula('H2O',density=1)
    H2O_coh,H2O_abs,H2O_inc = periodictable.neutron_sld(H2O)
    assert D2O_coh == D2O_coh2
    assert D2O_coh != H2O_coh

if __name__ == "__main__":
    test()
