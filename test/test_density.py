from periodictable import elements

def test():
    assert elements.Cm.density == 13.51
    assert elements.Cm.density_caveat == "calculated"
    assert elements.Cm.density_units == "g/cm^3"
    assert elements.D.density == elements.H.density * elements.D.mass/elements.H.mass

    #elements.list('symbol','density','density_caveat', format="%3s %10s   %s")

    #import mass
    #mass.init()
    #elements.list('symbol','interatomic_distance')

if __name__ == "__main__": test()
