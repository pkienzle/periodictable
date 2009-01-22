from elements import periodic_table

def test():
    assert periodic_table.Cm.density == 13.51
    assert periodic_table.Cm.density_caveat == "calculated"
    assert periodic_table.Cm.density_units == "g/cm**3"

    #periodic_table.list('symbol','density','density_caveat', format="%3s %10s   %s")
    
    #import mass
    #mass.init()
    #periodic_table.list('symbol','interatomic_distance')

if __name__ == "__main__": test()
