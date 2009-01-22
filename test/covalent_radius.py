from elements import periodic_table

def test():
    #assert periodic_table.Po.covalent_radius == 1.53
    assert periodic_table.Po.covalent_radius == 1.40
    assert periodic_table.Po.covalent_radius_uncertainty == 0.04
    
if __name__ == "__main__": test()
