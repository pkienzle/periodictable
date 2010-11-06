from periodictable import elements

def test():
    #assert elements.Po.covalent_radius == 1.53
    assert elements.Po.covalent_radius == 1.40
    assert elements.Po.covalent_radius_uncertainty == 0.04
    assert elements.Po.covalent_radius_units == 'angstrom'

if __name__ == "__main__": test()
