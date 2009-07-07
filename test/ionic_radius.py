from periodictable import elements

def test():
    assert elements.Cl.ionic_radius[7] == 0.26
    ion = elements.Cl.ion[7]
    assert ion.ionic_radius[ion.charge] == 0.26

if __name__ == "__main__": test()
