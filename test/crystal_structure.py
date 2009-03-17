from periodictable import elements
import periodictable.crystal_structure

def test():
    xtal = elements.Hg.crystal_structure
    assert xtal['symmetry'] == 'Rhombohedral'
    assert xtal['a'] == 2.99
    assert xtal['alpha'] == 70.45

if __name__ == "__main__": test()
