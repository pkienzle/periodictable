from elements import periodic_table
import elements.crystal_structure

def test():
    xtal = periodic_table.Hg.crystal_structure
    assert xtal['symmetry'] == 'Rhombohedral'
    assert xtal['a'] == 2.99
    assert xtal['alpha'] == 70.45

if __name__ == "__main__": test()
