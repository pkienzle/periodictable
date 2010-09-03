from periodictable import elements

def test():
    A,a,B,b,C,c,D = elements.Fe.magnetic_ff[2].j0
    assert b == 15.9435
    ion = elements.Fe.ion[2]
    assert ion.magnetic_ff[ion.charge].j0[3] == b

if __name__ == "__main__": test()
