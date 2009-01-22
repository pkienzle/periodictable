# -*- coding: latin-1 -*-
"""

Average covalent radii for the elements.

Fields::

    covalent_radius (A)
    covalent_radius_uncertainty (A)
    covalent_radius_units ("angstrom")

Data taken from Cordero 2008.[1]  

From the abstract::

    A new set of covalent atomic radii has been deduced from 
    crystallographic data for most of the elements with atomic 
    numbers up to 96. The proposed radii show a well behaved 
    periodic dependence that allows us to interpolate a few 
    radii for elements for which structural data is lacking, 
    notably the noble gases. The proposed set of radii therefore 
    fills most of the gaps and solves some inconsistencies in 
    currently used covalent radii. The transition metal and 
    lanthanide contractions as well as the differences in covalent 
    atomic radii between low spin and high spin configurations in 
    transition metals are illustrated by the proposed radii set.

Notes::

 #. Values are averages only.  The particular radius can be highly 
    dependent on oxidation state and chemical compound.

 #. The paper lists values for multiple spin states on select
    elements.  We are using sp3 for carbon and low spin for manganese,
    iron and cobalt.
    
 #. Elements with zero or one measurements of covalent radius are
    assigned an uncertainty of 0.00.  These are He, Ne, Pm, At, Rn,
    Fr, Ac, Pa

 #. Elements above 96 are assigned a covalent radius and uncertainty
    of None.
    
 #. Radii are measured from bonds to C, N or O.  The choice of which
    compound was used is element dependent.  Details are available in
    the references.


[1] Beatriz Cordero, Verónica Gómez, Ana E. Platero-Prats, Marc Revés, 
Jorge Echeverría, Eduard Cremades, Flavia Barragán and Santiago Alvarez. 
Covalent radii revisited. Dalton Trans., 2008, 2832-2838
doi:http://dx.doi.org/10.1039%2Fb801115j
"""

from core import periodic_table, Element

def _init():
    if 'covalent_radius' in periodic_table.properties: return
    periodic_table.properties.append('covalent_radius')

    hasattr(periodic_table[0],'covalent_radius') # TODO why is this needed?
    periodic_table[0].covalent_radius = 0.20
    Element.covalent_radius_units = 'angstrom'
    Element.covalent_radius = None
    Element.covalent_radius_uncertainty = None

    
    for line in Cordero.split('\n'):
        fields = line.split()
        if len(fields) == 3: fields += ['0','0']  # Fill in uncertainties
        if fields[0] == '-': continue  # Skip alternate spin states
        Z = int(fields[0])
        el = fields[1]
        r = float(fields[2])
        dr = float(fields[3])*0.01
        n = int(fields[4])
        
        #hasattr(periodic_table[Z],'covalent_radius') # TODO why is this needed
        periodic_table[Z].covalent_radius = r
        periodic_table[Z].covalent_radius_uncertainty = dr

# Table of radii from Cordero.  Note that in cases where there are
# multiple spin states (C,Mn,Fe,Co) only the first spin state is used.
#
#Z  Symbol radius(A) uncertainty number of measurements
Cordero = """\
1    H    0.31    5    129
2    He    0.28        
3    Li    1.28    7    5789
4    Be    0.96    3    310
5    B    0.84    3    1770
6    Csp3    0.76    1    10000
-    Csp2    0.73    2    10000
-    Csp    0.69    1    171
7    N    0.71    1    2200
8    O    0.66    2    10000
9    F    0.57    3    10000
10    Ne    0.58        
11    Na    1.66    9    1629
12    Mg    1.41    7    3234
13    Al    1.21    4    3206
14    Si    1.11    2    10000
15    P    1.07    3    10000
16    S    1.05    3    10000
17    Cl    1.02    4    1987
18    Ar    1.06    10    9
19    K    2.03    12    435
20    Ca    1.76    10    2647
21    Sc    1.70    7    32
22    Ti    1.60    8    231
23    V    1.53    8    389
24    Cr    1.39    5    916
25    Mnl.s.    1.39    5    321
-     Mnh.s.    1.61    8    929
26    Fel.s.    1.32    3    336
-     Feh.s.    1.52    6    1540
27    Col.s.    1.26    3    5733
-     Coh.s.    1.50    7    780
28    Ni    1.24    4    1030
29    Cu    1.32    4    1149
30    Zn    1.22    4    443
31    Ga    1.22    3    1330
32    Ge    1.20    4    1013
33    As    1.19    4    2015
34    Se    1.20    4    1717
35    Br    1.20    3    2140
36    Kr    1.16    4    5
37    Rb    2.20    9    23
38    Sr    1.95    10    1500
39    Y    1.90    7    30
40    Zr    1.75    7    93
41    Nb    1.64    6    18
42    Mo    1.54    5    97
43    Tc    1.47    7    96
44    Ru    1.46    7    1032
45    Rh    1.42    7    458
46    Pd    1.39    6    1892
47    Ag    1.45    5    1728
48    Cd    1.44    9    19
49    In    1.42    5    546
50    Sn    1.39    4    2999
51    Sb    1.39    5    609
52    Te    1.38    4    692
53    I    1.39    3    451
54    Xe    1.40    9    2
55    Cs    2.44    11    24
56    Ba    2.15    11    3076
57    La    2.07    8    190
58    Ce    2.04    9    47
59    Pr    2.03    7    58
60    Nd    2.01    6    96
61    Pm    1.99        
62    Sm    1.98    8    53
63    Eu    1.98    6    167
64    Gd    1.96    6    178
65    Tb    1.94    5    55
66    Dy    1.92    7    59
67    Ho    1.92    7    48
68    Er    1.89    6    66
69    Tm    1.90    10    15
70    Yb    1.87    8    122
71    Lu    1.87    8    61
72    Hf    1.75    10    53
73    Ta    1.70    8    88
74    W    1.62    7    219
75    Re    1.51    7    476
76    Os    1.44    4    99
77    Ir    1.41    6    131
78    Pt    1.36    5    1768
79    Au    1.36    6    114
80    Hg    1.32    5    137
81    Tl    1.45    7    291
82    Pb    1.46    5    112
83    Bi    1.48    4    51
84    Po    1.40    4    4
85    At    1.50        
86    Rn    1.50        
87    Fr    2.60        
88    Ra    2.21    2    3
89    Ac    2.15    0    1
90    Th    2.06    6    11
91    Pa    2.00    0    1
92    U    1.96    7    57
93    Np    1.90    1    22
94    Pu    1.87    1    9
95    Am    1.80    6    11
96    Cm    1.69    3    16\
"""

_init()
