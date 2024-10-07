# This program is public domain
# Author: Paul Kienzle
"""
Various fundamental constants.
"""

#: Avogadro's number (mol^-1)
avogadro_number = 6.02214179e23 #(30) mol^-1
#: Planck's constant (eV s)
plancks_constant = 4.13566733e-15 #(10) eV s
#: Electron volt (J/eV)
electron_volt = 1.602176487e-19 #(40) J / eV
#: speed of light c (m/s)
speed_of_light = 299792458 # m/s (exact)
#: electron radius r_e (m)
electron_radius = 2.8179402894e-15 #(58) m

# [CODATA] From NIST Reference on Constants, Units, and Uncertainty
#   http://physics.nist.gov/cuu/index.html
# [AME2020] "The Ame2020 atomic mass evaluation (II)"
#     by M.Wang, W.J.Huang, F.G.Kondev, G.Audi and S.Naimi
#     Chinese Physics C45, 030003, March 2021.
#: neutron mass (u)
neutron_mass = 1.00866491590 #(47) u [AME 2020]
neutron_mass_unc = 0.00000000047
#neutron_mass = 1.00866491597 #(43) u [CODATA 2010?]
#neutron_mass = 1.00866491595 #(49) u [CODATA 2018]
#: atomic mass constant (kg / u)
atomic_mass_constant = 1.660538782e-27 #(83) kg / u
#: electron mass (u)
electron_mass = 5.48577990946e-4 #(22) u
