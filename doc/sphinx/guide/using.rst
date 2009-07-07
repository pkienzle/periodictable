.. _using:

***********
Basic usage
***********

Access particular elements by name::

    from periodictable import hydrogen
    print "H mass", hydrogen.mass, hydrogen.mass_units

Access particular elements as symbols::

    from periodictable import H,B,Cu,Ni
    print "H mass", H.mass, H.mass_units
    print "B absorption", B.neutron.absorption, B.neutron.absorption_units
    print "Ni f1/f2 for Cu K-alpha X-rays", Ni.xray.scattering_factors(Cu.K_alpha)

Access isotopes using mass number subscripts::

    print "58-Ni vs 62-Ni scattering", \\
        Ni[58].neutron.coherent, Ni[62].neutron.coherent

Access elements indirectly::

    import periodictable
    print "Cd density", periodictable.Cd.density, periodictable.Cd.density_units

Import all elements::

    from periodictable import *

Deuterium and tritium are special isotopes named D and T
some neutron information is available as 'n'::

    print "D mass",D.mass
    print "neutron mass",n.mass

Process all the elements::

    for el in periodictable.elements: # if you used "import periodictable"
        print el.symbol,el.name

    for el in elements: # if you used "from periodictable import *"
        print el.symbol,el.number

Process all the isotopes for an element::

    for iso in periodictable.Fe:
        print iso,iso.mass

Look at ionic radii::

    for ion in periodictable.Cl.ions:
        print "ionic radius",ion,ion.ionic_radius[ion.charge]


Missing properties generally evaluate to None::

    print "Radon density",periodictable.Rn.density

Helper function for listing only the defined properties::

    elements.list('symbol','K_alpha',format="%s K-alpha=%s")

Work with molecules::

    SiO2 = periodictable.formula('SiO2')
    hydrated = SiO2 + periodictable.formula('3H2O')
    print hydrated,'mass',hydrated.mass
    rho,mu,inc = periodictable.neutron_sld('SiO2+3H2O',density=1.5,wavelength=4.75)
    print hydrated,'neutron sld','%.3g'%rho
    rho,mu = periodictable.xray_sld(hydrated,density=1.5,wavelength=Cu.K_alpha)
    print hydrated,'X-ray sld','%.3g'%rho

