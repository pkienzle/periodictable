.. _using:

***********
Basic usage
***********

This section provides various examples on how to get started with Periodic table.
These examples should work on all periodic table supported platforms.

Access particular elements by name:

.. doctest::

    >>> from periodictable import hydrogen
    >>> print "H mass", hydrogen.mass, hydrogen.mass_units
    H mass 1.00794 u

Access particular elements as symbols:

.. doctest::

    >>> from periodictable import H,B,Cu,Ni
    >>> print "B absorption", B.neutron.absorption
    B absorption 767.0
    >>> print "Ni f1/f2 for Cu K-alpha X-rays", Ni.xray.scattering_factors(wavelength=Cu.K_alpha)
    Ni f1/f2 for Cu K-alpha X-rays (25.022929905648375, 0.52493074546535157)

Access isotopes using mass number subscripts:

.. doctest::

    >>> print "58-Ni vs 62-Ni scattering", Ni[58].neutron.coherent, Ni[62].neutron.coherent
    58-Ni vs 62-Ni scattering 26.1 9.5

Access elements indirectly:

.. doctest::

    >>> import periodictable
    >>> print "Cd density", periodictable.Cd.density, periodictable.Cd.density_units
    Cd density 8.65 g/cm^3

Import all elements:

.. doctest::

    >>> from periodictable import *
    >>> print periodictable.H
    H
    >>> print periodictable.H.mass
    1.00794

Deuterium and tritium are special isotopes named D and T
some neutron information is available as 'n':

.. doctest::

    >>> print "D mass",D.mass
    D mass 2.014101778
    >>> print "neutron mass",n.mass
    neutron mass 1.00866491597

Process all the elements:

.. doctest::

    >>> import periodictable
    >>> for el in periodictable.elements: # doctest: +ELLIPSIS, +NORMALIZE_WHITESPACE
    ...     print el.symbol,el.name
    n neutron
    H hydrogen
    He helium
       ...
    Uuh ununhexium

Another example for processing all elements:

.. doctest::

    >>> from periodictable import elements
    >>> for el in elements: # doctest: +ELLIPSIS, +NORMALIZE_WHITESPACE
    ...     print el.symbol,el.number
    n 0
    H 1
    He 2
       ...

Process all the :class:`isotopes <periodictable.core.Isotope>` for an element:

.. doctest::

    >>> for iso in periodictable.H:
    ...     print iso,iso.mass
    1-H 1.0078250321
    D 2.014101778
    T 3.0160492675
    4-H 4.02783
    5-H 5.03954
    6-H 6.04494

You can create a unique handle to an individual ion.  In addition to storing
the ion charge, this can be used to reference the underlying properties of
the element or isotope:

.. doctest::

    >>> Ni58_2 = periodictable.Ni[58].ion[2]
    >>> Ni_2 = periodictable.Ni.ion[2]
    >>> print "charge for Ni2+",Ni_2.charge
    charge for Ni2+ 2
    >>> print "mass for Ni[58] and for natural abundance", Ni58_2.mass, Ni_2.mass
    mass for Ni[58] and for natural abundance 57.9353479 58.6934

The ion specific properties can be accessed from the ion using ion.charge
for the ion index:

.. doctest::

    >>> import pylab
    >>> import periodictable
    >>> Fe_2 = periodictable.Fe.ion[2]
    >>> print Fe_2.magnetic_ff[Fe_2.charge].M_Q([0,0.1,0.2])
    [ 1.          0.99935255  0.99741366]

The following is a plot of the magnetic form factor vs. Q:

    >>> Q = pylab.linspace(0,16,200)
    >>> M = Fe_2.magnetic_ff[Fe_2.charge].j0_Q(Q)
    >>> pylab.xlabel(r'Magnetic Form Factor for Fe') # doctest: +SKIP
    >>> pylab.ylabel(r'$\AA^{-1}$') # doctest: +SKIP
    >>> pylab.title('Ion specific property for Fe') # doctest: +SKIP
    >>> pylab.plot(Q,M) # doctest: +SKIP

.. plot:: plots/magnetic_ff.py

Missing properties generally evaluate to *None*:

.. doctest::

    >>> print "Radon density",periodictable.Rn.density
    Radon density None


Specific defined properties related to elements can be accessed in a table format as shown in following example :

.. doctest::

    >>> elements.list('symbol','K_alpha',format="%s K-alpha = %s") # doctest: +ELLIPSIS, +NORMALIZE_WHITESPACE
    Ti K-alpha = 2.7496
    Cr K-alpha = 2.2909
    Mn K-alpha = 2.1031
       ...
    Ag K-alpha = 0.5608

