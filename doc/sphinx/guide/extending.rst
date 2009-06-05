.. _private-table:

*************
Private table
*************

You can create your own private instance of PeriodicTable and
populate it with values.

Example::

	from periodictable import core, mass, density

    elements = core.PeriodicTable()
    mass.init(elements)
    density.init(elements)

    # Reset mass to H=1 rather than C=12
    scale = elements.H.mass
    for el in elements:
        el.mass /= scale
        el.density /= scale
        for iso in el.isotopes:
            iso.mass /= scale

You will need to add individual properties by hand for all additional
desired properties using module.init(elements).

Note that if you are using chemical formulas, you will need to
define your own parser using::

	from periodictable.formulas import formula_grammar
	parser = formula_grammar(table=elements)

Support for private tables could be made much smoother by delegating
all properties not defined in the private table back to the base table,
much like is currently done for Isotope and Ion.  That way you only
need to replace the properties of interest rather than defining all
new properties.

Instead of using private tables, you can replace a dataset in the
base table using e.g., mymass.init(elements, reload=True), but you
are strongly discouraged from doing so.

.. _extending:

**************
Extended table
**************

The periodic table is extensible.  Third party packages can
add attributes to the table, and they will appear in all of
the elements.  For example::

    import danse.gammaray  # loads gammaray data
    from periodictable import *  # loads symbols for H, He, ...

    for el in periodictable.elements: print el.symbol, el.gammadata

To implement this, you will need the following in gammaray/core.py::

    from periodictable.core import elements, Element

    def init(table, reload=False):
        if 'gammadata' in table.properties and not reload: return
        table.properties.append('gammadata')

        # Set the default, if any
        Element.gamma = None

        # Set the units
        Element.gamma_units = "Ev"

        # Load the data
        for s,data in gamma_table.iteritems():
            el = table.symbol(s)
            el.gammadata = data

    # Define the data
    gamma_table = dict(
        Si="Silicon gamma values",
        O="Oxygen gamma values",
        )

You can use similar tricks for isotope specific data.  For example,
in shelltable/core.py you might have::

    from periodictable.core import elements, Isotope

    def init(table, reload=False):
        if 'shells' in elements.properties and not reload: return
        elements.properties.append('shells')

        # Set the default.  This is required, even if it is only
        # setting it to None.  If the attribute is missing then the
        # isotope data reverts to the element to supply the value,
        # which is almost certainly not what you want.
        Isotope.shells = None

        # Load the data
        for symbol,data in shell_table.iteritems():
            el = table.symbol(symbol)
            for iso,isodata in data.iteritems():
                el[iso].shells = isodata

    # Define the data
    shell_table = dict(
        Fe={56: "56-Fe shell info",
            58: "58-Fe shell info",
            }
        )


Initializing the table
----------------------

Since your data table is in its own package, you will need some way to
load it.

The simplest for is to load it directly when your package is imported.
For the shelltable example above, you would use the following
in __init__.py::

     import periodictable
     from . import core
     core.init(periodictable.elements)

     del periodictable,core  # Clean up package namespace

In order to allow faster startup and smaller runtime, you may wish to
delay loading table attributes until they are referenced.  For example,
in order to load the isotope information from the shelltable package,
you can use delayed_load in __init__.py::


     from periodictable.core import delayed_load, elements

     # Delayed loading of shell info
     def _load_shell():
         '''
         Electron shell information for isotopes.

         T. Student, Tables of Shell Information
         '''
         from . import core
         core.init(elements)
     delayed_load(['shells'],_load_shell)

     del delayed_load, elements, _load_shell # Clean up package namespace

The first argument to delayed_load is the list of all attributes that will
be defined when the module is loaded.  The second argument is the loading
function, whose docstring will appear as the attribute description for
each attribute in the first list.
