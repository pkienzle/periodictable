.. _custom-table:

*************
Custom tables
*************

You can create your own private instance of PeriodicTable and
populate it with values.

Example:

.. doctest::

    >>> import periodictable
    >>> from periodictable import core, mass, density, elements, formula

    >>> mytable = core.PeriodicTable(table="H=1")
    >>> mass.init(mytable)
    >>> density.init(mytable)

    >>> # Reset mass to H=1 rather than C=12
    >>> scale = elements.H[1].mass
    >>> for el in mytable:
    ...    el._mass /= scale
    ...    if hasattr(el,'_density') and el._density != None: 
    ...        el._density /= scale
    ...    for iso in el:
    ...        iso._mass /= scale

    >>> print("%.10f %.10f"%(mytable.H[1].mass, mytable.C[12].mass))
    1.0000000000 11.9068286833
    >>> print("%.10f %.10f"%(periodictable.H[1].mass, periodictable.C[12].mass))
    1.0078250321 12.0000000000
    >>> print("%.10f"%formula('2H[1]', table=mytable).mass)
    2.0000000000

You will need to add individual properties by hand for all additional
desired properties using ``module.init(elements)``.

The table name (*H=1* above) must be unique within the session.  If you
are pickling elements from a custom table, you must create a custom
table of the same name before attempting to restore them. The default
table is just a custom table with the name *public*.

.. Note: If you are using chemical formulas, you will need to
         define your own parser using::

	from periodictable.formulas import formula_grammar
	parser = formula_grammar(table=elements)


Support for custom tables could be made much smoother by delegating
all properties not defined in the custom table back to the base table,
much like is currently done for :class:`Isotopes <periodictable.core.Isotope>`
and :class:`Ions <periodictable.core.Ion>`. That way you only
need to replace the properties of interest rather than defining all
new properties.

The alternative to using custom tables is to replace a dataset in the
base table using e.g., ``custom_mass.init(elements, reload = True)``, where
custom_mass is your own version of the periodic table values.  Be aware,
however, that this will have a global effect, changing the mass used by
all packages, so you are strongly discouraged from doing so.

