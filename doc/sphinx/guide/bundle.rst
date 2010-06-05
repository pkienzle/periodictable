
.. bundling:

********************
Bundling with py2exe
********************

When using periodictable as part of a bundled package, you need to be sure to
include the data associated with the tables.  This can be done by adding a
periodictable entry into the `package_data` property of the distutils setup file::

    import periodictable
    ...
    setup(..., package_data=periodictable.package_data, ...)

If you have a number of packages which add package data (for example, periodic
table extensions), then you can use the following::

    import periodictable

    package_data = {}
    ...
    package_data.update(periodictable.package_data)
    ...
    setup(..., package_data=package_data, ...)
