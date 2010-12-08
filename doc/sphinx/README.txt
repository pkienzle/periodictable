This is the top level build direcgtory for the periodictable
documentation.  All of the documentation is written using
sphinx, a python documentation system using restructured text.

See the User Guide section "Contributing Changes" for details
on how to build the documentation.

Manifest:

    Makefile    use "make clean html pdf" to build the docs
    guide       user documentation
    api         auto-generated interface documentation
    plots       programs to generate figures
    index.rst   top level include document
    genmods.py  creates the module index files
    conf.py     sphinx configuration
    rst_prolog  restructured text macros (for units)
    discoverer, 
    shelltable  sample extensions
    _static     images and style sheets
    _templates  page layouts
    _extensions sphinx extensions
    _build      output directory
