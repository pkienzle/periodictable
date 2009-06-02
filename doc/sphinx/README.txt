This is the top level build direcgtory for the periodictable
documentation.  All of the documentation is written using
sphinx, a python documentation system using restructured text.

To build the HTML documentation, install sphinx using easy_install
for example and type:

    make html pdf

Produces:

    _build/html/index.html
    _build/latex/PeriodicTable.pdf

Manifest:

    guide       user documentation
    api         auto-generated interface documentation
    index.rst   top level include document
    genmods.py  creates the module index files
    conf.py     sphinx configuration
    Makefile    make html to
    _static     images and style sheets
    _templates  page layouts
    _build      output directory
