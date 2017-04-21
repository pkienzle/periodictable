# This program is public domain
# Author: Paul Kienzle
"""

***Not implemented***

Chemical
    Class following the structure of Template::chembox new
    from wikipedia.

    The basic properties are

        * IUCname - standard name
        * CASid - standard identifier
        * all_names - other names
        * formula - chemical formula string

    A variety of numeric properties may be available such as

        * density - usual density
        * boiling_point - boiling temperature

    Numeric properties are rich structures with a number of fields

        * value: median property value
        * variance: range of values for the property (1-sigma)
        * units: SI units for property (string, eg |g/cm^3|)
        * notes: caveats on value to present to user (e.g., temperature
          or pressure for the measurement)
        * reference: source of data
        * value_as('units'): value in particular units
        * variance_as('units'): variance in particular units

lookup(chemical, properties=['formula', 'density'])
    Return information for the named chemical or formula,
    or *None* if the chemical is not found.  This command will
    search for the chemical in all available sources, including
    a number of proprietary sources if online access is available
    from the current IP address. This search is repeated for
    all desired properties.

wikipedia(chemical)
    Look up a chemical in wikipedia. Returns a list of
    chemicals which match. Wikipedia provides information on
    6000+ chemicals. Though it is not necessarily complete or
    accurate, it is freely available. This package will include a
    copy of the chemical information from the latest wikipedia database
    at the time of release.

    Properties available are density, melting_point,
    boiling_point, chemical_formula, CAS_id, and more.
    The initial version extracts only chemical_formula
    and density.

    wikipedia is an instance of ChemicalDatabase (see below),
    so it has methods for returning sorted lists of selected
    chemicals. The file is opened readonly for normal usage.

personal(chemical)
    Lookup a chemical in the personal database.  Returns a
    list of chemicals which match.

    The personal list is stored as a database in the user
    configuration directory. Chemicals found in other
    sources can be added as they are used, and entirely
    new entries can be created. Databases can be merged,
    allowing sharing between groups and between machines.

crc(chemical) [Not implemented]
    Search CRC online for the chemical information.  This
    returns a set of matching chemicals.

    While not an instance of ChemicalDatabase, as much as
    possible it follows the same interface.

    This resource is only available from institutions which
    have purchased a license to the CRC online database.

<other databases> [Not implemented]
    Interfaces to other databases may be provided depending
    on demand and resources.

register(source, rank=0)
    Add a new source to the chemical lookup function.  The
    source should be a class which follows the
    ChemicalDatabase interface as closely as possible.

    The source should have a properties attribute indicating
    which properties are available. When called with a
    chemical name, the source should return an object
    with attributes for the available properties.

    One of the attributes of ChemicalDatabase is a rank
    from 0-10, with 0 being the lowest. The rank determines
    the order in which the sources are searched.  Within a rank,
    sources are chosen arbitrarily. The search will stop
    when the required information is found.

ChemicalDatabase(filename, mode=['update'|'read'])
    A private list of chemicals associated with each user.
    This class supports a dictionary interface, with multiple
    keys per record and with an additional iterator methods
    for traversing sorted subsets of the database.

    The values stored in the table are of class Chemical.
    This list will be stored in a SQLite database, which
    is accessible from C and other languages.
"""
