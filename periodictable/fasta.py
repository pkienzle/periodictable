"""
Biomolecule support.

:class:`Molecule` lets you define biomolecules with labile hydrogen atoms
specified using tritium (T) in the chemical formula.  The biomolecule object
creates forms with natural isotope ratio, all hydrogen and all deuterium.
Density can be provided as natural density or cell volume.  A %D2O contrast
match value is computed for matching the molecule SLD in the presence of
labile hydrogens.  :meth:`Molecule.D2Osld` computes the neutron SLD for
the solvated molecule in a %D2O solvent.

:func:`D2Omatch` computes the %D2O constrast match value given the fully
hydrogenated and fully deuterated forms.

:class:`Sequence` lets you read amino acid and DNA/RNA sequences from FASTA
files.

Tables for common molecules are provided[1]:

    *AMINO_ACID_CODES* : amino acids indexed by FASTA code

    *RNA_CODES*, DNA_CODES* : nucleic bases indexed by FASTA code

    *RNA_BASES*, DNA_BASES* : individual nucleic acid bases

    *NUCLEIC_ACID_COMPONENTS*, *LIPIDS*, *CARBOHYDRATE_RESIDUES*

Neutron SLD for water at 20C is also provided as *H2O_SLD* and *D2O_SLD*.

For unmodified protein need to add 2*T and O for terminations.

Assumes that proteins were created in an environment with the usual H/D isotope
ratio on the non-swappable hydrogens.

[1] Perkins, Modern Physical Methods in Biochemistry Part B, 143-265 (1988)

"""
from __future__ import division

import periodictable as pt

class Molecule(object):
    """
    Specify a biomolecule by name, chemical formula, cell volume and charge.

    Labile hydrogen positions should be coded using tritium (T) rather than H.  That
    way the tritium can be changed to H[1] for solutions with pure water, H for solutions
    with a natural abundance of water or D for solutions with pure deuterium.

    **Attributes**

    *formula* is the original tritiated formula.  You can retrieve the hydrogenated or
    deuterated forms using :func:`isotope_substitution` with *formula*, periodictable.T
    and periodictable.H or periodictable.D.

    *D2Omatch* is the % D2O in H2O required to contrast match the molecule, including
    the the proton swapping effect.

    *sld*/*Hsld*/*Dsld* are the the scattering length densities of the molecule with tritium
    replaced by naturally occurring H/D ratios, pure H[1] and pure H[2] respectively.

    *mass*/*Hmass*/*Dmass* are the masses the three conditions.

    *charge* is the charge on the molecule

    *cell_volume* is the estimated cell volume for the molecule

    *density* is the estimated molecule density
    """
    def __init__(self, name, formula, cell_volume=None, density=None, charge=0):
        # Fill in density or cell_volume
        M = pt.formula(formula, natural_density=density)
        if cell_volume is not None:
            M.density = 1e24*M.molecular_mass/cell_volume if cell_volume > 0 else 0
            #print name, M.molecular_mass, cell_volume, M.density
        else:
            cell_volume = 1e24*M.molecular_mass/M.density

        Hnatural = isotope_substitution(M, pt.T, pt.H)
        H = isotope_substitution(M, pt.T, pt.H[1])
        D = isotope_substitution(M, pt.T, pt.D)

        self.name = name
        self.formula = M
        self.cell_volume = cell_volume
        self.sld = pt.neutron_sld(Hnatural, wavelength=5)[0]
        self.Hsld = pt.neutron_sld(H, wavelength=5)[0]
        self.Dsld = pt.neutron_sld(D, wavelength=5)[0]
        self.mass, self.Hmass, self.Dmass = Hnatural.mass, H.mass, D.mass
        self.D2Omatch = D2Omatch(self.Hsld, self.Dsld)
        self.charge = charge
        self.Hnatural = Hnatural

    def D2Osld(self, volume_fraction=1., D2O_fraction=0.):
        """
        Neutron SLD of the molecule in a %D2O solvent.
        """
        solvent_sld = D2O_fraction*D2O_SLD + (1-D2O_fraction)*H2O_SLD
        solute_sld = D2O_fraction*self.Dsld + (1-D2O_fraction)*self.Hsld
        return volume_fraction*solute_sld + (1-volume_fraction)*solvent_sld

class Sequence(Molecule):
    """
    Convert FASTA sequence into chemical formula.

    *name* sequence name

    *sequence* code string

    *type* is one of::

       aa: amino acid sequence
       dna: dna sequence
       rna: rna sequence

    Note: rna sequence files treat T as U and dna sequence files treat U as T.
    """
    @classmethod
    def loadall(self, filename, type=None):
        """
        Iterate over sequences in FASTA file, loading each in turn.

        Yields one FASTA sequence each cycle.
        """
        type = _guess_type_from_filename(filename, type)
        with open(filename, 'rt') as fh:
            for name, seq in read_fasta(fh):
                yield Sequence(name, seq, type=type)

    @classmethod
    def load(self, filename, type=None):
        """
        Load the first FASTA sequence from a file.
        """
        type = _guess_type_from_filename(filename, type)
        with open(filename, 'rt') as fh:
            name, seq = next(read_fasta(fh))
            return Sequence(name, seq, type=type)

    def __init__(self, name, sequence, type='aa'):
        codes = CODE_TABLES[type]
        sequence = sequence.split('*', 1)[0]  # stop at first '*'
        sequence = sequence.replace(' ', '')  # ignore spaces
        parts = tuple(codes[c] for c in sequence)
        cell_volume = sum(p.cell_volume for p in parts)
        charge = sum(p.charge for p in parts)
        structure = []
        for p in parts:
            structure.extend(list(p.formula.structure))
        formula = pt.formula(structure).hill

        Molecule.__init__(self, name, formula,
                          cell_volume=cell_volume, charge=charge)
        self.sequence = sequence

def _guess_type_from_filename(filename, type):
    if type is None:
        if filename.endswith('.fna'):
            type = 'dna'
        elif filename.endswith('.ffn'):
            type = 'dna'
        elif filename.endswith('.faa'):
            type = 'aa'
        elif filename.endswith('.frn'):
            type = 'rna'
        else:
            type = 'aa'
    return type

# Water density at 20C; neutron wavelength doesn't matter (use 5 A).
H2O_SLD = pt.neutron_sld(pt.formula("H2O@0.9982"), wavelength=5)[0]
D2O_SLD = pt.neutron_sld(pt.formula("D2O@0.9982"), wavelength=5)[0]
def D2Omatch(Hsld, Dsld):
    """
    Find the D2O% concentration of solvent such that neutron SLD of the
    material matches the neutron SLD of the solvent.

    *Hsld*, *Dsld* are the SLDs for the hydrogenated and deuterated forms
    of the material respectively, where *D* includes all the labile protons
    swapped for deuterons.  Water SLD is calculated at 20 C.

    Note that the resulting percentage is only meaningful between
    0% to 100%.  Beyond 100% you will need an additional constrast agent
    in the 100% D2O solvent to increase the SLD enough to match.
    """
    # SLD(%Dsample + (1-%)Hsample) = SLD(%D2O + (1-%)H2O)
    # %SLD(Dsample) + (1-%)SLD(Hsample) = %SLD(D2O) + (1-%)SLD(H2O)
    # %(SLD(Dsample) - SLD(Hsample) + SLD(H2O) - SLD(D2O)) = SLD(H2O) - SLD(Hsample)
    # % = 100*(SLD(H2O) - SLD(Hsample)) / (SLD(Dsample) - SLD(Hsample) + SLD(H2O) - SLD(D2O))
    return 100*(H2O_SLD - Hsld) / (Dsld - Hsld + H2O_SLD - D2O_SLD)

def read_fasta(fp):
    """
    Iterate over the sequences in a FASTA file.

    Each iteration is a pair (sequence name, sequence codes).
    """
    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name:
                yield (name, ''.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name:
        yield (name, ''.join(seq))


def isotope_substitution(formula, source, target, portion=1):
    """
    Substitute one atom/isotope in a formula with another in some proportion.

    *formula* is the formula being updated.

    *source* is the isotope/element to be substituted.

    *target* is the replacement isotope/element.

    *portion* is the proportion of source which is substituted for target.
    """
    atoms = formula.atoms
    if source in atoms:
        mass = formula.mass
        mass_reduction = atoms[source]*portion*(source.mass - target.mass)
        density = formula.density * (mass - mass_reduction)/mass
        atoms[target] = atoms.get(target, 0) + atoms[source]*portion
        if portion == 1:
            del atoms[source]
        else:
            atoms[source] *= 1-portion
    else:
        density = formula.density
    return pt.formula(atoms, density=density)

def _code_average(bases, code_table):
    """
    Compute average over possible nucleotides, assuming equal weight if
    precise nucleotide is not known
    """
    n = len(bases)
    formula, cell_volume, charge = pt.formula(), 0, 0
    for c in bases:
        base = code_table[c]
        formula += base.formula
        cell_volume += base.cell_volume
        charge += base.charge
    if n > 0:
        formula, cell_volume, charge = (1/n) * formula, cell_volume/n, charge/n
    return formula, cell_volume, charge

def _set_amino_acid_average(target, codes, name=None):
    formula, cell_volume, charge = _code_average(codes, AMINO_ACID_CODES)
    if name is None:
        name = "/".join(AMINO_ACID_CODES[c].name for c in codes)
    molecule = Molecule(name, formula, cell_volume=cell_volume, charge=charge)
    AMINO_ACID_CODES[target] = molecule

# FASTA code table
def _(code, V, formula, name):
    if formula[-1] == '-':
        charge = -1
        formula = formula[:-1]
    elif formula[-1] == '+':
        charge = +1
        formula = formula[:-1]
    else:
        charge = 0
    molecule = Molecule(name, formula, cell_volume=V, charge=charge)
    molecule.code = code  # Add code attribute so we can write as well as read
    return code, molecule

# pylint: disable=bad-whitespace
AMINO_ACID_CODES = dict((
    #code, volume, formula,        name
    _("A",  91.5, "C3H4TNO",    "alanine"),
    #B: D or N
    _("C", 105.6, "C3H3TNOS",   "cysteine"),
    _("D", 124.5, "C4H3TNO3-",  "aspartic acid"),
    _("E", 155.1, "C5H5TNO3-",  "glutamic acid"),
    _("F", 203.4, "C9H8TNO",    "phenylalanine"),
    _("G",  66.4, "C2H2TNO",    "glycine"),
    _("H", 167.3, "C6H5T3N3O+", "histidine"),
    _("I", 168.8, "C6H10TNO",   "isoleucine"),
    #J: L or I
    _("K", 171.3, "C6H9T4N2O+", "lysine"),
    _("L", 168.8, "C6H10TNO",   "leucine"),
    _("M", 170.8, "C5H8TNOS",   "methionine"),
    _("N", 135.2, "C4H3T3N2O2", "asparagine"),
    #O: _("O", ???.?, "C12H21N3O3", "pyrrolysine") -- update X below
    _("P", 129.3, "C5H7NO",     "proline"),
    _("Q", 161.1, "C5H5T3N2O2", "glutamine"),
    _("R", 202.1, "C6H7T6N4O+", "arginine"),
    _("S",  99.1, "C3H3T2NO2",  "serine"),
    _("T", 122.1, "C4H5T2NO2",  "threonine"),
    #U: selenocysteine -- update X below
    _("V", 141.7, "C5H8TNO",    "valine"),
    _("W", 237.6, "C11H8T2N2O", "tryptophan"),
    #X: any
    _("Y", 203.6, "C9H7T2NO2",  "tyrosine"),
    #Z: E or Q
    #-: gap
    ))
_set_amino_acid_average('B', 'DN')
_set_amino_acid_average('J', 'LI')
_set_amino_acid_average('Z', 'EQ')
_set_amino_acid_average('X', 'ACDEFGHIKLMNPQRSTVWY', name='any')
_set_amino_acid_average('-', '', name='gap')
__doc__ += "\n\n*AMINO_ACID_CODES*::\n\n    " + "\n    ".join("%s: %s"%(k, v.name) for k, v in sorted(AMINO_ACID_CODES.items()))

def _(formula, V, name):
    molecule = Molecule(name, formula, cell_volume=V)
    return name, molecule
NUCLEIC_ACID_COMPONENTS = dict((
    # formula, volume, name
    _("NaPO3",      60, "phosphate"),
    _("C5H6TO3",   125, "ribose"),
    _("C5H7O2",    115, "deoxyribose"),
    _("C5H2T2N5",  114, "adenine"),
    _("C4H2TN2O2",  99, "uracil"),
    _("C5H4TN2O2", 126, "thymine"),
    _("C5HT3N5O",  119, "guanine"),
    _("C4H2T2N3O", 103, "cytosine"),
    ))
__doc__ += "\n\n*NUCLEIC_ACID_COMPONENTS*::\n\n  " + "\n  ".join("%s: %s"%(k, v.formula) for k, v in sorted(NUCLEIC_ACID_COMPONENTS.items()))

CARBOHYDRATE_RESIDUES = dict((
    # formula, volume, name
    _("C6H7T3O5",    171.9, "Glc"),
    _("C6H7T3O5",    166.8, "Gal"),
    _("C6H7T3O5",    170.8, "Man"),
    _("C6H7T4O5",    170.8, "Man (terminal)"),
    _("C8H10T3NO5",  222.0, "GlcNAc"),
    _("C8H10T3NO5",  232.9, "GalNAc"),
    _("C6H7T3O4",    160.8, "Fuc (terminal)"),
    _("C11H11T5NO8", 326.3, "NeuNac (terminal)"),
    # Glycosaminoglycans
    _("C14H15T5NO11Na", 390.7, "hyaluronate"),  # GlcA.GlcNAc
    _("C14H17T5NO13SNa", 473.5, "keratan sulphate"), # Gal.GlcNAc.SO4
    _("C14H15T4NO14SNa", 443.5, "chondroitin sulphate"), # GlcA.GalNAc.SO4
    ))
__doc__ += "\n\n*CARBOHYDRATE_RESIDUES*::\n\n  " + "\n  ".join("%s: %s"%(k, v.formula) for k, v in sorted(CARBOHYDRATE_RESIDUES.items()))

LIPIDS = dict((
    # formula, volume, name
    _("CH2", 27, "methylene"),
    _("CD2", 27, "methylene-D"),
    _("C10H18NO8P", 350, "phospholipid headgroup"),
    _("C6H5O6", 240, "triglyceride headgroup"),
    _("C36H72NO8P", 1089, "DMPC"),
    _("C36H20D52NO8P", 1089, "DMPC-D52"),
    _("C29H55T3NO8P", 932, "DLPE"),
    _("C27H45TO", 636, "cholesteral"),
    _("C45H78O2", 1168, "oleate"),
    _("C57H104O6", 1617, "trioleate form"),
    _("C39H77T2N2O2P", 1166, "palmitate ester"),
    ))
__doc__ += "\n\n*LIPIDS*::\n\n  " + "\n  ".join("%s: %s"%(k, v.formula) for k, v in sorted(LIPIDS.items()))

def _(code, formula, V, name):
    molecule = Molecule(name, formula, cell_volume=V)
    molecule.code = code
    return code, molecule
RNA_BASES = dict((
    # code, formula, volume, name
    _("A",  "C10H8T3N5O6PNa", 299, "adenosine"),
    _("T",   "C9H8T2N2O8PNa", 284, "uridine"), # Use T for U in RNA
    _("G",  "C10H7T4N5O7PNa", 304, "guanosine"),
    _("C",   "C9H8T3N3O7PNa", 288, "cytidine"),
    ))
__doc__ += "\n\n*RNA_BASES*::\n\n  " + "\n  ".join("%s:%s"%(k, v.name) for k, v in sorted(RNA_BASES.items()))

DNA_BASES = dict((
    # code, formula, volume, %D2O matchpoint, name
    _("A",  "C10H9T2N5O5PNa", 289, "adenosine"),
    _("T", "C10H11T1N2O7PNa", 301, "thymidine"),
    _("G",  "C10H8T3N5O6PNa", 294, "guanosine"),
    _("C",   "C9H9T2N3O6PNa", 278, "cytidine"),
    ))
__doc__ += "\n\n*DNA_BASES*::\n\n  " + "\n  ".join("%s:%s"%(k, v.name) for k, v in sorted(DNA_BASES.items()))

def _(code, bases, name):
    D, V, _ = _code_average(bases, RNA_BASES)
    rna = Molecule(name, D.hill, cell_volume=V)
    rna.code = code
    D, V, _ = _code_average(bases, DNA_BASES)
    dna = Molecule(name, D.hill, cell_volume=V)
    rna.code = code
    return (code,rna), (code,dna)
RNA_CODES,DNA_CODES = [dict(v) for v in zip(
    #code, nucleotides,  name
    _("A", "A",    "adenosine"),
    _("C", "C",    "cytidine"),
    _("G", "G",    "guanosine"),
    _("T", "T",    "thymidine"),
    _("U", "T",    "uridine"), # RNA_BASES["T"] is uridine
    _("R", "AG",   "purine"),
    _("Y", "CT",   "pyrimidine"),
    _("K", "GT",   "ketone"),
    _("M", "AC",   "amino"),
    _("S", "CG",   "strong"),
    _("W", "AT",   "weak"),
    _("B", "CGT",  "not A"),
    _("D", "AGT",  "not C"),
    _("H", "ACT",  "not G"),
    _("V", "ACG",  "not T"),
    _("N", "ACGT", "any base"),
    _("X", "",     "masked"),
    _("-", "",     "gap"),
    )]
# pylint: enable=bad-whitespace


CODE_TABLES = {
    'aa': AMINO_ACID_CODES,
    'dna': DNA_CODES,
    'rna': RNA_CODES,
}

def fasta_table():
    rows = []
    rows += [v for k, v in sorted(AMINO_ACID_CODES.items())]
    rows += [v for k, v in sorted(NUCLEIC_ACID_COMPONENTS.items())]
    rows += [Sequence("beta casein", beta_casein)]

    print("%20s %7s %7s %7s %5s %5s %5s %5s %5s %5s"
          % ("name", "M(H2O)", "M(D2O)", "volume",
             "den", "#el", "xray", "nH2O", "nD2O", "%D2O match"))
    for v in rows:
        protons = sum(num*el.number for el, num in v.formula.atoms.items())
        electrons = protons - v.charge
        Xsld = pt.xray_sld(v.formula, wavelength=pt.Cu.K_alpha)
        print("%20s %7.1f %7.1f %7.1f %5.2f %5d %5.2f %5.2f %5.2f %5.1f"%(
            v.name, v.Hmass, v.Dmass, v.cell_volume, v.formula.density,
            electrons, Xsld[0], v.Hsld, v.Dsld, v.D2Omatch))

beta_casein = "RELEELNVPGEIVESLSSSEESITRINKKIEKFQSEEQQQTEDELQDKIHPFAQTQSLVYPFPGPIPNSLPQNIPPLTQTPVVVPPFLQPEVMGVSKVKEAMAPKHKEMPFPKYPVEPFTESQSLTLTDVENLHLPLPLLQSWMHQPHQPLPPTVMFPPQSVLSLSQSKVLPVPQKAVPYPQRDMPIQAFLLYQEPVLGPVRGPFPIIV"

def test():
    from periodictable.constants import avogadro_number
    # Beta casein results checked against Duncan McGillivray's spreadsheet
    # beta casein 23561.9 23880.9 30872.9  1.27 12614 11.55  1.68  2.75
    s = Sequence("beta casein", beta_casein)
    assert abs(s.Dmass-23880.9) < 0.1
    #print "density",s.mass/avogadro_number/s.cell_volume*1e24
    assert abs(s.mass/avogadro_number/s.cell_volume*1e24 - 1.267) < 0.01
    assert abs(s.Dsld-2.75) < 0.01

    # Check that X-ray sld is independent of isotope
    H = isotope_substitution(s.formula, pt.T, pt.H)
    D = isotope_substitution(s.formula, pt.T, pt.D)
    Hsld, Dsld = pt.xray_sld(H, wavelength=1.54), pt.xray_sld(D, wavelength=1.54)
    #print Hsld, Dsld
    assert abs(Hsld[0]-Dsld[0]) < 1e-10

if __name__=="__main__":
    fasta_table()
