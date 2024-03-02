# -*- coding: iso-8859-15 -*-
# This program is public domain
# Author: Paul Kienzle
# Based on spreadsheet by Les Slaback (1998).
r"""
Calculate expected neutron activation from time spent in beam line.

Notation information for activation product:

 m, m1, m2: indicate metastable states.  Decay may be to the ground state or to
 another nuclide.

 \+: indicates radioactive daughter production already included in daughter listing
 several parent t1/2's required to acheive calculated daughter activity.  All
 activity assigned at end of irradiation.  In most cases the added activity to the
 daughter is small.

 \*: indicates radioactive daughter production NOT calculated, approximately
 secular equilibrium

 s: indicates radioactive daughter of this nuclide in secular equilibrium after several
 daughter t1/2's

 t: indicates transient equilibrium via beta decay.  Accumulation of that nuclide
 during irradiation is separately calculated.

Reaction = b indicates production via decay from an activation produced parent

Accounts for burnup and 2n, g production.

This is only a gross estimate.  Many effects are not taken into account, such as
self-shielding in the sample and secondary activation from the decay products.

Example::

    >>> from periodictable import activation
    >>> env = activation.ActivationEnvironment(fluence=1e5, Cd_ratio=70, fast_ratio=50, location="BT-2")
    >>> sample = activation.Sample("Co30Fe70", 10)
    >>> sample.calculate_activation(env, exposure=10, rest_times=[0, 1, 24, 360])
    >>> sample.show_table()
                                          ----------------- activity (uCi) ------------------
    isotope  product  reaction  half-life        0 hrs        1 hrs       24 hrs      360 hrs
    -------- -------- -------- ---------- ------------ ------------ ------------ ------------
    Co-59    Co-60         act    5.272 y     0.000496     0.000496    0.0004958    0.0004933
    Co-59    Co-60m+       act     10.5 m        1.664       0.0317          ---          ---
    -------- -------- -------- ---------- ------------ ------------ ------------ ------------
                                    total        1.665      0.03221    0.0005084     0.000505
    -------- -------- -------- ---------- ------------ ------------ ------------ ------------

    >>> print("%.3f"%sample.decay_time(0.001)) # number of hours to reach 1 nCi
    2.053

The default rest times used above show the sample activity at the end of neutron
activation and after 1 hour, 1 day, and 15 days.

The neutron activation table, *activation.dat*,\ [#Shleien1998]_ contains
details about the individual isotopes, with interaction cross sections taken
from from IAEA-273\ [#IAEA1987]_.

Activation can be run from the command line using::

    $ python -m periodictable.activation FORMULA

where FORMULA is the chemical formula for the material.

.. [#Shleien1998] Shleien, B., Slaback, L.A., Birky, B.K., 1998.
   Handbook of health physics and radiological health.
   Williams & Wilkins, Baltimore.

.. [#IAEA1987] IAEA (1987)
   Handbook on Nuclear Activation Data.
   TR 273 (International Atomic Energy Agency, Vienna, Austria, 1987).
   http://cds.cern.ch/record/111089/files/IAEA-TR-273.pdf
"""

from __future__ import division, print_function

from math import exp, log, expm1
import os

from .formulas import formula as build_formula
from . import core

LN2 = log(2)

def NIST2001_isotopic_abundance(iso):
    """
    Isotopic abundance in % from the periodic table package.

    Böhlke, et al.
    Isotopic Compositions of the Elements, 2001.
    J. Phys. Chem. Ref. Data, Vol. 34, No. 1, 2005
    """
    return iso.abundance

def IAEA1987_isotopic_abundance(iso):
    """
    Isotopic abundance in % from the IAEA, as provided in the activation.dat table.

    Note: this will return an abundance of 0 if there is no neutron activation for
    the isotope even though for isotopes such as H[1], the natural abundance may in
    fact be rather large.

    IAEA 273: Handbook on Nuclear Activation Data, 1987.
    """
    try:
        return iso.neutron_activation[0].abundance
    except AttributeError:
        return 0

class Sample(object):
    """
    Sample properties.

    *formula* : chemical formula

        Chemical formula.  Any format accepted by
        :func:`.formulas.formula` can be used, including
        formula string.

    *mass* : float | g

        Sample mass.

    *name* : string

        Name of the sample (defaults to formula).
    """
    def __init__(self, formula, mass, name=None):
        self.formula = build_formula(formula)
        self.mass = mass               # cell F19
        self.name = name if name else str(self.formula) # cell F20
        self.activity = {}

        # The following are set in calculation_activation
        self.environment = None  # type: "ActivationEnvironment"
        self.exposure = 0.
        self.rest_times = ()

    def calculate_activation(self, environment, exposure=1,
                             rest_times=(0, 1, 24, 360),
                             abundance=NIST2001_isotopic_abundance):
        """
        Calculate sample activation after exposure to a neutron flux.

        *environment* is the exposure environment.

        *exposure* is the exposure time in hours (default is 1 h).

        *rest_times* is the list of deactivation times in hours (default is [0, 1, 24, 360]).

        *abundance* is a function that returns the relative abundance of an isotope.  By
        default it uses :func:`NIST2001_isotopic_abundance`, and there is the alternative
        :func:`IAEA1987_isotopic_abundance`.
        """
        self.activity = {}
        self.environment = environment
        self.exposure = exposure
        self.rest_times = rest_times
        for el, frac in self.formula.mass_fraction.items():
            if core.isisotope(el):
                A = activity(el, self.mass*frac, environment, exposure, rest_times)
                self._accumulate(A)
            else:
                for iso in el.isotopes:
                    iso_mass = self.mass*frac*abundance(el[iso])*0.01
                    if iso_mass:
                        A = activity(el[iso], iso_mass, environment, exposure, rest_times)
                        self._accumulate(A)

    def decay_time(self, target):
        """
        After determining the activation, compute the number of hours required to achieve
        a total activation level after decay.
        """
        if not self.rest_times or not self.activity:
            return 0

        # Find the small rest time (probably 0 hr)
        min_rest, To = min(enumerate(self.rest_times), key=lambda x: x[1])
        # Find the activity at that time, and the decay rate
        data = [
            (Ia[min_rest], LN2/a.Thalf_hrs)
            for a, Ia in self.activity.items()
            # TODO: not sure why Ia is zero, but it messes up the initial value guess if it is there
            if Ia[min_rest] > 0.0
            ]
        # Build functions for total activity at time T - target and its derivative
        # This will be zero when activity is at target
        f = lambda t: sum(Ia*exp(-La*(t-To)) for Ia, La in data) - target
        df = lambda t: sum(-La*Ia*exp(-La*(t-To)) for Ia, La in data)
        # Return target time, or 0 if target time is negative
        if f(0) < target:
            return 0
        # Need an initial guess near the answer otherwise find_root gets confused.
        # Small but significant activation with an extremely long half-life will
        # dominate at long times, but at short times they will not affect the
        # derivative. Choosing a time that satisfies the longest half-life seems
        # to work well enough.
        #print("data", data, [])
        initial = max(-log(target/Ia)/La + To for Ia, La in data)
        t, ft = find_root(initial, f, df)
        percent_error = 100*abs(ft)/target
        if percent_error > 0.1:
            #return 1e100*365*24 # Return 1e100 rather than raising an error
            msg = (
                "Failed to compute decay time correctly (%.1g error). Please"
                " report material, mass, flux and exposure.") % percent_error
            raise RuntimeError(msg)
        return t

    def _accumulate(self, activity):
        for el, activity_el in activity.items():
            el_total = self.activity.get(el, [0]*len(self.rest_times))
            self.activity[el] = [T+v for T, v in zip(el_total, activity_el)]

    def show_table(self, cutoff=0.0001, format="%.6e"):
        """
        Tabulate the daughter products.

        *cutoff=1* : float | uCi

              The minimum activation value to show.

        *format="%.1f"* : string

              The number format to use for the activation.
        """
        # TODO: need format="auto" which picks an appropriate precision based on
        # cutoff and/or activation level.

        # Track individual rows with more than 1 uCi of activation, and total activation
        # Replace any activation below the cutoff with '---'
        rows = []
        total = [0]*len(self.rest_times)
        for el, activity_el in sorted_activity(self.activity.items()):
            total = [t+a for t, a in zip(total, activity_el)]
            if all(a < cutoff for a in activity_el):
                continue
            activity_str = [format%a if a >= cutoff else "---" for a in activity_el]
            rows.append([el.isotope, el.daughter, el.reaction, el.Thalf_str]+activity_str)
        footer = ["", "", "", "total"] + [format%t if t >= cutoff else "---" for t in total]

        # If no significant total activation then don't print the table
        if all(t < cutoff for t in total):
            print("No significant activation")
            return

        # Print the table header, with an overbar covering the various rest times
        # Print a dashed separator above and below each column
        header = ["isotope", "product", "reaction", "half-life"] \
                 + ["%g hrs"%vi for vi in self.rest_times]
        separator = ["-"*8, "-"*9, "-"*8, "-"*10] + ["-"*12]*len(self.rest_times)
        cformat = "%-8s %-9s %8s %10s " + " ".join(["%12s"]*len(self.rest_times))

        width = sum(len(c)+1 for c in separator[4:]) - 1
        if width < 16:
            width = 16
        overbar = "-"*(width//2-8) + " activity (uCi) " + "-"*((width+1)//2-8)
        offset = sum(len(c)+1 for c in separator[:4]) - 1
        print(" "*(offset+1)+overbar)
        print(cformat%tuple(header))
        print(cformat%tuple(separator))

        # Print the significant table rows, or indicate that there were no
        # significant rows if the total is significant but none of the
        # individual isotopes
        if rows:
            for r in rows:
                print(cformat%tuple(r))
        else:
            print("No significant isotope activation")
        print(cformat%tuple(separator))

        # If there is more than one row, or if there is enough marginally
        # significant activation that the total is greater then the one row
        # print the total in the footer
        if len(rows) != 1 or any(c != t for c, t in zip(rows[0][4:], footer[4:])):
            print(cformat%tuple(footer))
            print(cformat%tuple(separator))

def find_root(x, f, df, max=20, tol=1e-10):
    r"""
    Find zero of a function.

    Returns when $|f(x)| < tol$ or when max iterations have been reached,
    so check that $|f(x)|$ is small enough for your purposes.

    Returns x, f(x).
    """
    fx = f(x)
    for _ in range(max):
        #print(f"step {_}: {x=} {fx=} f'={df(x)}")
        if abs(f(x)) < tol:
            break
        x -= fx / df(x)
        fx = f(x)
    return x, fx


def sorted_activity(activity_pair):
    """Interator over activity pairs sorted by isotope then daughter product."""
    return sorted(activity_pair, key=lambda x: (x[0].isotope, x[0].daughter))


class ActivationEnvironment(object):
    """
    Neutron activation environment.

    The activation environment provides details of the neutron flux at the
    sample position.

    *fluence* : float | n/cm^2/s

        Thermal neutron fluence on sample.  For COLD neutrons enter equivalent
        thermal neutron fluence.

    *Cd_ratio* : float

        Neutron cadmium ratio.  Use 0 to suppress epithermal contribution.

    *fast_ratio* : float

        Thermal/fast ratio needed for fast reactions. Use 0 to suppress fast contribution.
    """
    def __init__(self, fluence=1e5, Cd_ratio=0., fast_ratio=0., location=""):
        self.fluence = fluence     # cell F13
        self.Cd_ratio = Cd_ratio   # cell F15
        self.fast_ratio = fast_ratio # cell F17
        self.location = location   # cell F21

    # Cell Q1
    @property
    def epithermal_reduction_factor(self):
        """
        Used as a multiplier times the resonance cross section to add to the
        thermal cross section for all thermal induced reactions.
        """
        return 1./self.Cd_ratio if self.Cd_ratio >= 1 else 0

COLUMN_NAMES = [
    "_symbol",     # 0 AF
    "_index",      # 1 AG
    "Z",           # 2 AH
    "symbol",      # 3 AI
    "A",           # 4 AJ
    "isotope",     # 5 AK
    "abundance",   # 6 AL
    "daughter",    # 7 AM
    "_Thalf",      # 8 AN
    "_Thalf_unit", # 9 AO
    "isomer",      # 10 AP
    "percentIT",   # 11 AQ
    "reaction",    # 12 AR
    "fast",        # 13 AS
    "thermalXS",   # 14 AT
    "gT",          # 15 AU
    "resonance",   # 16 AV
    "Thalf_hrs",   # 17 AW
    "Thalf_str",   # 18 AX
    "Thalf_parent", # 19 AY
    "thermalXS_parent",  # 20 AZ
    "resonance_parent",  # 21 BA
    "comments",          # 22 BB
]
INT_COLUMNS = [1, 2, 4]
BOOL_COLUMNS = [13]
FLOAT_COLUMNS = [6, 11, 14, 15, 16, 17, 19, 20, 21]
UNITS_TO_HOURS = {'y': 8760, 'd': 24, 'h': 1, 'm': 1/60, 's': 1/3600}

def activity(isotope, mass, env, exposure, rest_times):
    """
    Compute isotope specific daughter products after the given exposure time and rest period.
    """
    result = {}
    if not hasattr(isotope, 'neutron_activation'):
        return result

    for ai in isotope.neutron_activation:
        # Ignore fast neutron interactions if not using fast ratio
        if ai.fast and env.fast_ratio == 0:
            continue
        # Column D: elemental % mass content of sample
        #    mass fraction and abundance already included in mass calculation, so not needed
        # Column E: target nuclide and comment
        #    str(isotope), ai.comment
        # Column F: Nuclide Produced
        #    ai.daughter
        # Column G: Half-life
        #    ai.Thalf_str
        # Column H: initial effective cross-section (b)
        initialXS = ai.thermalXS + env.epithermal_reduction_factor*ai.resonance
        # Column I: reaction
        #    ai.reaction
        # Column J: fast?
        #    ai.fast
        # Column K: effective reaction flux (n/cm^2/s)
        flux = env.fluence/env.fast_ratio if ai.fast else env.fluence
        # Column L: root part of activation calculation
        # Decay correction portion done in column M
        # The given mass is sample mass * sample fraction * isotope abundance
        root = flux * initialXS * 1e-24 * mass / isotope.isotope * 1.6278e19
        # Column M: 0.69/t1/2  (1/h) lambda of produced nuclide
        lam = LN2/ai.Thalf_hrs
        #print(ai.thermalXS, ai.resonance, env.epithermal_reduction_factor)
        #print(isotope, "D", mass, "F", ai.daughter, "G", ai.Thalf_str,
        #      "H", initialXS, "I", ai.reaction, "J", ai.fast, "K", flux,
        #      "L", root, "M", lam)

        # Column Y: activity at the end of irradiation (uCi)
        if ai.reaction == 'b':
            # Column N: 0.69/t1/2 (1/h) lambda of parent nuclide
            parent_lam = LN2 / ai.Thalf_parent
            # Column O: Activation if "b" mode production
            # 2022-05-18 PAK: addressed the following
            #    Note: problems resulting from precision limitation not addressed
            #    in "b" mode production
            # by rewriting equation to use expm1:
            #    activity = root*(1 - exp(-lam*exposure)/(1 - (lam/parent_lam))
            #                     + exp(-parent_lam*exposure)/((parent_lam/lam)-1))
            # Let x1=-lam*exposure x2=-parent_lam*exposure a=x1/x2=lam/parent_lam
            #    activity = root * (1 - exp(x1)/(1-a) + exp(x2)/(1/a - 1))
            #             = root * (1 - x2*exp(x1)/(x2-x1)) + x1*exp(x2)/(x2-x1))
            #             = root * ((x2-x1) - x2*exp(x1) + x1*exp(x2)) / (x2-x1)
            #             = root * (x2*(1-exp(x1)) + x1*(exp(x2)-1)) / (x2-x1)
            #             = root * (x1*expm1(x2) - x2*expm1(x1)) / (x2-x1)
            #             = root * (lam*expm1(x2) - parent_lam*expm1(x1)) / (parent_lam - lam)
            # Checked for each b-mode production that small halflife results are
            # unchanged to four digits and Eu[151] => Gd[152] no longer fails.
            activity = root/(parent_lam - lam) * (
                lam*expm1(-parent_lam*exposure) - parent_lam*expm1(-lam*exposure))
            #print("N", parent_lam, "O", activity)
        elif ai.reaction == '2n':
            # Column N: 0.69/t1/2 (1/h) lambda of parent nuclide
            parent_lam = LN2 / ai.Thalf_parent
            # Column P: effective cross-section 2n product and n, g burnup (b)
            # Note: This cross-section always uses the total thermal flux
            effectiveXS = ai.thermalXS_parent + env.epithermal_reduction_factor*ai.resonance_parent
            # Column Q: 2n mode effective lambda of stable target (1/h)
            lam_2n = flux*initialXS*1e-24*3600
            # Column R: radioactive parent (1/h)
            parent_activity = env.fluence*1e-24*3600*effectiveXS+parent_lam
            # Column S: resulting product (1/h)
            product_2n = lam if ai.reaction == '2n' else 0
            # Column T: activity if 2n mode
            activity = root*lam*(parent_activity-parent_lam)*(
                (exp(-lam_2n*exposure)
                 / ((parent_activity-lam_2n)*(product_2n-lam_2n)))
                + (exp(-parent_activity*exposure)
                   / ((lam_2n-parent_activity)*(product_2n-parent_activity)))
                + (exp(-product_2n*exposure)
                   / ((lam_2n-product_2n)*(parent_activity-product_2n)))
                )
            #print("N", parent_lam, "P", effectiveXS, "Q", lam_2n,
            #      "R", parent_activity, "S", product_2n, "T", activity)
        else:
            # Provide the fix for the limitied precision (15 digits) in the
            # floating point calculation.  For neutron fluence rates above
            # 1e16 the precision in certain cells needs to be improved to
            # avoid erroneous results.  Also, burnup for single capture
            # reactions (excluding 'b') is included here.
            # See README file for details.

            # Note: This cross-section always uses the total thermal flux
            effectiveXS = ai.thermalXS_parent + env.epithermal_reduction_factor*ai.resonance_parent
            # Column U: nv1s1t
            U = flux*initialXS*3600*1e-24*exposure
            # Column V: nv2s2t+L2*t
            V = (env.fluence*effectiveXS*3600*1e-24+lam)*exposure
            # Column W: L/(L-nvs1+nvs2)
            W = lam/(lam-flux*initialXS*3600*1e-24+env.fluence*effectiveXS*3600*1e-24)
            # Column X: W*(exp(-U)-exp(V)) if U,V > 1e-10 else W*(V-U+(V-U)*(V+U)/2)
            # [PAK 2024-02-28] Rewrite the precision correction using the
            # high accuracy expm1() function to avoid cancellation.
            # See _check_exp_diff() below for verification.
            X = W*exp(-V)*expm1(V-U) if U>V else -W*exp(-U)*expm1(U-V)
            # Column Y: O if "b" else T if "2n" else L*X
            activity = root*X
            #print(f"{ai.isotope}=>{ai.daughter} {U=} {V=} {W=} {X=} {activity=}")

            # Check the effects of changing the formula. To look at all elements
            # at once, enable the following and use:
            #    $ python -m periodictable.activation -m 1000 -f 1e8 -e 1000
            if 0:
                if abs(U) < 1e-10 and abs(V) < 1e-10:
                    X = W * (V-U+(V+U)/2)
                    trigger = '*'
                else:
                    X = exp(-U) - exp(-V)
                    trigger = ''
                delta = activity - root*X
                if activity > 0 and delta/activity > 1e-4:
                    print(f"{trigger}activity change {delta:10.4e} from {activity:10.4e} ({100*delta/activity:.2f}%) for {ai.isotope} => {ai.daughter} ({ai.reaction})")

            if activity < 0:
                msg = "activity %g less than zero for %g"%(activity, isotope)
                raise RuntimeError(msg)
            #print(ai.thermalXS_parent, ai.resonance_parent, exposure)
            #print("P", effectiveXS, "U", U, "V", V, "W", W, "X",
            #      precision_correction, "Y", activity)
            # columns: F32 H K L U V W X
            #data = env.fluence, initialXS, flux, root, U, V, W, precision_correction
            #print(" ".join("%.5e"%v for v in data))

        result[ai] = [activity*exp(-lam*Ti) for Ti in rest_times]
        #print([(Ti, Ai) for Ti, Ai in zip(rest_times, result[ai])])

    return result

def init(table, reload=False):
    """
    Add neutron activation levels to each isotope.
    """
    from math import floor
    if 'neutron_activation' in table.properties and not reload:
        return
    table.properties.append('neutron_activation')

    # Clear the existing activation table
    for el in table:
        for iso in el.isotopes:
            if hasattr(el[iso], 'neutron_activation'):
                del el[iso].neutron_activation

    # We are keeping the table as a simple export of the activation data
    # from the ncnr health physics excel spreadsheet so that it is easier
    # to validate that the table contains the same data. Unfortunately some
    # of the cells involved formulas, which need to be reproduced when loading
    # in order to match full double precision values.
    path = os.path.join(core.get_data_path('.'), 'activation.dat')
    with open(path, 'r') as fd:
        for row in fd:
            #print(row, end='')
            columns = row.split('\t')
            if columns[0].strip() in ('', 'xx'):
                continue
            columns = [c[1:-1] if c.startswith('"') else c
                    for c in columns]
            #print columns
            for c in INT_COLUMNS:
                columns[c] = int(columns[c])
            for c in BOOL_COLUMNS:
                columns[c] = (columns[c] == 'y')
            for c in FLOAT_COLUMNS:
                columns[c] = float(columns[c]) if columns[c].strip() else 0.
            # clean up comment column
            columns[-1] = columns[-1].replace('"', '').strip()
            kw = dict(zip(COLUMN_NAMES, columns))
            # Recreate Thalf_hrs column using double precision.
            # Note: spreadsheet is not converting half-life to hours in cell AW1462 (186-W => 188-W)
            kw['Thalf_hrs'] = float(kw['_Thalf']) * UNITS_TO_HOURS[kw['_Thalf_unit']]
            # Recreate Thalf_parent by fetching from the new Thalf_hrs
            # e.g., =IF(OR(AR1408="2n",AR1408="b"),IF(AR1407="b",AW1406,AW1407),"")
            # This requires that the parent is directly before the 'b' or 'nb'
            # with its activation list already entered into the isotope.
            # Note: 150-Nd has 'act' followed by two consecutive 'b' entries.
            if kw['reaction'] in ('b', '2n'):
                act = table[kw['Z']][kw['A']].neutron_activation
                parent = act[-2] if act[-1].reaction == 'b' else act[-1]
                kw['Thalf_parent'] = parent.Thalf_hrs
            else:
                #assert kw['Thalf_parent'] == 0
                kw['Thalf_parent'] = None
            # Half-lives use My, Gy, Ty, Py
            value, units = float(kw['_Thalf']), kw['_Thalf_unit']
            if units == 'y':
                if value >= 1e12:
                    value, units = value/1e12, 'Ty'
                elif value >= 1e9:
                    value, units = value/1e9, 'Gy'
                elif value >= 200e3: # above 200,000 years use My
                    value, units = value/1e6, 'My'
                elif value >= 20e3: # between 20,000 and 200,000 use ky
                    value, units = value/1e3, 'ky'
            formatted = f"{value:g} {units}"
            #if formatted.replace(' ', '') != kw['Thalf_str'].replace(' ', ''):
            #    print(f"{kw['_index']}: old {kw['Thalf_str']} != {formatted} new")
            kw['Thalf_str'] = formatted

            # Strip columns whose names start with underscore
            kw = dict((k, v) for k, v in kw.items() if not k.startswith('_'))

            # Create an Activation record and add it to the isotope
            iso = table[kw['Z']][kw['A']]
            activation = getattr(iso, 'neutron_activation', [])
            activation.append(ActivationResult(**kw))
            iso.neutron_activation = activation

            # Check abundance values
            #if abs(iso.abundance - kw['abundance']) > 0.001*kw['abundance']:
            #    percent = 100*abs(iso.abundance - kw['abundance'])/kw['abundance']
            #    print "Abundance of", iso, "is", iso.abundance, \
            #        "but activation.dat has", kw['abundance'], "(%.1f%%)"%percent

class ActivationResult(object):
    def __init__(self, **kw):
        self.__dict__ = kw
    def __str__(self):
        return str(self.__dict__)


def _check_exp_diff(p, delta):
    """
    Compare methods for computing exp(-U) - exp(-V).

    The output shows error relative as log10(|(f - m)/m|) where 
    m is computed with 500 bits of precision:

    * Δexp uses exp(-V) - exp(-U)
    * Δexpm1 uses exp(-V)*expm1(V-U)
    * Δtaylor uses (V-U) + (U**2-V**2)/2 = (V-U)*(1-(V+U)/2)
    * Δoriginal uses (V-U) + (V+U)/2

    Call with (p, δ) such that U = 10**p, V = 10**p + 10**(p-δ). That is,
    *p* gives the order of magnitude of the function and *delta* gives the
    relative magnitude of the difference. For example, (-3,2)=>(0.001,0.00101)

    Try with::

        for k in range(-20,3,3): _check_exp_diff(k,8)

    The spreadsheet uses the original expression when U, V are less than 1e-10
    but the approximation is poor even in that region.

    The Δexpm1 is strictly better than the direct calculatio almost everywhere.
    Where it is worse it is only slightly worse, and the calculated activity
    is small. It is at least as good as the Taylor approximation everywhere.

    The expm1 form gives us at least 0.1% precision relative to the high
    accuracy calculation. This occurs when U,V differ in the final digit
    of a double precision number. When U,V differ in the 10th digit we
    get 6-7 digits of precision. If they differ in the 3rd digit we get the
    full precision for floating point double values. The precision is roughly
    constant across the domain.
    """
    import mpmath as mp
    import math
    with mp.workprec(500):
        mp_x, mp_y = mp.mpf(10)**p, mp.mpf(10)**p + mp.mpf(10)**(p-delta)
        mp_diff = mp.exp(-mp_x) - mp.exp(-mp_y)
        #print(f"{mp_x} {mp_y} {mp_diff}")
    math_x, math_y = 10**p, 10**p + 10**(p-delta)
    #print(f"{math_y}-{math_x}={math_y - math_x}")
    math_diff = math.exp(-math_x) - math.exp(-math_y)
    math_alt = math.exp(-math_y)*math.expm1(math_y-math_x)
    math_orig = math_y-math_x + (math_y+math_x)/2
    math_tay = math_y-math_x + (math_x**2-math_y**2)/2
    math_tay = (math_y-math_x)*(1 - (math_x+math_y)/2)
    def err(v): return f"{math.log10(float(abs((v-mp_diff)/mp_diff))):4.1f}"
    #print(f"mp:{mp_diff}  direct:{math_diff}  alt:{math_alt}  taylor:{math_tay}")
    print(f"U:1e{p:<3d} mp: {mp.nstr(mp_diff,10):15s}  Δexp: {err(math_diff)}  Δexpm1: {err(math_alt)}  Δtaylor: {err(math_tay)}  Δoriginal: {err(math_orig)}")

def demo():  # pragma: nocover
    import sys
    import argparse
    import periodictable as pt

    # From chatGPT
    parser = argparse.ArgumentParser(description='Process some data with mass, flux, exposure, and decay options.')
    parser.add_argument('-m', '--mass', type=float, default=1, help='Specify the mass value (default: 1)')
    parser.add_argument('-f', '--flux', type=float, default=1e8, help='Specify the flux value (default: 1e8)')
    parser.add_argument('-e', '--exposure', type=float, default=10, help='Specify the exposure value (default: 10)')
    parser.add_argument('-d', '--decay', type=float, default=5e-4, help='Specify the decay value (default: 5e-4)')
    parser.add_argument('--cd-ratio', type=float, default=70, help='Specify the Cd ratio value (default: 70)')
    parser.add_argument('--fast-ratio', type=float, default=50, help='Specify the fast ratio value (default: 50)')

    parser.add_argument('formula', nargs='?', type=str, default=None, help='Specify the formula as a positional argument')
    args = parser.parse_args()

    formula = args.formula
    if formula is None:
        # Make sure all elements compute
        #formula = "".join(str(el) for el in pt.elements)[1:]
        formula = build_formula([(1/el.mass, el) for el in pt.elements][1:])
        # Use an enormous mass to force significant activation of rare isotopes
        mass, fluence = 1e15, 1e8
    env = ActivationEnvironment(fluence=args.flux, Cd_ratio=args.cd_ratio, fast_ratio=args.fast_ratio, location="BT-2")
    sample = Sample(formula, mass=args.mass)
    abundance = IAEA1987_isotopic_abundance
    #abundance=NIST2001_isotopic_abundance,
    sample.calculate_activation(
        env, exposure=args.exposure, rest_times=(0, 1, 24, 360),
        abundance=abundance,
        )
    decay_time = sample.decay_time(args.decay)
    print(f"{args.mass} g {formula} for {args.exposure} hours at {args.flux} n/cm^2/s")
    print(f"Time to decay to {args.decay} uCi is {decay_time} hours.")
    sample.calculate_activation(
        env, exposure=args.exposure, rest_times=(0, 1, 24, 360, decay_time),
        abundance=abundance,
        )
    sample.show_table(cutoff=0.0)

    ## Print a table of flux vs. activity so we can debug the
    ## precision_correction value in the activity() function.
    ## Note that you also need to uncomment the print statement
    ## at the end of activity() that shows the column values.
    #import numpy as np
    #sample = Sample('Co', mass=10)
    #for fluence in np.logspace(3, 20, 20-3+1):
    #    env = ActivationEnvironment(fluence=fluence)
    #    sample.calculate_activation(
    #        env, exposure=exposure, rest_times=[0],
    #        abundance=IAEA1987_isotopic_abundance)

if __name__ == "__main__":
    #for k in range(-20,3,3): _check_exp_diff(k,10)
    demo()  # pragma: nocover
