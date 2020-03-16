# This program is in the public domain
# Author: Paul Kienzle
"""

Adds magnetic_ff[charge].t for t in j0, j2, j4, j6, and J.
J should be the dipole approximation <j0> + (1 - 2/g) <j2>, according to the
documentation for CrysFML [#Brown]_ , but that does not seem to be the case in practice.


.. [#Brown] Brown. P. J. (Section 4.4.5) International Tables for Crystallography
        Volume C, Wilson. A. J. C.(ed).
"""
from __future__ import division

import numpy
from numpy import pi, exp

def formfactor_0(j0, q):
    """
    Returns the scattering potential for form factor *j0* at the given *q*.
    """
    q = numpy.asarray(q)
    s_sq = (q/(4*pi))**2
    A, a, B, b, C, c, D = j0
    return A * exp(-a*s_sq) + B * exp(-b*s_sq) + C * exp(-c*s_sq) + D

def formfactor_n(jn, q):
    """
    Returns the scattering potential for form factor *jn* at the given *q*.
    """
    q = numpy.asarray(q)
    s_sq = (q/(4*pi))**2
    A, a, B, b, C, c, D = jn
    return s_sq * (A * exp(-a*s_sq) + B * exp(-b*s_sq) + C * exp(-c*s_sq) + D)


class MagneticFormFactor(object):
    """
    Magnetic form factor for the ion.

    The available form factors are::

        M = <j0> form factor coefficients
        J = <j0> + C2 <j2> form factor coeffients
        jn = <jn> form factor coefficients for n = 0, 2, 4, 6

    Not all form factors are available for all ions.  Use the
    expression ``hasattr(ion.magnetic_ff, '<ff>')`` to test for the
    particular form factor <ff>.
    The form factor coefficients are a tuple (A, a, B, b, C, c, D).  The
    following expression computes the M/j0 and J form factors from the
    corresponding coefficients::

        s = q^2 / 16 pi^2
        ff = A exp(-a s^2) + B exp(-b s^2) + C exp(-c s^2) + D

    The remaining form factors *j2*, *j4* and *j6* are scalled by an additional s^2.
    The form factor calculation is performed by the <ff>_Q method for <ff>
    in *M*, *J*, *j0*, *j2*, *j4*, *j6*.  For example, here is the calculation for
    the *M* form factor for Fe^2+ computed at 0, 0.1 and 0.2:

    .. doctest::

        >>> import periodictable
        >>> ion = periodictable.Fe.ion[2]
        >>> print("[%.5f, %.5f, %.5f]"
        ...       % tuple(ion.magnetic_ff[ion.charge].M_Q([0, 0.1, 0.2])))
        [1.00000, 0.99935, 0.99741]

    """


    def _getM(self):
        return self.j0

    M = property(_getM, doc="j0")

    def j0_Q(self, Q):
        """Returns *j0* scattering potential at *Q* |1/Ang|"""
        return formfactor_0(self.j0, Q)

    def j2_Q(self, Q):
        """Returns *j2* scattering potential at *Q* |1/Ang|"""
        return formfactor_n(self.j2, Q)

    def j4_Q(self, Q):
        """Returns *j4* scattering potential at *Q* |1/Ang|"""
        return formfactor_n(self.j4, Q)

    def j6_Q(self, Q):
        """Returns j6 scattering potential at *Q* |1/Ang|"""
        return formfactor_n(self.j6, Q)

    def J_Q(self, Q):
        """Returns J scattering potential at *Q* |1/Ang|"""
        return formfactor_0(self.J, Q)

    M_Q = j0_Q


def init(table, reload=False):
    """Add magnetic form factor properties to the periodic table"""
    if 'magnetic_ff' in table.properties and not reload:
        return
    table.properties.append('magnetic_ff')

    # Function for interpreting ionization state and form factor tuple
    def Magnetic_Form_Type(state, values):
        return state, values

    # Remove fortran continuation character, which is & at the end of the line
    data = CFML_DATA.replace('&\n', '')

    # Parse each line
    for line in data.split('\n'):
        # Skip blank lines and comment lines
        line = line.strip()
        if '=' not in line:
            continue

        # Remaining lines contain one of:
        #   Magnetic_Form(n) = Magnetic_Form_Type("M<EL><ION>", (/FF/))
        #   Magnetic_Form(n) = Magnetic_Form_Type("J<EL><ION>", (/FF/))
        #   Magnetic_j2(n) = Magnetic_Form_Type("<EL><ION>", (/FF/))
        #   Magnetic_j4(n) = Magnetic_Form_Type("<EL><ION>", (/FF/))
        #   Magnetic_j6(n) = Magnetic_Form_Type("<EL><ION>", (/FF/))
        # The 'M' form is just j0.  The 'J' form should be j0 + (1-g/2)j2.
        # Split the two parts, remove the slash and figure out if we have
        # j0, j2, j4, j6 or J.
        a, b = line.split('=')
        b = b.replace('/', '')
        # The following eval is used for parsing the CFM_DATA table defined
        # below, so there is no risk that it can be invoked from a web
        # service with a user-defined malicious string.
        state, values = eval(b)  # eval checked for malicious usage
        if a.startswith('Magnetic_Form'):
            jn = "j0" if state[0] == 'M' else "J"
            state = state[1:]
        elif a.startswith('Magnetic_j2'):
            jn = "j2"
        elif a.startswith('Magnetic_j4'):
            jn = "j4"
        elif a.startswith('Magnetic_j6'):
            jn = "j6"

        # Parse <EL><ION> into element symbol and ion state
        if state[1].isdigit():
            symbol = state[0]
            charge = int(state[1])
        else:
            symbol = state[0:2].capitalize()
            charge = int(state[2])

        # Add the magnetic form factor info to the element
        el = table.symbol(symbol)
        if not hasattr(el, 'magnetic_ff'):
            el.magnetic_ff = {}
        if charge not in el.magnetic_ff:
            el.magnetic_ff[charge] = MagneticFormFactor()
        setattr(el.magnetic_ff[charge], jn, values)

CFML_DATA = """
       Magnetic_Form(  1) = Magnetic_Form_Type("MSC0", &
                                              (/  0.251200, 90.029602,  0.329000, 39.402100,  0.423500, 14.322200, -0.004300/) )
       Magnetic_Form(  2) = Magnetic_Form_Type("MSC1", &
                                              (/  0.488900, 51.160301,  0.520300, 14.076400, -0.028600,  0.179200,  0.018500/) )
       Magnetic_Form(  3) = Magnetic_Form_Type("MSC2", &
                                              (/  0.504800, 31.403500,  0.518600, 10.989700, -0.024100,  1.183100,  0.000000/) )
       Magnetic_Form(  4) = Magnetic_Form_Type("MTI0", &
                                              (/  0.465700, 33.589802,  0.549000,  9.879100, -0.029100,  0.323200,  0.012300/) )
       Magnetic_Form(  5) = Magnetic_Form_Type("MTI1", &
                                              (/  0.509300, 36.703300,  0.503200, 10.371300, -0.026300,  0.310600,  0.011600/) )
       Magnetic_Form(  6) = Magnetic_Form_Type("MTI2", &
                                              (/  0.509100, 24.976299,  0.516200,  8.756900, -0.028100,  0.916000,  0.001500/) )
       Magnetic_Form(  7) = Magnetic_Form_Type("MTI3", &
                                              (/  0.357100, 22.841299,  0.668800,  8.930600, -0.035400,  0.483300,  0.009900/) )
       Magnetic_Form(  8) = Magnetic_Form_Type("MV0 ", &
                                              (/  0.408600, 28.810900,  0.607700,  8.543700, -0.029500,  0.276800,  0.012300/) )
       Magnetic_Form(  9) = Magnetic_Form_Type("MV1 ", &
                                              (/  0.444400, 32.647900,  0.568300,  9.097100, -0.228500,  0.021800,  0.215000/) )
       Magnetic_Form( 10) = Magnetic_Form_Type("MV2 ", &
                                              (/  0.408500, 23.852600,  0.609100,  8.245600, -0.167600,  0.041500,  0.149600/) )
       Magnetic_Form( 11) = Magnetic_Form_Type("MV3 ", &
                                              (/  0.359800, 19.336399,  0.663200,  7.617200, -0.306400,  0.029600,  0.283500/) )
       Magnetic_Form( 12) = Magnetic_Form_Type("MV4 ", &
                                              (/  0.310600, 16.816000,  0.719800,  7.048700, -0.052100,  0.302000,  0.022100/) )
       Magnetic_Form( 13) = Magnetic_Form_Type("MCR0", &
                                              (/  0.113500, 45.199001,  0.348100, 19.493099,  0.547700,  7.354200, -0.009200/) )
       Magnetic_Form( 14) = Magnetic_Form_Type("MCR1", &
                                              (/ -0.097700,  0.047000,  0.454400, 26.005400,  0.557900,  7.489200,  0.083100/) )
       Magnetic_Form( 15) = Magnetic_Form_Type("MCR2", &
                                              (/  1.202400, -0.005500,  0.415800, 20.547501,  0.603200,  6.956000, -1.221800/) )
       Magnetic_Form( 16) = Magnetic_Form_Type("MCR3", &
                                              (/ -0.309400,  0.027400,  0.368000, 17.035500,  0.655900,  6.523600,  0.285600/) )
       Magnetic_Form( 17) = Magnetic_Form_Type("MCR4", &
                                              (/ -0.232000,  0.043300,  0.310100, 14.951800,  0.718200,  6.172600,  0.204200/) )
       Magnetic_Form( 18) = Magnetic_Form_Type("MMN0", &
                                              (/  0.243800, 24.962900,  0.147200, 15.672800,  0.618900,  6.540300, -0.010500/) )
       Magnetic_Form( 19) = Magnetic_Form_Type("MMN1", &
                                              (/ -0.013800,  0.421300,  0.423100, 24.667999,  0.590500,  6.654500, -0.001000/) )
       Magnetic_Form( 20) = Magnetic_Form_Type("MMN2", &
                                              (/  0.422000, 17.684000,  0.594800,  6.005000,  0.004300, -0.609000, -0.021900/) )
       Magnetic_Form( 21) = Magnetic_Form_Type("MMN3", &
                                              (/  0.419800, 14.282900,  0.605400,  5.468900,  0.924100, -0.008800, -0.949800/) )
       Magnetic_Form( 22) = Magnetic_Form_Type("MMN4", &
                                              (/  0.376000, 12.566100,  0.660200,  5.132900, -0.037200,  0.563000,  0.001100/) )
       Magnetic_Form( 23) = Magnetic_Form_Type("MFE0", &
                                              (/  0.070600, 35.008499,  0.358900, 15.358300,  0.581900,  5.560600, -0.011400/) )
       Magnetic_Form( 24) = Magnetic_Form_Type("MFE1", &
                                              (/  0.125100, 34.963299,  0.362900, 15.514400,  0.522300,  5.591400, -0.010500/) )
       Magnetic_Form( 25) = Magnetic_Form_Type("MFE2", &
                                              (/  0.026300, 34.959702,  0.366800, 15.943500,  0.618800,  5.593500, -0.011900/) )
       Magnetic_Form( 26) = Magnetic_Form_Type("MFE3", &
                                              (/  0.397200, 13.244200,  0.629500,  4.903400, -0.031400,  0.349600,  0.004400/) )
       Magnetic_Form( 27) = Magnetic_Form_Type("MFE4", &
                                              (/  0.378200, 11.380000,  0.655600,  4.592000, -0.034600,  0.483300,  0.000500/) )
       Magnetic_Form( 28) = Magnetic_Form_Type("MCO0", &
                                              (/  0.413900, 16.161600,  0.601300,  4.780500, -0.151800,  0.021000,  0.134500/) )
       Magnetic_Form( 29) = Magnetic_Form_Type("MCO1", &
                                              (/  0.099000, 33.125198,  0.364500, 15.176800,  0.547000,  5.008100, -0.010900/) )
       Magnetic_Form( 30) = Magnetic_Form_Type("MCO2", &
                                              (/  0.433200, 14.355300,  0.585700,  4.607700, -0.038200,  0.133800,  0.017900/) )
       Magnetic_Form( 31) = Magnetic_Form_Type("MCO3", &
                                              (/  0.390200, 12.507800,  0.632400,  4.457400, -0.150000,  0.034300,  0.127200/) )
       Magnetic_Form( 32) = Magnetic_Form_Type("MCO4", &
                                              (/  0.351500, 10.778500,  0.677800,  4.234300, -0.038900,  0.240900,  0.009800/) )
       Magnetic_Form( 33) = Magnetic_Form_Type("MNI0", &
                                              (/ -0.017200, 35.739201,  0.317400, 14.268900,  0.713600,  4.566100, -0.014300/) )
       Magnetic_Form( 34) = Magnetic_Form_Type("MNI1", &
                                              (/  0.070500, 35.856098,  0.398400, 13.804200,  0.542700,  4.396500, -0.011800/) )
       Magnetic_Form( 35) = Magnetic_Form_Type("MNI2", &
                                              (/  0.016300, 35.882599,  0.391600, 13.223300,  0.605200,  4.338800, -0.013300/) )
       Magnetic_Form( 36) = Magnetic_Form_Type("MNI3", &
                                              (/ -0.013400, 35.867699,  0.267800, 12.332600,  0.761400,  4.236900, -0.016200/) )
       Magnetic_Form( 37) = Magnetic_Form_Type("MNI4", &
                                              (/ -0.009000, 35.861401,  0.277600, 11.790400,  0.747400,  4.201100, -0.016300/) )
       Magnetic_Form( 38) = Magnetic_Form_Type("MCU0", &
                                              (/  0.090900, 34.983799,  0.408800, 11.443200,  0.512800,  3.824800, -0.012400/) )
       Magnetic_Form( 39) = Magnetic_Form_Type("MCU1", &
                                              (/  0.074900, 34.965599,  0.414700, 11.764200,  0.523800,  3.849700, -0.012700/) )
       Magnetic_Form( 40) = Magnetic_Form_Type("MCU2", &
                                              (/  0.023200, 34.968601,  0.402300, 11.564000,  0.588200,  3.842800, -0.013700/) )
       Magnetic_Form( 41) = Magnetic_Form_Type("MCU3", &
                                              (/  0.003100, 34.907398,  0.358200, 10.913800,  0.653100,  3.827900, -0.014700/) )
       Magnetic_Form( 42) = Magnetic_Form_Type("MCU4", &
                                              (/ -0.013200, 30.681700,  0.280100, 11.162600,  0.749000,  3.817200, -0.016500/) )
       Magnetic_Form( 43) = Magnetic_Form_Type("MY0 ", &
                                              (/  0.591500, 67.608101,  1.512300, 17.900400, -1.113000, 14.135900,  0.008000/) )
       Magnetic_Form( 44) = Magnetic_Form_Type("MZR0", &
                                              (/  0.410600, 59.996101,  1.054300, 18.647600, -0.475100, 10.540000,  0.010600/) )
       Magnetic_Form( 45) = Magnetic_Form_Type("MZR1", &
                                              (/  0.453200, 59.594799,  0.783400, 21.435699, -0.245100,  9.036000,  0.009800/) )
       Magnetic_Form( 46) = Magnetic_Form_Type("MNB0", &
                                              (/  0.394600, 49.229698,  1.319700, 14.821600, -0.726900,  9.615600,  0.012900/) )
       Magnetic_Form( 47) = Magnetic_Form_Type("MNB1", &
                                              (/  0.457200, 49.918201,  1.027400, 15.725600, -0.496200,  9.157300,  0.011800/) )
       Magnetic_Form( 48) = Magnetic_Form_Type("MMO0", &
                                              (/  0.180600, 49.056801,  1.230600, 14.785900, -0.426800,  6.986600,  0.017100/) )
       Magnetic_Form( 49) = Magnetic_Form_Type("MMO1", &
                                              (/  0.350000, 48.035400,  1.030500, 15.060400, -0.392900,  7.479000,  0.013900/) )
       Magnetic_Form( 50) = Magnetic_Form_Type("MTC0", &
                                              (/  0.129800, 49.661098,  1.165600, 14.130700, -0.313400,  5.512900,  0.019500/) )
       Magnetic_Form( 51) = Magnetic_Form_Type("MTC1", &
                                              (/  0.267400, 48.956600,  0.956900, 15.141300, -0.238700,  5.457800,  0.016000/) )
       Magnetic_Form( 52) = Magnetic_Form_Type("MRU0", &
                                              (/  0.106900, 49.423801,  1.191200, 12.741700, -0.317600,  4.912500,  0.021300/) )
       Magnetic_Form( 53) = Magnetic_Form_Type("MRU1", &
                                              (/  0.441000, 33.308601,  1.477500,  9.553100, -0.936100,  6.722000,  0.017600/) )
       Magnetic_Form( 54) = Magnetic_Form_Type("MRH0", &
                                              (/  0.097600, 49.882500,  1.160100, 11.830700, -0.278900,  4.126600,  0.023400/) )
       Magnetic_Form( 55) = Magnetic_Form_Type("MRH1", &
                                              (/  0.334200, 29.756399,  1.220900,  9.438400, -0.575500,  5.332000,  0.021000/) )
       Magnetic_Form( 56) = Magnetic_Form_Type("MPD0", &
                                              (/  0.200300, 29.363300,  1.144600,  9.599300, -0.368900,  4.042300,  0.025100/) )
       Magnetic_Form( 57) = Magnetic_Form_Type("MPD1", &
                                              (/  0.503300, 24.503700,  1.998200,  6.908200, -1.524000,  5.513300,  0.021300/) )
       Magnetic_Form( 58) = Magnetic_Form_Type("MCE2", &
                                              (/  0.295300, 17.684601,  0.292300,  6.732900,  0.431300,  5.382700, -0.019400/) )
       Magnetic_Form( 59) = Magnetic_Form_Type("MND2", &
                                              (/  0.164500, 25.045300,  0.252200, 11.978200,  0.601200,  4.946100, -0.018000/) )
       Magnetic_Form( 60) = Magnetic_Form_Type("MND3", &
                                              (/  0.054000, 25.029301,  0.310100, 12.102000,  0.657500,  4.722300, -0.021600/) )
       Magnetic_Form( 61) = Magnetic_Form_Type("MSM2", &
                                              (/  0.090900, 25.203199,  0.303700, 11.856200,  0.625000,  4.236600, -0.020000/) )
       Magnetic_Form( 62) = Magnetic_Form_Type("MSM3", &
                                              (/  0.028800, 25.206800,  0.297300, 11.831100,  0.695400,  4.211700, -0.021300/) )
       Magnetic_Form( 63) = Magnetic_Form_Type("MEU2", &
                                              (/  0.075500, 25.296000,  0.300100, 11.599300,  0.643800,  4.025200, -0.019600/) )
       Magnetic_Form( 64) = Magnetic_Form_Type("MEU3", &
                                              (/  0.020400, 25.307800,  0.301000, 11.474400,  0.700500,  3.942000, -0.022000/) )
       Magnetic_Form( 65) = Magnetic_Form_Type("MGD2", &
                                              (/  0.063600, 25.382299,  0.303300, 11.212500,  0.652800,  3.787700, -0.019900/) )
       Magnetic_Form( 66) = Magnetic_Form_Type("MGD3", &
                                              (/  0.018600, 25.386700,  0.289500, 11.142100,  0.713500,  3.752000, -0.021700/) )
       Magnetic_Form( 67) = Magnetic_Form_Type("MTB2", &
                                              (/  0.054700, 25.508600,  0.317100, 10.591100,  0.649000,  3.517100, -0.021200/) )
       Magnetic_Form( 68) = Magnetic_Form_Type("MTB3", &
                                              (/  0.017700, 25.509501,  0.292100, 10.576900,  0.713300,  3.512200, -0.023100/) )
       Magnetic_Form( 69) = Magnetic_Form_Type("MDY2", &
                                              (/  0.130800, 18.315500,  0.311800,  7.664500,  0.579500,  3.146900, -0.022600/) )
       Magnetic_Form( 70) = Magnetic_Form_Type("MDY3", &
                                              (/  0.115700, 15.073200,  0.327000,  6.799100,  0.582100,  3.020200, -0.024900/) )
       Magnetic_Form( 71) = Magnetic_Form_Type("MHO2", &
                                              (/  0.099500, 18.176100,  0.330500,  7.855600,  0.592100,  2.979900, -0.023000/) )
       Magnetic_Form( 72) = Magnetic_Form_Type("MHO3", &
                                              (/  0.056600, 18.317600,  0.336500,  7.688000,  0.631700,  2.942700, -0.024800/) )
       Magnetic_Form( 73) = Magnetic_Form_Type("MER2", &
                                              (/  0.112200, 18.122299,  0.346200,  6.910600,  0.564900,  2.761400, -0.023500/) )
       Magnetic_Form( 74) = Magnetic_Form_Type("MER3", &
                                              (/  0.058600, 17.980200,  0.354000,  7.096400,  0.612600,  2.748200, -0.025100/) )
       Magnetic_Form( 75) = Magnetic_Form_Type("MTM2", &
                                              (/  0.098300, 18.323601,  0.338000,  6.917800,  0.587500,  2.662200, -0.024100/) )
       Magnetic_Form( 76) = Magnetic_Form_Type("MTM3", &
                                              (/  0.058100, 15.092200,  0.278700,  7.801500,  0.685400,  2.793100, -0.022400/) )
       Magnetic_Form( 77) = Magnetic_Form_Type("MYB2", &
                                              (/  0.085500, 18.512300,  0.294300,  7.373400,  0.641200,  2.677700, -0.021300/) )
       Magnetic_Form( 78) = Magnetic_Form_Type("MYB3", &
                                              (/  0.041600, 16.094900,  0.284900,  7.834100,  0.696100,  2.672500, -0.022900/) )
       Magnetic_Form( 79) = Magnetic_Form_Type("MU3 ", &
                                              (/  0.505800, 23.288200,  1.346400,  7.002800, -0.872400,  4.868300,  0.019200/) )
       Magnetic_Form( 80) = Magnetic_Form_Type("MU4 ", &
                                              (/  0.329100, 23.547501,  1.083600,  8.454000, -0.434000,  4.119600,  0.021400/) )
       Magnetic_Form( 81) = Magnetic_Form_Type("MU5 ", &
                                              (/  0.365000, 19.803801,  3.219900,  6.281800, -2.607700,  5.301000,  0.023300/) )
       Magnetic_Form( 82) = Magnetic_Form_Type("MNP3", &
                                              (/  0.515700, 20.865400,  2.278400,  5.893000, -1.816300,  4.845700,  0.021100/) )
       Magnetic_Form( 83) = Magnetic_Form_Type("MNP4", &
                                              (/  0.420600, 19.804600,  2.800400,  5.978300, -2.243600,  4.984800,  0.022800/) )
       Magnetic_Form( 84) = Magnetic_Form_Type("MNP5", &
                                              (/  0.369200, 18.190001,  3.151000,  5.850000, -2.544600,  4.916400,  0.024800/) )
       Magnetic_Form( 85) = Magnetic_Form_Type("MNP6", &
                                              (/  0.292900, 17.561100,  3.486600,  5.784700, -2.806600,  4.870700,  0.026700/) )
       Magnetic_Form( 86) = Magnetic_Form_Type("MPU3", &
                                              (/  0.384000, 16.679300,  3.104900,  5.421000, -2.514800,  4.551200,  0.026300/) )
       Magnetic_Form( 87) = Magnetic_Form_Type("MPU4", &
                                              (/  0.493400, 16.835501,  1.639400,  5.638400, -1.158100,  4.139900,  0.024800/) )
       Magnetic_Form( 88) = Magnetic_Form_Type("MPU5", &
                                              (/  0.388800, 16.559200,  2.036200,  5.656700, -1.451500,  4.255200,  0.026700/) )
       Magnetic_Form( 89) = Magnetic_Form_Type("MPU6", &
                                              (/  0.317200, 16.050699,  3.465400,  5.350700, -2.810200,  4.513300,  0.028100/) )
       Magnetic_Form( 90) = Magnetic_Form_Type("MAM2", &
                                              (/  0.474300, 21.776100,  1.580000,  5.690200, -1.077900,  4.145100,  0.021800/) )
       Magnetic_Form( 91) = Magnetic_Form_Type("MAM3", &
                                              (/  0.423900, 19.573900,  1.457300,  5.872200, -0.905200,  3.968200,  0.023800/) )
       Magnetic_Form( 92) = Magnetic_Form_Type("MAM4", &
                                              (/  0.373700, 17.862499,  1.352100,  6.042600, -0.751400,  3.719900,  0.025800/) )
       Magnetic_Form( 93) = Magnetic_Form_Type("MAM5", &
                                              (/  0.295600, 17.372499,  1.452500,  6.073400, -0.775500,  3.661900,  0.027700/) )
       Magnetic_Form( 94) = Magnetic_Form_Type("MAM6", &
                                              (/  0.230200, 16.953300,  1.486400,  6.115900, -0.745700,  3.542600,  0.029400/) )
       Magnetic_Form( 95) = Magnetic_Form_Type("MAM7", &
                                              (/  0.360100, 12.729900,  1.964000,  5.120300, -1.356000,  3.714200,  0.031600/) )
       Magnetic_Form( 96) = Magnetic_Form_Type("MPR3", &
                                              (/  0.050400, 24.998900,  0.257200, 12.037700,  0.714200,  5.003900, -0.021900/) )
       Magnetic_Form( 97) = Magnetic_Form_Type("MO1", &
                                              (/  0.115285, 85.197300,  0.556229, 25.252200,  0.332476,  6.362070, -0.00460676/) )
       Magnetic_Form( 98) = Magnetic_Form_Type("JCE2", &
                                                 (/  0.031972,  8.926222,  0.265792,  7.678510,  0.682151,  2.329783,  0.020578/) )
       Magnetic_Form( 99) = Magnetic_Form_Type("JCE3", &
                                              (/  0.051183,  6.115375,  0.277738,  7.952485,  0.654079,  2.287000,  0.016355/) )
       Magnetic_Form(100) = Magnetic_Form_Type("JPR3", &
                                              (/  0.023288,  0.582954,  0.349391,  5.601756,  0.615363,  1.932779,  0.011454/) )
       Magnetic_Form(101) = Magnetic_Form_Type("JND2", &
                                              (/  0.089354,  2.282004,  0.206157,  1.708607,  0.669916,  2.297662,  0.048390/) )
       Magnetic_Form(102) = Magnetic_Form_Type("JND3", &
                                              (/  0.073287,  4.412361,  0.371485,  4.019648,  0.539459,  1.557985,  0.017335/) )
       Magnetic_Form(103) = Magnetic_Form_Type("JGD3", &
                                              (/  0.060537, 10.775218,  0.271475, 13.097898,  0.665241,  3.162837,  0.001566/) )
       Magnetic_Form(104) = Magnetic_Form_Type("JTB2", &
                                              (/  0.049801, 18.734161,  0.277437, 10.084129,  0.661194,  2.745624,  0.010774/) )
       Magnetic_Form(105) = Magnetic_Form_Type("JTB3", &
                                              (/  0.049792, 15.112189,  0.270644,  9.158312,  0.679388,  2.880260, -0.000131/) )
       Magnetic_Form(106) = Magnetic_Form_Type("JDY2", &
                                              (/  0.175586,  5.938148,  0.228867, 11.464046,  0.583298,  2.167554,  0.011186/) )
       Magnetic_Form(107) = Magnetic_Form_Type("JDY3", &
                                              (/  0.146536, 12.639305,  0.375822,  5.511785,  0.515731,  2.090789,  0.093576/) )
       Magnetic_Form(108) = Magnetic_Form_Type("JHO2", &
                                              (/  0.023234,  0.703240,  0.270745,  9.993475,  0.677581,  2.521403,  0.027101/) )
       Magnetic_Form(109) = Magnetic_Form_Type("JHO2", &
                                              (/  0.023234,  0.703240,  0.270745,  9.993475,  0.677581,  2.521403,  0.027101/) )
       Magnetic_Form(110) = Magnetic_Form_Type("JHO3", &
                                              (/  0.043204,  0.910121,  0.279392,  8.683387,  0.668537,  2.417518,  0.008207/) )
       Magnetic_Form(111) = Magnetic_Form_Type("JER2", &
                                              (/  0.037734,  6.081446,  0.256447,  9.598846,  0.679204,  2.139296,  0.025543/) )
       Magnetic_Form(112) = Magnetic_Form_Type("JER3", &
                                              (/  0.038871,  5.311772,  0.259781,  8.173226,  0.678414,  2.082836,  0.022169/) )
       Magnetic_Form(113) = Magnetic_Form_Type("JTM2", &
                                              (/  0.037670,  4.455198,  0.254184,  9.151058,  0.677308,  2.021746,  0.029718/) )
       Magnetic_Form(114) = Magnetic_Form_Type("JTM3", &
                                              (/  0.028279,  2.291633,  0.265583,  7.776700,  0.675720,  2.018924,  0.029883/) )
       Magnetic_Form(115) = Magnetic_Form_Type("JYB3", &
                                           (/  0.092380,  2.046342,  0.258408,  7.471918,  0.609716,  1.913869,  0.038824/) )
       Magnetic_Form(116) = Magnetic_Form_Type("JO1 ", &
                                              (/  0.115285, 85.197300,  0.556229, 25.252200,  0.332476,  6.362070,-0.00460676/) )

       !---- <j2> Coefficients ----!
       Magnetic_j2(  1) = Magnetic_Form_Type("SC0 ",(/10.8172,54.327, 4.7353,14.847, 0.6071, 4.218, 0.0011/))
       Magnetic_j2(  2) = Magnetic_Form_Type("SC1 ",(/ 8.5021,34.285, 3.2116,10.994, 0.4244, 3.605, 0.0009/))
       Magnetic_j2(  3) = Magnetic_Form_Type("SC2 ",(/ 4.3683,28.654, 3.7231,10.823, 0.6074, 3.668, 0.0014/))
       Magnetic_j2(  4) = Magnetic_Form_Type("TI0 ",(/ 4.3583,36.056, 3.8230,11.133, 0.6855, 3.469, 0.0020/))
       Magnetic_j2(  5) = Magnetic_Form_Type("TI1 ",(/ 6.1567,27.275, 2.6833, 8.983, 0.4070, 3.052, 0.0011/))
       Magnetic_j2(  6) = Magnetic_Form_Type("TI2 ",(/ 4.3107,18.348, 2.0960, 6.797, 0.2984, 2.548, 0.0007/))
       Magnetic_j2(  7) = Magnetic_Form_Type("TI3 ",(/ 3.3717,14.444, 1.8258, 5.713, 0.2470, 2.265, 0.0005/))
       Magnetic_j2(  8) = Magnetic_Form_Type("V0  ",(/ 3.8099,21.347, 2.3295, 7.409, 0.4333, 2.632, 0.0015/))
       Magnetic_j2(  9) = Magnetic_Form_Type("V1  ",(/ 4.7474,23.323, 2.3609, 7.808, 0.4105, 2.706, 0.0014/))
       Magnetic_j2( 10) = Magnetic_Form_Type("V2  ",(/ 3.4386,16.530, 1.9638, 6.141, 0.2997, 2.267, 0.0009/))
       Magnetic_j2( 11) = Magnetic_Form_Type("V3  ",(/ 2.3005,14.682, 2.0364, 6.130, 0.4099, 2.382, 0.0014/))
       Magnetic_j2( 12) = Magnetic_Form_Type("V4  ",(/ 1.8377,12.267, 1.8247, 5.458, 0.3979, 2.248, 0.0012/))
       Magnetic_j2( 13) = Magnetic_Form_Type("CR0 ",(/ 3.4085,20.127, 2.1006, 6.802, 0.4266, 2.394, 0.0019/))
       Magnetic_j2( 14) = Magnetic_Form_Type("CR1 ",(/ 3.7768,20.346, 2.1028, 6.893, 0.4010, 2.411, 0.0017/))
       Magnetic_j2( 15) = Magnetic_Form_Type("CR2 ",(/ 2.6422,16.060, 1.9198, 6.253, 0.4446, 2.372, 0.0020/))
       Magnetic_j2( 16) = Magnetic_Form_Type("CR3 ",(/ 1.6262,15.066, 2.0618, 6.284, 0.5281, 2.368, 0.0023/))
       Magnetic_j2( 17) = Magnetic_Form_Type("CR4 ",(/ 1.0293,13.950, 1.9933, 6.059, 0.5974, 2.346, 0.0027/))
       Magnetic_j2( 18) = Magnetic_Form_Type("MN0 ",(/ 2.6681,16.060, 1.7561, 5.640, 0.3675, 2.049, 0.0017/))
       Magnetic_j2( 19) = Magnetic_Form_Type("MN1 ",(/ 3.2953,18.695, 1.8792, 6.240, 0.3927, 2.201, 0.0022/))
       Magnetic_j2( 20) = Magnetic_Form_Type("MN2 ",(/ 2.0515,15.556, 1.8841, 6.063, 0.4787, 2.232, 0.0027/))
       Magnetic_j2( 21) = Magnetic_Form_Type("MN3 ",(/ 1.2427,14.997, 1.9567, 6.118, 0.5732, 2.258, 0.0031/))
       Magnetic_j2( 22) = Magnetic_Form_Type("MN4 ",(/ 0.7879,13.886, 1.8717, 5.743, 0.5981, 2.182, 0.0034/))
       Magnetic_j2( 23) = Magnetic_Form_Type("FE0 ",(/ 1.9405,18.473, 1.9566, 6.323, 0.5166, 2.161, 0.0036/))
       Magnetic_j2( 24) = Magnetic_Form_Type("FE1 ",(/ 2.6290,18.660, 1.8704, 6.331, 0.4690, 2.163, 0.0031/))
       Magnetic_j2( 25) = Magnetic_Form_Type("FE2 ",(/ 1.6490,16.559, 1.9064, 6.133, 0.5206, 2.137, 0.0035/))
       Magnetic_j2( 26) = Magnetic_Form_Type("FE3 ",(/ 1.3602,11.998, 1.5188, 5.003, 0.4705, 1.991, 0.0038/))
       Magnetic_j2( 27) = Magnetic_Form_Type("FE4 ",(/ 1.5582, 8.275, 1.1863, 3.279, 0.1366, 1.107,-0.0022/))
       Magnetic_j2( 28) = Magnetic_Form_Type("CO0 ",(/ 1.9678,14.170, 1.4911, 4.948, 0.3844, 1.797, 0.0027/))
       Magnetic_j2( 29) = Magnetic_Form_Type("CO1 ",(/ 2.4097,16.161, 1.5780, 5.460, 0.4095, 1.914, 0.0031/))
       Magnetic_j2( 30) = Magnetic_Form_Type("CO2 ",(/ 1.9049,11.644, 1.3159, 4.357, 0.3146, 1.645, 0.0017/))
       Magnetic_j2( 31) = Magnetic_Form_Type("CO3 ",(/ 1.7058, 8.859, 1.1409, 3.309, 0.1474, 1.090,-0.0025/))
       Magnetic_j2( 32) = Magnetic_Form_Type("CO4 ",(/ 1.3110, 8.025, 1.1551, 3.179, 0.1608, 1.130,-0.0011/))
       Magnetic_j2( 33) = Magnetic_Form_Type("NI0 ",(/ 1.0302,12.252, 1.4669, 4.745, 0.4521, 1.744, 0.0036/))
       Magnetic_j2( 34) = Magnetic_Form_Type("NI1 ",(/ 2.1040,14.866, 1.4302, 5.071, 0.4031, 1.778, 0.0034/))
       Magnetic_j2( 35) = Magnetic_Form_Type("NI2 ",(/ 1.7080,11.016, 1.2147, 4.103, 0.3150, 1.533, 0.0018/))
       Magnetic_j2( 36) = Magnetic_Form_Type("NI3 ",(/ 1.1612, 7.700, 1.0027, 3.263, 0.2719, 1.378, 0.0025/))
       Magnetic_j2( 37) = Magnetic_Form_Type("NI4 ",(/ 1.1612, 7.700, 1.0027, 3.263, 0.2719, 1.378, 0.0025/))
       Magnetic_j2( 38) = Magnetic_Form_Type("CU0 ",(/ 1.9182,14.490, 1.3329, 4.730, 0.3842, 1.639, 0.0035/))
       Magnetic_j2( 39) = Magnetic_Form_Type("CU1 ",(/ 1.8814,13.433, 1.2809, 4.545, 0.3646, 1.602, 0.0033/))
       Magnetic_j2( 40) = Magnetic_Form_Type("CU2 ",(/ 1.5189,10.478, 1.1512, 3.813, 0.2918, 1.398, 0.0017/))
       Magnetic_j2( 41) = Magnetic_Form_Type("CU3 ",(/ 1.2797, 8.450, 1.0315, 3.280, 0.2401, 1.250, 0.0015/))
       Magnetic_j2( 42) = Magnetic_Form_Type("CU4 ",(/ 0.9568, 7.448, 0.9099, 3.396, 0.3729, 1.494, 0.0049/))
       Magnetic_j2( 43) = Magnetic_Form_Type("Y0  ",(/14.4084,44.658, 5.1045,14.904,-0.0535, 3.319, 0.0028/))
       Magnetic_j2( 44) = Magnetic_Form_Type("ZR0 ",(/10.1378,35.337, 4.7734,12.545,-0.0489, 2.672, 0.0036/))
       Magnetic_j2( 45) = Magnetic_Form_Type("ZR1 ",(/11.8722,34.920, 4.0502,12.127,-0.0632, 2.828, 0.0034/))
       Magnetic_j2( 46) = Magnetic_Form_Type("NB0 ",(/ 7.4796,33.179, 5.0884,11.571,-0.0281, 1.564, 0.0047/))
       Magnetic_j2( 47) = Magnetic_Form_Type("NB1 ",(/ 8.7735,33.285, 4.6556,11.605,-0.0268, 1.539, 0.0044/))
       Magnetic_j2( 48) = Magnetic_Form_Type("MO0 ",(/ 5.1180,23.422, 4.1809, 9.208,-0.0505, 1.743, 0.0053/))
       Magnetic_j2( 49) = Magnetic_Form_Type("MO1 ",(/ 7.2367,28.128, 4.0705, 9.923,-0.0317, 1.455, 0.0049/))
       Magnetic_j2( 50) = Magnetic_Form_Type("TC0 ",(/ 4.2441,21.397, 3.9439, 8.375,-0.0371, 1.187, 0.0066/))
       Magnetic_j2( 51) = Magnetic_Form_Type("TC1 ",(/ 6.4056,24.824, 3.5400, 8.611,-0.0366, 1.485, 0.0044/))
       Magnetic_j2( 52) = Magnetic_Form_Type("RU0 ",(/ 3.7445,18.613, 3.4749, 7.420,-0.0363, 1.007, 0.0073/))
       Magnetic_j2( 53) = Magnetic_Form_Type("RU1 ",(/ 5.2826,23.683, 3.5813, 8.152,-0.0257, 0.426, 0.0131/))
       Magnetic_j2( 54) = Magnetic_Form_Type("RH0 ",(/ 3.3651,17.344, 3.2121, 6.804,-0.0350, 0.503, 0.0146/))
       Magnetic_j2( 55) = Magnetic_Form_Type("RH1 ",(/ 4.0260,18.950, 3.1663, 7.000,-0.0296, 0.486, 0.0127/))
       Magnetic_j2( 56) = Magnetic_Form_Type("PD0 ",(/ 3.3105,14.726, 2.6332, 5.862,-0.0437, 1.130, 0.0053/))
       Magnetic_j2( 57) = Magnetic_Form_Type("PD1 ",(/ 4.2749,17.900, 2.7021, 6.354,-0.0258, 0.700, 0.0071/))
       Magnetic_j2( 58) = Magnetic_Form_Type("CE2 ",(/ 0.9809,18.063, 1.8413, 7.769, 0.9905, 2.845, 0.0120/))
       Magnetic_j2( 59) = Magnetic_Form_Type("ND2 ",(/ 1.4530,18.340, 1.6196, 7.285, 0.8752, 2.622, 0.0126/))
       Magnetic_j2( 60) = Magnetic_Form_Type("ND3 ",(/ 0.6751,18.342, 1.6272, 7.260, 0.9644, 2.602, 0.0150/))
       Magnetic_j2( 61) = Magnetic_Form_Type("SM2 ",(/ 1.0360,18.425, 1.4769, 7.032, 0.8810, 2.437, 0.0152/))
       Magnetic_j2( 62) = Magnetic_Form_Type("SM3 ",(/ 0.4707,18.430, 1.4261, 7.034, 0.9574, 2.439, 0.0182/))
       Magnetic_j2( 63) = Magnetic_Form_Type("EU2 ",(/ 0.8970,18.443, 1.3769, 7.005, 0.9060, 2.421, 0.0190/))
       Magnetic_j2( 64) = Magnetic_Form_Type("EU3 ",(/ 0.3985,18.451, 1.3307, 6.956, 0.9603, 2.378, 0.0197/))
       Magnetic_j2( 65) = Magnetic_Form_Type("GD2 ",(/ 0.7756,18.469, 1.3124, 6.899, 0.8956, 2.338, 0.0199/))
       Magnetic_j2( 66) = Magnetic_Form_Type("GD3 ",(/ 0.3347,18.476, 1.2465, 6.877, 0.9537, 2.318, 0.0217/))
       Magnetic_j2( 67) = Magnetic_Form_Type("TB2 ",(/ 0.6688,18.491, 1.2487, 6.822, 0.8888, 2.275, 0.0215/))
       Magnetic_j2( 68) = Magnetic_Form_Type("TB3 ",(/ 0.2892,18.497, 1.1678, 6.797, 0.9437, 2.257, 0.0232/))
       Magnetic_j2( 69) = Magnetic_Form_Type("DY2 ",(/ 0.5917,18.511, 1.1828, 6.747, 0.8801, 2.214, 0.0229/))
       Magnetic_j2( 70) = Magnetic_Form_Type("DY3 ",(/ 0.2523,18.517, 1.0914, 6.736, 0.9345, 2.208, 0.0250/))
       Magnetic_j2( 71) = Magnetic_Form_Type("HO2 ",(/ 0.5094,18.515, 1.1234, 6.706, 0.8727, 2.159, 0.0242/))
       Magnetic_j2( 72) = Magnetic_Form_Type("HO3 ",(/ 0.2188,18.516, 1.0240, 6.707, 0.9251, 2.161, 0.0268/))
       Magnetic_j2( 73) = Magnetic_Form_Type("ER2 ",(/ 0.4693,18.528, 1.0545, 6.649, 0.8679, 2.120, 0.0261/))
       Magnetic_j2( 74) = Magnetic_Form_Type("ER3 ",(/ 0.1710,18.534, 0.9879, 6.625, 0.9044, 2.100, 0.0278/))
       Magnetic_j2( 75) = Magnetic_Form_Type("TM2 ",(/ 0.4198,18.542, 0.9959, 6.600, 0.8593, 2.082, 0.0284/))
       Magnetic_j2( 76) = Magnetic_Form_Type("TM3 ",(/ 0.1760,18.542, 0.9105, 6.579, 0.8970, 2.062, 0.0294/))
       Magnetic_j2( 77) = Magnetic_Form_Type("YB2 ",(/ 0.3852,18.550, 0.9415, 6.551, 0.8492, 2.043, 0.0301/))
       Magnetic_j2( 78) = Magnetic_Form_Type("YB3 ",(/ 0.1570,18.555, 0.8484, 6.540, 0.8880, 2.037, 0.0318/))
       Magnetic_j2( 79) = Magnetic_Form_Type("U3  ",(/ 4.1582,16.534, 2.4675, 5.952,-0.0252, 0.765, 0.0057/))
       Magnetic_j2( 80) = Magnetic_Form_Type("U4  ",(/ 3.7449,13.894, 2.6453, 4.863,-0.5218, 3.192, 0.0009/))
       Magnetic_j2( 81) = Magnetic_Form_Type("U5  ",(/ 3.0724,12.546, 2.3076, 5.231,-0.0644, 1.474, 0.0035/))
       Magnetic_j2( 82) = Magnetic_Form_Type("NP3 ",(/ 3.7170,15.133, 2.3216, 5.503,-0.0275, 0.800, 0.0052/))
       Magnetic_j2( 83) = Magnetic_Form_Type("NP4 ",(/ 2.9203,14.646, 2.5979, 5.559,-0.0301, 0.367, 0.0141/))
       Magnetic_j2( 84) = Magnetic_Form_Type("NP5 ",(/ 2.3308,13.654, 2.7219, 5.494,-0.1357, 0.049, 0.1224/))
       Magnetic_j2( 85) = Magnetic_Form_Type("NP6 ",(/ 1.8245,13.180, 2.8508, 5.407,-0.1579, 0.044, 0.1438/))
       Magnetic_j2( 86) = Magnetic_Form_Type("PU3 ",(/ 2.0885,12.871, 2.5961, 5.190,-0.1465, 0.039, 0.1343/))
       Magnetic_j2( 87) = Magnetic_Form_Type("PU4 ",(/ 2.7244,12.926, 2.3387, 5.163,-0.1300, 0.046, 0.1177/))
       Magnetic_j2( 88) = Magnetic_Form_Type("PU5 ",(/ 2.1409,12.832, 2.5664, 5.152,-0.1338, 0.046, 0.1210/))
       Magnetic_j2( 89) = Magnetic_Form_Type("PU6 ",(/ 1.7262,12.324, 2.6652, 5.066,-0.1695, 0.041, 0.1550/))
       Magnetic_j2( 90) = Magnetic_Form_Type("AM2 ",(/ 3.5237,15.955, 2.2855, 5.195,-0.0142, 0.585, 0.0033/))
       Magnetic_j2( 91) = Magnetic_Form_Type("AM3 ",(/ 2.8622,14.733, 2.4099, 5.144,-0.1326, 0.031, 0.1233/))
       Magnetic_j2( 92) = Magnetic_Form_Type("AM4 ",(/ 2.4141,12.948, 2.3687, 4.945,-0.2490, 0.022, 0.2371/))
       Magnetic_j2( 93) = Magnetic_Form_Type("AM5 ",(/ 2.0109,12.053, 2.4155, 4.836,-0.2264, 0.027, 0.2128/))
       Magnetic_j2( 94) = Magnetic_Form_Type("AM6 ",(/ 1.6778,11.337, 2.4531, 4.725,-0.2043, 0.034, 0.1892/))
       Magnetic_j2( 95) = Magnetic_Form_Type("AM7 ",(/ 1.8845, 9.161, 2.0746, 4.042,-0.1318, 1.723, 0.0020/))

       !---- <j4> Coefficients ----!
       Magnetic_j4(  1) = Magnetic_Form_Type("SC0 ",(/ 1.3420,10.200, 0.3837, 3.079, 0.0468, 0.118,-0.0328/))
       Magnetic_j4(  2) = Magnetic_Form_Type("SC1 ",(/ 7.1167,15.487,-6.6671,18.269, 0.4900, 2.992, 0.0047/))
       Magnetic_j4(  3) = Magnetic_Form_Type("SC2 ",(/-1.6684,15.648, 1.7742, 9.062, 0.4075, 2.412, 0.0042/))
       Magnetic_j4(  4) = Magnetic_Form_Type("TI0 ",(/-2.1515,11.271, 2.5149, 8.859, 0.3555, 2.149, 0.0045/))
       Magnetic_j4(  5) = Magnetic_Form_Type("TI1 ",(/-1.0383,16.190, 1.4699, 8.924, 0.3631, 2.283, 0.0044/))
       Magnetic_j4(  6) = Magnetic_Form_Type("TI2 ",(/-1.3242,15.310, 1.2042, 7.899, 0.3976, 2.156, 0.0051/))
       Magnetic_j4(  7) = Magnetic_Form_Type("TI3 ",(/-1.1117,14.635, 0.7689, 6.927, 0.4385, 2.089, 0.0060/))
       Magnetic_j4(  8) = Magnetic_Form_Type("V0  ",(/-0.9633,15.273, 0.9274, 7.732, 0.3891, 2.053, 0.0063/))
       Magnetic_j4(  9) = Magnetic_Form_Type("V1  ",(/-0.9606,15.545, 1.1278, 8.118, 0.3653, 2.097, 0.0056/))
       Magnetic_j4( 10) = Magnetic_Form_Type("V2  ",(/-1.1729,14.973, 0.9092, 7.613, 0.4105, 2.039, 0.0067/))
       Magnetic_j4( 11) = Magnetic_Form_Type("V3  ",(/-0.9417,14.205, 0.5284, 6.607, 0.4411, 1.967, 0.0076/))
       Magnetic_j4( 12) = Magnetic_Form_Type("V4  ",(/-0.7654,13.097, 0.3071, 5.674, 0.4476, 1.871, 0.0081/))
       Magnetic_j4( 13) = Magnetic_Form_Type("CR0 ",(/-0.6670,19.613, 0.5342, 6.478, 0.3641, 1.905, 0.0073/))
       Magnetic_j4( 14) = Magnetic_Form_Type("CR1 ",(/-0.8309,18.043, 0.7252, 7.531, 0.3828, 2.003, 0.0073/))
       Magnetic_j4( 15) = Magnetic_Form_Type("CR2 ",(/-0.8930,15.664, 0.5590, 7.033, 0.4093, 1.924, 0.0081/))
       Magnetic_j4( 16) = Magnetic_Form_Type("CR3 ",(/-0.7327,14.073, 0.3268, 5.674, 0.4114, 1.810, 0.0085/))
       Magnetic_j4( 17) = Magnetic_Form_Type("CR4 ",(/-0.6748,12.946, 0.1805, 6.753, 0.4526, 1.800, 0.0098/))
       Magnetic_j4( 18) = Magnetic_Form_Type("MN0 ",(/-0.5452,15.471, 0.4406, 4.902, 0.2884, 1.543, 0.0059/))
       Magnetic_j4( 19) = Magnetic_Form_Type("MN1 ",(/-0.7947,17.867, 0.6078, 7.704, 0.3798, 1.905, 0.0087/))
       Magnetic_j4( 20) = Magnetic_Form_Type("MN2 ",(/-0.7416,15.255, 0.3831, 6.469, 0.3935, 1.800, 0.0093/))
       Magnetic_j4( 21) = Magnetic_Form_Type("MN3 ",(/-0.6603,13.607, 0.2322, 6.218, 0.4104, 1.740, 0.0101/))
       Magnetic_j4( 22) = Magnetic_Form_Type("MN4 ",(/-0.5127,13.461, 0.0313, 7.763, 0.4282, 1.701, 0.0113/))
       Magnetic_j4( 23) = Magnetic_Form_Type("FE0 ",(/-0.5029,19.677, 0.2999, 3.776, 0.2576, 1.424, 0.0071/))
       Magnetic_j4( 24) = Magnetic_Form_Type("FE1 ",(/-0.5109,19.250, 0.3896, 4.891, 0.2810, 1.526, 0.0069/))
       Magnetic_j4( 25) = Magnetic_Form_Type("FE2 ",(/-0.5401,17.227, 0.2865, 3.742, 0.2658, 1.424, 0.0076/))
       Magnetic_j4( 26) = Magnetic_Form_Type("FE3 ",(/-0.5507,11.493, 0.2153, 4.906, 0.3468, 1.523, 0.0095/))
       Magnetic_j4( 27) = Magnetic_Form_Type("FE4 ",(/-0.5352, 9.507, 0.1783, 5.175, 0.3584, 1.469, 0.0097/))
       Magnetic_j4( 28) = Magnetic_Form_Type("CO0 ",(/-0.4221,14.195, 0.2900, 3.979, 0.2469, 1.286, 0.0063/))
       Magnetic_j4( 29) = Magnetic_Form_Type("CO1 ",(/-0.4115,14.561, 0.3580, 4.717, 0.2644, 1.418, 0.0074/))
       Magnetic_j4( 30) = Magnetic_Form_Type("CO2 ",(/-0.4759,14.046, 0.2747, 3.731, 0.2458, 1.250, 0.0057/))
       Magnetic_j4( 31) = Magnetic_Form_Type("CO3 ",(/-0.4466,13.391, 0.1419, 3.011, 0.2773, 1.335, 0.0093/))
       Magnetic_j4( 32) = Magnetic_Form_Type("CO4 ",(/-0.4091,13.194,-0.0194, 3.417, 0.3534, 1.421, 0.0112/))
       Magnetic_j4( 33) = Magnetic_Form_Type("NI0 ",(/-0.4428,14.485, 0.0870, 3.234, 0.2932, 1.331, 0.0096/))
       Magnetic_j4( 34) = Magnetic_Form_Type("NI1 ",(/-0.3836,13.425, 0.3116, 4.462, 0.2471, 1.309, 0.0079/))
       Magnetic_j4( 35) = Magnetic_Form_Type("NI2 ",(/-0.3803,10.403, 0.2838, 3.378, 0.2108, 1.104, 0.0050/))
       Magnetic_j4( 36) = Magnetic_Form_Type("NI3 ",(/-0.3715, 8.952, 0.1211, 2.940, 0.2526, 1.105, 0.0061/))
       Magnetic_j4( 37) = Magnetic_Form_Type("NI4 ",(/-0.3509, 8.157, 0.2220, 2.106, 0.1567, 0.925, 0.0065/))
       Magnetic_j4( 38) = Magnetic_Form_Type("CU0 ",(/-0.3204,15.132, 0.2335, 4.021, 0.2312, 1.196, 0.0068/))
       Magnetic_j4( 39) = Magnetic_Form_Type("CU1 ",(/-0.3572,15.125, 0.2336, 3.966, 0.2315, 1.197, 0.0070/))
       Magnetic_j4( 40) = Magnetic_Form_Type("CU2 ",(/-0.3914,14.740, 0.1275, 3.384, 0.2548, 1.255, 0.0103/))
       Magnetic_j4( 41) = Magnetic_Form_Type("CU3 ",(/-0.3671,14.082,-0.0078, 3.315, 0.3154, 1.377, 0.0132/))
       Magnetic_j4( 42) = Magnetic_Form_Type("CU4 ",(/-0.2915,14.124,-0.1065, 4.201, 0.3247, 1.352, 0.0148/))
       Magnetic_j4( 43) = Magnetic_Form_Type("Y0  ",(/-8.0767,32.201, 7.9197,25.156, 1.4067, 6.827,-0.0001/))
       Magnetic_j4( 44) = Magnetic_Form_Type("ZR0 ",(/-5.2697,32.868, 4.1930,24.183, 1.5202, 6.048,-0.0002/))
       Magnetic_j4( 45) = Magnetic_Form_Type("ZR1 ",(/-5.6384,33.607, 4.6729,22.338, 1.3258, 5.924,-0.0003/))
       Magnetic_j4( 46) = Magnetic_Form_Type("NB0 ",(/-3.1377,25.595, 2.3411,16.569, 1.2304, 4.990,-0.0005/))
       Magnetic_j4( 47) = Magnetic_Form_Type("NB1 ",(/-3.3598,25.820, 2.8297,16.427, 1.1203, 4.982,-0.0005/))
       Magnetic_j4( 48) = Magnetic_Form_Type("MO0 ",(/-2.8860,20.572, 1.8130,14.628, 1.1899, 4.264,-0.0008/))
       Magnetic_j4( 49) = Magnetic_Form_Type("MO1 ",(/-3.2618,25.486, 2.3596,16.462, 1.1164, 4.491,-0.0007/))
       Magnetic_j4( 50) = Magnetic_Form_Type("TC0 ",(/-2.7975,20.159, 1.6520,16.261, 1.1726, 3.943,-0.0008/))
       Magnetic_j4( 51) = Magnetic_Form_Type("TC1 ",(/-2.0470,19.683, 1.6306,11.592, 0.8698, 3.769,-0.0010/))
       Magnetic_j4( 52) = Magnetic_Form_Type("RU0 ",(/-1.5042,17.949, 0.6027, 9.961, 0.9700, 3.393,-0.0010/))
       Magnetic_j4( 53) = Magnetic_Form_Type("RU1 ",(/-1.6278,18.506, 1.1828,10.189, 0.8138, 3.418,-0.0009/))
       Magnetic_j4( 54) = Magnetic_Form_Type("RH0 ",(/-1.3492,17.577, 0.4527,10.507, 0.9285, 3.155,-0.0009/))
       Magnetic_j4( 55) = Magnetic_Form_Type("RH1 ",(/-1.4673,17.957, 0.7381, 9.944, 0.8485, 3.126,-0.0012/))
       Magnetic_j4( 56) = Magnetic_Form_Type("PD0 ",(/-1.1955,17.628, 0.3183,11.309, 0.8696, 2.909,-0.0006/))
       Magnetic_j4( 57) = Magnetic_Form_Type("PD1 ",(/-1.4098,17.765, 0.7927, 9.999, 0.7710, 2.930,-0.0006/))
       Magnetic_j4( 58) = Magnetic_Form_Type("CE2 ",(/-0.6468,10.533, 0.4052, 5.624, 0.3412, 1.535, 0.0080/))
       Magnetic_j4( 59) = Magnetic_Form_Type("ND2 ",(/-0.5416,12.204, 0.3571, 6.169, 0.3154, 1.485, 0.0098/))
       Magnetic_j4( 60) = Magnetic_Form_Type("ND3 ",(/-0.4053,14.014, 0.0329, 7.005, 0.3759, 1.707, 0.0209/))
       Magnetic_j4( 61) = Magnetic_Form_Type("SM2 ",(/-0.4150,14.057, 0.1368, 7.032, 0.3272, 1.582, 0.0192/))
       Magnetic_j4( 62) = Magnetic_Form_Type("SM3 ",(/-0.4288,10.052, 0.1782, 5.019, 0.2833, 1.236, 0.0088/))
       Magnetic_j4( 63) = Magnetic_Form_Type("EU2 ",(/-0.4145,10.193, 0.2447, 5.164, 0.2661, 1.205, 0.0065/))
       Magnetic_j4( 64) = Magnetic_Form_Type("EU3 ",(/-0.4095,10.211, 0.1485, 5.175, 0.2720, 1.237, 0.0131/))
       Magnetic_j4( 65) = Magnetic_Form_Type("GD2 ",(/-0.3824,10.344, 0.1955, 5.306, 0.2622, 1.203, 0.0097/))
       Magnetic_j4( 66) = Magnetic_Form_Type("GD3 ",(/-0.3621,10.353, 0.1016, 5.310, 0.2649, 1.219, 0.0147/))
       Magnetic_j4( 67) = Magnetic_Form_Type("TB2 ",(/-0.3443,10.469, 0.1481, 5.416, 0.2575, 1.182, 0.0104/))
       Magnetic_j4( 68) = Magnetic_Form_Type("TB3 ",(/-0.3228,10.476, 0.0638, 5.419, 0.2566, 1.196, 0.0159/))
       Magnetic_j4( 69) = Magnetic_Form_Type("DY2 ",(/-0.3206,12.071, 0.0904, 8.026, 0.2616, 1.230, 0.0143/))
       Magnetic_j4( 70) = Magnetic_Form_Type("DY3 ",(/-0.2829, 9.525, 0.0565, 4.429, 0.2437, 1.066, 0.0092/))
       Magnetic_j4( 71) = Magnetic_Form_Type("HO2 ",(/-0.2976, 9.719, 0.1224, 4.635, 0.2279, 1.005, 0.0063/))
       Magnetic_j4( 72) = Magnetic_Form_Type("HO3 ",(/-0.2717, 9.731, 0.0474, 4.638, 0.2292, 1.047, 0.0124/))
       Magnetic_j4( 73) = Magnetic_Form_Type("ER2 ",(/-0.2975, 9.829, 0.1189, 4.741, 0.2116, 1.004, 0.0117/))
       Magnetic_j4( 74) = Magnetic_Form_Type("ER3 ",(/-0.2568, 9.834, 0.0356, 4.741, 0.2172, 1.028, 0.0148/))
       Magnetic_j4( 75) = Magnetic_Form_Type("TM2 ",(/-0.2677, 9.888, 0.0925, 4.784, 0.2056, 0.990, 0.0124/))
       Magnetic_j4( 76) = Magnetic_Form_Type("TM3 ",(/-0.2292, 9.895, 0.0124, 4.785, 0.2108, 1.007, 0.0151/))
       Magnetic_j4( 77) = Magnetic_Form_Type("YB2 ",(/-0.2393, 9.947, 0.0663, 4.823, 0.2009, 0.965, 0.0122/))
       Magnetic_j4( 78) = Magnetic_Form_Type("YB3 ",(/-0.2121, 8.197, 0.0325, 3.153, 0.1975, 0.884, 0.0093/))
       Magnetic_j4( 79) = Magnetic_Form_Type("U3  ",(/-0.9859,16.601, 0.6116, 6.515, 0.6020, 2.597,-0.0010/))
       Magnetic_j4( 80) = Magnetic_Form_Type("U4  ",(/-1.0540,16.605, 0.4339, 6.512, 0.6746, 2.599,-0.0011/))
       Magnetic_j4( 81) = Magnetic_Form_Type("U5  ",(/-0.9588,16.485, 0.1576, 6.440, 0.7785, 2.640,-0.0010/))
       Magnetic_j4( 82) = Magnetic_Form_Type("NP3 ",(/-0.9029,16.586, 0.4006, 6.470, 0.6545, 2.563,-0.0004/))
       Magnetic_j4( 83) = Magnetic_Form_Type("NP4 ",(/-0.9887,12.441, 0.5918, 5.294, 0.5306, 2.263,-0.0021/))
       Magnetic_j4( 84) = Magnetic_Form_Type("NP5 ",(/-0.8146,16.581,-0.0055, 6.475, 0.7956, 2.562,-0.0004/))
       Magnetic_j4( 85) = Magnetic_Form_Type("NP6 ",(/-0.6738,16.553,-0.2297, 6.505, 0.8513, 2.553,-0.0003/))
       Magnetic_j4( 86) = Magnetic_Form_Type("PU3 ",(/-0.7014,16.369,-0.1162, 6.697, 0.7778, 2.450, 0.0000/))
       Magnetic_j4( 87) = Magnetic_Form_Type("PU4 ",(/-0.9160,12.203, 0.4891, 5.127, 0.5290, 2.149,-0.0022/))
       Magnetic_j4( 88) = Magnetic_Form_Type("PU5 ",(/-0.7035,16.360,-0.0979, 6.706, 0.7726, 2.447, 0.0000/))
       Magnetic_j4( 89) = Magnetic_Form_Type("PU6 ",(/-0.5560,16.322,-0.3046, 6.768, 0.8146, 2.426, 0.0001/))
       Magnetic_j4( 90) = Magnetic_Form_Type("AM2 ",(/-0.7433,16.416, 0.3481, 6.788, 0.6014, 2.346, 0.0000/))
       Magnetic_j4( 91) = Magnetic_Form_Type("AM3 ",(/-0.8092,12.854, 0.4161, 5.459, 0.5476, 2.172,-0.0011/))
       Magnetic_j4( 92) = Magnetic_Form_Type("AM4 ",(/-0.8548,12.226, 0.3037, 5.909, 0.6173, 2.188,-0.0016/))
       Magnetic_j4( 93) = Magnetic_Form_Type("AM5 ",(/-0.6538,15.462,-0.0948, 5.997, 0.7295, 2.297, 0.0000/))
       Magnetic_j4( 94) = Magnetic_Form_Type("AM6 ",(/-0.5390,15.449,-0.2689, 6.017, 0.7711, 2.297, 0.0002/))
       Magnetic_j4( 95) = Magnetic_Form_Type("AM7 ",(/-0.4688,12.019,-0.2692, 7.042, 0.7297, 2.164,-0.0011/))

       !---- <j6> Coefficients ----!
       Magnetic_j6(  1) = Magnetic_Form_Type("CE2 ",(/-0.1212, 7.994,-0.0639, 4.024, 0.1519, 1.096, 0.0078/))
       Magnetic_j6(  2) = Magnetic_Form_Type("ND2 ",(/-0.1600, 8.009, 0.0272, 4.028, 0.1104, 1.068, 0.0139/))
       Magnetic_j6(  3) = Magnetic_Form_Type("ND3 ",(/-0.0416, 8.014,-0.1261, 4.040, 0.1400, 1.087, 0.0102/))
       Magnetic_j6(  4) = Magnetic_Form_Type("SM2 ",(/-0.1428, 6.041, 0.0723, 2.033, 0.0550, 0.513, 0.0081/))
       Magnetic_j6(  5) = Magnetic_Form_Type("SM3 ",(/-0.0944, 6.030,-0.0498, 2.074, 0.1372, 0.645,-0.0132/))
       Magnetic_j6(  6) = Magnetic_Form_Type("EU2 ",(/-0.1252, 6.049, 0.0507, 2.085, 0.0572, 0.646, 0.0132/))
       Magnetic_j6(  7) = Magnetic_Form_Type("EU3 ",(/-0.0817, 6.039,-0.0596, 2.120, 0.1243, 0.764,-0.0001/))
       Magnetic_j6(  8) = Magnetic_Form_Type("GD2 ",(/-0.1351, 5.030, 0.0828, 2.025, 0.0315, 0.503, 0.0187/))
       Magnetic_j6(  9) = Magnetic_Form_Type("GD3 ",(/-0.0662, 6.031,-0.0850, 2.154, 0.1323, 0.891, 0.0048/))
       Magnetic_j6( 10) = Magnetic_Form_Type("TB2 ",(/-0.0758, 6.032,-0.0540, 2.158, 0.1199, 0.890, 0.0051/))
       Magnetic_j6( 11) = Magnetic_Form_Type("TB3 ",(/-0.0559, 6.031,-0.1020, 2.237, 0.1264, 1.107, 0.0167/))
       Magnetic_j6( 12) = Magnetic_Form_Type("DY2 ",(/-0.0568, 6.032,-0.1003, 2.240, 0.1401, 1.106, 0.0109/))
       Magnetic_j6( 13) = Magnetic_Form_Type("DY3 ",(/-0.0423, 6.038,-0.1248, 2.244, 0.1359, 1.200, 0.0188/))
       Magnetic_j6( 14) = Magnetic_Form_Type("HO2 ",(/-0.0725, 6.045,-0.0318, 2.243, 0.0738, 1.202, 0.0252/))
       Magnetic_j6( 15) = Magnetic_Form_Type("HO3 ",(/-0.0289, 6.050,-0.1545, 2.230, 0.1550, 1.260, 0.0177/))
       Magnetic_j6( 16) = Magnetic_Form_Type("ER2 ",(/-0.0648, 6.056,-0.0515, 2.230, 0.0825, 1.264, 0.0250/))
       Magnetic_j6( 17) = Magnetic_Form_Type("ER3 ",(/-0.0110, 6.061,-0.1954, 2.224, 0.1818, 1.296, 0.0149/))
       Magnetic_j6( 18) = Magnetic_Form_Type("TM2 ",(/-0.0842, 4.070, 0.0807, 0.849,-0.2087, 0.039, 0.2095/))
       Magnetic_j6( 19) = Magnetic_Form_Type("TM3 ",(/-0.0727, 4.073, 0.0243, 0.689, 3.9459, 0.002,-3.9076/))
       Magnetic_j6( 20) = Magnetic_Form_Type("YB2 ",(/-0.0739, 5.031, 0.0140, 2.030, 0.0351, 0.508, 0.0174/))
       Magnetic_j6( 21) = Magnetic_Form_Type("YB3 ",(/-0.0345, 5.007,-0.0677, 2.020, 0.0985, 0.549,-0.0076/))
       Magnetic_j6( 22) = Magnetic_Form_Type("U3  ",(/-0.3797, 9.953, 0.0459, 5.038, 0.2748, 1.607, 0.0016/))
       Magnetic_j6( 23) = Magnetic_Form_Type("U4  ",(/-0.1793,11.896,-0.2269, 5.428, 0.3291, 1.701, 0.0030/))
       Magnetic_j6( 24) = Magnetic_Form_Type("U5  ",(/-0.0399,11.891,-0.3458, 5.580, 0.3340, 1.645, 0.0029/))
       Magnetic_j6( 25) = Magnetic_Form_Type("NP3 ",(/-0.2427,11.844,-0.1129, 5.377, 0.2848, 1.568, 0.0022/))
       Magnetic_j6( 26) = Magnetic_Form_Type("NP4 ",(/-0.2436, 9.599,-0.1317, 4.101, 0.3029, 1.545, 0.0019/))
       Magnetic_j6( 27) = Magnetic_Form_Type("NP5 ",(/-0.1157, 9.565,-0.2654, 4.260, 0.3298, 1.549, 0.0025/))
       Magnetic_j6( 28) = Magnetic_Form_Type("NP6 ",(/-0.0128, 9.569,-0.3611, 4.304, 0.3419, 1.541, 0.0032/))
       Magnetic_j6( 29) = Magnetic_Form_Type("PU3 ",(/-0.0364, 9.572,-0.3181, 4.342, 0.3210, 1.523, 0.0041/))
       Magnetic_j6( 30) = Magnetic_Form_Type("PU4 ",(/-0.2394, 7.837,-0.0785, 4.024, 0.2643, 1.378, 0.0012/))
       Magnetic_j6( 31) = Magnetic_Form_Type("PU5 ",(/-0.1090, 7.819,-0.2243, 4.100, 0.2947, 1.404, 0.0015/))
       Magnetic_j6( 32) = Magnetic_Form_Type("PU6 ",(/-0.0001, 7.820,-0.3354, 4.144, 0.3097, 1.403, 0.0020/))
       Magnetic_j6( 33) = Magnetic_Form_Type("AM2 ",(/-0.3176, 7.864, 0.0771, 4.161, 0.2194, 1.339, 0.0018/))
       Magnetic_j6( 34) = Magnetic_Form_Type("AM3 ",(/-0.3159, 6.982, 0.0682, 3.995, 0.2141, 1.188,-0.0015/))
       Magnetic_j6( 35) = Magnetic_Form_Type("AM4 ",(/-0.1787, 7.880,-0.1274, 4.090, 0.2565, 1.315, 0.0017/))
       Magnetic_j6( 36) = Magnetic_Form_Type("AM5 ",(/-0.0927, 6.073,-0.2227, 3.784, 0.2916, 1.372, 0.0026/))
       Magnetic_j6( 37) = Magnetic_Form_Type("AM6 ",(/ 0.0152, 6.079,-0.3549, 3.861, 0.3125, 1.403, 0.0036/))
       Magnetic_j6( 38) = Magnetic_Form_Type("AM7 ",(/ 0.1292, 6.082,-0.4689, 3.879, 0.3234, 1.393, 0.0042/))
"""
