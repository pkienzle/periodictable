import numpy as np

import periodictable as pt
from periodictable.activation import Sample, ActivationEnvironment
from periodictable.activation import IAEA1987_isotopic_abundance, table_abundance

def test():
    # This is not a very complete test of the activation calculator.
    # Mostly just a smoke test to see that things run and produce the
    # same answers as before.  The target values have been checked
    # against the NCNR internal activation calculator spreadsheet.  The
    # values herein differ slightly from the spreadsheet since we are using
    # different tables for materials and different precision on our
    # constants. There has also been some small corrections to the
    # formulas over time.

    def _get_Au_activity(fluence=1e5):
        sample = Sample('Au', mass=1)
        env = ActivationEnvironment(fluence=fluence)
        sample.calculate_activation(env, rest_times=[0])
        for product, activity in sample.activity.items():
            if str(product.daughter) ==  'Au-198':
                return activity[0]
        else:
            raise RuntimeError("missing activity from Au-198")

    # Activity scales linearly with fluence and mass
    act1 = _get_Au_activity(fluence=1e5)
    act2 = _get_Au_activity(fluence=1e8)
    #print(f"Au: {act1} x 1000 = {act2}")
    assert (act2 - 1000*act1) < 1e-8

    # Smoke test: does every element run to completion?
    formula = "".join(str(el) for el in pt.elements)
    # Use an enormous mass to force significant activation of rare isotopes
    mass, fluence, exposure = 1e15, 1e8, 10
    env = ActivationEnvironment(fluence=fluence, Cd_ratio=70, fast_ratio=50, location="BT-2")
    sample = Sample(formula, mass=mass)
    sample.calculate_activation(
        env, exposure=exposure, rest_times=(0, 1, 24),
        abundance=IAEA1987_isotopic_abundance,
        #abundance=table_abundance,
        )
    #sample.show_table(cutoff=0)

    sample = Sample('Co', mass=10)
    env = ActivationEnvironment(fluence=1e8)
    sample.calculate_activation(
        env, rest_times=[0, 1, 24],
        abundance=IAEA1987_isotopic_abundance,
        #abundance=table_abundance,
        )
    #sample.show_table(cutoff=0)
    # ACT_2N_X.xls for 10g Co at Io=1e8 for t = [0, 1, 24]
    # There are duplicate entries for Co-61 because there are different
    # activation paths for 59Co + n -> 61Co.
    # Results are limited to 3 digits because formulas use ln(2) = 0.693. There
    # is also precision loss on the calculation of differences in exponentials,
    # with exp(a) - exp(b) converted to exp(b)(exp(a-b) - 1) = exp(b)expm1(a-b)
    # in activation.py.
    # results = {
    #     "C0-60m+":  [5.08800989419543E+03,	9.69933141298983E+01,	2.69898447909379E-38],
    #     "Co-61":    [7.30937687202485E-09,	4.80260282859366E-09,	3.06331725449002E-13],
    #     "Co-60":    [1.55042700053951E-01,	1.55040373560730E-01,	1.54986873850802E-01],
    #     "Co-61":    [1.36469792582999E-09,	8.96670432174800E-10,	5.71936948464344E-14],
    # }
    # Results from R2.0.0-pre
    results = {
        "Co-60m+":	[5.08746786932552e+03,	9.69014499692868e+01,	2.64477047705988e-38],
        "Co-61":	[7.30866736559643e-09,	4.80170831653643e-09,	3.05646958051663e-13],
        "Co-60":	[1.55056574647797e-01,	1.55054247452235e-01,	1.55000731593431e-01],
        "Co-61":	[1.36478223029061e-09,	8.96645839472098e-10,	5.70749106813738e-14],
    }
    #print(list(sample.activity.keys()))
    #print(dir(list(sample.activity.keys())[0]))
    for product, activity in sample.activity.items():
        #print(product)
        #print(dir(product))
        # Uncomment to show new table values
        activity_str = ",\t".join(f"{Ia:.14E}" for Ia in activity)
        #print(f'        "{product.daughter}":\t[{activity_str}],')
        assert product.daughter in results, f"Missing {product.daughter}"
        # TODO: include duplicate decay paths in test, or identical daughters
        # Test that results haven't changed since last update
        if product.daughter != "Co-61":
            assert np.allclose(results[product.daughter], activity, atol=0, rtol=1e-12)

    # 129-I has a long half-life (16 My) so any combination of exposure
    # fluence and mass that causes significant activation will yield a product
    # that is radioactive for a very long time.
    sample = Sample('Te', mass=1e13)
    env = ActivationEnvironment(fluence=1e8)
    sample.calculate_activation(env, rest_times=[1,10,100])
    #sample.show_table(cutoff=0)
    target = 1e-5
    t_decay = sample.decay_time(target)
    #print(f"{t_decay=}")
    sample.calculate_activation(env, rest_times=[t_decay])
    total = sum(v[-1] for v in sample.activity.values())
    assert abs(total - target) < 1e-10, f"total activity {total} != {target}"

    # Al and Si daughters have short half-lives
    sample = Sample('AlSi', mass=1e13)
    env = ActivationEnvironment(fluence=1e8)
    sample.calculate_activation(env, rest_times=[100,200])
    #sample.show_table(cutoff=0)
    target = 1e-5
    t_decay = sample.decay_time(target)
    #print(f"{t_decay=}")
    sample.calculate_activation(env, rest_times=[t_decay])
    total = sum(v[-1] for v in sample.activity.values())
    assert abs(total - target) < 1e-10, f"total activity {total} != {target}"


test()
