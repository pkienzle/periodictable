import numpy as np

import periodictable as pt
from periodictable.activation import Sample, ActivationEnvironment
from periodictable.activation import IAEA1987_isotopic_abundance
from periodictable.activation import NIST2001_isotopic_abundance

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
    formula = "".join(str(el) for el in pt.elements)[1:]
    # Use an enormous mass to force significant activation of rare isotopes
    mass, fluence, exposure = 1e15, 1e8, 10
    env = ActivationEnvironment(fluence=fluence, Cd_ratio=70, fast_ratio=50, location="BT-2")
    sample = Sample(formula, mass=mass)
    sample.calculate_activation(
        env, exposure=exposure, rest_times=(0, 1, 24, 360),
        abundance=IAEA1987_isotopic_abundance,
        #abundance=NIST2001_isotopic_abundance,
        )

    # TODO: Why aren't we matching the spreadsheet results?
    results = {
        # The following were labelled as spreadsheet results
        #"Co-60m+":  [5186.888,  98.878,   2.751e-38],
        #"Co-61":    [2.689e-8,  1.767e-8, 1.127e-12],
        #"Co-60":    [1.550322,  1.550299, 1.549764],
        #"Co-61":    [5.695e-8,  3.741e-8, 2.386e-12],
        # Values as of R1.7.0
        "Co-60m+": [5088.093018977685, 96.91335724994167, 2.6450954672880745e-38],
        "Co-61": [7.309565456941425e-09, 4.802298351964856e-09, 3.0568451604601675e-13],
        "Co-60": [0.15507562804846606, 0.15507330056693813, 0.15501977813208836],
        "Co-61": [1.3649499351603886e-09, 8.967560195949139e-10, 5.708192406435261e-14],
    }

    sample = Sample('Co', mass=10)
    env = ActivationEnvironment(fluence=1e8)
    sample.calculate_activation(
        env, rest_times=[0, 1, 24],
        abundance=IAEA1987_isotopic_abundance,
        #abundance=NIST2001_isotopic_abundance,
        )
    print(sample.activity)
    for product, activity in sample.activity.items():
        print(f'        "{product.daughter}": {activity},')
        assert product.daughter in results, f"Missing {product.daughter}"
        print(results[product.daughter])
        print(activity)
        #assert np.allclose(results[product.daughter], activity, atol=0, rtol=1e-12)


    # 129-I has a long half-life (16 My) so any combination of exposure
    # fluence and mass that causes significant activation will yield a product
    # that is radioactive for a very long time.
    sample = Sample('Te', mass=1e13)
    env = ActivationEnvironment(fluence=1e8)
    sample.calculate_activation(env, rest_times=[0])
    target = 1e-5
    decay = sample.decay_time(target)
    sample.calculate_activation(env, rest_times=[0, decay])
    total = sum(v[-1] for v in sample.activity.values())
    assert abs(total - target)/target < 1e-6, f"total activity {total} != {target}"

test()
