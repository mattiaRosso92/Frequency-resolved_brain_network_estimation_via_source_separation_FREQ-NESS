"""Seeded calibration of the production participant-contrast statistic.

Run with PYTHONPATH=FREQNESS_Toolbox/python python this_file.py.
Defaults: 2000 independent data sets per scenario, n=14, 499 Monte Carlo
draws, a 2x4 design, four adjacent voxels and two frequencies. These are
prepared contrast fields: symmetry is assumed AFTER map preparation, not
asserted for arbitrary raw patterns followed by nonlinear normalization.
The skewed scenario deliberately violates that assumption. No data files
are read or written; JSON summaries are printed. This is a calibration
experiment, not a proof of validity for other distributions.

For the larger exact covariance-stress replication, append:
--exact --scenario nonspherical --nsets 8000 --seed 9304
"""
import argparse
import json
import time

import numpy as np
from scipy import stats

from freqness._repeated_spatial import RepeatedSpatialModel, repeated_spatial_statistic


def maximum_cluster(statistic, p_values, nsets, threshold):
    statistic = statistic.reshape(nsets, 4, 2)
    selected = p_values.reshape(nsets, 4, 2) <= threshold
    running = np.zeros((nsets, 2))
    maximum = np.zeros((nsets, 2))
    # Four-voxel line adjacency; frequency planes are separate maps.
    for voxi in range(4):
        running = np.where(selected[:, voxi], running + statistic[:, voxi], 0)
        maximum = np.maximum(maximum, running)
    return maximum.max(axis=1)


def simulate(scenario, nsets, draws, seed, exact=False):
    started = time.perf_counter()
    rng = np.random.default_rng(seed)
    # participants x seven orthonormal effect coordinates x data sets x maps
    shape = (14, 7, nsets, 8)
    values = rng.normal(size=shape)
    if scenario != "independent":
        values = .6*values + .8*rng.normal(size=(14, 7, nsets, 1))
        covariance = .5 ** np.abs(np.arange(7)[:, None] - np.arange(7))
        values = np.einsum("qr,srdm->sqdm", np.linalg.cholesky(covariance), values)
    if scenario == "random_drug_slopes":
        values[:, 0] *= 8
    elif scenario == "nonspherical":
        values *= np.array([1, 1, 3, 9, 1, 3, 9])[None, :, None, None]
    elif scenario == "heavy_tails":
        values /= np.sqrt(rng.chisquare(3, size=(14, 1, nsets, 1))/3)
    elif scenario == "partial_null":
        values[:, 0] += 2.5  # Drug is nonnull; Condition and interaction are null.
    elif scenario == "skewed_stress":
        values = (rng.exponential(size=(14, 7, nsets, 1))-1
                  + .05*rng.normal(size=shape))
    models = []
    for indices in (slice(0, 1), slice(1, 4), slice(4, 7)):
        scores = values[:, indices].reshape(14, -1, nsets*8)
        models.append(RepeatedSpatialModel(
            scores, np.einsum("sqm,srm->qrm", scores, scores), np.empty((0, 0))))
    del values
    null_effects = (1, 2) if scenario == "partial_null" else (0, 1, 2)
    observed = np.zeros(nsets)
    for effecti in null_effects:
        result = repeated_spatial_statistic(models[effecti])
        observed = np.maximum(observed, maximum_cluster(
            result["statistic"], result["p_uncorrected"], nsets, .05))
    exceedances = np.zeros(nsets, dtype=int)
    if exact:
        # F and GG are invariant to reversing ALL signs; evaluate one of
        # each +/- pair. This is the full 2^14 distribution with duplicates
        # removed, including the observed configuration, not Monte Carlo.
        draws = 2**13
    for drawi in range(draws):
        # Shared signs across the full family, as in the public function.
        signs = (np.r_[1, 2*((drawi >> np.arange(13)) & 1)-1]
                 if exact else rng.choice((-1, 1), size=14))
        maximum = np.zeros(nsets)
        for model in models:
            result = repeated_spatial_statistic(model, signs)
            maximum = np.maximum(maximum, maximum_cluster(
                result["statistic"], result["p_uncorrected"], nsets, .05))
        exceedances += maximum >= observed
    p_values = exceedances/draws if exact else (1+exceedances)/(draws+1)
    rejected = np.count_nonzero(p_values <= .05)
    interval = stats.binomtest(rejected, nsets).proportion_ci(method="exact")
    return dict(scenario=scenario, nsets=nsets, draws=draws, exact=exact, seed=seed,
                false_rejections=int(rejected), family_error_rate=rejected/nsets,
                binomial_95_interval=[interval.low, interval.high],
                seconds=round(time.perf_counter()-started, 2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--nsets", type=int, default=2000)
    parser.add_argument("--draws", type=int, default=499)
    parser.add_argument("--exact", action="store_true")
    parser.add_argument("--scenario", default=None)
    parser.add_argument("--seed", type=int, default=None)
    args = parser.parse_args()
    for index, name in enumerate(("independent", "correlated", "random_drug_slopes",
                                  "nonspherical", "heavy_tails", "partial_null",
                                  "skewed_stress")):
        if args.scenario is None or args.scenario == name:
            seed = 8301+index if args.seed is None else args.seed
            print(json.dumps(simulate(name, args.nsets, args.draws, seed, args.exact)), flush=True)
