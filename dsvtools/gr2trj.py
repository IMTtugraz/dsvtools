#!/usr/bin/python3
import numpy as np

"""
Convert INF.dsv and GRX.dsv into trajectories (==sampled kspace position)

This does not incorporate RF effects. I.e., spin-/ stimulated echo is not possible,
perfect spoiling is assumed.

ADC sampling works as follows?
      |-*-|-*-|-*-|-*-|
      Î     Î
      |     \- sample
   adcstart
    --adcduration--

Gradient values are given for the *center* of the respective sample!
https://www.magnetom.net/t/code-for-dsvs/5152/8
"""

ADC_DT_SHIFT_DEFAULT = .5
GRADIENT_DT_SHIFT = 1

def get_index_safe(t, dt, r=0):
    i = t // dt
    if (diff := abs(i - t / dt)) > 1e-8:
        if r==1:
            i += 1
        elif r==-1:
            pass
        else:
            raise RuntimeError("t={} not an integer multiple of dt={}".format(t, dt))
    return i

def extract_measurements(x, dt, inf):
    maxtl = 0
    meas_indices = []
    adcs = []
    total_adc_count = 0;

    # For each block, find adc durations.
    # furthermore determine max block length.
    for i, block in enumerate(inf):
        if 'adcs' in block:
            meas_indices.append({})

            adc_count = len(block['adcs'])

            meas_indices[-1]['start'] = get_index_safe(block['start'], dt)  # Absolute start index
            meas_indices[-1]['stop'] = get_index_safe(block['adcs'][-1]['start'] + block['adcs'][-1]['duration'], dt, r=1)  # Absolute stop index

            meas_indices[-1]['adc_count'] = adc_count
            total_adc_count += adc_count

            ksp_start_ind = block['start'] + block['rf']['start'] + block['rf']['duration']//2

            meas_indices[-1]['ksp_start'] = get_index_safe(ksp_start_ind, dt, r=1)

            maxtl = max(maxtl, meas_indices[-1]['stop'] - meas_indices[-1]['ksp_start'])

            for adc in block['adcs']:
                adcs.append(dict(samples = int(adc['samples']), duration = int(adc['duration']),
                             start = int(adc['start'] - ksp_start_ind),
                             moment_idx = i))


    moments = np.zeros(shape=(int(maxtl), len(meas_indices)), dtype=float)

    for i, m in enumerate(meas_indices):
        gradients_here = x[m['ksp_start']:m['stop']]

        # moments = scipy.integrate.cumulative_trapezoid(moments, dx=args.dt, axis=0, initial=0)
        # cumulative sum (!!!) is correct approximation due to how sampling is chosen.!

        moments[:m['stop']-m['ksp_start'], i] = np.cumsum(gradients_here, axis = 0) * dt

    return moments, adcs

def resample(x, dt, adc_start, adc_dur, adc_n, adc_dt_shift = ADC_DT_SHIFT_DEFAULT, gradient_dt_shift = GRADIENT_DT_SHIFT):
    adc_dt = adc_dur / adc_n
    xt    = np.arange(len(x)) * dt                 + gradient_dt_shift * dt
    adc_t = np.arange(adc_n)  * adc_dt + adc_start + adc_dt_shift      * adc_dt

    trj = np.interp(adc_t, xt, x)
    return trj

if __name__ == '__main__':
    import argparse
    import json
    parser = argparse.ArgumentParser()

    parser.add_argument("g_i", nargs=3, help="Gradient fields")
    parser.add_argument("inf", help="Block information (json)")
    parser.add_argument("trajectory", help="Output")
    g = parser.add_mutually_exclusive_group()
    g.add_argument("--cfl", action="store_true", help="Use cfl format")

    parser.add_argument("--dt",    type=int, help="Gradient dt", default = 1)
    parser.add_argument("--unit",   type=float, help="Trajectory normalization", default = 1)
    parser.add_argument("--adcshift",   type=float, help="Fraction of dwelltime by which acquisition is shifted", default = ADC_DT_SHIFT_DEFAULT)
    args = parser.parse_args()

    if args.cfl:
        from .contrib import cfl

    with open(args.inf, 'r') as f:
        inf = json.load(f)

    trj = []
    for g in args.g_i:
        x = cfl.readcfl(g).real if args.cfl else np.load(g)
        moments, adcs = extract_measurements(x, args.dt, inf)
        lines = []
        for a in adcs:
            m = moments[:,a['moment_idx']]
            lines.append(resample(m, args.dt, a['start'], a['duration'], a['samples'], adc_dt_shift=args.adcshift))
        trj.append(np.stack(lines, axis=1))
    trj = np.stack(trj, axis=0)
    trj *= args.unit

    # simulation axis to real axis
    trj2 = np.zeros_like(trj)
    trj2[0] = -1 * trj[1]
    trj2[1] = -1 * trj[0]
    trj2[2] =  1 * trj[2]

    if args.cfl:
        cfl.writecfl(args.trajectory, trj2)
    else:
        np.save(args.trajectory, trj2)
