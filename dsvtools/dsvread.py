# Copyright 2024. Institute of Biomedical Imaging. TU Graz.
# All rights reserved. Use of this source code is governed by
# a BSD-style license which can be found in the LICENSE file.
#
# Authors:
# 2024 Philip Schaten <philip.schaten@tugraz.at>

import re
import numpy as np

state_names = ['start', 'FILETYPE', 'DEFINITIONS', 'VALUES']
states = { x: i for i,x in enumerate(state_names) }

def dsv_assert(expr, msg=None):
    if not expr:
        raise RuntimeError(msg)

def _statemachine(data, state, line):
    """Internally used function implementing the logic of .dsv files"""
    parsed = 0

    if 0 == parsed and line[:-1] == "[FILETYPE]":
        state = states['FILETYPE']
        parsed += 1

    if 0 == parsed and line[:-1] == "[DEFINITIONS]":
        state = states['DEFINITIONS']
        dsv_assert('definitions' not in data, msg="Two [DEFINITIONS]?")
        data['definitions'] = {}
        parsed += 1

    if 0 == parsed and line[:-1] == "[VALUES]":
        state = states['VALUES']
        dsv_assert('values' not in data, msg="Two [VALUES]?")
        data['values'] = []
        parsed += 1

    if 0 == parsed and state == states['FILETYPE'] and (match := re.match(r"^FORMAT=(?P<format>.*)$", line)):
        dsv_assert('format' not in data, msg="Two FORMAT=")
        data['format'] = match.groupdict()['format']
        dsv_assert(data['format'] == "DSV_V0100", msg="Format {} probably not understood by this script".format(data['format']))
        state = states['start']
        parsed += 1

    if 0 == parsed and state == states['DEFINITIONS']:
        if (match := re.match(r"^(?P<key>.*)=(?P<value>.*)$", line)):
            d = match.groupdict()
            dsv_assert(d['key'] not in data['definitions'], msg="Two keys {}".format(d['key']))
            data['definitions'][d['key']] = str(d['value'])
        else:
            state = states['start']
        parsed +=1

    if 0 == parsed and state == states['VALUES']:
        if len(line.strip()) > 0:
            data['values'].append(int(line[:-1]))
        else:
            state = states['start']
        parsed += 1

    return state

def dsvparse(fname):
    """Parse .dsv files"""
    data = {}
    state = states['start']
    with open(fname, 'r', encoding='latin1') as f:
        for line in f:
            state = _statemachine(data, state, line)

    return data

def dsvconvert(data):
    """Convert parsed dsv files into a plain format.
    DSV format:
    - if a value occurs twice, the following (third) value indicates the
      number of subsequent repetitions. I.e. the sequence
      1,1,1,1,2,1,1,2 becomes
      1,1,*2*,2,1,1,*0*,2
    - Additionally, only the first value is absolute, the others are deltas.
    - Data needs to be divided by "VERTFACTOR".
    - A de-dsv'ed total of "SAMPLES" is described in the file.
    - time distance between two uncompressed samples is HORIUNITNAME * HORIDELTA.
    """

    # Parse metadata

    meta = {}

    special_cases = ['VERTUNITNAME', 'HORIDELTA', 'HORIUNITNAME', 'VERTFACTOR']
    for k,v in data['definitions'].items():
        if k not in special_cases:
            meta[k] = v

    meta['yunitname'] = data['definitions']['VERTUNITNAME']
    meta['xunit']     = float(data['definitions']['HORIDELTA'])
    meta['xunitname'] = data['definitions']['HORIUNITNAME']
    if data['definitions']['HORIUNITNAME'] == 'µs':
        meta['xunit'] *= 1e-6
        meta['xunitname'] = 's'

    # Parse values

    values = np.zeros(int(data['definitions']['SAMPLES']), dtype=int)

    x = data['values'][0]
    values[0] = x
    last = x

    i = 1
    repeater = False
    for x in data['values'][1:]:
        if repeater:
            values[i:i + x] = values[i-1] + np.arange(1,x+1, dtype=int) * last
            i += x
            last = None
            repeater = False
        else:
            values[i] =       values[i-1] + x
            i += 1
            if x == last:
                repeater = True
            last = x

    dsv_assert(i == len(values), msg="decompression failed: Expected {} values, got {} values.".format(len(values), i))

    vf = float(data['definitions']['VERTFACTOR'])

    return meta, values/vf


if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser(description="Converter for .dsv files")
    parser.add_argument("DSV", help="Input .dsv")
    parser.add_argument("out", default=None, nargs="?", help="Output (default: Store metadata as json; if given: store samples")
    out_format = parser.add_mutually_exclusive_group()
    out_format.add_argument("--cfl", action="store_true", help="Save samples as .cfl instead of npy")
    out_format.add_argument("--txt", action="store_true", help="Save samples as txt, one per line.")
    out_format.add_argument("--meta", action="store_true", help="Save meta data as json")

    args = parser.parse_args()

    data = dsvparse(args.DSV)
    meta, data = dsvconvert(data)

    out_fname = args.out

    if args.cfl:
        assert(not out_fname is None)
        from .contrib import cfl
        cfl.writecfl(out_fname, data)

    elif args.txt:
        if out_fname:
            np.savetxt(out_fname, data)
        else:
            for x in data:
                print(x)

    elif args.meta or out_fname is None:
        import json
        f = open(out_fname, 'w') if out_fname else sys.stdout
        json.dump(meta, f)
        f.close()

    else:
        assert(not out_fname is None)
        np.save(out_fname, data)
