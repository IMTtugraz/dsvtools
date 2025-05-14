# Copyright 2024. Institute of Biomedical Imaging. TU Graz.
# All rights reserved. Use of this source code is governed by
# a BSD-style license which can be found in the LICENSE file.
#
# Authors:
# 2024 Philip Schaten <philip.schaten@tugraz.at>
# 2024 Christian Holme <holme@tugraz.at>

import re
import numpy as np


class LineParser:

    # States
    START = "start"
    INFO_BLOCK = "info_block"
    MDH = "mdh_block"

    # Types
    Block = "event block start and duration"
    Adc = "An ADC sample result"
    Mdh_value = "An MDH Key-Value Pair"
    RF = "RF pulse relative start and duration"
    Error = "parsing error"

    # Expressions
    block_start_ex = "^\[INFO\]$"
    info_block_start_ex = r"^\* EventBlock (?P<start>[0-9]*) (?P<duration>[0-9]*)$"
    info_block_end_ex = r"^#INFO-END$"
    adc_start_ex = r"^adc MeasHeader (?P<start>[0-9]*) (?P<duration>[0-9]*)$"
    adc_value_ex = r"^(?P<key>[^:]+) *: (?P<value>[^;]+)"
    rf_value_ex = r"^\|\s+(?P<start>\d+) \| \[\S+\]:\s+[\d.]+\/\s*(?P<duration>\d+)\s+\|.*$"


    def __init__(self):
        self.state = LineParser.START

    def parse(self, line):

        if self.state == LineParser.START:
            if match := re.match(LineParser.adc_start_ex, line):
                self.state = LineParser.MDH
                d = match.groupdict()
                return LineParser.Adc, dict(start = int(d['start']), duration = int(d['duration']))
            elif match := re.match(LineParser.info_block_start_ex, line):
                d = match.groupdict()
                self.state = LineParser.INFO_BLOCK
                return LineParser.Block, dict(start = int(d['start']), duration = int(d['duration']))

        if self.state == LineParser.INFO_BLOCK:
            if match := re.match(LineParser.rf_value_ex, line):
                d = match.groupdict()
                return LineParser.RF, dict(start = int(d['start']), duration = int(d['duration']))

        if self.state == LineParser.INFO_BLOCK:
            if re.match(LineParser.info_block_end_ex, line):
                self.state = LineParser.START

        if self.state == LineParser.MDH:
            if match := re.match(LineParser.adc_value_ex, line):
                d = match.groupdict()
                d['key'] = d['key'].strip()
                return LineParser.Mdh_value, d
            elif re.match(LineParser.block_start_ex, line):
                self.state = LineParser.START
            else:
                return LineParser.Error, "Found neither MDH Key-Value Pair nor End of MDH in MDH context."

        return None, None


def _join_adcs(in_data):
    data = in_data['blocks']
    last_index = 0
    for adc in in_data['adcs']:
        found = False
        i = last_index
        steps = 0
        while (not found) and (steps < len(data)):
            block = data[i]
            if adc['start'] > block['start'] and \
            adc['start'] + adc['duration'] < block['start'] + block['duration']:
                found = True
                last_index = i
                if not 'adcs' in block:
                    data[i]['adcs'] = []
                data[i]['adcs'].append(adc)
            steps += 1
            i = (i+1) % len(data)
        if not found:
            print("WARN: Discarding lonely ADC object with starttime {}!".format(adc['start']))
    return data


def dsvparse_inf(fname):
    """Parse .dsv files"""
    data = {'adcs': [], 'blocks': []}
    parser = LineParser()

    with open(fname, 'r', encoding='latin1') as f:
        for i,line in enumerate(f):
            typ, result = parser.parse(line)
            if typ == LineParser.Error:
                raise RuntimeError(f"In context {parser.state}, could not parse line {i}: {line}.\nCause: {result}")
            elif typ == LineParser.RF:
                data['blocks'][-1]['rf'] = result
            elif typ == LineParser.Block:
                data['blocks'].append(result)
            elif typ == LineParser.Adc:
                data['adcs'].append(result)
            elif typ == LineParser.Mdh_value:
                if result['key'] == "ushSamplesInScan":
                    data['adcs'][-1]['samples'] = int(result['value'])

                val = result['value'].strip()
                if re.match(r"^[0-9]+$", val):
                    val = int(val)
                elif re.match(r"^[0-9]+(\.[0-9]+)?$", val):
                    val = float(val)

                data['adcs'][-1][result['key']] = val



    return _join_adcs(data)


def main():
    import argparse
    import json
    parser = argparse.ArgumentParser(description="Converter for _INF.dsv files")
    parser.add_argument("DSV",  help="Input .dsv")
    parser.add_argument("out", default=None, nargs="?", help="Output (Store metadata as json; default stdout)")
    args = parser.parse_args()

    data = dsvparse_inf(args.DSV)

    if args.out:
        with open(args.out, 'w') as f:
            json.dump(data, f)
    else:
        print(json.dumps(data))

if __name__ == "__main__":
    main()
