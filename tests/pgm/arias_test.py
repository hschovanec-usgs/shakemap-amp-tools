#!/usr/bin/env python

# stdlib imports
import os.path

# local imports
from amptools.io.knet.core import read_knet
from amptools.stream import group_channels
from pgm.station_summary import StationSummary

from obspy import read


def test_arias():
    homedir = os.path.dirname(os.path.abspath(
        __file__))  # where is this script?
    data_dir = os.path.join(homedir, '..', 'data', 'knet')
    files = ['AIC0090010061330.EW', 'AIC0090010061330.NS',
            'AIC0090010061330.UD', 'SMN0020010061330.EW',
            'SMN0020010061330.UD', 'SMN0020010061330.NS']
    sts = []
    for file in files:
        sts += [read(os.path.join(data_dir, file))]
    sts = group_channels(sts)
    AIC009 = sts[0]
    SMN002 = sts[1]
    station_summary = StationSummary.from_stream(AIC009, ['channels'],
            ['arias'])
    station_dict = station_summary.pgms
    print(station_dict)
    station_summary = StationSummary.from_stream(SMN002, ['channels'],
            ['arias'])
    station_dict = station_summary.pgms
    print(station_dict)


if __name__ == '__main__':
    test_arias()
