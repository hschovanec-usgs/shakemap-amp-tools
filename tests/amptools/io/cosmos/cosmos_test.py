#!/usr/bin/env python

import os.path
import numpy as np
from amptools.io.cosmos.core import is_cosmos, read_cosmos

def test_cosmos():
    homedir = os.path.dirname(os.path.abspath(__file__)) #where is this script?
    datadir = os.path.join(homedir,'..','..','..','data','cosmos')
    one_channel = os.path.join(datadir,'Cosmos12TimeSeriesTest.v1')
    two_channels = os.path.join(datadir,'Cosmos12TimeSeriesTest2.v1')

    assert is_cosmos(one_channel)
    try:
        assert is_cosmos(os.path.abspath(__file__))
    except AssertionError as ae:
        assert 1==1

    # test a one channel cosmos file
    stream1 = read_cosmos(one_channel)

    stats = stream1[0].stats
    assert stats['station'] == 'J2236'
    assert stats['delta'] == .005000
    assert stats['location'] == '--'
    assert stats['network'] == 'CE'
    dt = '%Y-%m-%dT%H:%M:%SZ'
    assert stats['starttime'].strftime(dt) == '2005-06-16T20:53:04Z'
    assert stats.coordinates['latitude'] == 34.046
    assert stats.coordinates['longitude'] == -117.035
    assert stats.coordinates['elevation'] == 15
    assert stats.standard['station_name'] == 'Yucaipa - Bryant & Oak Glen'
    assert stats.standard['instrument'] == 'Kinemetrics FBA-11 accelerometer'
    assert stats.standard['sensor_serial_number'] == '1889'
    dt = '%Y-%m-%dT%H:%M:%SZ'
    assert stats.standard['process_time'].strftime(dt) == '2005-06-17T12:01:00Z'
    assert stats.format_specific['sensor_sensitivity'] == 220
    assert stats.standard['horizontal_orientation'] == 340
    assert stats.standard['instrument_period'] == 1.0 / 25
    assert stats.standard['instrument_damping'] == 0.20
    assert stats.standard['process_level'] == 'V2'
    assert stats.standard['source_format'] == 'cosmos'
    assert stats.standard['structure_type'] == 'Building'
    assert stats.standard['source'] == 'California Geological Survey'
    assert stats.format_specific['scaling_factor'] == 1
    assert stats.format_specific['v30'] == 120
    assert stats.format_specific['physical_units'] == 'cm/sec/sec'
    assert stats.format_specific['least_significant_bit'] == 123.45
    assert stats.format_specific['low_filter_type'] == 'Butterworth single direction'
    assert stats.format_specific['low_filter_corner'] == 4
    assert stats.format_specific['low_filter_decay'] == 3
    assert stats.format_specific['high_filter_type'] == 'Rectangular'
    assert stats.format_specific['high_filter_corner'] == 40
    assert stats.format_specific['high_filter_decay'] == 4
    assert stats.format_specific['maximum'] == -161.962
    assert stats.format_specific['maximum_time'] == 27.85
    assert stats.format_specific['station_code'] == 10
    assert stats.format_specific['record_flag'] == 'No problem'

    # test that one channel is created
    assert len(stream1) == 1

    # read the maximum from the text header check that the trace max
    # is the equivalent when rounded to the same number of decimal places
    with open(one_channel, 'rt') as f:
        file_line = f.readlines()[10].replace(' ', '').lower()
    file_max = file_line[file_line.find('max=') + 4: file_line.find('cm')]
    assert np.round(stream1[0].max(), 3) == float(file_max)

    # test a two channel cosmos file
    stream2 = read_cosmos(two_channels)

    # test that reading a file that is a valid station type returns a
    # stream with traces
    building_code = 10
    stream3 = read_cosmos(two_channels, valid_station_types=[building_code])
    assert stream3.count() == 1

    # test that reading a file that is not a valid station type returns an
    # empty stream
    stream4 = read_cosmos(two_channels, valid_station_types=[1, 2, 3, 4])
    assert stream4.count() == 0

    # test that two channels are created
    assert len(stream2) == 2

    # test that reading a file that is a valid station type returns a
    # stream with traces
    building_code = 10
    stream3 = read_cosmos(one_channel, valid_station_types=[building_code])
    assert stream3.count() == 1



if __name__ == '__main__':
    test_cosmos()
