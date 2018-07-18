#!/usr/bin/env python

# stdlib imports
import os.path
import warnings

# third party imports
import numpy as np

# local imports
from amptools.io.geonet.core import read_geonet
from pgm.station_summary import StationSummary


def test_stationsummary():
    homedir = os.path.dirname(os.path.abspath(
        __file__))  # where is this script?
    datafile = os.path.join(homedir, '..', 'data', 'geonet',
                            '20161113_110259_WTMC_20.V2A')
    target_imcs = np.sort(np.asarray(['GREATER_OF_TWO_HORIZONTALS',
                                      'HN1', 'HN2', 'HNZ']))
    target_imts = np.sort(np.asarray(['SA1.0', 'PGA', 'PGV']))
    stream = read_geonet(datafile)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        stream_summary = StationSummary.from_stream(stream,
                                                    ['greater_of_two_horizontals',
                                                     'channels', 'invalid'],
                                                    ['sa1.0', 'PGA', 'pgv', 'invalid'])
        original_stream = stream_summary.stream
        stream_summary.stream = []
        final_stream = stream_summary.stream
        assert original_stream == final_stream
        original_code = stream_summary.station_code
        stream_summary.station_code = ''
        final_code = stream_summary.station_code
        assert original_code == final_code
        original_oscillators = stream_summary.oscillators
        final_oscillators = stream_summary.oscillators
        assert original_oscillators == final_oscillators
        np.testing.assert_array_equal(np.sort(stream_summary.components),
                                      target_imcs)
        np.testing.assert_array_equal(np.sort(stream_summary.imts),
                                      target_imts)
        np.testing.assert_almost_equal(stream_summary.get_pgm('PGA', 'HN1'),
                                       99.3173469387755, decimal=1)
        target_available = np.sort(np.asarray([
            'calculate_greater_of_two_horizontals', 'calculate_channels',
            'calculate_gmrotd']))
        imcs = stream_summary.available_imcs
        np.testing.assert_array_equal(np.sort(imcs), target_available)
        target_available = np.sort(np.asarray(['calculate_pga',
                                               'calculate_pgv',
                                               'calculate_sa',
                                               'calculate_arias']))
        imts = stream_summary.available_imts
        np.testing.assert_array_equal(np.sort(imts), target_available)

    test_pgms = {'SA1.0':{
                     'HNZ': 27.41836158024183,
                     'HN1': 135.83462000177008,
                     'GREATER_OF_TWO_HORIZONTALS': 135.83462000177008,
                     'HN2': 84.17504097202713},
                 'PGA': {
                         'HNZ': 183.7722361866693,
                         'HN1': 99.24999872535474,
                         'GREATER_OF_TWO_HORIZONTALS': 99.24999872535474,
                         'HN2': 81.23467239067368},
                 'PGV': {'HNZ': 37.47740000000001,
                         'HN1': 100.81460000000004,
                         'GREATER_OF_TWO_HORIZONTALS': 100.81460000000004,
                         'HN2': 68.4354}}
    datafile = os.path.join(homedir, '..', 'data', 'geonet',
                            '20161113_110313_THZ_20.V2A')
    invalid_stream = read_geonet(datafile)
    station_code = 'WTMC'
    pgm_summary = StationSummary.from_pgms(station_code, test_pgms)

    # TODO: change to an almost equal assert
    assert pgm_summary.pgms == stream_summary.pgms

    # oscillators cannot be calculated without a stream
    try:
        pgm_summary.generate_oscillators(pgm_summary.imts, 0.05)
        success = True
    except Exception:
        success = False
    assert success == False
    # Invalid stream inputs should be rejected
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        pgm_summary.stream = []
        assert pgm_summary.stream is None
        pgm_summary.stream = invalid_stream
        assert pgm_summary.stream is None
        pgm_summary.stream = stream
        assert pgm_summary.stream == stream


if __name__ == '__main__':
    test_stationsummary()
