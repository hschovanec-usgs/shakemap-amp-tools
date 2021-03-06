#!/usr/bin/env python

# stdlib imports
import os.path
from collections import OrderedDict

# third party imports
import numpy as np

from amptools.io.smc.core import is_smc, read_smc
from amptools.stream import group_channels


def test_smc():
    homedir = os.path.dirname(os.path.abspath(
        __file__))  # where is this script?
    datadir = os.path.join(homedir, '..', '..', '..', 'data', 'smc')

    files = OrderedDict([('0111a.smc', (1.5057E+0, -2.8745E-1)),
                         ('0111b.smc', (-1.2518E+1, -1.6806E+0)),
                         ('0111c.smc', (-5.8486E+0, -1.1594E+0))])

    streams = []
    for tfilename, accvals in files.items():
        filename = os.path.join(datadir, tfilename)
        assert is_smc(filename)

        # test acceleration from the file
        stream = read_smc(filename)

        # test for one trace per file
        assert stream.count() == 1

        # test that the traces are acceleration
        for trace in stream:
            assert trace.stats.standard.units == 'acc'

        # compare the start/end points
        np.testing.assert_almost_equal(accvals[0], stream[0].data[0])
        np.testing.assert_almost_equal(accvals[1], stream[0].data[-1])

        # append to list of streams, so we can make sure these group together
        streams.append(stream)

    newstreams = group_channels(streams)
    assert len(newstreams) == 1


if __name__ == '__main__':
    test_smc()
