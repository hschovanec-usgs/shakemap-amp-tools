from obspy.core.stream import Stream
from obspy.core.trace import Trace
from obspy.signal.invsim import corn_freq_2_paz, simulate_seismometer
import warnings

from amptools.constants import GAL_TO_PCTG
from pgm.rotation import rotate


def get_acceleration(stream, units='%%g'):
    """
    Returns a stream of acceleration with specified units.

    Args:
        stream (obspy.core.stream.Stream): Strong motion timeseries
            for one station. With units of g (cm/s).
        units (str): Units of accelearation for output. Default is %g

    Returns:
        obpsy.core.stream.Stream: stream of acceleration.
    """

    accel_stream = Stream()
    for trace in stream:
        accel_trace = trace.copy()
        if units == '%%g':
            accel_trace.data = trace.data * GAL_TO_PCTG
            accel_trace.stats['units'] = '%%g'
        elif units == 'm/s/s':
            print('m/s/s')
            accel_trace.data = trace.data * 0.01
            accel_trace.stats['units'] = 'm/s/s'
        else:
            accel_trace.data = trace.data * 100
            accel_trace.stats['units'] = 'cm/s/s'
        accel_stream.append(accel_trace)

    return accel_stream

def get_spectral(period, stream, damping, rotation=None):
    """
    Returns a stream of spectral response with units of %%g.

    Args:
        period (float): Period for spectral response.
        stream (obspy.core.stream.Stream): Strong motion timeseries
            for one station.
        damping (float): Damping of oscillator.
        rotation (str): Wheter a rotation matrix should be return and the
            specific type or rotation. Default is None.

    Returns:
        obpsy.core.stream.Stream: stream of spectral response.
    """
    T = period
    freq = 1.0 / T
    omega = (2 * 3.14159 * freq) ** 2
    paz_sa = corn_freq_2_paz(freq, damp=damping)
    paz_sa['sensitivity'] = omega
    paz_sa['zeros'] = []
    spect_stream = Stream()

    horizontals = []
    for trace in stream:
        if trace.stats['channel'].upper().find('Z') < 0:
            horizontals += [trace.copy()]
    h1_stats = horizontals[0].stats

    if rotation is None:
        for trace in stream:
            samp_rate = trace.stats['sampling_rate']
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                dd = simulate_seismometer(trace.data, samp_rate,
                                          paz_remove=None,
                                          paz_simulate=paz_sa,
                                          taper=True,
                                          simulate_sensitivity=True,
                                          taper_fraction=0.05)
            period_str = 'T' + '{:04.2f}'.format(T)
            stats_out = trace.stats.copy()
            stats_out['period'] = period_str
            spect_trace = Trace(dd, stats_out)
            spect_trace.data = spect_trace.data * GAL_TO_PCTG
            spect_trace.stats['units'] = '%%g'
            spect_stream.append(spect_trace)
        return spect_stream
    elif rotation.lower() == 'nongm':
        if len(horizontals) != 2:
            warnings.warn('Spectral amplitude rotation could not be performed.')
            return
        rot = [rotate(horizontals[0], horizontals[1], combine=True)]
    elif rotation.lower() == 'gm':
        if len(horizontals) != 2:
            warnings.warn('Spectral amplitude rotation could not be performed.')
            return
        rot1, rot2 = rotate(horizontals[0], horizontals[1], combine=False)
        rot = [rot1, rot2]
    rotated = []
    for rot_matrix in rot:
        rotated_spectrals = np.zeros(rot_matrix.shape)
        for idx, row in enumerate(rot_matrix):
            samp_rate = h1_stats['sampling_rate']
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                dd = simulate_seismometer(row, samp_rate,
                                          paz_remove=None,
                                          paz_simulate=paz_sa,
                                          taper=True,
                                          simulate_sensitivity=True,
                                          taper_fraction=0.05)

            period_str = 'T' + '{:04.2f}'.format(T)
            stats_out = h1_stats.copy()
            stats_out['period'] = period_str
            spect_trace = Trace(dd, stats_out)
            spect_trace.data = spect_trace.data * GAL_TO_PCTG
            spect_trace.stats['units'] = '%%g'
            rotated_spectrals[idx] = spect_trace
        rotated += [rotated_spectrals]
    return rotated


def get_velocity(stream, rotation=None):
    """
    Returns a stream of velocity with units of cm/s.

    Args:
        stream (obspy.core.stream.Stream): Strong motion timeseries
            for one station.

    Returns:
        obpsy.core.stream.Stream: stream of velocity.
    """
    veloc_stream = Stream()
    for trace in stream:
        veloc_trace = trace.copy()
        veloc_trace.integrate()
        veloc_trace.stats['units'] = 'cm/s'
        veloc_stream.append(veloc_trace)
    return veloc_stream
