# stdlib imports
import warnings

# third party imports
import numpy as np
from obspy.core.stream import Stream
from obspy.core.trace import Trace
from obspy.signal.invsim import corn_freq_2_paz, simulate_seismometer

# local imports
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
            accel_trace.data = trace.data * 0.01
            accel_trace.stats['units'] = 'm/s/s'
        else:
            accel_trace.data = trace.data
            accel_trace.stats['units'] = 'cm/s/s'
        accel_stream.append(accel_trace)

    return accel_stream

def calculate_spectrals(times, acc, period, damping):
    cdef float dt = times[1] - times[0]
    cdef int npoints = len(acc)
    cdef int initialdisp = 0
    cdef int initialvel = 0
    cdef int damp = damping
    d = np.zeros(npoints)
    v = np.zeros(npoints)
    aa = np.zeros(npoints)

    d[0] = initialdisp
    v[0] = initialvel

    cdef float w = 2 * np.pi / period
    cdef float ww = w ** 2

    cdef float sqrtd  = np.sqrt(1 - damp**2)
    cdef float wd = w * sqrtd
    cdef float e = np.exp( - w * dt * damp)
    cdef float sine = np.sin(wd * dt)
    cdef float cosine = np.cos(wd * dt)
    aa[0] =  - ww * d[0] - (2. * damp * w) * v[0]

    cdef float a11 = e * ((damp / sqrtd) * sine + cosine)
    cdef float a12 = e * sine / wd
    cdef float cc = (e * (((1. - 2 * damp**2) / wd - (damp / sqrtd) * dt) * sine - ((2. * damp / w) + dt) * cosine) + (2. * damp / w)) * (-1. / (ww * dt))
    cdef float cd = (e * ( - (1. - 2 * damp**2) / wd * sine + (2. * damp / w) * cosine) + dt - (2. * damp / w)) * (-1. / (ww * dt))
    cdef float a21 =  - e * w * sine / sqrtd
    cdef float a22 = e * (cosine - (damp / sqrtd) * sine)
    cdef float ccp = (e * ((w * dt / sqrtd  + (damp / sqrtd)) * sine + cosine) - 1.) * (-1. / (ww * dt))
    cdef float cdp = (1. - a11) * (-1. / (ww * dt))

    for i in np.arange(1, npoints):
        d[i]  = a11 * d[i - 1]  +  \
                a12 * v[i - 1]  +  \
                cc * acc[i - 1]  +  \
                cd * acc[i]
        v[i]  = a21 * d[i - 1] + \
                a22 * v[i - 1] + \
                ccp * acc[i - 1] + \
                cdp * acc[i]
        aa[i] =  - ww * d[i] - (2. * damp * w) * v[i]
    return (aa, v, d)

def get_spectral(period, stream, damping=0.05, rotation=None):
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
    spect_stream = Stream()

    horizontals = []
    for idx, trace in enumerate(stream):
        # Group all of the max values from traces without
        # Z in the channel name
        if 'Z' not in trace.stats['channel'].upper():
            horizontals += [trace.copy()]
    h1_stats = horizontals[0].stats
    h1_times = horizontals[0].times()

    if rotation is None:
        for trace in stream:
            acc_sa = calculate_spectrals(trace.times(), trace.data,
                    period, damping)[0]
            stats = trace.stats.copy()
            stats['units'] = '%%g'
            acc_sa = acc_sa * GAL_TO_PCTG
            spect_trace = Trace(data=acc_sa, header=stats)
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
            acc_sa = calculate_spectrals(h1_times, row,
                    period, damping)[0]
            stats = h1_stats.copy()
            stats['units'] = '%%g'
            acc_sa = acc_sa * GAL_TO_PCTG
            spect_trace = Trace(data=acc_sa, header=stats)
            rotated_spectrals[idx] = spect_trace
        rotated += [rotated_spectrals]
    return rotated


def get_velocity(stream):
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
