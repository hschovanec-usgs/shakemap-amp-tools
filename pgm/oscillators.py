from obspy.core.stream import Stream


GAL_TO_PCTG = 1 / (9.8)
GAL_TO_MSS = 0.01


def get_acceleration(stream, units=None):
    """
    Returns a stream of acceleration with units of %%g.

    Args:
        stream (obspy.core.stream.Stream): Strong motion timeseries
            for one station.

    Returns:
        obpsy.core.stream.Stream: stream of acceleration.
    """
    accel_stream = Stream()
    for trace in stream:
        accel_trace = trace.copy()
        if units == '%%g':
            accel_trace.data = trace.data * GAL_TO_PCTG
            accel_trace.stats['units'] = '%%g'
        elif units == 'm/s':
            accel_trace.data = trace.data * GAL_TO_MSS
            accel_trace.stats['units'] = 'm/s'
        accel_stream.append(accel_trace)
    return accel_stream

def get_spectral(period, stream, damping):
    """
    Returns a stream of spectral response with units of %%g.

    Args:
        period (float): Period for spectral response.
        stream (obspy.core.stream.Stream): Strong motion timeseries
            for one station.
        damping (float): Damping of oscillator.

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
        self.period_str = period_str
        spect_trace = Trace(dd, stats_out)
        spect_trace.data = spect_trace.data * GAL_TO_PCTG
        spect_trace.stats['units'] = '%%g'
        spect_stream.append(spect_trace)
    return spect_stream

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
