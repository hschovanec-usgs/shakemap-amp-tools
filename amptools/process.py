"""Helper functions for processing strong ground motion data"""

# stdlib imports
import numpy as np
import warnings

# third party imports
from obspy.signal.util import next_pow_2
from obspy.signal.konnoohmachismoothing import konno_ohmachi_smoothing
from scipy.optimize import curve_fit

# local imports
from amptools.config import get_config


CONFIG = get_config()


def filter_detrend(trace, taper_type='cosine', taper_percentage=0.05,
                   filter_type='highpass', filter_frequency=0.02,
                   filter_zerophase=True, filter_corners=4):
    """
    Read files from a directory and return stream objects.

    Args:
        trace (obspy.core.trace.Trace): Trace of strong motion data.
        taper_type (str): Type of taper used for processing.
            Default is 'cosine'.
        taper_percentage (float): Maximum taper percentage.
            Default is 0.05.
        filter_type (str): Type of filter used for processing.
            Default is 'highpass'.
        filter_frequency (float): Filter corner frequency.
            Default is 0.02.
        filter_zerophase (bool): If True, applies a forward and backward
            filter. Results in a zero phase shift. Default is True.
        filter_corners (int): Filter corners / order.

    Returns:
        list : List of obspy.core.stream.Stream objects

    Notes:
        Depricated function. This will be removed.
    """
    trace.detrend('linear')
    trace.detrend('demean')
    trace.taper(max_percentage=taper_percentage, type=taper_type)
    trace.filter(filter_type, freq=filter_frequency,
                 zerophase=filter_zerophase, corners=filter_corners)
    trace.detrend('linear')
    trace.detrend('demean')
    return trace


def correct_baseline_mean(stream):
    """
    Performs a zero-order baseline correction by subtracting the mean from
    the entire waveform. This is step 1 as found in the Rennolet et. al
    paper (https://doi.org/10.1193/101916EQS175DP)

    Args:
        stream (obpsy.core.stream.Stream): Stream of raw waveform data.

    Returns:
        stream (obspy.core.stream.Stream): Stream of corrected waveform data.
    """
    for tr in stream:
        tr.detrend(type='demean')
    return stream


def remove_clipped(stream, max_count=2000000):
    """
    Identify clipped waveforms having a count greater than the specified
    max count. This will remove any clipped waveforms from the stream.
    This is step 2 as found in the Rennolet et. al paper
    (https://doi.org/10.1193/101916EQS175DP)

    Args:
        stream (obspy.core.stream.Stream): Stream of raw data in counts.
        max_count (int): Maximum count for clipping.
            Default is 2 million.

    Returns:
        stream (obspy.core.stream.Stream): Stream of raw data with clipped
        waveforms removed
    """
    for tr in stream:

        if abs(tr.max()) >= max_count:
            stream.traces.remove(tr)
            warnings.warn('Clipped trace was removed from the stream')
            warnings.warn(tr.get_id())
    return stream


def check_max_amplitude(trace, min_amp=10e-7, max_amp=5e3):
    """
    Checks that the maximum amplitude of the trace is within a defined
    range.

    Args:
        trace (obspy.core.trace.Trace): Trace of strong motion data.
        min_amp (float): Minimum amplitude for the acceptable range.
            Default is 10e-7.
        max_amp (float): Maximum amplitude for the acceptable range.
            Default is 5e3.

    Returns:
        bool: True if trace passes the check. False otherwise.
    """
    amplitude = {'min': min_amp, 'max': max_amp}
    trace = _update_params(trace, 'amplitude', amplitude)
    if (abs(trace.max()) >= min_amp and abs(trace.max()) <= max_amp):
        return True
    else:
        return False


def trim_total_window(trace, org_time, epi_dist, vmin=1.0):
    """
    Trims a trace to the window using the algorithm
    defined in the Rennolet data paper.

    Args:
        trace (obspy.core.trace.Trace): Trace of strong motion data.
        org_time (UTCDateTime): Event origin time.
        epi_dist (float): Distance from event epicenter to station.
        vmin (float): Minimum apparent velocity.
            Default is 1.0 km/s.

    Returns:
        trace (obspy.core.trace.Trace) after windowing.
        If trim failure occurs: -1
    """
    window = {'vmin': vmin}
    trace = _update_params(trace, 'window', window)

    end_time = org_time + max(120, epi_dist / vmin)
    # Check to make sure that the trimming end time is after
    # the start time of our trace
    if (end_time <= trace.stats.starttime):
        return -1
    else:
        trace.trim(endtime=end_time)
        return trace


def split_signal_and_noise(trace, event_time, epi_dist):
    """
    Identifies the noise and signal windows for the waveform.
    The noise window is defined from the start of the waveform through the
    arrival time of a 7 km/s phase. The signal window is defned from the
    arrival time of a 7 km/s phase through the end of the waveform.

    Args:
        trace (obspy.core.trace.Trace): Trace of strong motion data.
        event_time (UTCDateTime): Event origin time.
        epi_dist (float): Distance form event epicenter to station.

    Returns:
        tuple of two traces:
            1) Noise trace
            2) Signal trace
        If cannot separate noise from signal, returns -1.
    """

    # Calculate the arrival time of a 7 km/s phase
    phase_arrival_time = event_time + epi_dist / 7.0

    # Check if the arrival time is before or after the trace window
    if (phase_arrival_time <= trace.stats.starttime):
        return (-1, -1)
    elif (phase_arrival_time >= trace.stats.endtime):
        return (-1, -1)
    else:
        orig_trace_1 = trace.copy()
        orig_trace_2 = trace.copy()

    noise_trace = orig_trace_1.trim(endtime=phase_arrival_time)
    signal_trace = orig_trace_2.trim(starttime=phase_arrival_time)
    return (signal_trace, noise_trace)


def fft_smooth(trace, nfft):
    """
    Pads a trace to the nearest upper power of 2, takes the FFT, and
    smooths the amplitude spectra following the algorithm of
    Konno and Ohmachi.

    Args:
        trace (obspy.core.trace.Trace): Trace of strong motion data.
        nfft (int): Number of data points for the fourier transform.

    Returns:
        numpy.ndarray: Smoothed amplitude data and frequencies.
    """

    # Compute the FFT, normalizing by the number of data points
    spec = abs(np.fft.rfft(trace.data, n=nfft)) / nfft

    # Get the frequencies associated with the FFT
    freqs = np.fft.rfftfreq(nfft, 1 / trace.stats.sampling_rate)

    # Konno Omachi Smoothing using 20 for bandwidth parameter
    spec_smooth = konno_ohmachi_smoothing(spec.astype(float), freqs, 20)
    return spec_smooth, freqs


def get_corner_frequencies(trace, event_time, epi_dist, ratio=3.0,
                           max_low_corner=0.1, min_high_corner=5.0,
                           taper_type='hann', taper_percentage=0.05,
                           taper_side='both'):
    """
    Returns the corner frequencies for a trace.

    Args:
        trace (obspy.core.trace.Trace): Trace of strong motion data.
        event_time (UTCDateTime): Event origin time.
        epi_dist (float): Distance from event epicenter to station.
        ratio (float): Required signal-to-noise ratio.
            Default is 3.0.
        max_low_corner (float): Maxmimum low corner frequency allowed.
            Default is 0.1.
        min_high_corner(float): Minimum low corner frequency allowed.
            Default is 5.0.
        taper_type (str): Taper types allowed for filtering. Must be an Obspy
            supported tapering method.
            Default is 'hann' (Hanning taper).
        taper_percentage (float): Decimal percentage of taper.
            Default is 0.05 (5%).
        taper_side (str): Speicfy which sides should be tapered. Either 'left',
            'right', or 'both'.
            Default is 'both'.

    Returns:
        list : List of floats representing corner frequencies.
        Returns two -1 values if inadequate signal to noise ratio.
    """

    # Split the noise and signal into two separate traces
    signal, noise = split_signal_and_noise(trace, event_time, epi_dist)

    # Check if signal and noise splitting failed
    if (signal == -1 and noise == -1):
        return [-1, -1]

    # Taper the noise and signal traces
    noise = taper(noise, taper_type=taper_type,
                  max_percentage=taper_percentage, side=taper_side)

    signal = taper(signal, taper_type=taper_type,
                   max_percentage=taper_percentage, side=taper_side)

    # Find the number of points for the Fourier transform
    nfft = max(next_pow_2(signal.stats.npts), next_pow_2(noise.stats.npts))

    # Transform to frequency domain and smooth spectra using
    # konno-ohmachi smoothing
    sig_spec_smooth, freqs_signal = fft_smooth(signal, nfft)
    noise_spec_smooth, freqs_noise = fft_smooth(noise, nfft)

    # remove the noise level from the spectrum of the signal window
    sig_spec_smooth -= noise_spec_smooth

    # Loop through frequencies to find low corner and high corner
    corner_frequencies = []
    lows = []
    highs = []
    have_low = False
    for idx, freq in enumerate(freqs_signal):
        if have_low is False:
            if (sig_spec_smooth[idx] / noise_spec_smooth[idx]) >= ratio:
                lows.append(freq)
                have_low = True
            else:
                continue
        else:
            if (sig_spec_smooth[idx] / noise_spec_smooth[idx]) < ratio:
                highs.append(freq)
                have_low = False
            else:
                continue

    # If we didn't find any corners
    if not lows:
        return [-2, -2]

    # If we find an extra low, add another high for the maximum frequency
    if len(lows) > len(highs):
        highs.append(max(freqs_signal))

    # Check if any of the low/high pairs are valid
    found_valid = False
    for idx, val in enumerate(lows):
        if (val <= max_low_corner and highs[idx] > min_high_corner):
            low_corner = val
            high_corner = highs[idx]
            found_valid = True

    # Check if we found any valid pairs
    if not found_valid:
        return [-3, -3]
    else:
        corner_frequencies = [low_corner, high_corner]
        corners = {'get_dynamically': True,
                   'ratio': ratio,
                   'high_corner': corner_frequencies[1],
                   'low_corner': corner_frequencies[0]}
        trace = _update_params(trace, 'corners', corners)
        return corner_frequencies


def filter_waveform(trace, filter_type, freqmax=None, freqmin=None,
                    zerophase=False, corners=4):
    """
    Returns the filtered waveform using the provided corner frequencies.

    Args:
        trace (obspy.core.trace.Trace): Trace of strong motion data.
        filter_type (str): Type of filter to be applied.
        freqmax (float): Maximum frequency corner.
            Default is None.
        freqmin (float): Minimum frequency corner.
            Default is None.
        zerophase (bool): Whether to perform zerophase filtering.
            Default is False.
        corners (int): Number of corners (poles).

    Returns:
        trace (obspy.core.trace.Trace): Filtered trace
    """

    filter_type = filter_type.lower()
    two_freq = ['bandpass', 'bandstop']
    low_freq = ['lowpass', 'lowpass_cheby_2', 'lowpass_fir']
    high_freq = ['highpass']

    if filter_type in two_freq:
        trace.filter(filter_type, freqmin=freqmin, freqmax=freqmax,
                     corners=corners, zerophase=zerophase)
        filter_params = {'type': filter_type, 'corners': corners,
                         'zerophase': zerophase}
    elif filter_type in low_freq:
        # Only perform low pass frequency if corner is less than nyquist freq
        nyquist_freq = 0.5 * trace.stats.sampling_rate
        if (freqmax < nyquist_freq):
            trace.filter(filter_type, freq=freqmax,
                         corners=corners, zerophase=zerophase)
            filter_params = {'type': filter_type, 'corners': corners,
                             'zerophase': zerophase}
        else:
            warnings.warn('Low pass frequency is greater than or equal '
                          'to the Nyquist frequency. Low pass filter will '
                          'not be applied.')
            filter_params = None
    elif filter_type in high_freq:
        if freqmin != 0.0:
            trace.filter(filter_type, freq=freqmin,
                         corners=corners, zerophase=True)
            filter_params = {'type': filter_type, 'corners': corners,
                             'zerophase': zerophase}
        else:
            warnings.warn('High pass frequency is equal to zero. High pass '
                          'filter will not be applied')
            filter_params = None
    else:
        warnings.warn('Filter type not available %r. Available filters: '
                      '%r' % (filter_type, two_freq + low_freq + high_freq))
        filter_params = None
    try:
        if filter_params is not None:
            filters = trace.stats.processing_parameters['filters']
            filters += [filter_params]
            trace = _update_params(trace, 'filters', [filters])
    except KeyError:
        if filter_params is not None:
            trace = _update_params(trace, 'filters', [filter_params])
    return trace


def poly_func(x, a, b, c, d, e):
        """
        Model polynomial function for polynomial baseline correction.
        """
        return a*x**6 + b*x**5 + c*x**4 + d*x**3 + e*x**2


def correct_baseline(trace):
    """
    Performs a baseline correction following the method of Ancheta
    et al. (2013). This removes low-frequency, non-physical trends
    that remain in the time series following filtering.

    Args:
        trace (obspy.core.trace.Trace): Trace of strong motion data.

    Returns:
        trace (obspy.core.trace.Trace): Baseline-corrected trace.
    """

    # Make copies of the trace for our accleration data
    orig_trace = trace.copy()
    acc_trace = trace.copy()

    # Integrate twice to get the displacement time series
    disp_trace = (acc_trace.integrate()).integrate()

    # Fit a sixth order polynomial to displacement time series, requiring
    # that the 1st and 0th order coefficients are zero
    time_values = np.linspace(0, trace.stats.npts-1, trace.stats.npts)
    poly_cofs = list(curve_fit(poly_func, time_values, disp_trace.data)[0])
    poly_cofs += [0, 0]

    # Construct a polynomial from the coefficients and compute
    # the second derivative
    polynomial = np.poly1d(poly_cofs)
    polynomial_second_derivative = np.polyder(polynomial, 2)

    # Subtract the second derivative of the polynomial from the
    # acceleration trace
    for i in range(orig_trace.stats.npts):
        orig_trace.data[i] -= polynomial_second_derivative(i)
    orig_trace = _update_params(orig_trace, 'baseline_correct', True)
    return orig_trace


def process(stream, amp_min, amp_max, window_vmin, taper_type,
            taper_percentage, taper_side, get_corners, sn_ratio,
            max_low_freq, min_high_freq, default_low_frequency,
            default_high_frequency, filters, baseline_correct, event_time=None,
            epi_dist=None):
    """
    Processes an acceleration trace following the step-by-step process
    described in the Rennolet et al paper
    (https://doi.org/10.1193/101916EQS175DP)
    This function completes
    Step 4 through Step 11 from the paper.

    - Amplitude
    - Windowing
    - Identify corner frequencies for filtering
        - Split signal and noise windows
        - Taper, compute FFT, Konno Ohmachi smoothing
        - Signal-to-noise ratio
    - Filter
    - Polynomial baseline correction

    This processing should be performed on data with physical units (acc,
    vel, etc).

    Along with the stream, this function requires two pieces of metadata
    (origin time of the event and epicentral distance to the station) to
    complete all processing steps.

    To process with defaults use process_config with the amptools default
    config dictionary.

    Args:
        stream (obspy.core.stream.Stream): Stream for one station.
        amp_min (float): Lower amplitude limit for step 4.
        amp_max (float): Upper amplitude limit for step 4.
        window_vmin (float): Minimum velocity for step 5.
        taper_type (str): Type of taper for step 6.
        taper_percentage (float): Maximum taper percentage for step 6.
        taper_side (str): Sides to taper for step 6.
        get_corners (bool): Whether to complete step 8 or use defaults.
        sn_ratio (float): Signal to noise ratio.
        max_low_freq (float): Maxmium low corner frequency allowed.
        min_high_freq (float): Minimum high corner frequency allowed.
        default_low_frequency (float): Default minimum frequency used in
                place of corners calculated from step 8.
        default_high_frequency (float): Default maximum frequency used in
                place of corners calculated from step 8.
        filters (list): List of filters (dict) with type (str), corners (int),
                and zerophase (bool) defined.
        baseline_correct (bool): Whether or not to complete step 11.
        event_time (UTCDateTime): Origin time of the event. Default is None.
        epi_dist (float): Epicentral distance. Default is None.

    Returns:
        obspy.core.stream.Stream: Processed stream.
    """
    for idx, trace in enumerate(stream):
        trace_copy = trace.copy()
        # The stats need to be set even if the process checks fail
        trace_copy = _update_params(trace_copy, 'amplitude', {'min': amp_min,
                'max': amp_min})
        trace_copy = _update_params(trace_copy, 'window',
                {'vmin': window_vmin})
        trace_copy = _update_params(trace_copy, 'taper', {'type': taper_type,
                'side': taper_side, 'max_percentage': taper_percentage})
        trace_copy = _update_params(trace_copy, 'corners',
                {'get_dynamically': get_corners, 'sn_ratio': sn_ratio,
                'max_low_freq': max_low_freq, 'min_high_freq': min_high_freq,
                'default_low_frequency': default_low_frequency,
                'default_high_frequency': default_high_frequency})
        trace_copy = _update_params(trace_copy, 'filters', [])
        trace_copy = _update_params(trace_copy, 'baseline_correct',
                baseline_correct)
        trace_copy.stats['passed_tests'] = True

        # Check amplitude
        if not check_max_amplitude(trace_copy, amp_min, amp_max):
            trace_copy.stats['passed_tests'] = False
            err_msg = ('Processing: Trace maximum amplitude is not '
                       'within the acceptable range: %r to %r. Skipping '
                       'processing for trace: %r' % (amp_min,
                                                     amp_max, trace))
            trace_copy = _update_comments(trace_copy, err_msg)
            stream[idx] = trace_copy
            continue

        # Windowing
        if event_time is not None and epi_dist is not None:
            trace_trim = trim_total_window(trace_copy, event_time, epi_dist,
                                           vmin=window_vmin)
            windowed = True
            # Check if windowing failed
            if (trace_trim == -1):
                trace_copy.stats['passed_tests'] = False
                err_msg = ('Processing: Invalid time windowing. The start '
                           'time of the trace is after the calculated '
                           'end time. Skipping procesing for trace: %r', trace)
                trace_trim = _update_comments(trace_copy, err_msg)
                stream[idx] = trace_copy
                continue
        else:
            trace_copy.stats['passed_tests'] = False
            err_msg = ('Processing: No windowing test performed. Missing '
                       'event time and/or epicentral distance information to '
                       'perform calculation.')
            trace_copy = _update_comments(trace_copy, err_msg)
            trace_copy = _update_params(trace_copy, 'window',
                                        {'vmin': window_vmin})
            trace_trim = trace_copy
            # Corners cannot be calculated dynamically without windowing
            warnings.warn('Missing event information. Continuing processing '
                          'without windowing. Default frequencies will be '
                          'used for filtering.')
            windowed = False

        # Find corner frequencies
        if get_corners and windowed:
            corners = get_corner_frequencies(trace_trim, event_time,
                    epi_dist, sn_ratio, max_low_freq, min_high_freq,
                    taper_type, taper_percentage, taper_side)

            if (corners[0] < 0 or corners[1] < 0):
                trace_copy.stats['passed_tests'] = False
                dynamic = False
                low_freq = -9999
                high_freq = -9999

                if corners == [-1, -1]:
                    err_msg = ('Not enough pre-event noise to calculate '
                               'signal to noise ratio. Skipping processing '
                               'for trace: %r' % (trace))
                elif corners == [-2, -2]:
                    err_msg = ('Signal-to-noise ratio too low to find corner '
                               'frequencies, skipping processing for '
                               'trace: %r' % (trace))
                else:
                    err_msg = ('Did not find any corner frequencies within '
                               'the valid bandwidth. Skipping processing for '
                               'trace: %r' % (trace))
                trace_copy = _update_comments(trace_copy, err_msg)

            else:
                low_freq = corners[0]
                high_freq = corners[1]
                dynamic = True
        else:
            high_freq = default_high_frequency
            low_freq = default_low_frequency
            dynamic = False

        corner_params = {'get_dynamically': dynamic,
                         'default_high_frequency': default_high_frequency,
                         'sn_ratio': sn_ratio,
                         'default_low_frequency': default_low_frequency,
                         'max_low_freq': max_low_freq,
                         'min_high_freq': min_high_freq,
                         'low_corner': low_freq,
                         'high_corner': high_freq}

        if trace_copy.stats['passed_tests'] is False:
            trace_copy = _update_params(trace_copy, 'corners', corner_params)
            stream[idx] = trace_copy
            continue

        trace_trim = _update_params(trace_trim, 'corners', corner_params)

        # Filter
        trace_filt = trace_trim
        for filter_dict in filters:
            filter_type = filter_dict['type']
            corners = filter_dict['corners']
            zerophase = filter_dict['zerophase']
            trace_filt = filter_waveform(trace_filt, filter_type,
                                         high_freq, low_freq, zerophase,
                                         corners)

        # Correct baseline
        if baseline_correct:
            trace_cor = correct_baseline(trace_filt)
        else:
            trace_cor = _update_params(trace_filt, 'baseline_correct', False)
        stream[idx] = trace_cor
    return stream


def process_config(stream, config=None, event_time=None, epi_dist=None):
    """Implement the process method according to a config file.

    Args:
        stream (obspy.core.stream.Stream): Stream for one station.
        config (dictionary): Config dictionary. Default is None.
            When config is None, it will be set to the amptools config.
        event_time (UTCDateTime): Origin time of the event. Default is None.
        epi_dist (float): Epicentral distance. Default is None.

    Notes:
        This function looks for a config file ~/.amptools/config.yml, which
        has a processing_parameters subdictionary which should be formatted
        as follows:
        processing_parameters:
            amplitude:
                min: <float>
                max: <float>
            window:
                vmin: <float>
            taper:
                type: <str>
                max_percentage:<float>
                side: <str>
            corners:
                get_dynamically: <bool>
                default_low_frequency: <float>
                default_high_frequency: <float>
            filters:
                - type: <str>
                  corners: <int>
                  zerophase: <bool>
                    ...
            baseline_correct: <bool>
        If no config file is available, the default config is used. Other
        packages (such as strongmotion-database) that use amptools as a
        dependency have the option of overriding the config dictionary.
    """
    if config is None:
        config = CONFIG
    # Set all float/int params in the correct format
    params = config['processing_parameters']
    amp_min = float(params['amplitude']['min'])
    amp_max = float(params['amplitude']['max'])
    window_vmin = float(params['window']['vmin'])
    taper_type = params['taper']['type']
    taper_percentage = float(params['taper']['max_percentage'])
    taper_side = params['taper']['side']
    get_corners = params['corners']['get_dynamically']
    sn_ratio = params['corners']['sn_ratio']
    default_low_frequency = float(params['corners']['default_low_frequency'])
    default_high_frequency = float(params['corners']['default_high_frequency'])
    max_low_freq = float(params['corners']['max_low_freq'])
    min_high_freq = float(params['corners']['min_high_freq'])
    filters = params['filters']
    baseline_correct = params['baseline_correct']

    corrected_stream = process(stream, amp_min, amp_max, window_vmin,
            taper_type, taper_percentage, taper_side, get_corners, sn_ratio,
            max_low_freq, min_high_freq, default_low_frequency,
            default_high_frequency, filters, baseline_correct,
            event_time=event_time, epi_dist=epi_dist)
    return corrected_stream


def taper(trace, taper_type='hann', max_percentage=0.05, side='both'):
    """
    Taper a stream and record processing step.

    Args:
        trace (obspy.core.trace.Trace): Trace of strong motion data.
        taper_type (str): Type of taper. Default is 'hann'.
        max_percentage (float): Taper percentage at one end. Default is
            5% (0.05).
        side (str): Which side will be tapered. Default is 'both'.

    Returns:
        obspy.core.trace.Trace: Trace of tapered data.
    """
    trace.taper(type=taper_type,
                max_percentage=max_percentage, side=side)
    taper_params = {
            'type': taper_type,
            'max_percentage': max_percentage,
            'side': side
    }
    trace = _update_params(trace, 'taper', taper_params)
    return trace


def _update_comments(trace, comment):
    """
    Helper function to update a trace's stats.

    Args:
        trace (obspy.core.trace.Trace): Trace of strong motion dataself.
        comment (str): Comment to add to the trace's stats.

    Returns:
        obspy.core.trace.Trace: Trace with updated stats.
    """
    if 'standard' not in trace.stats:
        trace.stats['standard'] = {}
        trace.stats.standard['comments'] = ''
    if (isinstance(trace.stats.standard['comments'], str) and
            trace.stats.standard['comments'] != ''):
        trace.stats.standard['comments'] += ', ' + comment
    else:
        trace.stats.standard['comments'] = comment
    return trace


def _update_params(trace, process_type, parameters):
    """
    Helper function to update a trace's processing_parameters.

    Args:
        trace (obspy.core.trace.Trace): Trace of strong motion dataself.
        process_type (str): Key for processing_parameters subdictionary.
        parameters (dict or list): Parameters for the given key.

    Returns:
        obspy.core.trace.Trace: Trace with updated processing_parameters.
    """
    if 'processing_parameters' not in trace.stats:
        trace.stats['processing_parameters'] = {}
    trace.stats['processing_parameters'][process_type] = parameters
    return trace
