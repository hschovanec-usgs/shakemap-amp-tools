#!/usr/bin/env python

# stdlib imports
from datetime import datetime
import os
import re
import warnings

# third party
from obspy.core.trace import Trace
from obspy.core.stream import Stream
from obspy.core.trace import Stats
import numpy as np

# local imports
from amptools.exception import AmptoolsException
from amptools.io.usc.core import is_usc
from amptools.io.seedname import get_channel_name

V2_TEXT_HDR_ROWS = 25
V2_INT_HDR_ROWS = 7
V2_INT_FMT = [5] * 16
V2_REAL_HDR_ROWS = 13
V2_REAL_FMT = [10] * 8

V1_MARKER = 'UNCORRECTED ACCELEROGRAM DATA'
V2_MARKER = 'CORRECTED ACCELEROGRAM'
V3_MARKER = 'RESPONSE AND FOURIER AMPLITUDE SPECTRA'

homedir = os.path.dirname(os.path.abspath(__file__))
codedir = os.path.join(homedir, '..', 'fdsn_codes.csv')
CODES, SOURCES1, SOURCES2 = np.genfromtxt(codedir, skip_header=1, usecols=(0, 1, 2),
                                          unpack=True, dtype=bytes, delimiter=',')
CODES = CODES.astype(str)

UNITS = [
    'acc',
    'vel',
    'disp'
]


def is_dmg(filename):
    """Check to see if file is a DMG strong motion file.

    Notes:
        CSMIP is synonymous to as DMG in this reader.

    Args:
        filename (str): Path to possible COSMOS V0/V1 data file.
    Returns:
        bool: True if DMG , False otherwise.
    """
    try:
        f = open(filename, 'rt')
        first_line = f.readline().upper()
        second_line = f.readline().upper()
        third_line = f.readline().upper()
        f.close()

        # dmg/csmip both have the same markers so is_usc must be checked
        if first_line.find(V1_MARKER) >= 0 and not is_usc(filename):
            return True
        elif first_line.find(V2_MARKER) >= 0 and not is_usc(filename):
            if second_line.find(V1_MARKER) >= 0:
                return True
        elif first_line.find(V3_MARKER) >= 0 and not is_usc(filename):
            if second_line.find(V2_MARKER) >= 0 and third_line.find(V1_MARKER) >= 0:
                return True
        else:
            return False
    except UnicodeDecodeError:
        return False


def read_dmg(filename, **kwargs):
    """Read DMG strong motion file.

    Notes:
        CSMIP is synonymous to as DMG in this reader.

    Args:
        filename (str): Path to possible DMG data file.
        kwargs (ref):
            units (str): String determining which timeseries is return. Valid
                    options include 'acc', 'vel', 'disp'. Default is 'acc'.
            Other arguments will be ignored.
    Returns:
        Stream: Obspy Stream containing three channels of acceleration data (cm/s**2).
    """
    if not is_dmg(filename):
        raise Exception('Not a DMG file format.')

    # Check for type
    units = kwargs.get('units', 'acc')
    if units not in UNITS:
        raise Exception('Not a valid choice of units.')

    # Check for DMG format and determine volume type
    line = open(filename, 'rt').readline()
    if is_dmg(filename):
        if line.lower().find('uncorrected') >= 0:
            reader = 'V1'
        elif line.lower().find('corrected') >= 0:
            reader = 'V2'
        elif line.lower().find('response') >= 0:
            reader = 'V3'

    # Count the number of lines in the file
    with open(filename) as f:
        line_count = sum(1 for _ in f)

    # Read as many channels as are present in the file
    line_offset = 0
    trace_list = []
    while line_offset < line_count:
        if reader == 'V2':
            traces, line_offset = _read_volume_two(filename, line_offset)
            trace_list += traces
        else:
            line_offset = line_count

    stream = Stream([])
    for trace in trace_list:
        if trace.stats['standard']['units'] == units:
            stream.append(trace)
    return stream


def _read_volume_two(filename, line_offset):
    """Read channel data from DMG text file.

    Args:
        filename (str): Input DMG V2 filename.
        line_offset (int): Line offset to beginning of channel text block.
    Returns:
        tuple: (list of obspy Trace, int line offset)
    """
    # read station, location, and process level from text header
    with open(filename, 'rt') as f:
        for _ in range(line_offset):
            next(f)
        lines = [next(f) for x in range(V2_TEXT_HDR_ROWS)]

    # read in lines of integer data
    skip_rows = V2_TEXT_HDR_ROWS + line_offset
    int_data = _read_lines(skip_rows, V2_INT_HDR_ROWS, V2_INT_FMT, filename)
    int_data = int_data[0:100].astype(np.int32)

    # read in lines of float data
    skip_rows += V2_INT_HDR_ROWS
    flt_data = _read_lines(skip_rows, V2_REAL_HDR_ROWS, V2_REAL_FMT, filename)
    flt_data = flt_data[:100]
    skip_rows += V2_REAL_HDR_ROWS

    # according to the powers that defined the Network.Station.Channel.Location
    # "standard", Location is a two character field.  Most data providers,
    # including csmip/dmg here, don't always provide this.  We'll flag it as "--".
    hdr = _get_header_info(int_data, flt_data, lines, 'V2')

    traces = []
    # read acceleration data
    if hdr['npts'] > 0:
        acc_rows, acc_fmt = _get_data_format(filename, skip_rows, hdr['npts'])
        acc_data = _read_lines(skip_rows + 1, acc_rows, acc_fmt, filename)
        acc_data = acc_data[:hdr['npts']]
        acc_trace = Trace(acc_data.copy(), Stats(hdr.copy()))
        traces += [acc_trace]
        skip_rows += int(acc_rows) + 1

    # read acceleration data
    vel_hdr = hdr.copy()
    vel_hdr['standard']['units'] = 'vel'
    vel_hdr['npts'] = int_data[63]
    if vel_hdr['npts'] > 0:
        vel_rows, vel_fmt = _get_data_format(
            filename, skip_rows, vel_hdr['npts'])
        vel_data = _read_lines(skip_rows + 1, vel_rows, vel_fmt, filename)
        vel_data = vel_data[:vel_hdr['npts']]
        vel_trace = Trace(vel_data.copy(), Stats(vel_hdr.copy()))
        traces += [vel_trace]
        skip_rows += int(vel_rows) + 1

    # read displacement data
    disp_hdr = hdr.copy()
    disp_hdr['standard']['units'] = 'disp'
    disp_hdr['npts'] = int_data[65]
    if disp_hdr['npts'] > 0:
        disp_rows, disp_fmt = _get_data_format(
            filename, skip_rows, disp_hdr['npts'])
        disp_data = _read_lines(skip_rows + 1, disp_rows, disp_fmt, filename)
        disp_data = disp_data[:disp_hdr['npts']]
        disp_trace = Trace(disp_data.copy(), Stats(disp_hdr.copy()))
        traces += [disp_trace]
        skip_rows += int(disp_rows) + 1
    new_offset = skip_rows + 1  # there is an 'end of record' line after the data]
    return (traces, new_offset)


def _get_header_info(int_data, flt_data, lines, level):
    """Return stats structure from various headers.

    Output is a dictionary like this:
     - network (str): Default is 'ZZ'. Determined using COSMOS_NETWORKS
     - station (str)
     - channel (str)
     - location (str): Default is '--'
     - starttime (datetime)
     - sampling_rate (float)
     - delta (float)
     - coordinates:
       - latitude (float)
       - longitude (float)
       - elevation (float): Default is np.nan
    - standard (Defaults are either np.nan or '')
      - horizontal_orientation (float): Rotation from north (degrees)
      - instrument_period (float): Period of sensor (Hz)
      - instrument_damping (float): Fraction of critical
      - process_time (datetime): Reported date of processing
      - process_level: Either 'V0', 'V1', 'V2', or 'V3'
      - station_name (str): Long form station description
      - sensor_serial_number (str): Reported sensor serial
      - instrument (str)
      - comments (str): Processing comments
      - structure_type (str)
      - corner_frequency (float): Sensor corner frequency (Hz)
      - units (str)
      - source (str): Network source description
      - source_format (str): Always dmg
    - format_specific
      - sensor_sensitivity (float): Transducer sensitivity (cm/g)
      - time_sd (float): Standard deviaiton of time steop (millisecond)
      - fractional_unit (float): Units of digitized acceleration
            in file (fractions of g)
      - scaling_factor (float): Scaling used for converting acceleration
            from g/10 to cm/sec/sec
      - low_filter_corner (float): Filter corner for low frequency
            V2 filtering (Hz)
      - high_filter_corner (float): Filter corner for high frequency
            V2 filtering (Hz)

    Args:
        int_data (ndarray): Array of integer data
        flt_data (ndarray): Array of float data
        lines (list): List of text headers (str)
        level (str): Process level code (V0, V1, V2, V3)

    Returns:
        dictionary: Dictionary of header/metadata information
    """
    hdr = {}
    coordinates = {}
    standard = {}
    format_specific = {}

    # Required metadata
    name_length = int_data[29]
    station_name = re.sub(' +', ' ', lines[6][:name_length]).strip()
    code = re.sub(' +', ' ', lines[1][name_length:]).strip().split(' ')[-1][:2]
    if code.upper() in CODES:
        network = code.upper()
        idx = np.argwhere(CODES == network)[0][0]
        source = SOURCES1[idx].decode(
            'utf-8') + ', ' + SOURCES2[idx].decode('utf-8')
    else:
        network = 'ZZ'
        source = ''
    hdr['network'] = network
    station_line = lines[5]
    station = station_line[12:17].strip()
    hdr['station'] = station
    angle = int_data[26]

    hdr['delta'] = flt_data[60]
    hdr['sampling_rate'] = 1 / hdr['delta']

    if angle == 500 or angle == 600 or (angle >= 0 and angle <= 360):
        if angle == 500 or angle == 600:
            hdr['channel'] = get_channel_name(hdr['sampling_rate'],
                                              is_acceleration=True,
                                              is_vertical=True,
                                              is_north=False)
        elif angle > 315 or angle < 45 or (angle > 135 and angle < 225):
            hdr['channel'] = get_channel_name(hdr['sampling_rate'],
                                              is_acceleration=True,
                                              is_vertical=False,
                                              is_north=True)
        else:
            hdr['channel'] = get_channel_name(hdr['sampling_rate'],
                                              is_acceleration=True,
                                              is_vertical=False,
                                              is_north=False)
    else:
        errstr = ('Not enough information to distinguish horizontal from '
                  'vertical channels.')
        raise AmptoolsException(errstr)

    hdr['location'] = '--'
    trigger_line = lines[4][35:77]
    starttime_message = "Start time must be specified with " + \
        "full date on header line 5, columns 36=77 " + \
        "(Example: 'Start time: 3/29/2014, 04:09:34.0 UTC')"
    if trigger_line.find('-') >= 0 or trigger_line.find('/') >= 0:
        if trigger_line.find('-') >= 0:
            delimeter = '-'
        elif trigger_line.find('/') >= 0:
            delimeter = '/'
        date = trigger_line.split(delimeter)
        try:
            month = int(date[0][-2:])
            day = int(date[1])
            time = trigger_line.split(':')
            hour = int(time[1][-2:])
            minute = int(time[2])
            second = float(time[3][:2])
            microsecond = int((second - int(second)) * 1e6)
            year = int(date[2][:4])
            hdr['starttime'] = datetime(year, month, day, hour, minute,
                                        int(second), microsecond)
        except ValueError:
            # Looking for full year in integer header
            try:
                month = int(date[0][-2:])
                day = int(date[1])
                time = trigger_line.split(':')
                hour = int(time[1][-2:])
                minute = int(time[2])
                second = float(time[3][:2])
                microsecond = int((second - int(second)) * 1e6)
                year = int_data[23]
                hdr['starttime'] = datetime(year, month, day, hour, minute,
                                            int(second), microsecond)
            except ValueError:
                raise AmptoolsException(starttime_message)
    else:
        raise AmptoolsException(starttime_message)

    hdr['npts'] = int_data[52]

    # Coordinates
    latitude_str = station_line[20:27].strip()
    longitude_str = station_line[29:37].strip()
    try:
        latitude = float(latitude_str[:-1])
        if latitude_str.upper().find('S') >= 0:
            latitude = -1 * latitude
    except Exception:
        warnings.warn('No latitude or invalid latitude format provided. '
                      'Setting to np.nan.', Warning)
        latitude = np.nan
    try:
        longitude = float(longitude_str[:-1])
        if longitude_str.upper().find('W') >= 0:
            longitude = -1 * longitude
    except:
        warnings.warn('No longitude or invalid longitude format provided.',
                      'Setting to np.nan.', Warning)
        longitude = np.nan
    coordinates['latitude'] = latitude
    coordinates['longitude'] = longitude
    coordinates['elevation'] = np.nan

    # Standard metadata
    standard['horizontal_orientation'] = angle
    standard['instrument_period'] = flt_data[0]
    standard['instrument_damping'] = flt_data[1]
    process_line = lines[1].lower()
    if process_line.find('processed:') >= 0:
        process_info = process_line[process_line.find('processed:'):]

        try:
            process_info = process_line[process_line.find('processed:'):]
            date = process_info.split('/')
            month = int(date[0][-2:])
            day = int(date[1])
            try:
                process_year = int(date[2][:4])
            except:
                process_year = date[2][:2]
            if len(process_year) == 2 and str(process_year) == str(year)[-2:]:
                process_year = year
            standard['process_time'] = datetime(process_year, month, day)
        except:
            standard['process_time'] = ''
    else:
        standard['process_time'] = ''
    standard['process_level'] = level
    standard['comments'] = ''
    standard['structure_type'] = lines[7][46:80].strip()
    standard['instrument'] = station_line[39:47].strip()
    standard['sensor_serial_number'] = station_line[53:57].strip()
    standard['corner_frequency'] = ''
    standard['units'] = 'acc'
    standard['source'] = source
    standard['source_format'] = 'dmg'
    standard['station_name'] = station_name

    # Format specific metadata
    format_specific['fractional_unit'] = flt_data[4]
    format_specific['sensor_sensitivity'] = flt_data[5]
    if flt_data[13] == -999:
        format_specific['time_sd'] = np.nan
    else:
        format_specific['time_sd'] = flt_data[13]
    format_specific['scaling_factor'] = flt_data[51]
    format_specific['low_filter_corner'] = flt_data[61]
    format_specific['high_filter_corner'] = flt_data[72]
    # Set dictionary
    hdr['coordinates'] = coordinates
    hdr['standard'] = standard
    hdr['format_specific'] = format_specific
    return hdr


def _read_lines(skip_rows, max_rows, widths, filename):
    """Read lines of headers and.

    Args:
        skip_rows (int): Number of rows to skip.
        filename (str): Path to possible DMG data file.
    Returns:
        array-like: List of comments or array of data.
    """
    data_arr = np.genfromtxt(filename, skip_header=skip_rows,
                             max_rows=max_rows, dtype=np.float64,
                             delimiter=widths).flatten()
    return data_arr


def _get_data_format(filename, skip_rows, npts):
    """Read data header and return the format.

    Args:
        skip_rows (int): Number of rows to skip.
        filename (str): Path to possible DMG data file.
        npts (int): Number of data points.
    Returns:
        tuple: (int number of rows, list list of widths).
    """
    fmt = np.genfromtxt(filename, skip_header=skip_rows,
                        max_rows=1, dtype=str)[-1]

    # Check for a format in header or use default
    if fmt.find('f') >= 0 and fmt.find('(') >= 0 and fmt.find(')') >= 0:
        fmt = fmt.replace('(', '').replace(')', '')
        cols = int(fmt.split('f')[0])
        widths = int(fmt.split('f')[-1].split('.')[0])
    else:
        cols = 8
        widths = 10
    fmt = [widths] * cols
    rows = np.ceil(npts/cols)
    return (rows, fmt)
