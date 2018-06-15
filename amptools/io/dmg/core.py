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

V2_TEXT_HDR_ROWS = 25
V2_INT_HDR_ROWS = 7
V2_INT_FMT = [5] * 16
V2_REAL_HDR_ROWS = 13
V2_REAL_FMT = [10] * 8

VALID_MARKERS = [
    'UNCORRECTED ACCELEROGRAM',
    'CORRECTED ACCELEROGRAM',
    'RESPONSE AND FOURIER AMPLITUDE SPECTRA'
]

homedir = os.path.dirname(os.path.abspath(__file__))
codedir = os.path.join(homedir, '..', 'fdsn_codes.txt')
VALID_CODES = np.loadtxt(codedir, skiprows=1, usecols=0, dtype=str)

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
    line = open(filename, 'rt').readline()
    for marker in VALID_MARKERS:
        if line.lower().find(marker.lower()) >= 0:
            return True
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

    # Parse name and code information
    name_length = int_data[29]
    name = re.sub(' +', ' ', lines[6][:name_length]).strip().replace(' ', '_')
    hdr['name'] = name
    code = re.sub(' +', ' ', lines[1][name_length:]).strip().split(' ')[-1][:2]
    if code.upper() in VALID_CODES:
        hdr['network'] = code.upper()
    else:
        hdr['network'] = 'UNK'

    # set statistics
    hdr['units'] = 'acc'
    lat = lines[5][21:27].replace(' ', '')
    if lat[-1].upper() == 'S':
        lat = -1 * float(lat[0:-1])
    lon = lines[5][30:37].replace(' ', '')
    if lon[-1].upper() == 'W':
        lon = -1 * float(lon[0:-1])
    hdr['lat'] = lat
    hdr['lon'] = lon
    hdr['location'] = '--'
    try:
        datestr = lines[4][36:80].lower().replace(' ', '')
        date = datestr.split(',')[0].split(':')[-1]
        time = datestr.split(',')[-1]
        time = time[:time.find('.') + 2]
        month = int(date.split('/')[0])
        day = int(date.split('/')[1])
        year = int(date.split('/')[2])
        if int_data[23] != 0:
            year = int_data[23]
        hour = int(time.split(':')[0])
        minute = int(time.split(':')[1])
        seconds = float(time.split(':')[2])
        microseconds = int((seconds - int(seconds)) * 1e6)
    except ValueError:
        message = "Incorrectformat for trigger time.\n" +\
            "Correct format examle:\n" + \
            "TRIGGER TIME: 05/27/80, 14:51:00.9 GMT\n" +\
            "Attempting to retreive time from integer header."
        warnings.warn(message, Warning)
        month = int_data[21]
        day = int_data[22]
        year = int_data[23]
        hour = int_data[16]
        minute = int_data[17]
        seconds = int_data[18]
        microseconds = 0
        if (
            year == 0 and month == 0 and day == 0 and hour == 0 and
            minute == 0 and seconds == 0
        ):
            raise Exception('Missing start time data from integer header.')
    hdr['starttime'] = datetime(
        year, month, day, hour, minute, int(seconds), microseconds)
    hdr['delta'] = flt_data[60]
    hdr['sampling_rate'] = 1 / hdr['delta']
    hdr['npts'] = int_data[52]
    hdr['source'] = hdr['network']
    angle = int_data[26]
    if angle == 500 or angle == 600:
        hdr['channel'] = 'HHZ'
    elif angle > 315 or angle < 45 or (angle > 135 and angle < 225):
        hdr['channel'] = 'HHN'
    else:
        hdr['channel'] = 'HHE'

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
     - network (str): Default is '--'. Determined using COSMOS_NETWORKS
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
    station_name = re.sub(' +',' ',lines[6][:name_length]).strip()
    code = re.sub(' +',' ',lines[1][name_length:]).strip().split(' ')[-1][:2]
    if code.upper() in CODES:
        network = code.upper()
        idx = np.argwhere(CODES == network)[0][0]
        source = SOURCES1[idx].decode('utf-8') +  ', ' + SOURCES2[idx].decode('utf-8')
    else:
        network = '--'
        source = ''
    hdr['network'] = network
    station_line = lines[5]
    station = station_line[12:17].strip()
    hdr['station'] = station
    angle = int_data[26]
    if angle == 500 or angle == 600 or (angle >= 0 and angle <= 360):
        if angle == 500 or angle == 600:
            hdr['channel'] = 'Z'
        elif angle > 315 or angle < 45 or (angle > 135 and angle < 225):
            hdr['channel'] = 'H1'
        else:
            hdr['channel'] = 'H2'
    else:
        warnings.warn('Missing channel orientation. Setting channel to '
                      'channel number. Setting horizontal_orientation to '
                      'np.nan.', Warning)
        angle = np.nan
        # Recorder channel number
        hdr['channel'] = int_data[0]
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
            # Looking for full year on earthquake line
            try:
                month = int(date[0][-2:])
                day = int(date[1])
                time = trigger_line.split(':')
                hour = int(time[1][-2:])
                minute = int(time[2])
                second = float(time[3][:2])
                microsecond = int((second - int(second)) * 1e6)
                earthquake_year_line = lines[8]
                date = earthquake_year_line.split(',')[1].split(':')[0].strip()
                year = int(date.split(' ')[0])
                hdr['starttime'] = datetime(year, month, day, hour, minute,
                       int(second), microsecond)
            except ValueError:
                 raise AmptoolsException(starttime_message)
    else:
        raise AmptoolsException(starttime_message)
    hdr['delta'] = flt_data[60]
    hdr['sampling_rate'] = 1 / hdr['delta']
    hdr['npts'] = int_data[52]

    # Coordinates
    latitude_str = station_line[20:27].strip()
    longitude_str = station_line[29:37].strip()
    try:
        latitude = float(latitude_str[:-1])
        if latitude_str.upper().find('S') >= 0:
            latitude = -1 * latitude
    except:
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
    format_specific['fractional_unit'] =  flt_data[4]
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
