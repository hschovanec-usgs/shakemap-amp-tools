#!/usr/bin/env python

# stdlib imports
from datetime import datetime

# third party
from obspy.core.trace import Trace
from obspy.core.stream import Stream
import numpy as np

# local imports
from amptools.io.seedname import get_channel_name

ASCII_HEADER_LINES = 11
INTEGER_HEADER_LINES = 6
FLOAT_HEADER_LINES = 10
INT_HEADER_WIDTHS = 10
FLOAT_HEADER_WIDTHS = 15
DATA_COLUMNS = 8
FLOAT_DATA_WIDTHS = 10

VALID_HEADERS = {'1 UNCORRECTED ACCELEROGRAM': 'V1',
                 '2 CORRECTED ACCELEROGRAM': 'V2'}

INSTRUMENTS = {2: 'Sprengnether SA-3000 3-component fba',
               30: 'Kinemetrics FBA-13 3-component fba',
               31: 'Kinemetrics FBA-11 1-component fba',
               101: 'SMA-1',
               125: 'Kinemetrics FBA-23',
               102: 'C&GS Standard',
               126: 'Kinemetrics Episensor',
               103: 'AR-240',
               127: 'Kinemetrics FBA-4g',
               104: 'RFT-250',
               128: 'Kinemetrics FBA-2g',
               105: 'RFT-350',
               129: 'Kinemetrics FBA-1g',
               106: 'MO-2',
               130: 'Kinemetrics FBA-0.5g',
               107: 'RMT-280',
               131: 'Kinemetrics FBA-0.25g',
               108: 'SMA-2/3',
               132: 'Kinemetrics FBA-0.1g',
               109: 'DSA-1/DSA-3	133 WR1',
               110: 'DCA-300',
               134: 'S6000',
               111: 'DCA-333',
               135: 'Mark Products L22',
               112: 'A-700',
               136: 'Products L4C',
               113: 'SSA-1',
               137: 'CMG3',
               114: 'CRA-1',
               138: 'CMG3T',
               115: 'MO-2',
               139: 'CMG40T',
               116: 'FBA-3',
               140: 'CMG5',
               117: 'SMA-2',
               141: 'KS-2000',
               118: 'DCA-310',
               900: 'custom instrument',
               119: 'FBA-13',
               1302: 'Reftek Model 130-ANSS/02',
               120: 'SSA-2',
               121: 'SSR-1',
               122: 'BIDRA',
               123: 'CR-1',
               124: 'PDR-1'}

CONSTRUCTION_TYPES = {1: 'Reinforced concrete gravity',
                      2: 'Reinforced concrete arch',
                      3: 'earth fill',
                      4: 'other'}

STRUCTURES = {1: 'building',
              2: 'bridge',
              3: 'dam',
              4: 'other',
              np.nan: 'not a structure'}

BRIDGE_LOCATIONS = {0: 'free field',
                    1: 'at the base of a pier or abutment',
                    2: 'on an abutment',
                    3: 'on the deck at the top of a pier',
                    4: 'on the deck between piers or between an abutment and a pier'}

DAM_LOCATIONS = {0: 'upstream or downstream free field',
                 1: 'at the base of the dam',
                 2: 'on the crest of the dam',
                 3: 'on the abutment of the dam'}


def is_smc(filename):
    """Check to see if file is a SMC (corrected, in acc.) strong motion file.

    Args:
        filename (str): Path to possible SMC corrected data file.
    Returns:
        bool: True if SMC, False otherwise.
    """
    try:
        line = open(filename, 'rt').readline().strip()
        if line in VALID_HEADERS:
            return True
        return False
    except UnicodeDecodeError:
        return False


def read_smc(filename, **kwargs):
    """Read SMC strong motion file.

    Args:
        filename (str): Path to possible SMC data file.
        kwargs (ref):
            any_structure (bool): Read data from any type of structure,
                raise Exception if False and structure type is not free-field.
            Other arguments will be ignored.
    Returns:
        Stream: Obspy Stream containing one channel of acceleration data (cm/s**2).
    """
    any_structure = kwargs.get('any_structure', False)
    if not is_smc(filename):
        raise Exception('Not an SMC file.')

    stats, num_comments = _get_header_info(
        filename, any_structure=any_structure)

    skip = ASCII_HEADER_LINES + INTEGER_HEADER_LINES + \
        num_comments + FLOAT_HEADER_LINES

    # read float data (8 columns per line)
    nrows = int(np.floor(stats['npts']/DATA_COLUMNS))
    data = np.genfromtxt(filename,
                         max_rows=nrows,
                         skip_header=skip,
                         delimiter=FLOAT_DATA_WIDTHS)
    data = data.flatten()
    if stats['npts'] % DATA_COLUMNS:
        lastrow = np.genfromtxt(filename, max_rows=1,
                                skip_header=skip+nrows,
                                delimiter=FLOAT_DATA_WIDTHS)
        data = np.append(data, lastrow)

    trace = Trace(data, header=stats)
    stream = Stream([trace])
    return stream


def _get_header_info(filename, any_structure=False):
    """Return stats structure from various headers.

    Output is a dictionary like this:
     - network
     - station
     - channel
     - location
     - starttime
     - sampling_rate
     - npts
     - coordinates:
       - latitude
       - longitude
       - elevation
    - standard
      - horizontal_orientation
      - instrument_period
      - instrument_damping
      - process_level
      - station_name
      - sensor_serial_number
      - instrument
      - comments
      - structure_type
      - corner_frequency
      - units
      - source
      - source_format
    - format_specific
      - vertical_orientation
      - building_floor (0=basement, 1=floor above basement, -1=1st sub-basement, etc.
      - bridge_number_spans
      - bridge_transducer_location ("free field",
                                    "at the base of a pier or abutment",
                                    "on an abutment",
                                    "on the deck at the top of a pier"
                                    "on the deck between piers or between an abutment and a pier."
        dam_transducer_location ("upstream or downstream free field",
                                 "at the base of the dam",
                                 "on the crest of the dam",
                                 on the abutment of the dam")
        construction_type ("Reinforced concrete gravity",
                           "Reinforced concrete arch",
                           "earth fill",
                           "other")

        filter_poles
        data_source
    """
    stats = {}
    standard = {}
    format_specific = {}
    coordinates = {}
    # read the ascii header lines
    with open(filename) as f:
        ascheader = [next(f).strip() for x in range(ASCII_HEADER_LINES)]

    standard['process_level'] = VALID_HEADERS[ascheader[0]]

    # station code is in the third line
    stats['station'] = ''
    if len(ascheader[2]) >= 4:
        stats['station'] = ascheader[2][0:4]

    standard['process_time'] = ''
    standard['station_name'] = ascheader[5][10:40].strip()
    # sometimes the data source has nothing in it,
    # most of the time it seems has has USGS in it
    # sometimes it's something like JPL/USGS, CDOT/USGS, etc.
    # if it's got USGS in it, let's just say network=US, otherwise "--"
    stats['network'] = 'ZZ'
    if ascheader[7].find('USGS') > -1:
        stats['network'] = 'US'

    standard['source'] = ascheader[7].split('=')[2].strip()
    standard['source_format'] = 'smc'

    # read integer header data
    intheader = np.genfromtxt(filename, dtype=np.int32,
                              max_rows=INTEGER_HEADER_LINES,
                              skip_header=ASCII_HEADER_LINES,
                              delimiter=INT_HEADER_WIDTHS)
    # 8 columns per line
    # first line is start time information, and then inst. serial number
    missing_data = intheader[0, 0]
    year = intheader[0, 1]
    jday = intheader[0, 2]
    hour = intheader[0, 3]
    minute = intheader[0, 4]
    second = intheader[0, 5]
    microsecond = 0
    if not intheader[0, 6] == missing_data:
        # convert milliseconds to microseconds
        microsecond = intheader[0, 6] / 1e3
    datestr = '%i %00i %i %i %i %i' % (
        year, jday, hour, minute, second, microsecond)
    stats['starttime'] = datetime.strptime(datestr, '%Y %j %H %M %S %f')

    standard['sensor_serial_number'] = ''
    if intheader[1, 3] != missing_data:
        standard['sensor_serial_number'] = intheader[1, 3]

    # second line is information about number of channels, orientations
    # we care about orientations
    format_specific['vertical_orientation'] = np.nan
    if intheader[1, 4] != missing_data:
        format_specific['vertical_orientation'] = intheader[1, 4]

    standard['horizontal_orientation'] = np.nan
    if intheader[1, 5] != missing_data:
        standard['horizontal_orientation'] = intheader[1, 5]

    if intheader[1, 6] == missing_data:
        standard['instrument'] = ''
    else:
        standard['instrument'] = INSTRUMENTS[intheader[1, 6]]

    num_comments = intheader[1, 7]

    # third line contains number of data points
    stats['npts'] = intheader[2, 0]
    problem_flag = intheader[2, 1]
    if problem_flag == 1:
        fmt = 'Record found in file %s has a problem flag!'
        raise Exception(fmt % filename)
    stype = intheader[2, 2]
    if stype == missing_data:
        stype = np.nan
    standard['structure_type'] = STRUCTURES[stype]
    if standard['structure_type'] == 'building' and not any_structure:
        fmt = 'Record found in file %s is not a free-field sensor!'
        raise Exception(fmt % filename)

    format_specific['building_floor'] = np.nan
    if intheader[3, 0] != missing_data:
        format_specific['building_floor'] = intheader[3, 0]

    format_specific['bridge_number_spans'] = np.nan
    if intheader[3, 1] != missing_data:
        format_specific['bridge_number_spans'] = intheader[3, 1]

    format_specific['bridge_transducer_location'] = BRIDGE_LOCATIONS[0]
    if intheader[3, 2] != missing_data:
        format_specific['bridge_transducer_location'] = intheader[3, 2]

    format_specific['dam_transducer_location'] = DAM_LOCATIONS[0]
    if intheader[3, 3] != missing_data:
        format_specific['dam_transducer_location'] = intheader[3, 3]

    c1 = format_specific['bridge_transducer_location'].find('free field') == -1
    c2 = format_specific['dam_transducer_location'].find('free field') == -1
    if c1 and c2 and not any_structure:
        raise Exception(fmt % filename)

    format_specific['construction_type'] = CONSTRUCTION_TYPES[4]
    if intheader[3, 4] != missing_data:
        format_specific['construction_type'] = CONSTRUCTION_TYPES[intheader[3, 4]]

    # station is repeated here if all numeric
    if not len(stats['station']):
        stats['station'] = '%i' % intheader[3, 5]

    # read float header data
    skip = ASCII_HEADER_LINES + INTEGER_HEADER_LINES
    floatheader = np.genfromtxt(filename,
                                max_rows=FLOAT_HEADER_LINES,
                                skip_header=skip,
                                delimiter=FLOAT_HEADER_WIDTHS)

    # float headers are 10 lines of 5 floats each
    missing_data = floatheader[0, 0]
    stats['sampling_rate'] = floatheader[0, 1]
    coordinates['latitude'] = floatheader[2, 0]
    coordinates['longitude'] = floatheader[2, 1]
    coordinates['elevation'] = 0
    if floatheader[2, 2] != missing_data:
        coordinates['elevation'] = floatheader[2, 2]

    # figure out the channel code
    if format_specific['vertical_orientation'] in [0, 180]:
        stats['channel'] = get_channel_name(stats['sampling_rate'],
                                            is_acceleration=True,
                                            is_vertical=True,
                                            is_north=False)
    else:
        ho = standard['horizontal_orientation']
        quad1 = ho > 315 and ho <= 360
        quad2 = ho > 0 and ho <= 45
        quad3 = ho > 135 and ho <= 225
        if quad1 or quad2 or quad3:
            stats['channel'] = get_channel_name(stats['sampling_rate'],
                                                is_acceleration=True,
                                                is_vertical=False,
                                                is_north=True)
        else:
            stats['channel'] = get_channel_name(stats['sampling_rate'],
                                                is_acceleration=True,
                                                is_vertical=False,
                                                is_north=False)

    # we never get a two character location code
    stats['location'] = '--'

    sensor_frequency = floatheader[4, 1]
    standard['instrument_period'] = 1/sensor_frequency
    standard['instrument_damping'] = floatheader[4, 2]

    standard['corner_frequency'] = floatheader[3, 4]
    format_specific['filter_poles'] = floatheader[4, 0]
    standard['units'] = 'acc'

    # read in the comment lines
    with open(filename) as f:
        skip = ASCII_HEADER_LINES + INTEGER_HEADER_LINES + FLOAT_HEADER_LINES
        _ = [next(f) for x in range(skip)]
        standard['comments'] = [next(f).strip().lstrip('|')
                                for x in range(num_comments)]

    standard['comments'] = ' '.join(standard['comments'])
    stats['coordinates'] = coordinates
    stats['standard'] = standard
    stats['format_specific'] = format_specific

    return (stats, num_comments)
