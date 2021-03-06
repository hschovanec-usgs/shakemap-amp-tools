import numpy as np
from amptools.io.seedname import get_channel_name
from obspy import UTCDateTime
from obspy.clients.fdsn import Client
from obspy.core.util.attribdict import AttribDict

TIMEFMT = '%Y-%m-%dT%H:%M:%S'


def request_raw_waveforms(fdsn_client, org_time, lat, lon,
                          before_time=120, after_time=600,
                          dist_min=0, dist_max=1.5, networks=['*'],
                          stations=['*'], channels=['*']):
    """
    Requests raw waveform data from an FDSN client.
    The requested data can be constrained by time, distance from event,
    network, station, and channels.

    Args:
        fdsn_client (str): A string for the FDSN client. Valid choices are:
            BGR, EMSC, ETH, GEONET, GFZ, ICGC, INGV, IPGP, IRIS, ISC,
            KOERI, LMU, NCEDC, NIEP, NOA, ODC, ORFEUS, RESIF, SCEDC, TEXTNET,
            USGS, and USP
        org_time (str): Origin time of event. Must be able to convert to
            UTCDateTime.
        lat (float): Event latitude.
        lon (float): Event longitude.
        before_time (int): Number of seconds before origin time to request.
            Default is 120.
        after_time (int): Number of seconds after origin time to request.
            Default is 600.
        dist_min (float): Minimum distance from event epicenter.
            Default is 0.
        dist_max (float): Maximum distance from event epicenter.
            Default is 1.5
        network (list): List of strings for desired networks.
            Default is ['*'].
        channels (list): List of strings for desired channels.
            Default is ['*'].

    Returns:
        stream (obspy.core.trace.Trace): Stream of requested, raw data.
        inventory (obspy.core.inventory): Inventory object for the event.
    """

    client = Client(fdsn_client)

    # Time information
    origin_time = UTCDateTime(org_time)
    t1 = origin_time - before_time
    t2 = origin_time + after_time

    # Convert lists to comma-separated strings
    networks = ','.join(networks)
    stations = ','.join(stations)
    channels = ','.join(channels)

    # Get an inventory of all stations for the event
    inventory = client.get_stations(startbefore=t1, endafter=t2,
                                    latitude=lat, longitude=lon,
                                    minradius=dist_min, maxradius=dist_max,
                                    network=networks,
                                    channel=channels,
                                    station=stations,
                                    level='response', includerestricted=False)

    # Get the list of channels from the inventory
    channels = inventory.get_contents()['channels']
    print('Found {0} channels.'.format(len(channels)))

    # Set up the bulk data for the bulk data request
    bulk = [chan.split('.') for chan in channels]
    for b in bulk:
        b.append(t1)
        b.append(t2)

    # Perform the bulk data request
    print('Requesting waveforms for {0} channels.'.format(len(channels)))
    st = client.get_waveforms_bulk(bulk, attach_response=True)
    return st, inventory


def add_channel_metadata(tr, inv, client):
    """
    Adds the channel metadata for each channel in the stream.

    Args:
        tr (obspy.core.trace.Trace): Trace of requested data.
        inv (obspy.core.inventory): Inventory object corresponding to
            to the stream.
        client (str): FDSN client indicator.

    Returns:
        trace (obspy.core.trace.Trace): Trace with metadata added.
    """

    time = tr.stats.starttime
    id_string = tr.stats.network + '.' + tr.stats.station + '.'
    id_string += tr.stats.location + '.' + tr.stats.channel

    metadata = inv.get_channel_metadata(id_string, time)

    coordinates = {'latitude': metadata['latitude'],
                   'longitude': metadata['longitude'],
                   'elevation': metadata['elevation']}

    standard = {'horizontal_orientation': metadata['azimuth'],
                'instrument_period': np.nan,
                'instrument_damping': np.nan,
                'process_level': 'V0',
                'station_name': tr.stats.station,
                'sensor_serial_number': '',
                'instrument': '',
                'comments': '',
                'structure_type': '',
                'corner_frequency': np.nan,
                'units': 'raw',
                'source': client,
                'source_format': 'fdsn'}

    tr.stats['coordinates'] = coordinates
    tr.stats['standard'] = standard

    if metadata['dip'] in [90, -90, 180, -180]:
        tr.stats['channel'] = get_channel_name(tr.stats['sampling_rate'],
                                               is_acceleration=True,
                                               is_vertical=True,
                                               is_north=False)
    else:
        ho = metadata['azimuth']
        quad1 = ho > 315 and ho <= 360
        quad2 = ho >= 0 and ho <= 45
        quad3 = ho > 135 and ho <= 225
        if quad1 or quad2 or quad3:
            tr.stats['channel'] = get_channel_name(tr.stats['sampling_rate'],
                                                   is_acceleration=True,
                                                   is_vertical=False,
                                                   is_north=True)
        else:
            tr.stats['channel'] = get_channel_name(tr.stats['sampling_rate'],
                                                   is_acceleration=True,
                                                   is_vertical=False,
                                                   is_north=False)
    return tr


def clean_stats(my_stats):
    """
    Function for making dictionary json serializable.

    Args:
        stats (dict): Dictionary of stats.

    Returns:
        dictionary: Dictionary of cleaned stats.
    """
    stats = dict()
    for key, value in my_stats.items():
        stats[key] = value

    for key, value in stats.items():
        if isinstance(value, (dict, AttribDict)):
            stats[key] = dict(clean_stats(value))
        elif isinstance(value, UTCDateTime):
            stats[key] = value.strftime(TIMEFMT)
        elif isinstance(value, float) and np.isnan(value) or value == '':
            stats[key] = 'null'
    return stats


def remove_response(stream, output='ACC'):
    """
    Remove the instrument response from a stream and converts the waveforms
    to either acceleration, velocity, or displacement. This requires that
    the response data has been attached to the trace through the
    bulk waveform request. This is step 3 as found in the Rennolet
    et. al paper (https://doi.org/10.1193/101916EQS175DP)

    Args:
        stream (obspy.core.stream.Stream): Stream of raw data.
        output (str): Output units. Must be 'DISP, 'VEL', or 'ACC'.
            Default is 'ACC'.

    Returns:
        stream (obspy.core.stream.Stream): Stream of data with instrument
            responses removed.
    """
    stream.remove_response(output=output)
    return stream
