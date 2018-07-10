# stdlib imports
import warnings

# local imports
from pgm.exception import PGMException
from pgm.gather import get_pgm_classes, group_imcs


def calculate_arias(stream, imcs):
    """
    Calculate the peak ground acceleration.

    Args:
        stream (obspy.core.stream.Stream): streams of strong ground motion.
            Traces in stream must be in units of m/s/s.
        imcs (list): list of imcs.

    Returns:
        dictionary: Dictionary of arias for different components.
    """
    arias_dict = {}
    # check units and add channel arias
    for idx, trace in enumerate(stream):
        if trace.stats['units'] == '%%g':
            stream[idx].data = trace.data * 0.01
        elif trace.stats['units'] != 'm/s/s'
            raise PGMException('Invalid units for arias: %r. '
                               'Units must be m/s/s' % trace.stats['units'])
    # sort imcs
    grouped_imcs = group_imcs(imcs)
    # gather imc classes
    pgm_classes = get_pgm_classes('imc')
    # store arias for imcs
    for imc in grouped_imcs:
        if 'calculate_' + imc in pgm_classes:
            arias_func = pgm_classes['calculate_' + imc]
            arias = arias_func(stream, percentiles=grouped_imcs[imc])
            if imc.find('rot') >= 0:
                for percentile in arias:
                    arias_dict[imc.upper() + str(percentile)] = arias[percentile]
            elif imc.find('channels') >= 0:
                for channel in arias:
                    arias_dict[channel] = arias[channel]
            else:
                arias_dict[imc.upper()] = arias
        else:
            warnings.warn('Not a valid IMC: %r. Skipping...' % imc)
    return arias_dict
