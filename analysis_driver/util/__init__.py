import os
import itertools
from egcg_core.app_logging import logging_default as log_cfg
from egcg_core.exceptions import EGCGError
from analysis_driver.reader.run_info import Reads
from . import bash_commands, helper_functions

app_logger = log_cfg.get_logger(__name__)


def find_all_fastq_pairs_for_lane(location, lane):
    fastqs = []
    for name, dirs, files in os.walk(location):
        fastqs.extend(os.path.join(name, f) for f in files if f.endswith('.fastq.gz') and '_L00%s_' % lane in f)
    if len(fastqs) % 2 != 0:
        raise EGCGError('Expected even number of fastq files in %s, found %s' % (location, len(fastqs)))
    fastqs.sort()
    app_logger.info('Found %s fastqs in %s for lane %s', len(fastqs), location, lane)
    return list(zip(*[iter(fastqs)] * 2))


def get_ranges(list_int):
    for a, b in itertools.groupby(enumerate(sorted(list_int)), lambda x: x[1] - x[0]):
        b = list(b)
        yield b[0][1], b[-1][1]


def get_trim_values_for_bad_cycles(bad_cycle_list, run_info):
    if not bad_cycle_list:
        return None, None

    read1_length = Reads.num_cycles(run_info.reads.upstream_read)
    read2_length = Reads.num_cycles(run_info.reads.upstream_read)
    index_length = sum(run_info.reads.index_lengths)

    # split bad cycles between read1 and read2
    bad_cycle_list1 = [c for c in bad_cycle_list if c <= read1_length]
    bad_cycle_list2 = [c - (read1_length + index_length) for c in bad_cycle_list if c > read1_length + index_length]

    ranges1 = list(get_ranges(bad_cycle_list1))
    ranges2 = list(get_ranges(bad_cycle_list2))

    if ranges1 and ranges1[-1][1] == read1_length and ranges1[-1][0] < read1_length:
        trim_r1 = ranges1[-1][0]-1
    else:
        trim_r1 = None
    if ranges2 and ranges2[-1][1] == read2_length and ranges2[-1][0] < read2_length:
        trim_r2 = ranges2[-1][0]-1
    else:
        trim_r2 = None

    return trim_r1, trim_r2
