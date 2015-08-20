__author__ = 'mwham'
import csv
import os.path
from logging import getLogger

app_logger = getLogger(__name__)


def write_bcbio_csv(run_dir, sample_id, fastqs):
    csv_file = os.path.join(run_dir, 'samples_' + sample_id + '.csv')
    app_logger.info('Writing BCBio sample csv ' + csv_file)

    with open(csv_file, 'w') as f:
        writer = csv.writer(f)

        writer.writerow(['samplename', 'description'])
        for fq in fastqs:
            writer.writerow([fq, sample_id])

    return csv_file
