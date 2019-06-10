import os
import csv
from datetime import date
from collections import defaultdict
from egcg_core.config import cfg
from egcg_core import executor, util, rest_communication
from analysis_driver import segmentation, transfer_data
from analysis_driver.util import bash_commands
from analysis_driver.reader.demultiplexing_parsers import parse_interop_summary


class RapidStage(segmentation.Stage):
    def rapid_output_dir(self, lane):
        return os.path.join(self.job_dir, 'rapid_analysis_%s' % lane)


class Dragen(RapidStage):
    def _run(self): 
        cmds = []
        for lane, sample in self.dataset.rapid_samples_by_lane.items():
            # clean up any previous run on the Dragen server's scratch space
            cmd = 'if [ -d {tmp_dir} ]; then rm -r {tmp_dir}; fi && mkdir -p {tmp_dir} && ' \
                  'echo "Starting at $(date +%Y-%m-%d\ %H:%M:%S)" && '

            # run Dragen on the lane with bam/vcf outputs, dbSNP and duplicate marking enabled
            cmd += '{dragen} -r {ref} --bcl-input-dir {run_dir} --bcl-only-lane {lane} --output-directory {tmp_dir} ' \
                   '--output-file-prefix {user_sample_id} --enable-variant-caller true ' \
                   '--vc-sample-name {user_sample_id} --bcl-sample-sheet {sample_sheet} ' \
                   '--enable-map-align-output true --enable-duplicate-marking true --dbsnp {dbsnp}; ' \
                   'exit_status=$?; '

            # somehow, writing the output data to scratch and then rsyncing to the main file system is faster than
            # outputting directly
            cmd += 'echo "Finished setup/Dragen with exit status $exit_status at $(date +%Y-%m-%d\ %H:%M:%S)"; ' \
                   'echo "Outputting data and cleaning up"; ' \
                   'rsync -rLD {tmp_dir}/ {out_dir}; rm -r {tmp_dir}; ' \
                   'echo "Done at $(date +%Y-%m-%d\ %H:%M:%S)"; ' \
                   'exit $exit_status'

            cmds.append(
                cmd.format(
                    dragen=cfg['dragen']['executable'], ref=cfg['dragen']['reference'], run_dir=self.input_dir,
                    lane=lane, out_dir=self.rapid_output_dir(lane), user_sample_id=sample['User Sample Name'],
                    sample_sheet=self.dataset.sample_sheet_file, dbsnp=cfg['dragen']['dbsnp'],
                    tmp_dir=os.path.join(cfg['dragen']['staging'], self.dataset.name + '_' + lane)
                )
            )

        return executor.execute(
            *cmds,
            # mimic the ulimits set on the Dragen server
            prelim_cmds=['ulimit -n 65535', 'ulimit -c 0', 'ulimit -s 8192', 'ulimit -u 16384', 'ulimit -i 1029522'],
            job_name='dragen',
            job_queue='dragen',
            working_dir=self.job_dir,
            exclusive=True,
            log_commands=False
        ).join()


class DragenMetrics(RapidStage):
    @staticmethod
    def parse_metrics_file(metrics_file):
        data = defaultdict(dict)
        with open(metrics_file) as f:
            reader = csv.reader(f)
            for line in reader:
                if not line:
                    continue

                heading, sample, key, *values = line
                data[heading][key] = values

        return data

    @staticmethod
    def _interop_metrics(interop_metrics):
        yield_r1 = interop_metrics.get('yield_r1')
        yield_r2 = interop_metrics.get('yield_r2')
        if not yield_r1 or not yield_r2:
            return {}

        pc_q30_r1 = interop_metrics['pc_q30_r1']
        pc_q30_r2 = interop_metrics['pc_q30_r2']

        projected_yield_q30_r1 = (yield_r1 * pc_q30_r1) / 100
        projected_yield_q30_r2 = (yield_r2 * pc_q30_r2) / 100

        projected_yield_q30 = projected_yield_q30_r1 + projected_yield_q30_r2
        total_yield = yield_r1 + yield_r2
        pc_q30 = (projected_yield_q30 / total_yield) * 100
        return {'yield': total_yield, 'pc_q30': pc_q30}

    def _run(self):
        interop_metrics = parse_interop_summary(util.find_file(self.job_dir, 'fastq', 'interop_summary.txt'))
        data = []
        for lane in sorted(self.dataset.rapid_samples_by_lane):
            sample = self.dataset.rapid_samples_by_lane[lane]
            user_sample_id = sample['User Sample Name']
            output_dir = self.rapid_output_dir(lane)
            map_metrics = self.parse_metrics_file(
                util.find_file(output_dir, '%s.mapping_metrics.csv' % user_sample_id)
            )['MAPPING/ALIGNING SUMMARY']

            vc_metrics = self.parse_metrics_file(
                util.find_file(output_dir, '%s.vc_metrics.csv' % user_sample_id)
            )['VARIANT CALLER POSTFILTER']

            payload = {
                'sample_id': sample['sample_id'],
                'rapid_analysis': {
                    'data_source': [self.dataset.name + '_' + lane],
                    'var_calling': {
                         'ti_tv_ratio': float(vc_metrics['Ti/Tv ratio'][0]),
                         'het_hom_ratio': float(vc_metrics['Het/Hom ratio'][0]),
                    },
                    'mapping': {
                        'total_reads': int(map_metrics['Total input reads'][0]),
                        'duplicate_reads': int(map_metrics['Number of duplicate reads (marked)'][0]),
                        'pc_duplicates': float(map_metrics['Number of duplicate reads (marked)'][1]),
                        'mapped_reads': int(map_metrics['Mapped reads'][0]),
                        'pc_mapped': float(map_metrics['Mapped reads'][1]),
                        'unique_mapped_reads': int(map_metrics['Number of unique & mapped reads (excl. dups)'][0]),
                        'pc_unique_mapped': float(map_metrics['Number of unique & mapped reads (excl. dups)'][1]),
                        'mean_coverage': float(map_metrics['Average autosomal coverage'][0]),
                        'median_coverage': float(map_metrics['Median autosomal coverage'][0]),
                        'genome_size': int(map_metrics['Bases in reference genome'][0]),
                        'pc_genome_with_coverage': {
                            '100+': float(map_metrics['PCT of genome with coverage [100x:inf)'][0]),
                            '50-100': float(map_metrics['PCT of genome with coverage [50x:100x)'][0]),
                            '40-50': float(map_metrics['PCT of genome with coverage [40x:50x)'][0]),
                            '30-40': float(map_metrics['PCT of genome with coverage [30x:40x)'][0]),
                            '20-30': float(map_metrics['PCT of genome with coverage [20x:30x)'][0]),
                            '10-20': float(map_metrics['PCT of genome with coverage [10x:20x)'][0]),
                            '5-10': float(map_metrics['PCT of genome with coverage [ 5x:10x)'][0]),
                            '2-5': float(map_metrics['PCT of genome with coverage [ 2x: 5x)'][0]),
                            '1-2': float(map_metrics['PCT of genome with coverage [ 1x: 2x)'][0]),
                            '0-1': float(map_metrics['PCT of genome with coverage [ 0x: 1x)'][0]),
                        }
                    }
                }
            }
            payload['rapid_analysis'].update(self._interop_metrics(interop_metrics[lane]))
            data.append(payload)

        rest_communication.post_or_patch('samples', data, id_field='sample_id')
        return 0


class DragenOutput(RapidStage):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._datestamp = None

    @property
    def datestamp(self):
        if self._datestamp is None:
            self._datestamp = date.today().isoformat()
        return self._datestamp

    def delivery_folder(self, sample):
        return os.path.join(
            cfg['delivery']['dest'],
            sample.project.name,
            self.datestamp + '_rapid',
            sample.udf['User Sample Name']
        )

    def _run(self):
        exit_status = 0

        for lane in sorted(self.dataset.rapid_samples_by_lane):
            sample = self.dataset.rapid_samples_by_lane[lane]
            sample_id = sample['sample_id']

            qc_passed = self.review(sample_id)
            if not qc_passed:
                continue

            output_dir = os.path.join(cfg['sample']['output_dir'], sample['project_id'], sample_id, 'rapid_analysis')
            output_status = self.output_data(lane, sample, output_dir)
            if output_status != 0:
                self.warning('Data output failed for sample %s', sample_id)
                continue

            exit_status += self.deliver_data(sample_id, sample['User Sample Name'], sample['project_id'], output_dir)

        return exit_status

    def review(self, sample_id):
        rest_communication.post_entry(
            'actions',
            {'action_type': 'automatic_rapid_review', 'sample_id': sample_id},
            use_data=True
        )

        rapid_metrics = rest_communication.get_document('samples', where={'sample_id': sample_id})['rapid_analysis']
        if rapid_metrics['reviewed'] == 'pass':
            return True
        else:
            self.warning(
                'Sample %s has status %s: %s',
                sample_id,
                rapid_metrics['reviewed'],
                rapid_metrics.get('review_comments', '(no comments)')
            )
            return False

    def output_data(self, lane, sample, output_dir):
        user_sample_id = sample['User Sample Name']
        staging_dir = os.path.join(self.job_dir, 'linked_output_files_%s' % lane)

        os.makedirs(staging_dir, exist_ok=True)
        transfer_data.create_output_links(
            self.job_dir,
            'rapid_analysis',
            staging_dir,
            user_sample_id=user_sample_id,
            lane=lane
        )

        # TODO: #274 - implement generic MD5Sum stage
        self.dataset.start_stage('md5sum')
        md5sum_exit_status = executor.execute(
            *[
                bash_commands.md5sum(os.path.join(staging_dir, f))
                for f in os.listdir(staging_dir)
                if not f.endswith('md5')
            ],
            job_name='md5sum',
            working_dir=self.job_dir,
            cpus=1,
            mem=2,
            log_commands=False
        ).join()
        self.dataset.end_stage('md5sum')

        transfer_exit_status = transfer_data.output_data_and_archive(
            staging_dir.rstrip('/') + '/',
            output_dir.rstrip('/')
        )

        return md5sum_exit_status + transfer_exit_status

    def deliver_data(self, sample_id, user_sample_id, project_id, output_dir):
        today = self.datestamp + '_rapid'
        delivery_folder = os.path.join(
            cfg['delivery']['dest'],
            project_id,
            today,
            user_sample_id
        )
        os.makedirs(delivery_folder, exist_ok=False)

        exit_status = 0
        for f in util.find_files(output_dir, '*'):  # for now, deliver everything that has been output
            new_link = os.path.join(delivery_folder, os.path.basename(f))
            exit_status += executor.local_execute('ln %s %s' % (f, new_link)).join()

        self.dataset.ntf.notify_rapid_completion(sample_id)
        return exit_status


def build_pipeline(dataset, setup):
    """
    Not currently intended for use on its own - used by demultiplexing.build_pipeline.
    :param dataset.RunDataset dataset:
    :param analysis_driver.segmentation.Stage setup: The Setup stage from demultiplexing
    """
    def stage(cls, **kwargs):
        return cls(dataset=dataset, **kwargs)

    dragen = stage(Dragen, previous_stages=[setup])
    dragen_metrics = stage(DragenMetrics, previous_stages=[dragen])
    dragen_output = stage(DragenOutput, previous_stages=[dragen_metrics])
    return dragen_output
