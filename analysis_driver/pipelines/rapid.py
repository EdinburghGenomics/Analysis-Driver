import os
import csv
from collections import defaultdict
from egcg_core.config import cfg
from egcg_core import executor, util, rest_communication
from analysis_driver import segmentation, transfer_data
from analysis_driver.util import bash_commands


class RapidStage(segmentation.Stage):
    def rapid_output_dir(self, lane):
        return os.path.join(self.job_dir, 'rapid_analysis_%s' % lane)


class Dragen(RapidStage):
    def _run(self): 
        cmds = []
        for lane, sample in self.dataset.rapid_samples_by_lane.items():
            cmds.append(
                'if [ -d {tmp_dir} ]; then rm -r {tmp_dir}; fi; mkdir -p {tmp_dir}; {dragen} -r {ref} '
                '--bcl-input-dir {run_dir} --bcl-only-lane {lane} --output-directory {tmp_dir} '
                '--output-file-prefix {user_sample_id} --enable-variant-caller true --vc-sample-name {sample_id} '
                '--bcl-sample-sheet {sample_sheet} --enable-map-align-output true --enable-duplicate-marking true '
                '--dbsnp {dbsnp}; rsync -rLD {tmp_dir}/ {out_dir}; rm -r {tmp_dir}'.format(
                    dragen=cfg['dragen']['executable'],
                    ref=cfg['dragen']['reference'], run_dir=self.input_dir, lane=lane,
                    out_dir=self.rapid_output_dir(lane), user_sample_id=sample.udf['User Sample Name'],
                    sample_id=sample.name, sample_sheet=self.dataset.sample_sheet_file,
                    tmp_dir=os.path.join(cfg['dragen']['staging'], self.dataset.name + '_' + lane),
                    dbsnp=cfg['dragen']['dbsnp']
                )
            )

        return executor.execute(
            *cmds,
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

                data[line[0]][line[2]] = line[3:]

        return data

    def _run(self):
        data = []
        for lane in sorted(self.dataset.rapid_samples_by_lane):
            sample = self.dataset.rapid_samples_by_lane[lane]
            sample_id = sample.name
            user_sample_id = sample.udf['User Sample Name']
            output_dir = self.rapid_output_dir(lane)
            map_metrics = self.parse_metrics_file(
                util.find_file(output_dir, '%s.mapping_metrics.csv' % user_sample_id)
            )['MAPPING/ALIGNING SUMMARY']

            vc_metrics = self.parse_metrics_file(
                util.find_file(output_dir, '%s.vc_metrics.csv' % user_sample_id)
            )['VARIANT CALLER POSTFILTER']

            data.append(
                {
                    'sample_id': sample_id,
                    'rapid_metrics': {
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
            )

        rest_communication.post_or_patch('samples', data, id_field='sample_id')
        return 0


class DragenOutput(RapidStage):
    def _run(self):
        exit_status = 0
        for lane, sample in self.dataset.rapid_samples_by_lane.items():
            exit_status += self.output_data(lane, sample)

        return exit_status

    def output_data(self, lane, sample):
        user_sample_id = sample.udf['User Sample Name']
        staging_dir = os.path.join(self.job_dir, 'linked_output_files_%s' % lane)

        os.makedirs(staging_dir, exist_ok=True)
        transfer_data.create_output_links(
            self.job_dir,
            'rapid_analysis',
            staging_dir,
            user_sample_id=user_sample_id,
            lane=lane
        )

        # TODO: implement generic MD5Sum stage
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

        output_dir = os.path.join(cfg['sample']['output_dir'], sample.project.name, sample.name, 'rapid_analysis')
        transfer_exit_status = transfer_data.output_data_and_archive(
            staging_dir.rstrip('/') + '/',
            output_dir.rstrip('/')
        )
        return md5sum_exit_status + transfer_exit_status
