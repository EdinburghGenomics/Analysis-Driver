__author__ = 'tcezard'
import sys
import os
import argparse
import logging
import datetime
from collections import defaultdict

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from analysis_driver.config import default as cfg, logging_default as log_cfg
from analysis_driver.app_logging import AppLogger
from analysis_driver import executor, clarity, rest_communication
from analysis_driver.exceptions import AnalysisDriverError



hs_list_files = [
    '{ext_sample_id}.g.vcf.gz',
    '{ext_sample_id}.g.vcf.gz.tbi',
    '{ext_sample_id}.vcf.gz',
    '{ext_sample_id}.vcf.gz.tbi',
    '{ext_sample_id}.bam',
    '{ext_sample_id}.bam.bai',
    '{ext_sample_id}_R1.fastq.gz',
    '{ext_sample_id}_R2.fastq.gz',
    '{ext_sample_id}_R1_fastqc.html',
    '{ext_sample_id}_R2_fastqc.html',
]

other_list_files = [
    '{ext_sample_id}_R1.fastq.gz',
    '{ext_sample_id}_R2.fastq.gz',
    '{ext_sample_id}_R1_fastqc.html',
    '{ext_sample_id}_R2_fastqc.html',
]


class DataDelivery(AppLogger):

    def __init__(self, dry_run):
        self.all_commands = []
        self.dry_run = dry_run


    def get_deliverable_projects_samples(self, project_id=None, sample_id=None):
        project_to_samples = defaultdict(list)
        #Get sample that have been review from REST API but not marked as delivered
        where_clause={"useable":"yes"}
        if project_id:
            where_clause["project_id"] = project_id
        if sample_id:
            where_clause["sample_id"] = sample_id

        samples = rest_communication.get_documents(
                'samples',
                depaginate=True,
                embedded={"analysis_driver_procs":1, "run_elements":1},
                where=where_clause
        )
        self.all_samples_values =samples
        for sample in samples:
            processes = sample.get('analysis_driver_procs',[{}])
            processes.sort(key=lambda x: datetime.datetime.strptime(x.get('_created', '01_01_1970_00:00:00'), '%d_%m_%Y_%H:%M:%S'))
            if not processes[-1].get('status','new') == 'finished':
                raise AnalysisDriverError("Reviewed sample %s not marked as finished"%sample.get('sample_id'))
            if sample.get('delivered','no') == 'no':
                project_to_samples[sample.get('project_id')].append(sample.get('sample_id'))
        return project_to_samples

    def summarise_metrics_per_sample(self, project_id):
        headers=['Project', 'Sample Id', 'User sample id', 'Read pair sequenced', 'Yield', 'Yield Q30', 'Nb reads in bam', 'mapping rate', 'properly mapped reads rate', 'duplicate rate', 'Mean coverage', 'Callable bases rate']
        lines = []
        for sample in self.all_samples_values:
            if sample.get('project_id') == project_id:
                res = [sample.get('project_id'),sample.get('sample_id'), sample.get('user_sample_id')]
                clean_reads = sum([int(e.get('clean_reads', '0')) for e in sample.get('run_elements') if e.get('useable')=='yes'])
                clean_bases_r1 = sum([int(e.get('clean_bases_r1', '0')) for e in sample.get('run_elements') if e.get('useable')=='yes'])
                clean_bases_r2 = sum([int(e.get('clean_bases_r2', '0')) for e in sample.get('run_elements') if e.get('useable')=='yes'])
                clean_q30_bases_r1 = sum(int(e.get('clean_q30_bases_r1', '0')) for e in sample.get('run_elements') if e.get('useable')=='yes')
                clean_q30_bases_r2 = sum(int(e.get('clean_q30_bases_r2', '0')) for e in sample.get('run_elements') if e.get('useable')=='yes')
                res.append(str(clean_reads))
                res.append(str((clean_bases_r1+clean_bases_r2)/1000000000))
                res.append(str((clean_q30_bases_r1+clean_q30_bases_r2)/1000000000))
                theoritical_cov = (clean_bases_r1 + clean_bases_r2)/3200000000.0
                tr = sample.get('bam_file_reads', 0)
                mr = sample.get('mapped_reads', 0)
                dr = sample.get('duplicate_reads', 0)
                pmr = sample.get('properly_mapped_reads', 0)
                if not tr:
                    raise AnalysisDriverError('Sample %s has no total number of reads'%sample.get('sample_id'))
                res.append(str(tr))
                res.append(str(float(mr)/float(tr)*100))
                res.append(str(float(pmr)/float(tr)*100))
                res.append(str(float(dr)/float(tr)*100))
                theoritical_cov = theoritical_cov*(float(mr)/float(tr))*(1-(float(dr)/float(tr)))
                res.append(str(theoritical_cov))
                #res.append(str(sample.get('median_coverage', '')))
                res.append(str(sample.get('pc_callable', 0)*100))
                lines.append('\t'.join(res))
        return headers, lines

    def append_create_batch_delivery_folder(self, deliver_dir, project):
        batch = datetime.date.today().isoformat()
        batch_delivery_folder = os.path.join(deliver_dir, project, batch)
        command = 'mkdir -p %s'%batch_delivery_folder
        self.all_commands.append(command)
        return batch_delivery_folder

    def append_create_sample_folder(self, batch_delivery_folder, sample):
        sample_folder = os.path.join(batch_delivery_folder, sample)
        command = 'mkdir -p %s'%sample_folder
        self.all_commands.append(command)
        return sample_folder

    def append_move_file_to_sample_folder(self, file_to_move, sample_folder):
        name = os.path.basename(file_to_move)
        command = 'mv %s %s'%(file_to_move, os.path.join(sample_folder, name))
        self.all_commands.append(command)

    def get_list_of_file_to_move(self, sample_name):
        species = clarity.get_species_from_sample(sample_name)
        external_sample_name = clarity.get_user_sample_name(sample_name)
        if species is None:
            raise AnalysisDriverError('No species information found in the LIMS for ' + sample_name)
        elif species == 'Homo sapiens':
            list_of_file = hs_list_files
        else:
            list_of_file = other_list_files
        final_list = []
        for f in list_of_file:
            final_list.append(f.format(ext_sample_id=external_sample_name))
            final_list.append(f.format(ext_sample_id=external_sample_name) + '.md5')
        return final_list

    def mark_samples_as_released(self, samples):
        for sample_name in samples:
            rest_communication.patch_entry('samples',payload={'delivered': 'yes'}, id_field='sample_id', element_id=sample_name)


    def mark_only(self, project_id=None, sample_id=None):
        project_to_samples = self.get_deliverable_projects_samples(project_id, sample_id)
        all_samples=[]
        for project in project_to_samples:
            samples = project_to_samples.get(project)
            for sample in samples:
                all_samples.append(sample)
        if self.dry_run:
            self.info('Mark %s samples as delivered'%len(all_samples))
        else:
            self.mark_samples_as_released(all_samples)

    def deliver_data(self, project_id=None, sample_id=None):
        delivery_source = cfg.query('delivery_source')
        delivery_dest = cfg.query('delivery_dest')
        project_to_samples = self.get_deliverable_projects_samples(project_id, sample_id)
        all_samples = []
        project_to_delivery_folder = {}

        for project in project_to_samples:
            samples = project_to_samples.get(project)
            self.info('deliver %s samples from %s'%(len(samples), project))
            batch_delivery_folder = self.append_create_batch_delivery_folder(delivery_dest, project)
            project_to_delivery_folder[project] = batch_delivery_folder
            for sample in samples:
                origin_sample_dir = os.path.join(delivery_source, project, sample)
                if not os.path.isdir(origin_sample_dir):
                    raise AnalysisDriverError('Directory for sample %s in project %s does not exist'%(sample, project))
                sample_folder = self.append_create_sample_folder(batch_delivery_folder, sample)
                list_of_file_to_move = self.get_list_of_file_to_move(sample_name=sample)

                for file_to_move in list_of_file_to_move:
                    origin_file = os.path.join(origin_sample_dir, file_to_move)
                    if not os.path.isfile(origin_file):
                        raise AnalysisDriverError('File %s for sample %s in project %s does not exist'%(file_to_move, sample, project))
                    self.append_move_file_to_sample_folder(origin_file, sample_folder)
                all_samples.append(sample)

        if self.dry_run:
            print('====== Commands ======')
            print('\n'.join(self.all_commands))
            print('====== Metrics ======')
            for project in project_to_samples:
                header, lines = self.summarise_metrics_per_sample(project)
                print('\t'.join(header))
                print('\n'.join(lines))
        else:
            self.mark_samples_as_released(all_samples)
            for command in self.all_commands:
                exit_status = executor.execute([command], env='local').join()
                if exit_status != 0:
                    raise AnalysisDriverError('command %s exited with status %s'%(command, exit_status))

            for project in project_to_delivery_folder:
                header, lines = self.summarise_metrics_per_sample(project)
                summary_metrics_file = os.path.join(project_to_delivery_folder.get(project),'summary_metrics.csv')
                if os.path.isfile(summary_metrics_file):
                    with open(summary_metrics_file, 'a') as open_file:
                        open_file.write('\n'.join(lines) + '\n')
                else:
                    with open(summary_metrics_file, 'w') as open_file:
                        open_file.write('\t'.join(header) + '\n')
                        open_file.write('\n'.join(lines) + '\n')

            # TODO: Generate project report


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--dry_run', action='store_true')
    p.add_argument('--debug', action='store_true')
    p.add_argument('--work_dir', type=str, required=True)
    p.add_argument('--mark_only', action='store_true')
    p.add_argument('--project_id', type=str)
    p.add_argument('--sample_id', type=str)
    args = p.parse_args()

    if args.debug:
        log_cfg.default_level = logging.DEBUG
        log_cfg.add_handler('stdout', logging.StreamHandler(stream=sys.stdout), logging.DEBUG)

    cfg.merge(cfg['sample'])
    dd = DataDelivery(args.dry_run)
    if args.mark_only:
        dd.mark_only(project_id=args.project_id, sample_id=args.sample_id)
    else:
        dd.deliver_data(project_id=args.project_id, sample_id=args.sample_id)

if __name__ == '__main__':
    main()
