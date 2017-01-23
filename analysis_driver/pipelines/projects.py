import os
from egcg_core.util import find_file
from egcg_core import rest_communication, clarity
from analysis_driver.config import default as cfg
from analysis_driver.quality_control import Relatedness
from analysis_driver.exceptions import PipelineError
from analysis_driver.transfer_data import output_project_data

def project_pipeline(dataset):

    project_id = dataset.name
    samples_for_project = rest_communication.get_documents('aggregate/samples', match={'project_id': project_id, 'status': 'finished'})
    species_in_project = []
    for sample in samples_for_project:
        species = sample.get('species_name')
        if not species:
            species = clarity.get_species_from_sample(sample.get('sample_id'))
        species_in_project.append(species)
    if len(set(species_in_project)) != 1:
        raise PipelineError('Wrong number of species in this project: expected 1')
    species = list(set(species_in_project))[0]

    working_dir = os.path.join(cfg['jobs_dir'], project_id)
    reference = cfg['references'][species]['fasta']
    project_source = os.path.join(cfg.query('sample', 'delivery_source'), project_id)
    gvcf_files = []
    for sample in samples_for_project:
        gvcf_file = find_file(project_source, sample['sample_id'], sample['user_sample_id'] + '.g.vcf.gz')
        if gvcf_file:
            gvcf_files.append(gvcf_file)
    if not len(gvcf_files) > 1:
        raise PipelineError('Incorrect number of gVCF files: require at least two')

    dataset.start_stage('relatedness')
    r = Relatedness(dataset, working_dir, gvcf_files, reference, project_id)
    r.start()
    vcftools_relatedness_expected_outfile, exit_status = r.join()
    dataset.end_stage('relatedness')

    exit_status += output_project_data(working_dir, project_id)

    return exit_status
