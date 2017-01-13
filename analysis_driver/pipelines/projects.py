import glob
import os
from egcg_core import rest_communication, clarity
from analysis_driver.config import default as cfg
from analysis_driver.quality_control import Relatedness
from analysis_driver.exceptions import PipelineError

def project_pipeline(dataset):

    project_id = dataset.name
    samples_for_project = rest_communication.get_documents('samples', where={'project_id': project_id})
    species_in_project = []
    for sample in samples_for_project:
        sample_id = sample.get('sample_id')
        species_in_lims = clarity.get_species_from_sample(sample_id)
        species_in_project.append(species_in_lims)
    if len(set(species_in_project)) != 1:
        raise PipelineError('Wrong number of species in this project: expected 1')
    species = list(set(species_in_project))[0]

    working_dir = os.path.join(cfg['jobs_dir'], project_id)
    reference = cfg['references'][species]['fasta']
    project_source = os.path.join(cfg.query('sample', 'delivery_source'), project_id)
    gvcf_paths = project_source + '/*/*g.vcf.gz'
    gvcf_files = glob.glob(gvcf_paths)
    if not len(gvcf_files) > 1:
        raise PipelineError('Incorrect number of gVCF files: require at least two')

    dataset.start_stage('relatedness')
    r = Relatedness(dataset, working_dir, gvcf_files, reference, project_id)
    r.start()
    vcftools_relatedness_expected_outfile, exit_status = r.join()
    dataset.end_stage('relatedness')
    return exit_status
