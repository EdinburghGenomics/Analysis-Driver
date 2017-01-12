import glob
import os
from egcg_core import rest_communication, clarity
from analysis_driver.config import default as cfg
from analysis_driver.quality_control import Relatedness

def project_pipeline(dataset):
    species = clarity.get_species_from_sample(dataset.name)
    project_id = dataset.name
    working_dir = os.path.join(cfg['jobs_dir'], project_id)
    reference = cfg['references'][species]
    project = rest_communication.get_document('projects', where={'project_id': project_id})
    project_source = os.path.join(cfg.query('sample', 'delivery_source'), project_id)
    gvcf_paths = project_source + '/*/*g.vcf.gz'
    gvcf_files = glob.glob(gvcf_paths)

    dataset.start_stage('relatedness')
    r = Relatedness(dataset, working_dir, gvcf_files, reference, project_id)
    r.start()
    vcftools_relatedness_expected_outfile, exit_status = r.join()
    dataset.end_stage('relatedness')
    return exit_status
