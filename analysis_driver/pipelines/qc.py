from analysis_driver.pipelines import common

toolset_type = 'non_human_sample_processing'
name = 'qc'


def build_pipeline(dataset):
    bam_file_production = common.build_bam_file_production(dataset)
    data_output = common.SampleDataOutput(dataset=dataset, previous_stages=bam_file_production, output_fileset='non_human_qc')
    cleanup = common.Cleanup(dataset=dataset, previous_stages=[data_output])
    return cleanup
