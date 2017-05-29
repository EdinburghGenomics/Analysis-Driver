from analysis_driver.pipelines import common


def build_pipeline(dataset):
    bam_file_production = common.build_bam_file_production(dataset)
    data_output = common.DataOutput(dataset=dataset, previous_stages=bam_file_production, output_fileset='non_human_qc')
    cleanup = common.Cleanup(dataset=dataset, previous_stages=[data_output])
    return cleanup