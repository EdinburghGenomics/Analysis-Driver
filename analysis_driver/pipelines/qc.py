from analysis_driver.pipelines import common, variant_calling as v
from analysis_driver import quality_control as qc

toolset_type = 'non_human_sample_processing'
name = 'qc'


def build_pipeline(dataset):

    def stage(cls, **params):
        return cls(dataset=dataset, **params)

    g = v.GATKStage()
    bam_file_production = common.build_bam_file_production(dataset)
    haplotype = stage(v.HaplotypeCaller, input_bam=g.sorted_bam, previous_stages=bam_file_production)
    genotype = stage(v.GenotypeGVCFs, previous_stages=[haplotype])
    select_snp = stage(v.SelectVariants, command=g.select_snp_command, previous_stages=[genotype])
    filter_snp = stage(v.VariantFiltration, command=g.filter_snp_command, previous_stages=[select_snp])
    vcfstats = stage(qc.VCFStats, vcf_file=g.filter_snp_vcf, previous_stages=[filter_snp])
    bgzip = stage(v.BGZip, previous_stages=[vcfstats])
    tabix = stage(v.Tabix, previous_stages=[bgzip])
    data_output = common.SampleDataOutput(dataset=dataset, previous_stages=[tabix], output_fileset='non_human_qc')
    cleanup = common.Cleanup(dataset=dataset, previous_stages=[data_output])
    review = common.SampleReview(dataset=dataset, previous_stages=[cleanup])
    return review
