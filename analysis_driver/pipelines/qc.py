from analysis_driver.pipelines import common, variant_calling as v
from analysis_driver import quality_control as qc

toolset_type = 'non_human_sample_processing'
name = 'qc'


def build_pipeline(dataset):

    def stage(cls, **params):
        return cls(dataset=dataset, **params)

    g = v.GATKStage(dataset=dataset)
    bam_file_production = common.build_bam_file_production(dataset)
    haplotype = stage(v.HaplotypeCaller, input_bam=g.sorted_bam, previous_stages=bam_file_production)
    bgzip_gvcf = stage(v.BGZipGvcf, previous_stages=[haplotype])
    tabix_gvcf = stage(v.TabixGvcf, previous_stages=[bgzip_gvcf])
    genotype = stage(v.GenotypeGVCFs, previous_stages=[tabix_gvcf])
    select_snp = stage(v.SelectVariants, previous_stages=[genotype])
    filter_snp = stage(v.VariantFiltration, previous_stages=[select_snp])
    vcfstats = stage(qc.VCFStats, vcf_file=g.filter_snp_vcf, previous_stages=[filter_snp])
    bgzip_vcf = stage(v.BGZipVcf, previous_stages=[vcfstats])
    tabix_vcf = stage(v.TabixVcf, previous_stages=[bgzip_vcf])
    data_output = stage(common.SampleDataOutput, previous_stages=[tabix_vcf], output_fileset='non_human_qc')
    cleanup = stage(common.Cleanup, previous_stages=[data_output])
    review = stage(common.SampleReview, previous_stages=[cleanup])
    return review
