from analysis_driver.pipelines import Pipeline, common, variant_calling as v
from analysis_driver import quality_control as qc


class QC(Pipeline):
    toolset_type = 'non_human_sample_processing'

    def build(self):
        g = self.stage(v.GATKStage)
        bam_file_production = common.build_bam_file_production(self.dataset)
        haplotype = self.stage(v.HaplotypeCaller, input_bam=g.sorted_bam, previous_stages=bam_file_production)
        genotype = self.stage(v.GenotypeGVCFs, previous_stages=[haplotype])
        select_snp = self.stage(v.SelectVariants, previous_stages=[genotype])
        filter_snp = self.stage(v.VariantFiltration, previous_stages=[select_snp])
        vcfstats = self.stage(qc.VCFStats, vcf_file=g.filter_snp_vcf + '.gz', previous_stages=[filter_snp])
        data_output = self.stage(common.SampleDataOutput, previous_stages=[vcfstats], output_fileset='non_human_qc')
        cleanup = self.stage(common.Cleanup, previous_stages=[data_output])
        review = self.stage(common.SampleReview, previous_stages=[cleanup])
        return review
