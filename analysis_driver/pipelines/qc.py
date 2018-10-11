import os

from analysis_driver.pipelines import common, variant_calling as v
from analysis_driver import quality_control as qc
from analysis_driver.util import create_regions_of_genome

toolset_type = 'non_human_sample_processing'
name = 'qc'


def build_pipeline(dataset):

    def stage(cls, **params):
        return cls(dataset=dataset, **params)

    g = v.GATKStage(dataset=dataset)

    # limit variant call to a subset of contigs to speed up QC
    region_bed = os.path.join(g.job_dir, 'regions.bed')
    # creates the bed file from the genome
    create_regions_of_genome(dataset.reference_genome, region_bed)

    bam_file_production = common.build_bam_file_production(dataset)
    haplotype = stage(v.HaplotypeCaller, input_bam=g.sorted_bam, region_bed=region_bed, previous_stages=bam_file_production)
    genotype = stage(v.GenotypeGVCFs, region_bed=region_bed, previous_stages=[haplotype])
    select_snp = stage(v.SelectVariants, region_bed=region_bed, previous_stages=[genotype])
    filter_snp = stage(v.VariantFiltration, region_bed=region_bed, previous_stages=[select_snp])
    vcfstats = stage(qc.VCFStats, vcf_file=g.filter_snp_vcf + '.gz', previous_stages=[filter_snp])
    data_output = stage(common.SampleDataOutput, previous_stages=[vcfstats], output_fileset='non_human_qc')
    cleanup = stage(common.Cleanup, previous_stages=[data_output])
    review = stage(common.SampleReview, previous_stages=[cleanup])
    return review
