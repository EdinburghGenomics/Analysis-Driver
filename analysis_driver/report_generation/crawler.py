from egcg_core import clarity
from egcg_core import constants as c
from egcg_core.app_logging import AppLogger
from analysis_driver.quality_control.sex_check import sex_alias


class Crawler(AppLogger):
    @classmethod
    def get_sample_information_from_lims(cls, sample_name):
        lims_sample = clarity.get_sample(sample_name)
        sample_info = {
            c.ELEMENT_SAMPLE_EXTERNAL_ID: clarity.get_user_sample_name(sample_name, lenient=True),
            c.ELEMENT_SAMPLE_PLATE: clarity.get_plate_id_and_well(sample_name)[0],  # returns [plate_id, well]
            c.ELEMENT_SEX_VALIDATION: {c.ELEMENT_PROVIDED_SEX: sex_alias(clarity.get_sample_sex(sample_name))},
            c.ELEMENT_SAMPLE_SPECIES: clarity.get_species_from_sample(sample_name)
        }
        if 'Yield for Quoted Coverage (Gb)' in lims_sample.udf:
            sample_info[c.ELEMENT_SAMPLE_REQUIRED_YIELD_Q30] = lims_sample.udf.get('Yield for Quoted Coverage (Gb)') * 1000000000
        if 'Required Yield (Gb)' in lims_sample.udf:
            sample_info[c.ELEMENT_SAMPLE_REQUIRED_YIELD] = lims_sample.udf.get('Required Yield (Gb)') * 1000000000
        if 'Coverage (X)' in lims_sample.udf:
            sample_info[c.ELEMENT_SAMPLE_REQUIRED_COVERAGE] = lims_sample.udf.get('Coverage (X)')

        return sample_info
