from egcg_core import clarity, rest_communication
from egcg_core import constants as c
from egcg_core.app_logging import AppLogger
from analysis_driver.quality_control.sex_validation import sex_alias


class Crawler(AppLogger):
    @classmethod
    def get_sample_information_from_lims(cls, sample_name):
        rest_data = rest_communication.get_document('lims/sample_info', match={'sample_id': sample_name})
        sample_info = {
            c.ELEMENT_SAMPLE_EXTERNAL_ID: clarity.get_user_sample_name(sample_name),
            c.ELEMENT_SAMPLE_PLATE: clarity.get_plate_id_and_well(sample_name)[0],  # returns [plate_id, well]
            c.ELEMENT_SEX_VALIDATION: {c.ELEMENT_PROVIDED_SEX: sex_alias(clarity.get_sample_sex(sample_name))},
            c.ELEMENT_SAMPLE_SPECIES: clarity.get_species_from_sample(sample_name)
        }
        if 'Yield for Quoted Coverage (Gb)' in rest_data:
            sample_info[c.ELEMENT_SAMPLE_REQUIRED_YIELD_Q30] = int(rest_data.get('Yield for Quoted Coverage (Gb)')) * 1000000000
        if 'Required Yield (Gb)' in rest_data:
            sample_info[c.ELEMENT_SAMPLE_REQUIRED_YIELD] = int(rest_data.get('Required Yield (Gb)')) * 1000000000
        if 'Coverage (X)' in rest_data:
            sample_info[c.ELEMENT_SAMPLE_REQUIRED_COVERAGE] = int(rest_data.get('Coverage (X)'))

        return sample_info
