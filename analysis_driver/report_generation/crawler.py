from egcg_core import clarity
from egcg_core import constants as c
from egcg_core.app_logging import AppLogger


class Crawler(AppLogger):
    _gender_aliases = {'female': ['f', 'female', 'girl', 'woman'], 'male': ['m', 'male', 'boy', 'man']}

    @classmethod
    def gender_alias(cls, gender):
        for key in cls._gender_aliases:
            if str(gender).lower() in cls._gender_aliases[key]:
                return key
        return 'unknown'

    @classmethod
    def get_sample_information_from_lims(cls, sample_name):
        sample_info = {
            c.ELEMENT_SAMPLE_EXTERNAL_ID: clarity.get_user_sample_name(sample_name, lenient=True),
            c.ELEMENT_SAMPLE_PLATE: clarity.get_plate_id_and_well(sample_name)[0],  # returns [plate_id, well]
            c.ELEMENT_PROVIDED_GENDER: cls.gender_alias(clarity.get_sample_gender(sample_name)),
            c.ELEMENT_SAMPLE_SPECIES: clarity.get_species_from_sample(sample_name),
            c.ELEMENT_SAMPLE_EXPECTED_YIELD: clarity.get_expected_yield_for_sample(sample_name)
        }
        lims_sample = clarity.get_sample(sample_name)
        coverage = lims_sample.udf.get('Coverage')
        if coverage:
            sample_info[c.ELEMENT_SAMPLE_EXPECTED_COVERAGE] = coverage

        return sample_info