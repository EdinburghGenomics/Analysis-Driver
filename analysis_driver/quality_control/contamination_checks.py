from analysis_driver.config import default as cfg


class ContaminationCheck():

    def __init__(self, fastq_files, sample_id):
        self.fastq_files = fastq_files
        self.sample_id = sample_id
        self.contamination_cfg = cfg.get('contamination-check')


    def _kontaminant_command(self, referenceSpecies):
        """
        :param referenceSpecies: a string passed from _loop_through_references
        takes the name of a reference species (referenceSpecies) and returns a kontnaminant command for
        read filtering using your reads and the prepared referenceSpecies kmer database
        :return str: the command used to run Kontaminant
        """
        if len(self.fastq_files) == 1:
            kontaminant_command = ("This is my single end kontaminant command for %s" % referenceSpecies)
            return kontaminant_command
        elif len(self.fastq_files) == 2:
            kontaminant_command = ("This is my paired end kontaminant command for %s" % referenceSpecies)
            return kontaminant_command
        else:
           raise ValueError('Bad number of fastqs: ' + str(self.fastq_files))



    def _loop_through_references(self):
        """
        :return: dict: a dict of species:species specific kontaminant commands
        """
        kontaminantReferences = self.contamination_cfg.get('kontaminant_references')
        kontaminantCommands = {}
        for species in kontaminantReferences:
            speciesCommand = self._kontaminant_command(species)
            kontaminantCommands[species] = (speciesCommand)
        return(kontaminantCommands)









