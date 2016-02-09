import os
from analysis_driver.config import default as cfg
from analysis_driver import executor
from analysis_driver.notification import default as ntf


class ContaminationCheck():

    def __init__(self, fastq_files, sample_id):
        self.fastq_files = fastq_files
        self.sample_id = sample_id
        self.contamination_cfg = cfg.get('contamination-check')
        self.workDir = self.contamination_cfg.get('kontaminant_running')


    def createOutputDirectories(self,species):
        if os.path.exists(os.path.join(self.workDir, self.sample_id, species)) == True:
            os.rmdir(os.path.join(self.workDir, self.sample_id, species))
            outputDirectory = (os.path.join(self.workDir, self.sample_id, species))
        else:
            outputDirectory = (os.path.join(self.workDir, self.sample_id, species))
        return outputDirectory




    def _kontaminant_command(self, refSpecies, outputDirectory):
        """
        :param referenceSpecies: a string passed from _loop_through_references
        takes the name of a reference species (referenceSpecies) and returns a kontaminant command for
        read filtering using your reads and the prepared referenceSpecies kmer database
        :return str: the command used to run Kontaminant
        """
        if len(self.fastq_files) == 1:

            kontaminantBin = (self.contamination_cfg.get('kontaminant_bin').rstrip('/') + '/')
            referenceDirectory = (self.contamination_cfg.get('kontaminant_reference_directory').rstrip('/') + '/')
            referencePath = os.path.join(referenceDirectory,refSpecies)
            fastq1 = self.fastq_files[0]

            kontaminant_command = ("{}kmer_filter_31"
                                   " --reference {} "
                                   "--read1 {} "
                                   "--output_prefix filtered_K21_Thr1 "
                                   "--output_folder {} "
                                   "--threashold 1".format(kontaminantBin,
                                                           referencePath,
                                                           fastq1,
                                                           outputDirectory))
            return kontaminant_command

        elif len(self.fastq_files) == 2:

            kontaminantBin = (self.contamination_cfg.get('kontaminant_bin').rstrip('/') + '/')
            referenceDirectory = (self.contamination_cfg.get('kontaminant_reference_directory').rstrip('/') + '/')
            referencePath = os.path.join(referenceDirectory,refSpecies)
            fastq1 = self.fastq_files[0]
            fastq2 = self.fastq_files[1]

            kontaminant_command = ("{}kmer_filter_31"
                                   " --reference {} "
                                   "--read1 {} "
                                   "--read2 {} "
                                   "--output_prefix filtered_K21_Thr1 "
                                   "--output_folder {} "
                                   "--threashold 1".format(kontaminantBin,
                                                           referencePath,
                                                           fastq1,
                                                           fastq2,
                                                           outputDirectory))
            return kontaminant_command
        else:
           raise ValueError('Bad number of fastqs: ' + str(self.fastq_files))



    def _loop_through_references(self):
        """
        :return: list: a list of species specific kontaminant commands
        """
        kontaminantReferences = self.contamination_cfg.get('kontaminant_references')
        kontaminantCommands = {}
        for species in kontaminantReferences:
            outputDirectory = self.createOutputDirectories(species)
            speciesCommand = self._kontaminant_command(species, outputDirectory)
            kontaminantCommands[species] = (speciesCommand)
        return(kontaminantCommands)






    def _runKontaminant(self):
        kontaminantCommands = self._loop_through_references()
        expectedOutfiles = {}
        commands = []
        for species in kontaminantCommands:
            expectedOutfile = (os.path.join(self.workDir, self.sample_id, species, 'filtered_K21_Thr1' + species + '.FASTQ',))
            expectedOutfiles[species] = expectedOutfile
            commands.append(kontaminantCommands[species])
        ntf.start_stage('kontaminant_contamination_check')
        kontaminant_executor = executor.execute(
            commands,
            job_name='kontaminant',
            sample_id=self.sample_id,
            cpus=4,
            mem=25
        )
        exit_status = kontaminant_executor.join()
        ntf.end_stage('kontaminant_contamination_check', exit_status)

        return expectedOutfiles



