import os.path
from .pbs_writer import PBSWriter


class BCL2FastqPBSWriter(PBSWriter):
    def __init__(self, pbs_name, job_name, log_file, walltime='24', cpus='12', mem='24'):
        super().__init__(pbs_name, walltime, cpus, mem, job_name, log_file)

    def _bcl2fastq(self, mask, input_dir, fastq_path):
        self.write_line(
            'bcl2fastq -l INFO --runfolder-dir %s --output-dir %s --sample-sheet %s --use-bases-mask %s\n' % (
                input_dir, fastq_path, os.path.join(input_dir, 'SampleSheet.csv'), mask
            )
        )

    @staticmethod
    def _generate_mask(mask):  # mask = [this, that other]
        chain = []

        # mask is a list of triplets, so we reverse step in threes
        current_read_number = 0
        for i in range(len(mask)-3, -1, -3):

            read_number = mask[i]
            num_cycles = mask[i+1]
            is_indexed_read = mask[i+2]

            # build the mask string in reverse
            end_char = ''
            if read_number != current_read_number:
                # New read
                current_read_number = read_number
                end_char = 'n'

            mask_part = num_cycles + end_char

            if is_indexed_read.lower() == 'n':
                mask_part = 'y' + mask_part
            else:
                mask_part = 'i' + mask_part

            chain.append(mask_part)
        print(chain)

        return ','.join(chain[::-1])

    def write(self, mask, input_dir, fastq_path):
        self._bcl2fastq(mask, input_dir, fastq_path)
        self.save()


# unit test
if __name__ == '__main__':
    import sys
    w = BCL2FastqPBSWriter(sys.stdout, '1', '1', '1', 'test', sys.stdout)
    # generate fake mask list
    mask_list = ['1', '128', 'Y', '2', '8', 'N', '3', '128', 'Y']
    # call to generate string
    mask_string = w._generate_mask(mask_list)
    print(mask_string)
    # test PBS script generation for BCL2FASTQ
    w.write()
