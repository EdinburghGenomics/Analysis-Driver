import os
import yaml

home = os.getenv('HOME')
user_config_file = os.path.join(home, '.analysisdriver.yaml')
local_config_file = os.path.join(os.path.dirname(__file__), 'example_analysisdriver.yaml')

if os.path.exists(user_config_file):
    config_file = user_config_file
else:
    config_file = local_config_file

run_env = os.getenv('ANALYSISDRIVERENV')
if not run_env:
    run_env = 'testing'

config = yaml.load(open(config_file, 'r'))[run_env]

work_home = config['shared']['work_home']
fastq = config['shared']['fastq_dir']
jobs = config['shared']['jobs_dir']

job_execution = config['analysisdriver']['job_execution']

