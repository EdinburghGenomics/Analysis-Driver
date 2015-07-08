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

fastq = config['shared']['fastq']
jobs = config['shared']['jobs']

job_execution = config['analysisdriver']['job_execution']

