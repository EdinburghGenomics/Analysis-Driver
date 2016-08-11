import luigi
from os.path import join
from time import sleep, ctime
from analysis_driver.config import default as cfg, load_config
import analysis_driver.pipeline
from analysis_driver.dataset import RunDataset


class ExampleStage(analysis_driver.pipeline.Stage):
    def _run(self):
        print(self.stage_name)
        print('starting')
        sleep(5)
        with open(join(self.job_dir, self.stage_name + '.txt'), 'a') as f:
            f.write('Finished %s on %s' % (self.stage_name, ctime()))
        print('done')
        return 0


class Stage1(ExampleStage):
    pass


class Stage2(ExampleStage):
    previous_stages = Stage1


class Stage3(ExampleStage):
    previous_stages = Stage1


class Stage4(ExampleStage):
    previous_stages = (Stage2, Stage3)


def pipeline(dataset):
    analysis_driver.pipeline._dataset = dataset
    final_stage = Stage4()
    luigi.build(
        [final_stage],
        local_scheduler=True
    )
    return final_stage.exit_status

if __name__ == '__main__':
    load_config()
    cfg.merge(cfg['run'])
    print(
        pipeline(
            RunDataset(
                '150723_E00306_0025_BHCHK3CCXX',
                cfg['input_dir'] + '/150424_E00307_0017_AH3KGTCCXX',
                False,
                {'proc_id': 'run_150723_E00306_0025_BHCHK3CCXX_15_02_2016_12:14:31'}
            )
        )
    )
