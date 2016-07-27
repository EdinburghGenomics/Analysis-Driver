import luigi
from os.path import join
from time import sleep, ctime
from analysis_driver.pipeline import Stage
from analysis_driver.dataset import NoCommunicationDataset


class ExampleStage(Stage):
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
    previous_stages = (Stage1,)


class Stage3(ExampleStage):
    previous_stages = (Stage1,)


class Stage4(ExampleStage):
    previous_stages = (Stage2, Stage3)


def pipeline(dataset):
    final_stage = Stage4(dataset=dataset)
    luigi.build(
        [final_stage],
        local_scheduler=True
    )
    return final_stage.exit_status

if __name__ == '__main__':
    print(pipeline(NoCommunicationDataset('150424_E00307_0017_AH3KGTCCXX', {'things': 'thangs'})))
