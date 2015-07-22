__author__ = 'mwham'
import pytest
import sys


def main(profile, coverage):
    if profile:
        import cProfile
        pr = cProfile.Profile()
    else:
        pr = Null()

    if coverage:
        import coverage
        cov = coverage.coverage(omit='*site-packages*', branch=True)
    else:
        cov = Null()

    cov.start()
    pr.enable()
    pytest.main()
    pr.disable()
    cov.stop()

    if type(pr) is not Null:
        import pstats
        ps = pstats.Stats(pr, stream=sys.stdout).sort_stats('cumulative')
        print('\n================== reporting profiling ==================')
        ps.print_stats(20)

    if type(cov) is not Null:
        print('\n=================== reporting coverage ===================')
        cov.report(file=sys.stdout)


class Null:
    @staticmethod
    def start():
        pass

    @staticmethod
    def stop():
        pass

    @staticmethod
    def enable():
        pass

    @staticmethod
    def disable():
        pass


if __name__ == '__main__':
    main(True, True)
