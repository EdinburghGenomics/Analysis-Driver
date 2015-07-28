__author__ = 'mwham'
import pytest
import sys


def main(profile=False, coverage=False, pr_ignore=('site-packages', 'lib', 'tests')):
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
        ps.print_title()

        ignore_func = lambda path: any([ignorable for ignorable in pr_ignore if ignorable in path])
        counter = 0

        width, l = ps.get_print_list(tuple())
        for func in l:
            if ignore_func(func[0]):
                pass
            else:
                ps.print_line(func)
                counter += 1

            if counter > 20:
                break

    if type(cov) is not Null:
        print('\n=================== reporting coverage ===================')
        cov.report(file=sys.stdout, omit='tests/*')


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
    main(
        profile=True,
        coverage=True
    )
