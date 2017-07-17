import os
from egcg_core import util, archive_management
from egcg_core.app_logging import logging_default as log_cfg

app_logger = log_cfg.get_logger(__name__)


def create_output_links(input_dir, output_cfg, link_dir, **kwargs):
    exit_status = 0
    links = []

    for output_record in output_cfg.content.values():
        src_pattern = os.path.join(
            input_dir,
            os.path.join(*output_record.get('location', [''])),
            output_record['basename']
        ).format(**kwargs)

        source = util.find_file(src_pattern)
        if source:
            link_file = os.path.join(
                link_dir,
                output_record.get('new_name', os.path.basename(source))
            ).format(**kwargs)
            if os.path.islink(link_file):
                os.unlink(link_file)
            os.symlink(source, link_file)
            links.append(link_file)
        else:
            app_logger.warning('No file found for pattern ' + src_pattern)
            if output_record.get('required', True):
                exit_status += 1
    if exit_status == 0:
        return links
    else:
        app_logger.error('link creation failed with exit status ' + str(exit_status))


def output_data_and_archive(source_dir, output_dir):
    exit_status = util.move_dir(source_dir, output_dir)
    if exit_status == 0:
        if archive_management.archive_directory(output_dir):
            return 0
        else:
            return 1
    return exit_status
