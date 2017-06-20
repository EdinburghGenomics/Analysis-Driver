import os
from egcg_core import clarity, util, archive_management
from egcg_core.app_logging import logging_default as log_cfg

app_logger = log_cfg.get_logger(__name__)


def create_output_links(sample_id, input_dir, output_cfg, link_dir):
    exit_status = 0
    user_sample_id = clarity.get_user_sample_name(sample_id, lenient=True)

    links = []

    for output_record in output_cfg.content.values():
        src_pattern = os.path.join(
            input_dir,
            os.path.join(*output_record.get('location', [''])),
            output_record['basename']
        ).format(sample_id=sample_id, user_sample_id=user_sample_id)

        sources = util.find_files(src_pattern)
        if sources:
            source = sources[-1]
            link_file = os.path.join(
                link_dir,
                output_record.get('new_name', os.path.basename(source))
            ).format(user_sample_id=user_sample_id)
            if not os.path.islink(link_file):
                os.symlink(source, link_file)
            links.append(link_file)
        else:
            app_logger.warning('No files found for pattern ' + src_pattern)
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
