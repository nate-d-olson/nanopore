import attr
import click
import datetime
import getpass
import logging
import os
import shutil
import subprocess
import time

from pathlib import Path

FORMAT = '%(asctime)s - %(name)-25s - %(levelname)-5s - %(message)s'
DATEFMT = "%Y-%m-%d %H:%M:%S"
logging.basicConfig(format=FORMAT, level=logging.DEBUG, datefmt=DATEFMT)
logger = logging.getLogger(__name__)

@attr.s
class ReadBatch:
    sample_name = attr.ib()
    run_name = attr.ib()
    read_dir = attr.ib()

@attr.s
class Options:
    data_dir = attr.ib(converter=Path)
    staging_dir = attr.ib(converter=Path)
    server_addr = attr.ib()
    remote_dir = attr.ib(converter=Path)
    user = attr.ib()


def get_run_directories(data_dir, sample_name):
    yield from data_dir.glob(f"*_{sample_name}/fast5")

def get_read_batches(data_dir, sample_name):
    for run_dir in get_run_directories(data_dir, sample_name):
        for sub_dir in run_dir.iterdir():
            if sub_dir.is_dir():
                # could do more checks here, but at least at the moment we don't 
                # expect any other crap here
                # yield sub_dir

                read_batch = ReadBatch(sample_name, run_dir.parts[-2], sub_dir)

                yield read_batch


def should_archive_reads(read_batch):
    reads = list(read_batch.read_dir.glob('*.fast5'))
    num_reads = len(reads)

    if num_reads == 0:
        return False

    if num_reads >= 1000:
        return True

    most_recent = max(read.stat().st_mtime for read in reads)
    most_recent = datetime.datetime.fromtimestamp(most_recent)
    since = datetime.datetime.now() - most_recent

    if since > datetime.timedelta(minutes=5) and num_reads >= 50:
        return True

    if since > datetime.timedelta(minutes=30) and num_reads >= 10:
        return True

    return False

def get_staging_area(staging_dir, read_batch):
    read_chunk = read_batch.read_dir.parts[-1]
    i = 0
    while True:
        cur_staging_dir = \
            staging_dir / read_batch.sample_name / f"{read_batch.run_name}_{read_chunk}_{i}"

        if not os.path.exists(cur_staging_dir):
            break
        i += 1

    os.makedirs(cur_staging_dir)

    return cur_staging_dir


def do_archiving(options, read_batch):
    cur_staging_dir = get_staging_area(options.staging_dir, read_batch)

    logger.info("moving...")
    count = 0
    for fast5 in read_batch.read_dir.glob("*.fast5"):
        # logger.info("renaming:", fast5, "to", cur_staging_dir / fast5.name)
        # break
        fast5.rename(cur_staging_dir / fast5.name)
        count += 1
        if count >= 1000:
            break

    logger.info("tar'ing...")
    cur_tar_file = cur_staging_dir.with_name(cur_staging_dir.name + ".tar")

    # this adjust the behavior of tar on macOS X to ignore
    # the extended attributes, which by default are added to 
    # the tarfile as a ._x file for each file x; see
    # https://superuser.com/questions/61185/why-do-i-get-files-like-foo-in-my-tarball-on-os-x
    tar_cmd = f"COPYFILE_DISABLE=1 tar -cf {cur_tar_file} {cur_staging_dir}"

    logger.info(tar_cmd)
    subprocess.check_call(tar_cmd, shell=True)
    shutil.rmtree(cur_staging_dir)

    logger.info("copying to server...")

    # this should run separately, so we rsync everything in case there were
    # any connection errors previously, or data lost on server, etc

    cur_remote_dir = options.remote_dir / read_batch.sample_name / "fast5"
    mkdir_and_rsync_cmd = f"mkdir -p {cur_remote_dir} && rsync"
    rsync_cmd = f'rsync -aP --rsync-path "{mkdir_and_rsync_cmd}" ' \
                f'{cur_tar_file} {options.user}@{options.server_addr}:{cur_remote_dir}'

    logger.info(rsync_cmd)
    subprocess.check_call(rsync_cmd, shell=True)
    
    logger.info("...done")




def run_copy(options, sample_name):
    read_batches = sorted(get_read_batches(options.data_dir, sample_name))
    for read_batch in read_batches:
        if should_archive_reads(read_batch):
            print(read_batch.read_dir)
            do_archiving(options, read_batch)

    return len(read_batches)




@click.command()
@click.argument("data-dir")#, help="the reads directory, eg /Library/MinKNOW/data/reads")
@click.argument("sample-name")
@click.argument("staging-dir")
@click.option("--event-loop-interval", default=5, type=float)
@click.option("--server", default="oak-dtn.stanford.edu")
@click.option("--remote-dir", default="/oak/stanford/groups/msalit/nspies/nanopore/raw")
@click.option("--user", default=getpass.getuser())
def main(data_dir, sample_name, staging_dir, event_loop_interval, server, remote_dir, user):

    # data_dir = Path(data_dir)
    # staging_dir = Path(staging_dir)
    options = Options(data_dir, staging_dir, server, remote_dir, user)
    # poor man's main loop
    t0 = 0

    while True:
        t1 = time.time()
        if t1 - t0 < event_loop_interval:
            time.sleep(event_loop_interval - (t1-t0))
            t1 = time.time()

        logger.info("looping")
        count = run_copy(options, sample_name)

        if count == 0:
            logger.info(f"No read directories found for sample name {sample_name} in "
                        f"data directory {data_dir}; did you type these in correctly?")

        t0 = t1

if __name__ == '__main__':
    main()