import zipfile
import bz2
import gzip
import io
import os
from datetime import date
import logging
import subprocess
import tempfile
from .vcf_helper import VCF_FIELDS, CHROM_ORDER
logger = logging.getLogger(__name__)


def sort_vcf(input_file):
    outputfile = tempfile.TemporaryFile()
    sortingfile = tempfile.TemporaryFile()
    next_line = input_file.readline()
    while next_line and next_line.startswith('#'):
        outputfile.write(next_line.encode())
        try:
            next_line = input_file.readline()
        except StopIteration:
            next_line = None
            break
    while next_line:
        for key in CHROM_ORDER:
            if next_line.startswith(key + '\t'):
                out_line = "{}\t{}".format(CHROM_ORDER[key], next_line)
                sortingfile.write(out_line.encode())
                break
        try:
            next_line = input_file.readline()
        except StopIteration:
            next_line = None
            break
    sortingfile.seek(0)
    sort_proc = subprocess.Popen(['sort', '-k', '1n,1', '-k', '3n,3'],
                                 stdin=sortingfile,
                                 stdout=subprocess.PIPE)
    cut_proc = subprocess.Popen(['cut', '-f', '2-'],
                                stdin=sort_proc.stdout,
                                stdout=subprocess.PIPE)
    for line in cut_proc.stdout:
        outputfile.write(line)
    outputfile.seek(0)
    return outputfile


def vcf_header(source=None, reference=None, format_info=None):
    """Generate a VCF header."""
    header = []
    today = date.today()
    header.append('##fileformat=VCFv4.1')
    header.append('##fileDate=%s%s%s' % (str(today.year),
                                         str(today.month).zfill(2),
                                         str(today.day).zfill(2)))
    if source:
        header.append('##source=' + source)
    if reference:
        header.append('##reference=%s' % reference)
    for item in format_info:
        header.append('##FORMAT=' + item)
    header.append('#' + '\t'.join(VCF_FIELDS))
    return header


def filter_archive(zip_file):
    return [f for f in zip_file.namelist()
            if not f.startswith('__MACOSX/')]


def open_archive(input_file):
    error_message = ("Input file is expected to be either '.txt', "
                     "'.gz', '.bz2', or a single '.txt' file in a "
                     "'.zip' ZIP archive.")
    if input_file.name.endswith('.zip'):
        zip_file = zipfile.ZipFile(input_file)
        zip_files = filter_archive(zip_file)

        if len(zip_files) != 1:
            logger.warn(error_message)
            raise ValueError(error_message)

        return io.TextIOWrapper(zip_file.open(zip_files[0]))
    elif input_file.name.endswith('.gz'):
        return io.TextIOWrapper(gzip.open(input_file.name))
    elif input_file.name.endswith('.bz2'):
        return io.TextIOWrapper(bz2.BZ2File(input_file.name))
    elif input_file.name.endswith('.txt'):
        return open(input_file.name)

    logger.warn(error_message)
    raise ValueError(error_message)


def temp_join(tmp_directory, path):
    return os.path.join(tmp_directory, path)
