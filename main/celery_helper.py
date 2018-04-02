import zipfile
import bz2
import gzip
import io
import os
from datetime import date
import logging
import subprocess
import tempfile

logger = logging.getLogger(__name__)

VCF_FIELDS = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER',
              'INFO', 'FORMAT', 'ANCESTRYDNA_DATA']

CHROM_ORDER = {
    'chr1': '1',
    'chr2': '2',
    'chr3': '3',
    'chr4': '4',
    'chr5': '5',
    'chr6': '6',
    'chr7': '7',
    'chr8': '8',
    'chr9': '9',
    'chr10': '10',
    'chr11': '11',
    'chr12': '12',
    'chr13': '13',
    'chr14': '14',
    'chr15': '15',
    'chr16': '16',
    'chr17': '17',
    'chr18': '18',
    'chr19': '19',
    'chr20': '20',
    'chr21': '21',
    'chr22': '22',
    'chrX': '23',
    'chrY': '24',
    'chrM': '25',
    'chrMT': '25',
    'Chr1': '1',
    'Chr2': '2',
    'Chr3': '3',
    'Chr4': '4',
    'Chr5': '5',
    'Chr6': '6',
    'Chr7': '7',
    'Chr8': '8',
    'Chr9': '9',
    'Chr10': '10',
    'Chr11': '11',
    'Chr12': '12',
    'Chr13': '13',
    'Chr14': '14',
    'Chr15': '15',
    'Chr16': '16',
    'Chr17': '17',
    'Chr18': '18',
    'Chr19': '19',
    'Chr20': '20',
    'Chr21': '21',
    'Chr22': '22',
    'ChrX': '23',
    'ChrY': '24',
    'ChrM': '25',
    'ChrMT': '25',
    '1': '1',
    '2': '2',
    '3': '3',
    '4': '4',
    '5': '5',
    '6': '6',
    '7': '7',
    '8': '8',
    '9': '9',
    '10': '10',
    '11': '11',
    '12': '12',
    '13': '13',
    '14': '14',
    '15': '15',
    '16': '16',
    '17': '17',
    '18': '18',
    '19': '19',
    '20': '20',
    '21': '21',
    '22': '22',
    'X': '23',
    'Y': '24',
    'M': '25',
    'MT': '25',
}


def sort_vcf(input_file):
    outputfile = tempfile.TemporaryFile()
    sortingfile = tempfile.TemporaryFile()
    next_line = input_file.read()
    while next_line and next_line.startswith('#'):
        outputfile.write(next_line.encode())
        try:
            next_line = input_file.read()
        except StopIteration:
            next_line = None
            break
    while next_line:
        for key in CHROM_ORDER:
            if next_line.startswith(key + '\t'):
                sortingfile.write(CHROM_ORDER[key] + '\t' + next_line)
                break
        try:
            next_line = input_file.read()
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
        outputfile.write(line.encode())
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
