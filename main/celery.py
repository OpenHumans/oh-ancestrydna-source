from django.conf import settings
import os
from celery import Celery
import tempfile
import json
from ohapi import api
import requests

import bz2
import logging
import re
import shutil

from io import StringIO
from datetime import datetime

import arrow
from .celery_helper import vcf_header, temp_join, open_archive, sort_vcf
from .celery_helper import HEADER_V1, HEADER_V2, HEADER_V3

# set the default Django settings module for the 'celery' program.
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'oh_data_uploader.settings')

OH_BASE_URL = settings.OPENHUMANS_OH_BASE_URL

REF_ANCESTRYDNA_FILE = os.path.join(
    os.path.dirname(__file__), 'references/reference_b37.txt')

# Was used to generate reference genotypes in the previous file.
REFERENCE_GENOME_URL = ('http://hgdownload-test.cse.ucsc.edu/' +
                        'goldenPath/hg19/bigZips/hg19.2bit')

VCF_FIELDS = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER',
              'INFO', 'FORMAT', 'ANCESTRYDNA_DATA']


# The only non-commented-out header line. We want to ignore it.
EXPECTED_COLUMNS_HEADER = 'rsid\tchromosome\tposition\tallele1\tallele2'

CHROM_MAP = {
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
    '23': 'X',
    '24': 'Y',
    '25': 'X',
}

logger = logging.getLogger(__name__)

app = Celery('proj')

# Using a string here means the worker doesn't have to serialize
# the configuration object to child processes.
# - namespace='CELERY' means all celery-related configuration keys
#   should have a `CELERY_` prefix.
app.config_from_object('django.conf:settings', namespace='CELERY')
app.conf.update(CELERY_BROKER_URL=os.environ['REDIS_URL'],
                CELERY_RESULT_BACKEND=os.environ['REDIS_URL'])

# Load task modules from all registered Django app configs.
app.autodiscover_tasks()
# app.autodiscover_tasks(lambda: settings.INSTALLED_APPS)


def read_reference(ref_file):
    reference = dict()

    with open(ref_file) as f:
        for line in f:
            data = line.rstrip().split('\t')

            if data[0] not in reference:
                reference[data[0]] = dict()

            reference[data[0]][data[1]] = data[2]
    return reference


def check_header_lines(input_lines, header_lines, header_name):
    if not len(input_lines) == len(header_lines):
        if header_name:
            logger.debug("Header line count != {}".format(header_name))
        return False
    matched_lines = [
        header_lines[i] == input_lines[i] for i in range(len(header_lines))]
    if header_name:
        logger.debug("Header line matching for {}: {}".format(
            header_name, matched_lines))
    return all(matched_lines)


def vcf_from_raw_ancestrydna(raw_ancestrydna, genome_sex):
    output = StringIO()

    reference = read_reference(REF_ANCESTRYDNA_FILE)

    header = vcf_header(
        source='open_humans_data_importer.ancestry_dna',
        reference=REFERENCE_GENOME_URL,
        format_info=['<ID=GT,Number=1,Type=String,Description="Genotype">'])
    for line in header:
        output.write(line + '\n')
    for line in raw_ancestrydna:
        # Skip header
        if line.startswith('#'):
            continue
        if line == EXPECTED_COLUMNS_HEADER:
            continue

        data = line.rstrip().split('\t')

        # Skip uncalled and genotyping without explicit base calls
        if not re.match(r'^[ACGT]$', data[3]):
            continue
        if not re.match(r'^[ACGT]$', data[4]):
            continue
        vcf_data = {x: '.' for x in VCF_FIELDS}

        # Chromosome. Determine correct reporting according to genome_sex.
        try:
            vcf_data['REF'] = reference[data[1]][data[2]]
        except KeyError:
            continue
        vcf_data['CHROM'] = CHROM_MAP[data[1]]
        if data[1] == '24' and genome_sex == 'Female':
            continue
        if data[1] in ['23', '24'] and genome_sex == 'Male':
            alleles = data[3]
        else:
            alleles = data[3] + data[4]

        # Position, dbSNP ID, reference. Skip if we don't have ref.
        vcf_data['POS'] = data[2]
        if data[0].startswith('rs'):
            vcf_data['ID'] = data[0]

        # Figure out the alternate alleles.
        alt_alleles = []
        for alle in alleles:
            if alle != vcf_data['REF'] and alle not in alt_alleles:
                alt_alleles.append(alle)
        if alt_alleles:
            vcf_data['ALT'] = ','.join(alt_alleles)
        else:
            vcf_data['ALT'] = '.'
            vcf_data['INFO'] = 'END=' + vcf_data['POS']

        # Get allele-indexed genotype.
        vcf_data['FORMAT'] = 'GT'
        all_alleles = [vcf_data['REF']] + alt_alleles
        genotype_indexed = '/'.join([str(all_alleles.index(x))
                                     for x in alleles])
        vcf_data['ANCESTRYDNA_DATA'] = genotype_indexed
        output_line = '\t'.join([vcf_data[x] for x in VCF_FIELDS])
        output.write(output_line + '\n')

    return output


def clean_raw_ancestrydna(closed_input_file):
    """
    Create clean file in AncestryDNA format from downloaded version
    Obsessively careful processing that ensures AncestryDNA file format changes
    won't inadvertantly result in unexpected information, e.g. names.
    """
    inputfile = open_archive(closed_input_file)

    output = StringIO()

    header_l1 = inputfile.readline()
    logger.warn(header_l1)
    expected_header_l1 = '#AncestryDNA raw data download'
    if header_l1.rstrip() == expected_header_l1:
        output.write(header_l1)
    dateline = inputfile.readline()
    re_datetime_string = (r'([0-1][0-9]/[0-3][0-9]/20[1-9][0-9] ' +
                          r'[0-9][0-9]:[0-9][0-9]:[0-9][0-9]) MDT')
    if re.search(re_datetime_string, dateline.rstrip()):
        datetime_string = re.search(re_datetime_string, dateline.rstrip()).groups()[0]

        datetime_ancestrydna = datetime.strptime(datetime_string,
                                                 '%m/%d/%Y %H:%M:%S')

        output.write(
            '#This file was generated by AncestryDNA at: {}'.format(
                datetime_ancestrydna.strftime('%a %b %d %H:%M:%S %Y MDT')))

    re_array_version = (
        r'#Data was collected using AncestryDNA array version: V\d\.\d')

    header_array_version = inputfile.readline()

    if re.match(re_array_version, header_array_version):
        output.write(header_array_version)

    re_converter_version = (
        r'#Data is formatted using AncestryDNA converter version: V\d\.\d')

    header_converter_version = inputfile.readline()

    if re.match(re_converter_version, header_converter_version):
        output.write(header_converter_version)

    next_line = inputfile.readline()
    header_p_lines = []

    while next_line.startswith('#'):
        header_p_lines.append(next_line.rstrip())
        next_line = inputfile.readline()

    if check_header_lines(header_p_lines, HEADER_V1, 'HEADER_V1'):
        for line in HEADER_V1:
            output.write(line+'\r\n')
    elif check_header_lines(header_p_lines, HEADER_V2, 'HEADER_V2'):
        for line in HEADER_V2:
            output.write(line+'\r\n')
    elif check_header_lines(header_p_lines, HEADER_V3, 'HEADER_V3'):
        for line in HEADER_V3:
            output.write(line+'\r\n')
    else:
        logger.warn("AncestryDNA header didn't match expected formats")

    data_header = next_line
    if data_header.rstrip() == EXPECTED_COLUMNS_HEADER:
        output.write(EXPECTED_COLUMNS_HEADER)

    next_line = inputfile.readline()
    bad_format = False
    # AncestryDNA always reports two alleles for all X and Y positions.
    # For XY individuals, haplozygous positions are redundantly reported.
    # For XX individuals this means Y positions are "0".
    # Note the above two statements are not ALWAYS true! The raw data
    # ocassionally reports 'heterozygous' calls for X and Y in XY individuals,
    # and Y calls in XX individuals. So our test is forgiving of these.
    genome_sex = 'Female'
    called_Y = 0
    reported_Y = 0

    LINE_RE = re.compile(
        r'(rs|VGXS)[0-9]+\t[1-9][0-9]?\t[0-9]+\t[ACGTDI0]\t[ACGTDI0]')
    REPORTED_Y = re.compile(r'(rs|VGXS)[0-9]+\t24\t[0-9]+\t[ACGTDI0]\t[ACGTDI0]')
    CALLED_Y = re.compile(r'(rs|VGXS)[0-9]+\t24\t[0-9]+\t[ACGTDI]\t[ACGTDI]')

    while next_line:
        if LINE_RE.match(next_line):
            if REPORTED_Y.match(next_line):
                reported_Y += 1

                if CALLED_Y.match(next_line):
                    called_Y += 1

            output.write(next_line)
        else:
            # Only report this type of format issue once.
            if not bad_format:
                bad_format = True
                logger.warn('AncestryDNA body did not conform to expected format.')
                logger.warn('Bad format: "%s"', next_line)

        try:
            next_line = inputfile.readline()
        except StopIteration:
            next_line = None

    if reported_Y == 0:
        genome_sex = 'Female'
    elif called_Y / reported_Y > 0.5:
        genome_sex = 'Male'

    return output, genome_sex


def process_file(dfile, access_token, member, metadata):
    try:
        infile_suffix = dfile['basename'].split(".")[-1]
        tf_in = tempfile.NamedTemporaryFile(suffix="."+infile_suffix)
        tf_in.write(requests.get(dfile['download_url']).content)
        tf_in.flush()
        tmp_directory = tempfile.mkdtemp()
        filename_base = 'AncestryDNA-genotyping'
        raw_ancestry, chr_sex = clean_raw_ancestrydna(tf_in)
        raw_ancestry.seek(0)
        vcf_ancestry_unsort = vcf_from_raw_ancestrydna(raw_ancestry, chr_sex)

        # Save raw Ancestry genotyping to temp file.
        raw_filename = filename_base + '.txt'
        raw_filename = temp_join(tmp_directory, raw_filename)
        metadata = {
                    'description':
                    'AncestryDNA full genotyping data, original format',
                    'tags': ['AncestryDNA', 'genotyping'],
                    'creation_date': arrow.get().format(),
            }
        with open(raw_filename, 'w') as raw_file:
            raw_ancestry.seek(0)
            shutil.copyfileobj(raw_ancestry, raw_file)
            raw_file.flush()

        api.upload_aws(raw_filename, metadata,
                       access_token, base_url=OH_BASE_URL,
                       project_member_id=str(member['project_member_id']))

        # Save VCF Ancestry genotyping to temp file.
        vcf_filename = filename_base + '.vcf.bz2'
        vcf_filename = temp_join(tmp_directory, vcf_filename)

        metadata = {
            'description': 'AncestryDNA full genotyping data, VCF format',
            'tags': ['AncestryDNA', 'genotyping', 'vcf'],
            'creation_date': arrow.get().format()
        }

        vcf_ancestry_unsort.seek(0)
        vcf_ancestry = sort_vcf(vcf_ancestry_unsort)

        with bz2.BZ2File(vcf_filename, 'w') as vcf_file:
            vcf_ancestry.seek(0)
            for i in vcf_ancestry:
                vcf_file.write(i)

        api.upload_aws(vcf_filename, metadata,
                       access_token, base_url=OH_BASE_URL,
                       project_member_id=str(member['project_member_id']))

    except:
        api.message("AncestryDNA integration: A broken file was deleted",
                    "While processing your AncestryDNA file "
                    "we noticed that your file does not conform "
                    "to the expected specifications and it was "
                    "thus deleted. Please make sure you upload "
                    "the right file:\nWe expect the file to be a "
                    "single txt file (either unzipped, bz2 zipped or gzipped) "
                    "or a .zip file that contains a single txt file (this is "
                    " what you can download from Ancestry right away) Please "
                    "do not alter the original txt file, as unexpected "
                    "additions can invalidate the file.",
                    access_token, base_url=OH_BASE_URL)
        raise

    finally:
        api.delete_file(access_token,
                        str(member['project_member_id']),
                        file_id=str(dfile['id']),
                        base_url=OH_BASE_URL)


@app.task(bind=True)
def clean_uploaded_file(self, access_token, file_id):
    member = api.exchange_oauth2_member(access_token, base_url=OH_BASE_URL)
    for dfile in member['data']:
        if dfile['id'] == file_id:
            process_file(dfile, access_token, member, dfile['metadata'])
    pass
