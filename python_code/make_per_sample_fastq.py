#!/usr/bin/env python

from collections import defaultdict
from os.path import join, exists
from os import mkdir

import click
from skbio.parse.sequences.fastq import parse_fastq

def write_fastq(output_file, fastq_data, ascii_increment=33):
    """Writes tuples of (defline, seq, qual) to an output file

    Parameters
    ----------
    output_file : file
        Data in `fastq_data` will be written to this file
    fastq_data : iterable of tuples of (str, str, str)
        Each item in `fastq_data` is a FASTQ entry, each of which is
        represented by a tuple of (defline, sequence, quality).
    """
    for defline, seq, qual in fastq_data:
        qual = ''.join([chr(x+ascii_increment) for x in qual])
        output_file.write('@%s\n%s\n+\n%s\n' % (defline, seq, qual))

def split_helper(input_fastq, output_directory, sequence_buffer_size=1000,
                 ascii_increment=33):
    """Splits a demultiplexed FASTQ file into per-sample FASTQ files

    Parameters
    ----------
    input_fastq : file
        The input demultiplexed FASTQ file.
    output_directory : str
        Path to the output directory. It will be created if it does not already
        exist.
    sequence_buffer_size : int
        The number of sequences to hold in memory for each sample before
        writing them to disk.

    Notes
    -----
    The sequence identifiers in `input_fastq` should be of the form output by
    QIIME's demultiplexing scripts; namely, they should be:
    ``SampleID_SequenceNumber And Additional Notes if Applicable``
    """
    if not exists(output_directory):
        mkdir(output_directory)

    per_sample_seqs = defaultdict(list)
    per_sample_counts = defaultdict(int)
    for defline, seq, qual in parse_fastq(input_fastq):
        label = defline.split()[0]
        sample_name, sequence_number = label.rsplit('_', 1)

        per_sample_seqs[sample_name].append((defline, seq, qual))
        per_sample_counts[sample_name] += 1

        if per_sample_counts[sample_name] > sequence_buffer_size:
            with open(join(output_directory, sample_name+'.fastq'), 'a') \
                    as outfile:
                write_fastq(outfile, per_sample_seqs[sample_name])

            per_sample_seqs[sample_name] = []
            per_sample_counts[sample_name] = 0
    
    for sample_name, entries in per_sample_seqs.iteritems():
        if not entries:
            continue

        with open(join(output_directory, sample_name+'.fastq'), 'a') \
                as outfile:
            write_fastq(outfile, per_sample_seqs[sample_name])

@click.group()
def cli():
    pass

@cli.command()
@click.option('--input-fastq', required=True, type=click.File('r'))
@click.option('--sequence-buffer-size', default=1000)
@click.option('--output-directory', required=True, type=str)
@click.option('--ascii-increment', default=33)
def split(input_fastq, output_directory, sequence_buffer_size,
          ascii_increment):
    split_helper(input_fastq, output_directory, sequence_buffer_size,
                 ascii_increment)

if __name__ == '__main__':
    cli()
