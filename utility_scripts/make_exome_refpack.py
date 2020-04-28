#!/usr/bin/env python3
"""WARNING: Your targets/baits files may need pre-processing!
The intervals here should be compatible with the reference you are planning
to use. That means the sequence names (first column) must match the refdict
sequence names (SN field). AND The start/stop coordinates must fall
inside the sequence lengths described by the refdict (LN field).
Additionally, padding is added to the targets when creating the filter
files. It is not included in the Picard metrics files. The targets intervals
do NOT need to be padded before starting this process.
Usage:
See argument help with:
    make_exome_refpack.py -h
Typical use case:
    make_exome_refpack.py -t <targets.bed> [-b <baits>] -r <refdict> -g <gtf>
    
Summary:
The goal of this script is to generate the files needed to add a new capture
kit to the TGen pipeline. Three types of files will need to be created:
    Picard .interval_list
    These are used for Picard metrics modules and describe the locations
    targetted by the capture kit. Coverage metrics generated by the pipeline
    will be relative to the regions defined in these files. There is one
    '.interval_list' for the targets and one for the baits. If baits is
    not given, both files will be generated from the targets.
    https://gatkforums.broadinstitute.org/gatk/discussion/1319/collected-faqs-about-interval-lists
    http://broadinstitute.github.io/picard/command-line-overview.html#BedToIntervalList
    Extended .bed
    This file is generated by extending the regions in the targets intervals
    file to create a bed file that can be used for filtering results. This file
    is a union of the following intervals:
        A) Target intervals file padded by 100bp
        B) All exonic intervals (taken from the reference GTF) that
           intersect any interval in "A"
    CNA index .bed
    This file is an index for tCoNut, in bed format, which describes the
    targeted regions in 100bp intervals. This is generated by intersecting
    a template bed file with the extended .bed. The result is an
    index file that describes the locations covered by the kit in 100bp
    increments.
Glossary:
    Targets
    Also called "regions" or "target regions". This is a file that describes the
    regions that the exome kit aims to capture. This is usually distributed by
    the manufacturer of the kit as a .bed file. It may need some pre-processing
    before starting this application.
    https://www.ensembl.org/info/website/upload/bed.html
    Baits
    Also called "probes". This is a file that describes the hybrid capture
    probes in the kit. This is usually distributed by the manufacturer of the
    kit as a .bed file. It may need some pre-processing before starting this
    application.
"""
import os
import argparse
import logging
import subprocess
import tempfile
from collections import OrderedDict
from formats import intervals
from formats.refdict import refdict_to_bedtools_genome


log = logging.getLogger(__name__)


class TempStore(object):
    """A simple container for managing many tempfiles"""
    def __init__(self, return_objects=False):
        self.objs = OrderedDict()

        # Control whether the tempfile object is returned when creating a new 
        # item, by default only the path (tempfile.name) is returned.
        self.return_objects = return_objects

    @property
    def paths(self):
        """Returns an OrderedDict that includes the paths to each object
        :return: OrderedDict
        """
        return OrderedDict([(n, obj.name) for n, obj in self.objs.items()])

    def cleanup(self):
        """Cleanup all objects in this TempStore
        :return: None
        """
        for obj in self.objs.values():
            try:
                obj.cleanup()
            except AttributeError:
                obj.close()
            except FileNotFoundError:
                pass
        self.objs = OrderedDict()

    def save(self, item, outpath=None):
        """Save an item from the tempstore"""
        # Instead of copying the file, this writes the tempfile out
        # to a new file handle. This avoids transferring the strange file
        # permissions of tempfile to the out path"""
        
        if outpath is None:
            outpath = item

        temp = self.objs[item]
        temp.flush()

        log.critical(f'Saving: {item} to {outpath}')
        with open(temp.name, 'rb') as ifp, open(outpath, 'wb') as ofp:
            for b in ifp:
                ofp.write(b)

    def save_all(self, prefix=None):
        """Save all items in the tempstore"""
        for item in self.objs:
            if prefix is None:
                self.save(item)
            else:
                self.save(item, prefix+item)

    def create(self, name):
        """Creates a new tempfile and returns the path
        :param name: Name of the file to create
        :return: str
        """
        if name in self.objs:
            raise ValueError('{} already exists'.format(name, self))

        temp = tempfile.NamedTemporaryFile(suffix=name)
        self.objs[name] = temp

        if self.return_objects:
            return temp
        else:
            return temp.name

    def get(self, item, fallback=None):
        return self.objs.get(item, fallback)


TEMPFILES = TempStore()
RESULTSFILES = TempStore()
CNA_TEMPLATE = '${REF}/binaries/capture_specific_jetstream_file_creatio'\
               'n/Copy_Number_100bp_Interval_Template.bed'
KNOWN_SV = '${REF}/homo_sapiens/grch37_hg19/hs37d5/tool_specific_re' \
           'sources/tconut/nonmatched_tumor_normal/NA12878_specific/Merged' \
           '_SV_DGV_1kg.bed'
GATK_PATH = 'gatk'
BEDTOOLS_PATH = 'bedtools'


def arg_parser():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description='Generate a reference pack for a new exome kit. See epilog '
                    'for detailed usage instructions.',
        epilog=__doc__
    )

    parser.add_argument(
        '-t', '--targets',
        required=True,
        help='Path to the targets bed file.'
    )

    parser.add_argument(
        '-b', '--baits',
        help='Path to a baits bed file. This file is optional.'
    )

    parser.add_argument(
        '-g', '--gtf',
        required=True,
        help='Path to a reference GTF gene model used for creating the '
             'extended bed file'
    )

    parser.add_argument(
        '-r', '--ref-prefix',
        required=True,
        help='Prefix to the reference files, it is expected that a .fai '
             'and .dict will exist at this prefix'
    )

    parser.add_argument(
        '-o', '--out-prefix',
        required=True,
        help='Out file prefix, several files will be created during this process'
    )


    parser.add_argument(
        '-c', '--cna-template',
        default=CNA_TEMPLATE,
        help='Path to the CNA template bed file. [%(default)s]'
        # TODO: This should be replaced in the future with code that
        # that generates the template from a given refdict
    )

    parser.add_argument(
        '-k', '--known-sv',
        default=KNOWN_SV,
        help='Path to the known SVs bed file. [%(default)s]'
    )

    parser.add_argument(
        '--log-file',
        default='log.txt',
        help='Path to store command logs [%(default)s]'
    )

    parser.add_argument(
        '--save-temp', 
        action='store_true', 
        default=False,
        help='Also saves all temp files created during the process'
    )

    parser.add_argument(
        '--cna',
        action='store_true',
        help='Enable cna ref file creation'
    )

    parser.add_argument(
        '--style',
        action='store_true',
        help='Style of reference - chr vs no chr'
    )

    parser.add_argument(
        '--parent_dir',
        action='store_true',
        help='Location of where capture kits are stored'
    )

    parser.add_argument(
        '--exome_code',
        action='store_true',
        help='Code for the exome (helps with naming)'
    )

    parser.add_argument(
        '--exome_path',
        action='store_true',
        help='Path to exome capture kit subfolder'
    )


    return parser


def run_ext_process(cmd_args, check=True, **kwargs):
    """Run an external process and log when starting"""
    if isinstance(cmd_args, str):
        cmd = cmd_args
    else:
        cmd = ' '.join(cmd_args)

    log.info(f'Starting process: {cmd}')
    subprocess.run(cmd_args, check=check, **kwargs)


def picard_bedtointervallist(bed, refdict, out_path, no_header_out_path):
    """Starts a Picard BedToIntervalList process that writes to out_path"""
    cmd = f'{GATK_PATH} BedToIntervalList -SD "{refdict}" --INPUT "{bed}" --OUTPUT "{out_path}"'
    run_ext_process(cmd, shell=True)

    # Some GATK-Picard modules do not currently report an accurate exit 
    # status this is indended to catch errors with this GATK command
    if not os.path.exists(out_path):
        msg = f'Expected output ({out_path}) not found for command: {cmd}'
        raise ChildProcessError(msg)

    cmd = f'grep -v "@" "{out_path}" > {no_header_out_path}'
    run_ext_process(cmd, shell=True)

    if not os.path.exists(no_header_out_path):
        msg = f'Expected output ({no_header_out_path}) not found for command: {cmd}'
        raise ChildProcessError(msg)
	

def bedtools_intersect(a, b, out_path):
    """Writes bedfile to out_path that includes all intervals in "a" which
    intersect any interval in "b". """
    cmd = f'{BEDTOOLS_PATH} intersect -wa -a "{a}" -b "{b}" > "{out_path}"'
    run_ext_process(cmd, shell=True)


def bedtools_intersect_v(a, b, out_path):
    """Writes bedfile to out_path that includes all intervals in "a" which
    intersect any interval in "b". Uses -v option: Only report those entries 
    in A that have no overlap in B."""
    cmd = f'{BEDTOOLS_PATH} intersect -v -a "{a}" -b "{b}" > "{out_path}"'
    run_ext_process(cmd, shell=True)


def bedtools_slop(bed, genome, out_path, b=100):
    """bed can be path to a BED/GFF/VCF"""
    cmd = f'{BEDTOOLS_PATH} slop -b "{b}" -g "{genome}" -i "{bed}" > "{out_path}"'
    run_ext_process(cmd, shell=True)


def bedtools_sort_faidx(bed, faidx, out_path):
    cmd = f'{BEDTOOLS_PATH} sort -faidx "{faidx}" -i "{bed}" > "{out_path}"'
    run_ext_process(cmd, shell=True)


def bedtools_smerge(bed, faidx, out_path):
    cmd = f'set -o pipefail && ' \
          f'{BEDTOOLS_PATH} sort -faidx "{faidx}" -i "{bed}" | ' \
          f'"{BEDTOOLS_PATH}" merge -i stdin > "{out_path}"'
    run_ext_process(cmd, shell=True)


def union(a, b, out_path):
    with open(out_path, 'w') as fp:
        for interval in intervals.read_bed(a):
            print(intervals.bed.write(interval), file=fp)

        for interval in intervals.read_bed(b):
            print(intervals.bed.write(interval), file=fp)


def gtf_to_bed(gtf_path, out_path):
    """Using the libraries from formats/intervals/ this function
    converts a GTF to a bed file by selecting all the 'exon' features"""
    log.info(f'Converting gtf to bed: {gtf_path} > {out_path}')
    ints = intervals.read_gffv2(gtf_path)
    ints = ints.filter(lambda i: i['feature'] == 'exon')  # Select only exons

    with open(out_path, 'w') as fp:
        fp.write(intervals.to_bed(ints))


def make_extended_bed(targets, refdict, reffai, gtf, out_path):
    log.info('Making extended .bed')
    exons_bed = TEMPFILES.create('.exons.bed')
    refdict_genome = TEMPFILES.create('.refdict.genome')
    padded100_targets = TEMPFILES.create('.padded100_targets.bed')
    captured_exons = TEMPFILES.create('.captured_exons.bed')
    all_intervals = TEMPFILES.create('.all_intervals.bed')

    # A: Make an exon bed file from the GTF
    gtf_to_bed(gtf, out_path=exons_bed)

    # B: Make padded bed for generating filter files
    # First create a "bedtools genome" file from the refdict
    refdict_to_bedtools_genome(refdict, out_path=refdict_genome)

    # Now pad the targets file
    bedtools_slop(
        bed=targets,
        genome=refdict_genome,
        b=100,
        out_path=padded100_targets
    )

    # C: Intersect A and B, keep A
    # This keeps only those exons which intersect a target in the
    # padded targets file
    bedtools_intersect(
        a=exons_bed,
        b=padded100_targets,
        out_path=captured_exons
    )

    # D: Generate Union of A and C
    union(
        a=padded100_targets,
        b=captured_exons,
        out_path=all_intervals
    )

    # Sort and merge all the intervals in D
    bedtools_smerge(
        bed=all_intervals,
        faidx=reffai,
        out_path=out_path
    )


def make_cna_index(bed_path, cna_template_bed, out_path, x='23', y='24'):
    """Generates a CNA index from template and filter bed files."""
    bed = intervals.read_bed(bed_path)

    # First we need to chang X Y to integers in order to match template
    for i in bed:
        if i['seqname'] == 'X':
            i['seqname'] = x
        elif i['seqname'] == 'Y':
            i['seqname'] = y

    tmp_bed_path = TEMPFILES.create('.temp_23_24.bed')
    with open(tmp_bed_path, 'w') as fp:
        print(intervals.to_bed(bed), file=fp)

    cmd_args = [
        BEDTOOLS_PATH,
        'intersect',
        '-c',
        '-a', cna_template_bed,
        '-b', tmp_bed_path
    ]

    with open(out_path, 'w') as fp:
        subprocess.run(cmd_args, stdout=fp, check=True)


def make_unmatched_cna_index(intervals, known_sv, cna_template_bed, out_path):
    ints_not_in_sv = TEMPFILES.create('.ints_not_in_sv.bed')

    bedtools_intersect_v(
        a=intervals,
        b=known_sv,
        out_path=ints_not_in_sv
    )

    make_cna_index(
        bed_path=ints_not_in_sv,
        cna_template_bed=cna_template_bed,
        out_path=out_path
    )


def main(args=None):
    parser = arg_parser()
    args = parser.parse_args(args)
    logging.basicConfig(
        level=logging.INFO,
        format="[%(asctime)s] %(message)s",
        handlers=[logging.FileHandler("log.txt"), logging.StreamHandler()]
    )

    msg = 'Command options:\n'
    for k, v in vars(args).items():
        msg += f'{k}: {v}'
    log.critical(msg)

    refdict = args.ref_prefix + '.dict'
    reffai = args.ref_prefix + '.fa.fai'

    if not os.path.exists(refdict):
        raise FileNotFoundError(refdict)

    if not os.path.exists(reffai):
        raise FileNotFoundError(reffai)

    try:
        targets_intervallist = RESULTSFILES.create('.targets.interval_list')
        no_header_targets_intervallist = RESULTSFILES.create('.no.header.targets.interval_list')
        extended_bed = RESULTSFILES.create('.extended.bed')

        picard_bedtointervallist(
            bed=args.targets,
            refdict=refdict,
            out_path=targets_intervallist,
            no_header_out_path=no_header_targets_intervallist
        )

        if args.baits:
            baits_intervallist = RESULTSFILES.create('.baits.interval_list')
            no_header_baits_intervallist = RESULTSFILES.create('.no.header.baits.interval_list')

            picard_bedtointervallist(
                bed=args.baits,
                refdict=refdict,
                out_path=baits_intervallist,
                no_header_out_path=no_header_baits_intervallist
            )

        make_extended_bed(
            targets=args.targets,
            refdict=refdict,
            reffai=reffai,
            gtf=args.gtf,
            out_path=extended_bed
        )
        
        if args.cna:
            cna_index = RESULTSFILES.create('.cna_index.bed')
            ucna_index = RESULTSFILES.create('.unmatched_cna_index.bed')
            make_cna_index(
                bed_path=extended_bed,
                cna_template_bed=args.cna_template,
                out_path=cna_index
            )

            make_unmatched_cna_index(
                intervals=extended_bed,
                known_sv=args.known_sv,
                cna_template_bed=args.cna_template,
                out_path=ucna_index
            )

        log.critical("Saving results...")
        RESULTSFILES.save_all(prefix=args.out_prefix)
    finally:
        if args.save_temp:
            log.critical("Saving temp files...")
            TEMPFILES.save_all(prefix='temp')


if __name__ == "__main__":
    main()
