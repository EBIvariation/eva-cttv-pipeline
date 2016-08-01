import argparse
import gzip
import subprocess
import sys
from types import SimpleNamespace


class ArgParser:
    def __init__(self, argv):
        parser = argparse.ArgumentParser()

        parser.add_argument("-i", dest="infilepath", required=True)
        parser.add_argument("-o", dest="outfilepath", required=True)
        parser.add_argument('-s', dest="structural", action='store_true')

        args = parser.parse_args(args=argv[1:])

        self.infilepath = args.infilepath
        self.outfilepath = args.outfilepath
        self.structural = args.structural


def skip_line(record, structural):
    if record.assembly != "GRCh38":
        return True
    if all(x not in record.clin_sig for x in ["Pathogenic", "Likely Pathogenic"]):
        return True
    if structural:
        if record.nsv == "-":
            return True
    else:
        if not ((record.ref != "na" and record.alt != "na") or record.rs != "-1"):
            return True
    return False


def make_output_lines(record):
    output_lines = []
    for rcv in record.rcvs:
        output_line = '\t'.join([record.chrom, record.start, record.end, record.allele, "+",
                                 record.rs, rcv, record.gene_id]) + "\n"
        output_lines.append(output_line)
    return output_lines


def make_output_lines_structural(record):
    output_lines = []
    for rcv in record.rcvs:
        output_line = '\t'.join([record.chrom, record.start, record.end, record.allele, "+", rcv,
                                 record.nsv, record.gene_id, record.pheno_ids]) + "\n"
        output_lines.append(output_line)
    return output_lines


def get_variant_summary_record(line):
    line = line.rstrip()
    line_list = line.split("\t")
    record = SimpleNamespace()

    record.assembly = record[12]
    record.chrom = line_list[13]
    record.start = line_list[14]
    record.end = line_list[15]
    record.ref = line_list[25]
    record.alt = line_list[26]
    record.nsv = line_list[7]
    record.gene_id = line_list[3]
    record.pheno_ids = line_list[10]

    record.rcvs = line_list[8].split(";")
    record.rs = "rs{}".format(line_list[6])
    record.allele = "{}/{}".format(record.ref, record.alt)

    return record


def process_file(infilepath, outfilepath, structural):
    if structural:
        make_output_lines_func = make_output_lines_structural
    else:
        make_output_lines_func = make_output_lines

    with gzip.open(infilepath, "rt") as infile:
        with gzip.open(outfilepath, "wt") as outfile:
            infile.readline()
            for line in infile:
                record = get_variant_summary_record(line)
                if skip_line(record, structural):
                    continue
                output_lines = make_output_lines_func(line)
                if output_lines:
                    for output_line in output_lines:
                        outfile.write(output_line)


def post_process(outfilepath):
    print("running post_process")
    temppath = "{}.temp".format(outfilepath)
    commands = ["zcat {} | sort | uniq | gzip -c > {}".format(outfilepath, temppath),
                            "mv {} {}".format(temppath, outfilepath)]
    for command in commands:
        subprocess.check_call(command, shell=True)


def main():
    args = ArgParser(sys.argv)
    process_file(args.infilepath, args.outfilepath, args.structural)
    post_process(args.outfilepath)


if __name__ == '__main__':
    main()
