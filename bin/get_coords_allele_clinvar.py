import argparse
import gzip
import subprocess
import sys


COUNTER = 0


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
    global COUNTER
    assembly = record[12]
    if assembly != "GRCh38":
        return True
    clin_sig = record[5]
    if all(x not in clin_sig for x in ["Pathogenic", "Likely Pathogenic"]):
        return True
    if structural:
        nsv = record[7]
        if nsv == "-":
            return True
        COUNTER += 1
    else:
        rsid = record[6]
        ref = record[25]
        alt = record[26]
        if not ((ref != "na" and alt != "na") or rsid != "-1"):
            return True
    return False


def make_output_lines(record):
    chrom = record[13]
    start = record[14]
    end = record[15]
    ref = record[25]
    alt = record[26]
    if record[6] != "-1":
        var_id = "rs{}".format(record[6])
    else:
        var_id = "{}_{}_{}_{}_{}".format(chrom, start, end, ref, alt)
    clinvar_id = record[8]
    clinvar_gene = record[3]
    allele = "{}/{}".format(ref, alt)
    output_lines = []
    output_line = '\t'.join([chrom, start, end, allele, "+", var_id, clinvar_id, clinvar_gene]) + "\n"
    output_lines.append(output_line)
    return output_lines


def make_output_lines_struct(record):
    chrom = record[13]
    start = record[14]
    end = record[15]
    ref = record[25]
    alt = record[26]
    rcvs = record[8].split(";")
    nsv = record[7]
    gene_id = record[3]
    pheno_ids = record[10]

    allele = "{}/{}".format(ref, alt)

    output_lines = []
    for rcv in rcvs:
        output_line = '\t'.join([chrom, start, end, allele, "+", rcv, nsv, gene_id, pheno_ids]) + "\n"
        output_lines.append(output_line)
    return output_lines


def process_line(line, outfile, structural):
    line = line.rstrip()
    record = line.split("\t")
    if skip_line(record, structural):
        return
    if structural:
        output_lines = make_output_lines_struct(record)
    else:
        output_lines = make_output_lines(record)
    for output_line in output_lines:
        outfile.write(output_line)


def process_file(infilepath, outfilepath, structural):
    with gzip.open(infilepath, "rt") as infile:
        with gzip.open(outfilepath, "wt") as outfile:
            infile.readline()
            for line in infile:
                process_line(line, outfile, structural)


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

    print(COUNTER)


if __name__ == '__main__':
    main()
