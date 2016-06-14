import argparse
import gzip
import subprocess
import sys


COUNTER_BAD_RS_OR_AL = 0


class ArgParser:
    def __init__(self, argv):
        parser = argparse.ArgumentParser()

        parser.add_argument("-i", dest="infilepath", required=True)
        parser.add_argument("-o", dest="outfilepath", required=True)

        args = parser.parse_args(args=argv[1:])

        self.infilepath = args.infilepath
        self.outfilepath = args.outfilepath


def skip_line(record):
    global COUNTER_BAD_RS_OR_AL
    assembly = record[12]
    if assembly != "GRCh38":
        return True
    clin_sig = record[5]
    if all(x not in clin_sig for x in ["Pathogenic", "Likely Pathogenic"]):
        return True
    rsid = record[6]
    ref = record[25]
    alt = record[26]
    if not ((ref != "na" and alt != "na") or rsid != "-1"):
        COUNTER_BAD_RS_OR_AL += 1
        return True
    return False


def make_output_line(record):
    chrom = record[13]
    start = record[14]
    end = record[15]
    ref = record[25]
    alt = record[26]
    if record[6] != "-1":
        var_id = "rs{}".format(record[6])
    else:
        var_id = "{}_{}_{}_{}_{}".format(chrom, start, end, ref, alt)
    allele = "{}/{}".format(ref, alt)
    output_line = '\t'.join([chrom, start, end, allele, "+", var_id, "\n"])
    return output_line


def process_line(line, outfile):
    line = line.rstrip()
    record = line.split("\t")
    if skip_line(record):
        return
    output_line = make_output_line(record)
    outfile.write(output_line)


def process_file(infilepath, outfilepath):
    with gzip.open(infilepath, "rt") as infile:
        with gzip.open(outfilepath, "wt") as outfile:
            infile.readline()
            for line in infile:
                process_line(line, outfile)


def post_process(outfilepath):
    print("running post_process")
    temppath = "{}.temp".format(outfilepath)
    commands = ["zcat {} | sort | uniq | gzip -c > {}".format(outfilepath, temppath),
                            "mv {} {}".format(temppath, outfilepath)]
    for command in commands:
        subprocess.check_call(command, shell=True)


def main():
    global COUNTER_BAD_RS_OR_AL

    args = ArgParser(sys.argv)

    process_file(args.infilepath, args.outfilepath)

    post_process(args.outfilepath)

    print(COUNTER_BAD_RS_OR_AL)


if __name__ == '__main__':
    main()
