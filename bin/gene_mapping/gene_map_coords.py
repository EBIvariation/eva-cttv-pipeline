import argparse
import gzip
import sys


def main():
    parser = ArgParser(sys.argv)

    infile_path = parser.infile_path
    outfile_path = parser.outfile_path

    with open_file(infile_path, "rt") as varsum_file:
        with open_file(outfile_path, "wt") as outfile:
            for line in varsum_file:
                line_list = line.rstrip().split("\t")
                if skip_varsum_line(line_list):
                    continue

                output_lines = get_output_lines(line_list)
                for output_line in output_lines:
                    outfile.write(output_line + "\n")


def get_output_lines(line_list):
    chrom = line_list[13]
    start = line_list[14]
    stop = line_list[15]
    ref = line_list[25]
    alt = line_list[26]
    strand = "+"
    rcvs = line_list[8].split(";")
    rs = "rs" + line_list[6] if line_list[6] != "-1" else "-1"
    nsv = line_list[7] if line_list[7] != "-" else "-1"
    ncbi_geneid = line_list[3]

    output_lines = []
    for rcv in rcvs:
        output_line = build_output_line(chrom, start, stop, ref, alt, strand, rcv, rs, nsv, ncbi_geneid)
        output_lines.append(output_line)

    return output_lines


def build_output_line(chrom, start, stop, ref, alt, strand, rcv, rs, nsv, ncbi_geneid):
    output_line_list = [chrom, start, stop, ref, alt, strand, rs, rcv, ncbi_geneid, nsv]
    output_line = "\t".join(output_line_list)

    return output_line


def skip_varsum_line(line_list):
    if line_list[0].startswith("#"):
        return True

    clin_sig = line_list[5]
    if "pathogenic" not in clin_sig.lower():
        return True

    assembly = line_list[12]
    if assembly.lower() != "grch38":
        return True

    return False


def open_file(file_path, mode):
    if file_path.endswith(".gz"):
        return gzip.open(file_path, mode)
    else:
        return open(file_path, mode)


class ArgParser:
    def __init__(self, argv):
        description = """
                Script for extracting the coordinates, variant IDs, and NCBI gene ID for ClinVar
                records from a ClinVar variant_summary file.Output is a tab separated file with one
                ClinVar record per line, and containing the columns: chromosome, start pos,
                stop pos, reference and alternate alleles, strand, RS ID, RCV ID, NCBI gene ID, NSV ID.
                If a column doesn't have a value then a value of '-1' is used.
                """
        parser = argparse.ArgumentParser(description=description)

        parser.add_argument("-i", dest="infile_path", required=True,
                            help="path to variant summary file from ClinVar")
        parser.add_argument("-o", dest="outfile_path", required=True,
                            help="Tab separated output file, with one ClinVar record data per " +
                                 "line. '-1' represents an empty value.")

        args = parser.parse_args(args=argv[1:])

        self.infile_path = args.infile_path
        self.outfile_path = args.outfile_path


if __name__ == '__main__':
    main()
