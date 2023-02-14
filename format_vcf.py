import argparse
import gzip
import os

from pysam import VariantFile
from Bio import SeqIO


class VariationsFile:
    file_path = ""
    variations_file_id = ""

    def __init__(self, vcf_file_path, id=""):
        self.file_path = vcf_file_path
        if id != "":
            self.variations_file_id = id
        else:
            self.variations_file_id = str(hash(vcf_file_path))

    def print(self):
        print(self.variations_file_id)

    def split(self, output_dir, contigs_list):
        variations_map = dict()
        out_dir_with_id = os.path.join(output_dir, self.variations_file_id)
        os.mkdir(out_dir_with_id)

        bcf_in = VariantFile(self.file_path)
        contigs = [contig for contig in bcf_in.header.contigs]

        non_empty_vcf = dict()
        for contig in contigs:
            non_empty_vcf[contig] = False

        

        with gzip.open(self.file_path, 'rt') as handle:
            for record in SeqIO.parse(handle, "fasta"):
                out_file_name = os.path.join(out_dir_with_id, record.id + ".fasta.gz")
                variations_map[record.id] = os.path.relpath(out_file_name)
                with gzip.open(out_file_name, "wt") as output_handle:
                    SeqIO.write(record, output_handle, "fasta")

        return variations_map


class ContigsFile:
    file_path = ""
    contigs_file_id = ""

    def __init__(self, fasta_file_path, id=""):
        self.file_path = fasta_file_path
        if id != "":
            self.contigs_file_id = id
        else:
            self.contigs_file_id = str(hash(fasta_file_path))

    def print(self):
        print(self.contigs_file_id)

    def split(self, output_dir):
        contigs_map = dict()
        out_dir_with_id = os.path.join(output_dir, self.contigs_file_id)
        os.mkdir(out_dir_with_id)
        with gzip.open(self.file_path, 'rt') as handle:
            for record in SeqIO.parse(handle, "fasta"):
                out_file_name = os.path.join(out_dir_with_id, record.id + ".fasta.gz")
                contigs_map[record.id] = os.path.relpath(out_file_name)
                with gzip.open(out_file_name, "wt") as output_handle:
                    SeqIO.write(record, output_handle, "fasta")

        return contigs_map


# ------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description="", formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-v', '--variations', help='List of VCF files', nargs="+", dest="variations_files",
                        required=True)
    parser.add_argument('-c', '--contigs', help='List of contigs FASTA files', nargs="+", dest="contigs_files",
                        required=True)
    parser.add_argument('-o', '--out-dir', help='Output directory.', type=str, dest="output", required=True)
    args = parser.parse_args()

    # read all contigs files in
    contigs_files = list()
    for contigs_file_path in args.contigs_files:
        contigs_files = ContigsFile(contigs_file_path)

    # split all contigs in their directories
    contigs_maps = list()
    for contigs_file in contigs_files:
        contigs_maps.append(contigs_file.split("."))

    # check for multiple contigs with same id
    contigs_seen = set()
    for contig_map in contigs_maps:
        for contig_id in contig_map.keys():
            if contig_id in contigs_seen:
                print("Error: contig id already seen somewhere else. abort")
                exit(1)
            else:
                contigs_seen.add(contig_id)

    # read in all variations files
    variations_files = list()
    for variations_file_path in args.variations_files:
        variations_files.append(VariationsFile(variations_file_path))

    # split all variations in their directories
    for contig in contigs_seen:



if __name__ == '__main__':
    main()
