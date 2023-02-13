import argparse
import gzip
import os

from pysam import VariantFile
from Bio import SeqIO


class VCF_file:
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


class Contigs_file:
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
    parser.add_argument('-v', '--variations', help='List of VCF files', nargs="+", dest="variation_files",
                        required=True)
    parser.add_argument('-c', '--contigs', help='List of contigs FASTA files', nargs="+", dest="contigs_files",
                        required=True)
    parser.add_argument('-o', '--out-dir', help='Output directory.', type=str, dest="output", required=True)
    args = parser.parse_args()

    vcf_test = VCF_file(args.variation_files[0])
    vcf_test.print()

    contigs_test = Contigs_file(args.contigs_files[0])
    contigs_test.print()


if __name__ == '__main__':
    main()
