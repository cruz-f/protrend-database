import pandas as pd
from Bio import SeqIO
from neo4j.exceptions import Neo4jError, DriverError

from protrend.io.utils import read_promoters
from protrend.transform.functional_tfbs.base import FunctionalTFBSTransformer


def read_fasta(file_name):
    fasta_file = SeqIO.parse(open(file_name), 'fasta')
    fasta_seq = [seq for seq in fasta_file][0]
    return fasta_seq

def read_genbank(file_name):
    gb_file = SeqIO.read(file_name, "gb")
    return gb_file

class Promoter(FunctionalTFBSTransformer,
                      source='functional_tfbs',
                      version='0.0.0',
                      node=TFBS,
                      order=100,
                      register=True):

    def fetch_nodes(self):
        try:
            return self.node.node_to_df()
        except (Neo4jError, DriverError) as e:
            print(e)
            return pd.DataFrame()

    def retrieve_promoter_seq(self, in_gb, in_fasta, file_output, file2_output):
        prom_output = ""
        seq_output = ""

        genome = read_fasta(in_fasta)
        reverse_genome = genome[::-1]
        GBrecord = read_genbank(in_gb)

        prom_len = 150

        for feature in GBrecord.features:

            if feature.type == 'gene':
                my_start = feature.location.start.position  # Identifies the start position of the gene on the sense strand (5'-3')
                my_end = feature.location.end.position  # Identifies the end position of the gene on the sense strand (5'-3')
                if (my_start - prom_len) < 0:
                    start_500 = 0
                else:
                    start_500 = my_start - prom_len
                if (my_end + prom_len) > len(genome):
                    end_500 = len(genome)
                else:
                    end_500 = my_end + prom_len

                locus_tag = feature.qualifiers['locus_tag'][0]
                feat_loc = str(feature.location)

                if feature.strand == -1:
                    prom = genome[my_end:end_500].reverse_complement()
                    gene_seq = reverse_genome[my_start:my_end].reverse_complement()

                elif feature.strand == 1:
                    prom = genome[start_500:my_start]
                    gene_seq = genome[my_start:my_end]

                else:
                    continue

                prom_output += f">Promoter ___ {feat_loc} ___ {locus_tag} \n"
                seq_output += f">Gene sequence ___ {feat_loc} ___ {locus_tag} \n"

                if len(prom.seq) == prom_len:
                    prom_output += str(prom.seq) + "\n\n"

                else:
                    diff = prom_len - len(prom.seq)
                    complement = 'A' * diff
                    prom_output += complement + str(prom.seq) + "\n\n"

                seq_output += str(gene_seq.seq) + "\n\n"

        file = open(file_output, 'w')
        file.write(prom_output)
        file2 = open(file2_output, 'w')
        file2.write(seq_output)
        file.close()
        file2.close()

    def transform(self) -> pd.DataFrame:
        tfbs = self.fetch_nodes()
        promoters = read_promoters(source=self.source, version=self.version, columns=[])

        aligned_tfbs = self.align_tfbs(tfbs, promoters)

        final_tfbs = self.calculate_descriptors(aligned_tfbs)

        self.stack_transformed_nodes(final_tfbs)
        return final_tfbs