# TODO: implement the promoter transformer

class Promoter(FunctionalTFBSTransformer,
                      source='functional_tfbs',
                      version='0.0.0',
                      node=TFBS,
                      order=100, # confirmar se a order é suposto ser maior que a dos tfbs (90)
                      register=True):

    def fetch_nodes(self):
        # TODO: method to fetch promoter sequences from NCBI or another DB using biopython
        # se não houver maneira, temos de ir buscar ao genoma x pb para trás e para a frente do gene
        # ver se há ferramentas caso o biopython não dê
        pass

    def align_tfbs(self, tfbs: pd.DataFrame, promoters: pd.DataFrame) -> pd.DataFrame:
        # TODO: method to align tfbs against prometers
        pass

    def calculate_descriptors(self, aligned_tfbs: pd.DataFrame) -> pd.DataFrame:
        # TODO: method to calculate descriptors
        pass

    def transform(self) -> pd.DataFrame:
        tfbs = self.fetch_nodes()
        promoters = read_promoters(source=self.source, version=self.version, columns=[])

        aligned_tfbs = self.align_tfbs(tfbs, promoters)

        final_tfbs = self.calculate_descriptors(aligned_tfbs)

        self.stack_transformed_nodes(final_tfbs)
        return final_tfbs