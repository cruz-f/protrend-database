import pandas as pd

from protrend.io import read_from_stack, read_csv
from protrend.model.model import Evidence
from protrend.transform.coryneregnet.base import CoryneRegNetTransformer
from protrend.utils import SetList


class EvidenceTransformer(CoryneRegNetTransformer):
    default_node = Evidence
    default_transform_stack = {'bsub': 'bsub_regulation.csv',
                               'cglu': 'cglu_regulation.csv',
                               'cglu_rna': 'cglu_rna.csv',
                               'ecol': 'ecol_regulation.csv',
                               'mtub': 'mtub_regulation.csv'}
    default_order = 100
    columns = SetList(['protrend_id',
                       'name', 'description',
                       'TF_locusTag', 'TF_altLocusTag', 'TF_name', 'TF_role',
                       'TG_locusTag', 'TG_altLocusTag', 'TG_name', 'Operon',
                       'Binding_site', 'Role', 'Is_sigma_factor', 'Evidence', 'PMID', 'Source', 'taxonomy',
                       'locus_tag', 'evidence', 'srna_class', 'start_position', 'end_position',
                       'genome', 'orientation', 'evidence_functional', 'functional_rna', 'sequence'])
    cglu_rna_columns = SetList(['locus_tag', 'evidence', 'srna_class', 'start_position', 'end_position',
                                'genome', 'orientation', 'evidence_functional', 'functional_rna', 'sequence'])

    def _transform_evidence(self, regulation: pd.DataFrame, rna: pd.DataFrame) -> pd.DataFrame:
        regulation = self.drop_duplicates(df=regulation, subset=['Evidence'], perfect_match=True, preserve_nan=True)
        regulation = regulation.dropna(subset=['Evidence'])
        regulation['name'] = regulation['Evidence']
        regulation['description'] = None

        rna = self.drop_duplicates(df=rna, subset=['evidence_functional'], perfect_match=True, preserve_nan=True)
        rna = rna.dropna(subset=['evidence_functional'])
        rna['name'] = rna['evidence_functional']
        rna['description'] = None

        df = pd.concat([regulation, rna], axis=0)

        df = self.drop_duplicates(df=df, subset=['name'], perfect_match=True, preserve_nan=True)
        df = df.dropna(subset=['name'])

        return df

    def transform(self):
        regulation = self._build_regulations()
        rna = read_from_stack(self.transform_stack, file='cglu_rna', default_columns=self.cglu_rna_columns,
                              reader=read_csv, sep='\t')

        evidence = self._transform_evidence(regulation, rna)

        self._stack_transformed_nodes(evidence)
        return evidence
