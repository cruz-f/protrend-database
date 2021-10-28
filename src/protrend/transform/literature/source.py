import pandas as pd

from protrend.model.model import Source
from protrend.transform.literature.base import LiteratureTransformer
from protrend.utils import SetList


class SourceTransformer(LiteratureTransformer):
    name = ('bsub_faria_et_al_2017', 'ecol_fang_et_al_2017', 'mtub_turkarslan_et_al_2015', 'paer_vasquez_et_al_2011')
    type = ('literature', 'literature', 'literature', 'literature')
    url = ('https://www.frontiersin.org/articles/10.3389/fmicb.2016.00275',
           'https://www.pnas.org/content/114/38/10286',
           'https://www.nature.com/articles/sdata201510',
           'https://microbialinformaticsj.biomedcentral.com/articles/10.1186/2042-5783-1-3')
    doi = ('10.3389/fmicb.2016.00275',
           '10.1073/pnas.1702581114',
           '10.1038/sdata.2015.10',
           '10.1186/2042-5783-1-3')
    authors = (['José P. Faria', 'Ross Overbeek', 'Ronald C. Taylor', 'Neal Conrad', 'Veronika Vonstein',
                'Anne Goelzer', 'Vincent Fromion', 'Miguel Rocha', 'Isabel Rocha', 'Christopher S. Henry'],
               ['Xin Fang', 'Anand Sastry', 'Nathan Mih', 'Donghyuk Kim', 'Justin Tan', 'James T. Yurkovich',
                'Colton J. Lloyd', 'Ye Gao', 'Laurence Yang', 'Bernhard O. Palsson'],
               ['Serdar Turkarslan', 'Eliza J R Peterson', 'Tige R Rustad', 'Kyle J Minch', 'David J Reiss',
               'Robert Morrison', 'Shuyi Ma', 'Nathan D Price', 'David R Sherman', 'Nitin S Baliga'],
               ['Edgardo Galán-Vásquez', 'Beatriz Luna', 'Agustino Martínez-Antonio'])
    description = ('Reconstruction of the Regulatory Network for Bacillus subtilis and Reconciliation with Gene Expression Data',
                   'Global transcriptional regulatory network for Escherichia coli robustly connects gene expression to transcription factor activities',
                   'A comprehensive map of genome-wide gene regulation in Mycobacterium tuberculosis',
                   'The Regulatory Network of Pseudomonas aeruginosa')

    default_node = Source
    default_order = 100
    columns = SetList(['protrend_id', 'name', 'type', 'url', 'doi', 'authors', 'description'])

    def transform(self):
        db = dict(name=self.name,
                  type=self.type,
                  url=self.url,
                  doi=self.doi,
                  authors=self.authors,
                  description=self.description)

        df = pd.DataFrame(db)

        self._stack_transformed_nodes(df)

        return df