import pandas as pd

from protrend.model.model import Source
from protrend.transform.regulondb.base import RegulondbTransformer


class SourceTransformer(RegulondbTransformer):
    name = 'regulondb'
    type = 'database'
    url = 'http://regulondb.ccg.unam.mx/'
    doi = '10.1093/nar/gky1077'
    authors = ['Alberto Santos-Zavaleta', 'Heladia Salgado', 'Socorro Gama-Castro', 'Mishael Sánchez-Pérez',
               'Laura Gómez-Romero', 'Daniela Ledezma-Tejeida', 'Jair Santiago García-Sotelo',
               'Kevin Alquicira-Hernández', 'Luis José Muñiz-Rascado', 'Pablo Peña-Loredo',
               'Cecilia Ishida-Gutiérrez', 'David A Velázquez-Ramírez', 'Víctor Del Moral-Chávez',
               'César Bonavides-Martínez', 'Carlos-Francisco Méndez-Cruz', 'James Galagan', 'Julio Collado-Vides']
    description = 'RegulonDB v 10.5: tackling challenges to unify classic and high throughput knowledge of gene regulation in E. coli K-12'

    default_node = Source
    default_node_factors = ('name',)
    default_order = 100
    columns = {'protrend_id', 'name', 'type', 'url', 'doi', 'authors', 'description'}

    def transform(self):
        collectf = dict(name=[self.name],
                        type=[self.type],
                        url=[self.url],
                        doi=[self.doi],
                        authors=[self.authors],
                        description=[self.description])

        df = pd.DataFrame(collectf, index=[0])

        self._stack_transformed_nodes(df)

        return df