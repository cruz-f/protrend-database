from itemloaders import ItemLoader
from scrapy import Request
from scrapy.http import FormRequest, Response
from scrapy.spiders import CSVFeedSpider
from scrapy.utils.iterators import csviter
from scrapy.utils.spider import iterate_spider_output

from protrend.extract.items.collectf import TaxonomyItem, OrganismItem, RegulonItem, OperonItem, GeneItem, TFBSItem, \
    ExperimentalEvidenceItem, TranscriptionFactorItem, CollecTFItem


class CollecTFSpider(CSVFeedSpider):
    name = "collectf"

    regulon_url = 'http://www.collectf.org/uniprot/'
    export_url = 'http://www.collectf.org/browse/export/'
    motif_reports_url = 'http://www.collectf.org/browse/view_motif_reports_by_taxonomy'
    motif_reports_tf_url = 'http://www.collectf.org/browse/view_motif_reports_by_TF_family'

    custom_settings = {
        'ITEM_PIPELINES': {'extract.pipelines.collectf.CollecTFPipeline': 800}
    }

    start_urls = ("http://www.collectf.org/browse/browse_by_taxonomy/",
                  "http://www.collectf.org/browse/list_all_TFs/")

    allowed_domains = ["collectf.org"]

    delimiter = '\t'
    headers = ['TF', 'TF_accession', 'genome_accession', 'organism',
               'site_start', 'site_end', 'site_strand', 'sequence',
               'mode', 'experimental_evidence', 'regulated genes (locus_tags)']

    def start_requests(self):

        for url in self.start_urls:

            if url.endswith("browse_by_taxonomy/"):
                yield Request(url, self.parse_browse_by_taxonomy)

            elif url.endswith("list_all_TFs/"):
                yield Request(url, self.parse_all_tfs)

    def parse_browse_by_taxonomy(self, response: Response):

        taxa_xpath = "/html/body/div/div/div[1]/div/ul/li/div/a"
        taxa = response.xpath(taxa_xpath)

        for taxon in taxa:

            taxonomy_loader = ItemLoader(item=TaxonomyItem(), selector=taxon)

            taxonomy_loader.add_xpath("taxonomy_id", ".//@onclick")
            taxonomy_loader.add_xpath("name", ".//text()")
            taxonomy_item = taxonomy_loader.load_item()

            tax_id = taxonomy_item.get("taxonomy_id")
            if tax_id is not None:
                url = f'{self.motif_reports_url}/{tax_id}/'

                taxonomy_loader = ItemLoader(item=taxonomy_item)
                taxonomy_loader.add_value("url", url)

                taxonomy_item = taxonomy_loader.load_item()
                cb_kwargs = dict(taxonomy_item=taxonomy_item)

                yield response.follow(url=url, callback=self.parse_motif_report, cb_kwargs=cb_kwargs)

    def parse_motif_report(self, response: Response, taxonomy_item: TaxonomyItem):

        form_data = [('tsv', 'Download spreadsheet (TSV)')]

        form_inputs_xpath = "//*[@id='export']/div/form/input"
        form_inputs = response.xpath(form_inputs_xpath)

        for form_input in form_inputs:

            name = form_input.attrib.get('name')
            value = form_input.attrib.get('value')

            if name and value:
                form_data.append((name, value))

        cookies = {}
        cb_kwargs = {}

        for key, val in form_data:
            if 'csrfmiddlewaretoken' == key:
                cookies['csrftoken'] = val
                cb_kwargs['taxonomy_item'] = taxonomy_item
                break

        if cookies and cb_kwargs:

            # self.parse calls self.parse_row
            yield FormRequest(self.export_url,
                              formdata=form_data,
                              cookies=cookies,
                              cb_kwargs=cb_kwargs)

    def parse_row(self, response: Response, row: dict, taxonomy_item: TaxonomyItem = None):

        if row['TF'] == 'TF':
            return

        if taxonomy_item is None:
            return

        # organism
        organism_item = self.parse_organism(row=row)

        organism_loader = ItemLoader(item=organism_item)
        organism_loader.add_value("taxonomy", taxonomy_item["taxonomy_id"])
        organism_item = organism_loader.load_item()

        taxonomy_loader = ItemLoader(item=taxonomy_item)
        taxonomy_loader.add_value("organism", organism_item["name"])
        taxonomy_item = taxonomy_loader.load_item()

        # regulon
        regulon_item = self.parse_regulon(row=row)

        regulon_loader = ItemLoader(item=regulon_item)
        regulon_loader.add_value('organism', organism_item['name'])
        regulon_item = regulon_loader.load_item()

        organism_loader = ItemLoader(item=organism_item)
        organism_loader.add_value('regulon', regulon_item['uniprot_accession'])
        organism_item = organism_loader.load_item()

        # operon
        operon_item = self.parse_operon(row=row)

        # if there is an operon there is also genes
        if operon_item.get('operon_id'):

            operon_loader = ItemLoader(item=operon_item)
            operon_loader.add_value('regulon', regulon_item['uniprot_accession'])
            operon_item = operon_loader.load_item()

            regulon_loader = ItemLoader(item=regulon_item)
            regulon_loader.add_value('operon', operon_item['operon_id'])
            regulon_item = regulon_loader.load_item()

        # genes
        genes_items = self.parse_genes(row=row)

        operon_loader = ItemLoader(item=operon_item)
        regulon_loader = ItemLoader(item=regulon_item)

        # if there is an operon there is also genes
        for gene_item in genes_items:

            if gene_item.get('locus_tag'):

                operon_loader.add_value('gene', gene_item['locus_tag'])
                regulon_loader.add_value('gene', gene_item['locus_tag'])

                gene_loader = ItemLoader(item=gene_item)
                gene_loader.add_value('regulon', regulon_item['uniprot_accession'])
                gene_loader.add_value('operon', operon_item['operon_id'])
                gene_loader.load_item()

        operon_loader.load_item()
        regulon_loader.load_item()

        # tfbs
        tfbs_item = self.parse_tfbs(row=row, organism_name=organism_item['name'])

        organism_loader = ItemLoader(item=organism_item)
        organism_loader.add_value('tfbs', tfbs_item['tfbs_id'])
        organism_loader.load_item()

        regulon_loader = ItemLoader(item=regulon_item)
        regulon_loader.add_value('tfbs', tfbs_item['tfbs_id'])
        regulon_loader.load_item()

        tfbs_loader = ItemLoader(item=tfbs_item)
        tfbs_loader.add_value('organism', organism_item['name'])
        tfbs_loader.add_value('regulon', regulon_item['uniprot_accession'])

        if operon_item.get('operon_id'):

            operon_loader = ItemLoader(item=operon_item)
            operon_loader.add_value('tfbs', tfbs_item['tfbs_id'])
            operon_loader.load_item()

            tfbs_loader.add_value('operon', operon_item['operon_id'])

            for gene_item in genes_items:

                if gene_item.get('locus_tag'):

                    gene_loader = ItemLoader(item=gene_item)
                    gene_loader.add_value('tfbs', tfbs_item['tfbs_id'])
                    gene_loader.load_item()

                    tfbs_loader.add_value('gene', gene_item['locus_tag'])

        tfbs_loader.load_item()

        # experimental evidence
        exp_items = self.parse_experimental_evidences(row=row)

        regulon_loader = ItemLoader(item=regulon_item)
        tfbs_loader = ItemLoader(item=tfbs_item)
        for exp_item in exp_items:

            if exp_item.get('exp_id'):

                regulon_loader.add_value('experimental_evidence', exp_item['exp_id'])
                tfbs_loader.add_value('experimental_evidence', exp_item['exp_id'])

                exp_loader = ItemLoader(item=exp_item)
                exp_loader.add_value('regulon', regulon_item['uniprot_accession'])
                exp_loader.add_value('tfbs', tfbs_item['tfbs_id'])
                exp_loader.load_item()

        regulon_loader.load_item()
        tfbs_loader.load_item()

        collectf_loader = ItemLoader(item=CollecTFItem())

        collectf_loader.add_value('taxonomy_item', taxonomy_item)
        collectf_loader.add_value('organism_item', organism_item)
        collectf_loader.add_value('regulon_item', regulon_item)
        collectf_loader.add_value('operon_item', operon_item)
        collectf_loader.add_value('genes_items', genes_items)
        collectf_loader.add_value('tfbs_item', tfbs_item)
        collectf_loader.add_value('exp_items', exp_items)

        return collectf_loader.load_item()

    @staticmethod
    def parse_organism(row: dict):

        organism_loader = ItemLoader(item=OrganismItem())

        organism_name = row.get("organism")
        organism_loader.add_value("name", organism_name)

        genome_accession = row.get("genome_accession")
        organism_loader.add_value("genome_accession", genome_accession)

        return organism_loader.load_item()

    def parse_regulon(self, row: dict):

        regulon_loader = ItemLoader(item=RegulonItem())

        uniprot_accession = row.get('TF_accession')
        regulon_loader.add_value('uniprot_accession', uniprot_accession)

        name = row.get('TF')
        regulon_loader.add_value('name', name)

        url = f'{self.regulon_url}{uniprot_accession}'
        regulon_loader.add_value('url', url)

        return regulon_loader.load_item()

    @staticmethod
    def parse_operon(row: dict):

        operon_loader = ItemLoader(item=OperonItem())

        genes = row.get('regulated genes (locus_tags)', '').split(', ')
        operon_loader.add_value('operon_id', genes)

        return operon_loader.load_item()

    @staticmethod
    def parse_genes(row: dict):

        genes = row.get('regulated genes (locus_tags)', '').split(', ')

        genes_items = []

        for gene in genes:
            gene_loader = ItemLoader(item=GeneItem())
            gene = gene.replace(' ', '')
            if gene:
                gene_loader.add_value('locus_tag', gene)
                gene_item = gene_loader.load_item()

                genes_items.append(gene_item)

        return genes_items

    @staticmethod
    def parse_tfbs(row: dict, organism_name: str):

        tfbs_loader = ItemLoader(item=TFBSItem())

        tfbs_loader.add_value('tfbs_id', organism_name)

        site_start = row.get('site_start')
        tfbs_loader.add_value('tfbs_id', site_start)
        tfbs_loader.add_value('site_start', site_start)

        site_end = row.get('site_end')
        tfbs_loader.add_value('tfbs_id', site_end)
        tfbs_loader.add_value('site_end', site_end)

        site_strand = row.get('site_strand')
        tfbs_loader.add_value('tfbs_id', site_strand)
        tfbs_loader.add_value('site_strand', site_strand)

        mode = row.get('mode')
        tfbs_loader.add_value('tfbs_id', mode)
        tfbs_loader.add_value('mode', mode)

        sequence = row.get('sequence')
        tfbs_loader.add_value('sequence', sequence)

        pubmed = row.get('experimental_evidence')
        tfbs_loader.add_value('pubmed', pubmed)

        return tfbs_loader.load_item()

    @staticmethod
    def parse_experimental_evidences(row: dict):

        evidences = row.get('experimental_evidence', '').split(', ')

        evidences_items = []

        for evidence in evidences:
            exp_loader = ItemLoader(item=ExperimentalEvidenceItem())
            exp_loader.add_value('exp_id', evidence)

            exp_item = exp_loader.load_item()
            evidences_items.append(exp_item)

        return evidences_items

    def parse_rows(self, response, **kwargs):
        for row in csviter(response, self.delimiter, self.headers, self.quotechar):
            ret = iterate_spider_output(self.parse_row(response, row, **kwargs))
            for result_item in self.process_results(response, ret):
                yield result_item

    def _parse(self, response, **kwargs):
        response = self.adapt_response(response)
        return self.parse_rows(response, **kwargs)

    @staticmethod
    def parse_all_tfs(response: Response):

        tfs_xpath = "//*[@id='tf-list-table']/tbody/tr"
        tfs = response.xpath(tfs_xpath)

        for tf in tfs:
            tf_loader = ItemLoader(item=TranscriptionFactorItem(), selector=tf)

            name_xpath = './/td[1]/text()'
            tf_loader.add_xpath('name', name_xpath)

            family_xpath = './/td[2]/text()'
            tf_loader.add_xpath('family', family_xpath)

            description_xpath = './/td[3]/text()'
            tf_loader.add_xpath('description', description_xpath)

            pubmed_xpath = './/td[3]/a/text()'
            tf_loader.add_xpath('pubmed', pubmed_xpath)

            regulons_xpath = './/td[4]/a/text()'
            tf_loader.add_xpath('regulon', regulons_xpath)

            tf_item = tf_loader.load_item()

            yield tf_item
