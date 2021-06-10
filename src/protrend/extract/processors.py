from typing import Union


class RegPreciseProcessors:

    @staticmethod
    def process_href(href: str) -> Union[int, None]:
        # href="collection_tax.jsp?collection_id=1001"

        if "=" in href:
            *_, identifier = href.split("=")
            return int(identifier)

        return

    @staticmethod
    def process_operon_name(gene_name: str) -> Union[str, None]:

        if gene_name:

            operon_name = ''.join(letter for letter in gene_name if letter.islower())

            if operon_name:
                return operon_name
        return

    @staticmethod
    def process_span(span_text: str, prefix: str, to: type) -> Union[str, int, None]:

        # " Position: -213"
        # ' Locus tag: BSU03290'
        # " Position: -213"
        # " Score: 8"
        # " Sequence: CGCAAAAT-(1)-AATTTGCG"

        if prefix.lower() in span_text.lower():
            striped = span_text.strip()
            no_ws = striped.replace(" ", "")
            prefix_no_ws = prefix.replace(" ", "")
            no_prefix = no_ws.replace(prefix_no_ws, "")

            return to(no_prefix)

        return

    @staticmethod
    def process_locus_tag(locus_tag: str) -> Union[str, None]:

        # ' Locus tag: BSU03290'

        return RegPreciseProcessors.process_span(locus_tag, "Locus tag:", to=str)

    @staticmethod
    def process_name(name: str) -> Union[str, None]:

        # 'Name: nasE'

        return RegPreciseProcessors.process_span(name, "Name:", to=str)

    @staticmethod
    def process_function(func: str) -> Union[str, None]:

        # 'Funciton: Nitrite reductase [NAD(P)H] small subunit (EC 1.7.1.4)'

        if "Funciton:" in func:
            return func.replace("Funciton: ", "")

        return

    @staticmethod
    def process_position(position: str) -> Union[int, None]:

        # " Position: -213"

        return RegPreciseProcessors.process_span(position, "Position:", to=int)

    @staticmethod
    def process_score(score: str) -> Union[int, None]:

        # " Score: 8"

        return RegPreciseProcessors.process_span(score, "Score:", to=str)

    @staticmethod
    def process_score_str(score: int) -> Union[str, None]:

        # " Score: 8"

        return str(score)

    @staticmethod
    def process_sequence(sequence: str) -> Union[str, None]:

        # " Sequence: CGCAAAAT-(1)-AATTTGCG"

        return RegPreciseProcessors.process_span(sequence, "Sequence:", to=str)
