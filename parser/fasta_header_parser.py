from fasta_headers.uniprot_fasta_header import UniprotFastaHeader
from fasta_headers.ups_fasta_header import UPSFastaHeader
from fasta_headers.con_fasta_header import CONFastaHeader
from fasta_headers.not_found_fasta_header import NotFoundFastaHeader
from fasta_headers.human_proteome_resource_recombinant_fasta_header import HumanProteomeResourceRecombinantFastaHeader

class FastaHeaderParser:
  def parse(self, headerString):
    headers = [
      UniprotFastaHeader(),
      HumanProteomeResourceRecombinantFastaHeader(),
      UPSFastaHeader(),
      CONFastaHeader(),
      NotFoundFastaHeader()
    ]

    for header in headers:
      if header.doesMatch(headerString):
        return header.fromString(headerString)