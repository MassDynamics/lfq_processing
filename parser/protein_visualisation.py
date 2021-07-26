from fasta_header_parser import FastaHeaderParser

class ProteinVisualisation:
  def __init__(self, protein_data):
    self.protein_data = protein_data
    self.rewritten_protein_data = []

  def rewriteFastaHeaders(self):
    for item in self.protein_data:
      parsed = FastaHeaderParser().parse(item['FastaHeaders'])

      new_item = {}
      new_item.update(item)
      new_item.update({
        'ProteinId': parsed['proteinId'],
        'GeneName': parsed['gene'],
        'ProteinDescription': parsed['proteinName'],
      })

      self.rewritten_protein_data.append(new_item)

  def toJSON(self):
    return self.rewritten_protein_data