import json

from protein_visualisation import ProteinVisualisation

class VisualisationRunner:
  def __init__(self, protein_viz_path, protein_viz_output_path, protein_counts_path, protein_counts_output_path):
    self.protein_viz_path = protein_viz_path
    self.protein_viz_output_path = protein_viz_output_path
    self.protein_counts_path = protein_counts_path
    self.protein_counts_output_path = protein_counts_output_path

  def run(self):
    self.rewrite_protein_viz()
    self.rewrite_protein_counts()

  def rewrite_protein_viz(self):
    rewritten = []

    with open(self.protein_viz_path, 'r') as infile:
      protein_visualisation_array = json.load(infile)

    for condition_comparison in protein_visualisation_array:
      protein_visualisation = ProteinVisualisation(condition_comparison['data'])
      protein_visualisation.rewriteFastaHeaders()

      rewritten_condition_comparison = {
        'conditionComparison': condition_comparison['conditionComparison'],
        'up.condition': condition_comparison['up.condition'],
        'down.condition': condition_comparison['down.condition'],
        'fdrLimit': condition_comparison['fdrLimit'],
        'data': protein_visualisation.toJSON()
      }

      rewritten.append(rewritten_condition_comparison)

    with open(self.protein_viz_output_path, 'w') as outfile:
      json.dump(rewritten, outfile)

  def rewrite_protein_counts(self):
    with open(self.protein_counts_path, 'r') as infile:
      protein_counts_array = json.load(infile)

    protein_visualisation = ProteinVisualisation(protein_counts_array)
    protein_visualisation.rewriteFastaHeaders()
    rewritten = protein_visualisation.toJSON()

    with open(self.protein_counts_output_path, 'w') as outfile:
      json.dump(rewritten, outfile)
