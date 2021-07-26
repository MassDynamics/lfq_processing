import argparse

from visualisation_runner import VisualisationRunner

if __name__ == "__main__":
  parser = argparse.ArgumentParser(prog='FASTA Parser')
  parser.add_argument("--protein-viz-path", dest='protein_viz_path')
  parser.add_argument("--protein-viz-output-path", dest='protein_viz_output_path')
  parser.add_argument("--protein-counts-path", dest='protein_counts_path')
  parser.add_argument("--protein-counts-output-path", dest='protein_counts_output_path')

  args = parser.parse_args()

  VisualisationRunner(
    protein_viz_path=args.protein_viz_path,
    protein_viz_output_path=args.protein_viz_output_path,
    protein_counts_path=args.protein_counts_path,
    protein_counts_output_path=args.protein_counts_output_path,
  ).run()
