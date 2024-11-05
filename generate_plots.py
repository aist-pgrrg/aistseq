# generate_plots.py
import plotly.express as px
import pandas as pd

def generate_coverage_plot(input_coverage_file, output_html_file):
    # Load coverage data
    coverage_data = pd.read_csv(input_coverage_file, sep='\t', names=["chrom", "pos", "coverage"])
    coverage_data["pos"] = coverage_data["pos"].astype(int)

    # Create coverage plot
    coverage_fig = px.line(coverage_data, x="pos", y="coverage", title="Coverage Plot")

    # Save coverage plot as HTML
    coverage_fig.write_html(output_html_file)

def generate_mutations_plot(input_vcf_file, output_html_file):
    # Load mutations data
    mutations_data = pd.read_csv(input_vcf_file, sep='\t', comment="#", header=None,
                                 names=["chrom", "pos", "id", "ref", "alt", "qual", "filter", "info", "format", "sample"])

    # Extract information about each mutation
    mutations_data["mutation"] = mutations_data.apply(
        lambda row: f"{row['ref']} to {row['alt']}", axis=1)

    # Create mutations plot
    mutations_fig = px.scatter(mutations_data, x="pos", y="qual", color="mutation",
                               title="Mutations Plot", hover_data=["chrom", "mutation", "qual"])

    # Customize hover labels
    #mutations_fig.update_traces(hovertemplate="Chromosome: %{customdata[0]}<br>Position: %{x}<br>Mutation: %{customdata[1]}<br>Quality: %{customdata[2]}")
    mutations_fig.update_traces(hovertemplate="Chromosome: %{customdata[0]}<br>Position: %{x}<br>Mutation: %{customdata[1]}<br>Quality: %{y}")


    # Save mutations plot as HTML
    mutations_fig.write_html(output_html_file)

if __name__ == "__main__":
    import sys

    if len(sys.argv) != 3:
        print("Usage: python generate_plots.py input_file output_html_file")
        sys.exit(1)

    input_file = sys.argv[1]
    output_html_file = sys.argv[2]

    if "coverage" in input_file:
        generate_coverage_plot(input_file, output_html_file)
    elif "mutations" in input_file:
        generate_mutations_plot(input_file, output_html_file)
    else:
        print("Unknown input file type")
