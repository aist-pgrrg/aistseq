import sys
import matplotlib.pyplot as plt

def main(report_file, output_dir, reference_success_count, reference_unsuccess_count, denovo_success_count, denovo_unsuccess_count):
    labels = ['Reference Success', 'Reference Unsuccess', 'Denovo Success', 'Denovo Unsuccess']
    counts = [int(reference_success_count), int(reference_unsuccess_count), int(denovo_success_count), int(denovo_unsuccess_count)]
    colors = ['#ADF802', '#FFCE44', '#1589FF', '#797979']  # Color codes for green, red, black, yellow

    plt.figure(figsize=(10, 6))  # Adjust the figure size as needed
    plt.bar(labels, counts, color=colors)
    plt.xlabel('Run Types')
    plt.ylabel('Count')
    plt.title('Assembly Run Results')
    
    # Rotate x-axis labels for better visibility
    plt.xticks(rotation=45, ha="right")
    
    # Add space at the bottom to prevent label cutoff
    plt.tight_layout()

    plt.savefig(f'{output_dir}/success_bar_graph.png')

if __name__ == "__main__":
    if len(sys.argv) != 7:
        print("Usage: python generate_summary_plot.py <report_file> <output_dir> <reference_success_count> <reference_unsuccess_count> <denovo_success_count> <denovo_unsuccess_count>")
        sys.exit(1)
    
    report_file = sys.argv[1]
    output_dir = sys.argv[2]
    reference_success_count = sys.argv[3]
    reference_unsuccess_count = sys.argv[4]
    denovo_success_count = sys.argv[5]
    denovo_unsuccess_count = sys.argv[6]
    
    main(report_file, output_dir, reference_success_count, reference_unsuccess_count, denovo_success_count, denovo_unsuccess_count)
