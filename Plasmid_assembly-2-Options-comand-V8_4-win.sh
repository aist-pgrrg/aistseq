#!/bin/bash

run_plasmid_assembly_with_reference() {
# Activate the plasmid_analysis Conda environment
source activate plasmid_analysis

# Set the reference genome and bcftools path
local reference="$1"
local bcftools_path="$2"
export PATH="$bcftools_path:$PATH"

# Set the data directory specific to the sample
data_directory="$final_output_subdir"

# Run the pipeline
output_prefix="$data_directory"
echo "Running pipeline for sample: ${plasmid_name}"

# Create a directory for this plasmid
mkdir -p "${output_prefix}"

# Index the reference genome if it hasn't been indexed already
index_file="${reference%.*}-index"
if [[ ! -e "${index_file}.1.bt2" ]]; then
    echo "Building Bowtie2 index for reference: $reference"
    bowtie2-build "$reference" "$index_file"
fi

# Run Fastp for quality control, trimming, and filtering
echo "Running Fastp for quality control, trimming, and filtering"
fastp -w 4 -i "$read1_path" -I "$read2_path" -o "${output_prefix}/trimmed_1P.fastq.gz" -O "${output_prefix}/trimmed_2P.fastq.gz" --adapter_fasta "$adapter_fasta" --cut_front --cut_tail --cut_right --length_required 50 --cut_mean_quality 20 -h "${output_prefix}/fastp.html" -j "${output_prefix}/fastp.json" > "${output_prefix}/fastp_stdout.txt"

# Align reads to the reference genome using Bowtie 2
bowtie2 -x "$index_file" -1 "${output_prefix}/trimmed_1P.fastq.gz" -2 "${output_prefix}/trimmed_2P.fastq.gz" -p 4 -S "${output_prefix}/aligned.sam"

# Convert SAM to BAM and sort the BAM file
samtools view -Sb "${output_prefix}/aligned.sam" | samtools sort -o "${output_prefix}/aligned.sorted.bam"

# Index the sorted BAM file
samtools index "${output_prefix}/aligned.sorted.bam"

# Generate coverage information using samtools depth
samtools depth "${output_prefix}/aligned.sorted.bam" > "${output_prefix}/coverage.txt"

# Obtain low coverage positions using Bedtools
bedtools genomecov -ibam "${output_prefix}/aligned.sorted.bam" -bga | awk '$4 < 10' > "${output_prefix}/low_coverage.bed"

# Function to fetch nucleotide using pyfaidx
    fetch_nucleotide_pyfaidx() {
        local id=$1
        local start=$2
        local end=$3
        python -c "from pyfaidx import Fasta; reference = Fasta('$reference'); print(reference['$id'][int($start)-1:int($end)].seq)"
    }
 
    # Print positions and nucleotides to the console and append to all_positions.txt
    while read id start end dp; do
        # Adjust start if it is 0 to be 1
        if [[ "$start" -eq 0 ]]; then
            start=1
        fi
 
        if [[ "$USE_PYFAIDX" -eq 1 ]]; then
            # Fetch nucleotide using pyfaidx
            nucleotide=$(fetch_nucleotide_pyfaidx "$id" "$start" "$end")
        else
            # Fetch nucleotide information using samtools faidx
            nucleotide=$(samtools faidx "$reference" "${id}:${start}-${end}" | tail -n 1)
        fi
 
        # Print the required information to the console and append to all_positions.txt
        echo -e "${id}\t${start}\t${nucleotide}\t${dp}"
 
    done < "${output_prefix}/low_coverage.bed" >> "${output_prefix}/all_positions.txt"

# Check if there is at least one position with coverage >= 10
if [[ "$(awk '$3 >= 10' "${output_prefix}/coverage.txt" | wc -l)" -ge 1 ]]; then
    # Generate VCF file of mutations using bcftools mpileup and call
    bcftools mpileup -Ou -f "$reference" "${output_prefix}/aligned.sorted.bam" | \
        bcftools call -mv -O v --pval-threshold 0.05 | \
        bcftools filter -e 'DP<10 || AF<0.005 || QUAL<20' -o "${output_prefix}/mutations.vcf"

    # Create a new VCF file for low coverage positions
    awk -v OFS='\t' '$4 < 10 {print $1, $2, ".", $3,"N", "0", ".", "INDEL;DP=" $4 "\tGT\t0/0"}' "${output_prefix}/all_positions.txt" > "${output_prefix}/mutations_low_coverage.vcf"
    # Concatenate the low coverage VCF file with the original mutations VCF file
    cat "${output_prefix}/mutations.vcf" "${output_prefix}/mutations_low_coverage.vcf" > "${output_prefix}/mutations_combined.vcf"
    bgzip "${output_prefix}/mutations_combined.vcf"

    # Sort and index the combined VCF file
    bcftools sort -Oz -o "${output_prefix}/mutations_combined_sorted.vcf.gz" "${output_prefix}/mutations_combined.vcf.gz"
    bcftools index "${output_prefix}/mutations_combined_sorted.vcf.gz"

    # Generate consensus sequence using bcftools consensus
    echo Generating consensus sequence using bcftools consensus
    bcftools consensus -f "$reference" "${output_prefix}/mutations_combined_sorted.vcf.gz" --mark-ins "uc" -m "${output_prefix}/low_coverage.bed" > "${output_prefix}/${plasmid_name}_consensus.fa"

    # Set the status message for consensus generation
    status_message="Pipeline plasmid assembly for ${plasmid_name} was completed successfully."
    
    # Create status message HTML
    echo "<html><body><p>Status: ${status_message}</p></body></html>" > "${output_prefix}/status_message.html"

    # Activate the plannotate environment
    conda activate plannotate

    # Run pLannotate on the consensus sequence
    echo "Running pLannotate on the consensus sequence"
    plannotate batch -i "${output_prefix}/${plasmid_name}_consensus.fa" --html -o "$output_prefix" --file_name "assembly" -s ''

    # Rename and move the output files to the appropriate directories
    echo "Rename and move the output files to the appropriate directories"
    mv "$output_prefix/assembly.gbk" "$output_prefix/${plasmid_name}.gbk"
    mv "$output_prefix/assembly.html" "$output_prefix/${plasmid_name}.html"

    # Generate coverage plot
    python generate_plots.py "${output_prefix}/coverage.txt" "${output_prefix}/coverage_plot.html"

    # Check if the mutations_combined_sorted.vcf.gz file exists
    if [[ -e "${output_prefix}/mutations_combined_sorted.vcf.gz" ]]; then
        # Generate mutations plot
        python generate_plots.py "${output_prefix}/mutations_combined_sorted.vcf.gz" "${output_prefix}/mutations_plot.html"
    fi

    # Combine the status message, coverage plot, mutations plot, and pLannotate HTML files into the report
    cat "${output_prefix}/status_message.html" "${output_prefix}/fastp.html" "${output_prefix}/coverage_plot.html" "${output_prefix}/mutations_plot.html" > "${output_prefix}/temp.html"

    # Add the message before the pLannotate result
    echo "<p>The annotation of engineered plasmids by pLannotate</p>" >> "${output_prefix}/temp.html"

    # Concatenate pLannotate HTML
    cat "${output_prefix}/${plasmid_name}.html" >> "${output_prefix}/temp.html"

    # Combine everything into the final report
    cat "${output_prefix}/temp.html" > "${output_prefix}/report_of_${plasmid_name}.html"

    # Remove the intermediate HTML files
    rm "${output_prefix}/status_message.html" "${output_prefix}/fastp.html" "${output_prefix}/coverage_plot.html" "${output_prefix}/mutations_plot.html" "${output_prefix}/temp.html" "${output_prefix}/${plasmid_name}.html"

else
    # Set the status message indicating that consensus sequence won't be generated
    status_message="Pipeline plasmid assembly for ${plasmid_name} was completed, but consensus sequence was not generated due to low coverage."

    # Still generate coverage plot and show "Mutations plot was not generated" message
    python generate_plots.py "${output_prefix}/coverage.txt" "${output_prefix}/coverage_plot.html"
    
    # Check if the mutations.vcf.gz file exists
    if [[ -e "${output_prefix}/mutations.vcf.gz" ]]; then
        mutations_message="<h2>Mutations Plot</h2><iframe src='${output_prefix}/mutations_plot.html'></iframe>"
    else
        mutations_message="<p>Mutations plot was not generated due to low coverage.</p>"
    fi

    # Create status message HTML
    echo "<html><body><p>Status: ${status_message}</p></body></html>" > "${output_prefix}/status_message.html"
    echo "${mutations_message}" > "${output_prefix}/mutations_plot.html"

    # Combine Fastp results with the coverage and mutations plots
    combine_reports="${output_prefix}/combine_reports.html"
    cat "${output_prefix}/fastp.html" "${output_prefix}/coverage_plot.html" "${output_prefix}/mutations_plot.html" > "$combine_reports"

    # Remove the intermediate HTML files
    rm "${output_prefix}/fastp.html" "${output_prefix}/coverage_plot.html" "${output_prefix}/mutations_plot.html"

    # Combine status message, combined report, and mutations plot (if any)
    cat "${output_prefix}/status_message.html" "$combine_reports" > "${output_prefix}/report_of_${plasmid_name}.html"

    # Remove intermediate files
    rm "$combine_reports" "${output_prefix}/status_message.html"
fi

echo "Sample processing complete: ${plasmid_name}"
    # Deactivate the Conda environment
    conda deactivate
    conda activate base

}

#Function to run plasmid assembly without reference genome
function run_plasmid_assembly_without_reference {
    # Activate the plasmid_analysis environment
    source activate plasmid_analysis
    echo "Start running plasmid assembly pipeline for $plasmid_name"
    #local spades_path_PY="$1"

    # Define paths for outputs
    plannotate_output="$final_output_subdir/plannotate_output"
    fastp_output="$final_output_subdir/fastp_output"
    trim_output="$final_output_subdir/trim_output"
    qc_output="$final_output_subdir/qc_output"

    # Create output directories
    mkdir -p "$plannotate_output" "$fastp_output" "$trim_output" "$qc_output"
    # Run Fastp for quality control, trimming, and filtering
    echo "Running Fastp for quality control, trimming, and filtering"
    fastp -w 8 -i "$output_dir/$dir_name/$read1" -I "$output_dir/$dir_name/$read2" -o "$fastp_output/trimmed_1P.fastq.gz" -O "$fastp_output/trimmed_2P.fastq.gz" --adapter_fasta "$adapter_fasta" --cut_front --cut_tail --cut_right --cut_mean_quality 20 --length_required 50 -h "$fastp_output/fastp.html" -j "$fastp_output/fastp.json" > "$fastp_output/stdout.txt"

    # Run Unicycler in conservative mode
    unicycler_output="$final_output_subdir/unicycler_output"
    mkdir -p "$unicycler_output"  # Create the directory only when needed
    echo "Running Unicycler in conservative mode on the trimmed input files"
    #unicycler -1 "$fastp_output/trimmed_1P.fastq.gz" -2 "$fastp_output/trimmed_2P.fastq.gz" \
    #-o "$unicycler_output" --keep 0 --mode conservative --no_rotate
    unicycler -1 "$fastp_output/trimmed_1P.fastq.gz" -2 "$fastp_output/trimmed_2P.fastq.gz" \
    -o "$unicycler_output" --keep 0 --mode conservative --no_rotate

    

    # Check if the output file contains a single FASTA sequence
    num_sequences=$(grep -c '^>' "$unicycler_output/assembly.fasta")

    if [ "$num_sequences" -eq 1 ]; then
        unicycler_final="$unicycler_output"
    else
        # If not one sequence, re-run Unicycler in normal mode
        unicycler_normal_option="$final_output_subdir/unicycler_normal-option"
        mkdir -p "$unicycler_normal_option"  # Create the directory only when needed
        echo "re-running unicycler with normal-option"
        unicycler -1 "$fastp_output/trimmed_1P.fastq.gz" -2 "$fastp_output/trimmed_2P.fastq.gz" \
        -o "$unicycler_normal_option" --keep 0 --mode normal --no_rotate
        
        # Check if the normal mode output contains a single FASTA sequence
        num_sequences=$(grep -c '^>' "$unicycler_normal_option/assembly.fasta")
        
        if [ "$num_sequences" -eq 1 ]; then
            unicycler_final="$unicycler_normal_option"
        else
            # If still not one sequence, re-run Unicycler in bold mode
            unicycler_bold_option="$final_output_subdir/unicycler_bold-option"
            mkdir -p "$unicycler_bold_option"  # Create the directory only when needed
            echo "re-running unicycler with bold-option"
            unicycler -1 "$fastp_output/trimmed_1P.fastq.gz" -2 "$fastp_output/trimmed_2P.fastq.gz" \
            -o "$unicycler_bold_option" --keep 0 --mode bold --no_rotate
            
            unicycler_final="$unicycler_bold_option"
        fi
    fi

    # Check if the output file "$unicycler_final/${plasmid_name}.fasta" contains a single FASTA sequence
         num_sequences_FINAL=$(grep -c '^>' "$unicycler_final/assembly.fasta")

        if [ "$num_sequences_FINAL" -eq 1 ]; then
            # Align the trimmed reads to the Unicycler assembly using BWA
            echo "Running Align the trimmed reads to the Unicycler assembly using BWA"
            bwa index "$unicycler_final/assembly.fasta" > "$unicycler_final/bwa_index_stdout.txt"
            bwa mem -t 8 "$unicycler_final/assembly.fasta" "$fastp_output/trimmed_1P.fastq.gz" "$fastp_output/trimmed_2P.fastq.gz" | samtools sort -o "$unicycler_final/aligned.bam" > "$unicycler_final/bwa_Samtools_aligned_stdout.txt"

            # Index the aligned BAM file using Samtools
            echo "Running Index the aligned BAM file using Samtools"
            samtools index "$unicycler_final/aligned.bam" > "$unicycler_final/Samtools_index_stdout.txt"

            # Generate coverage information using samtools depth
            echo "Generating coverage information"
            samtools depth "$unicycler_final/aligned.bam" > "$unicycler_final/coverage.txt"

            # Invoke Python script to generate the coverage plot
            echo "Generating coverage plot"
            python generate_plots.py "$unicycler_final/coverage.txt" "$unicycler_final/coverage_plot.html"

            # Set the status message for consensus generation
            status_message="Pipeline plasmid assembly for ${plasmid_name} was completed successfully."

            # Create status message HTML
            echo "<html><body><p>Status: ${status_message}</p></body></html>" > "$unicycler_final/status_message.html"

            # Activate the plannotate environment
            conda activate plannotate

            # Run pLannotate on the Unicycler assembly
            echo "Running pLannotate on the Unicycler assembly"
            plannotate batch -i "$unicycler_final/assembly.fasta" --html -o "$plannotate_output" --file_name "assembly" -s '' > "$plannotate_output/stdout.txt"

            # Rename and move the output files to the appropriate directories
            echo "Rename and move the output files to the appropriate directories"
            mv "$plannotate_output/assembly.gbk" "${output_dir}/${output_subdir}/${plasmid_name}.gbk"
            mv "$plannotate_output/assembly.html" "${output_dir}/${output_subdir}/${plasmid_name}.html"
            mv "${output_dir}/${output_subdir}/${plasmid_name}.gbk" "$plannotate_output"
            mv "${output_dir}/${output_subdir}/${plasmid_name}.html" "$plannotate_output"

            # Combine the status message, coverage plot, and pLannotate HTML files into the report
            cat "$unicycler_final/status_message.html" "$fastp_output/fastp.html" "$unicycler_final/coverage_plot.html" "$plannotate_output/${plasmid_name}.html" > "$final_output_subdir/report_of_${plasmid_name}.html"

            # Remove the intermediate HTML files
            rm "$unicycler_final/status_message.html" "$fastp_output/fastp.html" "$unicycler_final/coverage_plot.html" "$plannotate_output/${plasmid_name}.html"

        else
            # Set the status message for consensus generation
            status_message="Pipeline plasmid assembly for ${plasmid_name} was completed but the output file does not contain a single FASTA sequence"
            
            # Align the trimmed reads to the Unicycler assembly using BWA
            echo "Running Align the trimmed reads to the Unicycler assembly using BWA"
            bwa index "$unicycler_final/assembly.fasta" > "$unicycler_final/bwa_index_stdout.txt"
            bwa mem -t 8 "$unicycler_final/assembly.fasta" "$fastp_output/trimmed_1P.fastq.gz" "$fastp_output/trimmed_2P.fastq.gz" | samtools sort -o "$unicycler_final/aligned.bam" > "$unicycler_final/bwa_Samtools_aligned_stdout.txt"

            # Index the aligned BAM file using Samtools
            echo "Running Index the aligned BAM file using Samtools"
            samtools index "$unicycler_final/aligned.bam" > "$unicycler_final/Samtools_index_stdout.txt"

            # Generate coverage information using samtools depth
            echo "Generating coverage information"
            samtools depth "$unicycler_final/aligned.bam" > "$unicycler_final/coverage.txt"

            # Invoke Python script to generate the coverage plot
            echo "Generating coverage plot"
            python generate_plots.py "$unicycler_final/coverage.txt" "$unicycler_final/coverage_plot.html"

            # Create status message HTML
            echo "<html><body><p>Status: ${status_message}</p></body></html>" > "$unicycler_final/status_message.html"

            # Combine the status message, coverage plot, and pLannotate HTML files into the report
            cat "$unicycler_final/status_message.html" "$fastp_output/fastp.html" "$unicycler_final/coverage_plot.html" > "$final_output_subdir/report_of_${plasmid_name}.html"

            # Remove the intermediate HTML files
            rm "$unicycler_final/status_message.html" "$fastp_output/fastp.html" "$unicycler_final/coverage_plot.html"
        fi              

        # Deactivate the plasmid_analysis environment
        conda deactivate
        conda activate base
    echo "Sample processing complete: ${plasmid_name}"
    }

# Default reference database path
default_reference_database="/usr/src/app/Adapter_and_reference_database/reference_database"
reference_database="$default_reference_database"
adapter_fasta="/usr/src/app/Adapter_and_reference_database/Adapter/NexteraPE-PE.fa"
bcftools_PATH="/usr/src/app/bcftools-1.20/bin"

# Default output directory
output_dir="."

print_help() {
    echo "Usage: $0 [OPTIONS]"
    echo "Options:"
    echo "  -i, --input-file <file>     Path to the input zip file."
    echo "  -t, --text-file <file>      Path to the text file containing analysis information."
    echo "  -o, --output-dir <dir>      Path to the output directory. (default: current directory)"
    echo "  -DB, --reference-database <dir>  Path to the reference database directory. (default: "$default_reference_database")"
    echo "  -h, --help                  Display this help message."
    echo ""
    echo "Example usage:"
    echo "  $0 -i input.zip -t analysis.txt"
    echo ""
    echo "Example of input file format:"
    echo "Primary key   read1   read2   expected_reads  plasmid_name    output_directory    Type_of_run    Reference"
    echo "17    PS-01A03_S3_L001_R1_001.fastq.gz    PS-01A03_S3_L001_R2_001.fastq.gz    10000   pTWIST_OsPYL1   Game    y    pENTR.fasta"
    echo "18    PS-01A03_S3_L001_R1_001.fastq.gz    PS-01A03_S3_L001_R2_001.fastq.gz    10000   pTWIST_OsPYL1   Game    n"
    echo ""
    echo "Options Details:"
    echo "  -i, --input-file <file>: Path to the input zip file containing sequencing data."
    echo "  -t, --text-file <file>: Path to the text file containing analysis information."
    echo "  -o, --output-dir <dir>: Path to the output directory where analysis results will be stored. (default: current directory)"
    echo "  -DB, --reference-database <dir>: Path to the reference database directory. (default: "$default_reference_database")"
    echo "  -h, --help: Display this help message."
    echo ""
    echo "Note: Specify the type of run using 'y' for reference-guided assembly, 'n' for de novo assembly, or 's' for a special case with a vector backbone."
}

# Default output directory
output_dir="."

# Function to convert CRLF to LF
convert_line_endings() {
    local file=$1
    # Check if the file has CRLF line endings and convert them to LF
    if grep -q $'\r$' "$file"; then
        echo "Converting CRLF to LF in $file..."
        sed -i 's/\r$//' "$file"
    else
        echo "File $file is already in LF format."
    fi
}

# Read the input options
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        -i|--input-file)
            input_file="$2"
            shift 2
            ;;
        -t|--text-file)
            text_file="$2"
            shift 2
            ;;
        -o|--output-dir)
            output_dir="$2"
            shift 2
            ;;
        -DB|--reference-database)
            reference_database="$2"
            shift 2
            ;;
        -h|--help)
            print_help
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            print_help
            exit 1
            ;;
    esac
done

# Convert line endings for the text file
convert_line_endings "$text_file"

# If no reference database path is provided, use the default path
if [[ -z "$reference_database" ]]; then
    reference_database="$default_reference_database"
fi

# Check if the input file exists and is a zip file
if [[ ! -f "$input_file" || ! "$input_file" =~ \.zip$ ]]; then
    echo "Error: Input file does not exist or is not a zip file."
    exit 1
fi

# Extract the directory name without the ".zip" extension
dir_name=$(basename "${input_file%.zip}")

# Check if the input file exists and is a plain text file
if [[ ! -f "$text_file" || ! "$text_file" =~ \.txt$ ]]; then
    echo "Error: Text file does not exist or is not a plain text file."
    exit 1
fi

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

# Check if the output directory exists or create it
output_dir="$output_dir/$dir_name"
mkdir -p "$output_dir"

# Extract the contents of the zip file to the output directory
unzip -q "$input_file" -d "$output_dir"

# Define the header
header="Plasmid Name\tType of Run Label\tStart Time\tEnd Time\tEstimated_Size(kb)\tAssembly Report\tCode Report"

# Append the run information to the report file 
echo -e "$header" >> "$output_dir/Plasmid-Assemble-report.txt"

start_time_pro=$(date +"%d-%b-%Y %H:%M:%S")

# Loop through each row in the text file
tail -n +2 "$text_file" | while IFS=$'\t' read -r line; do
    # Extract the relevant columns from the line
    read1=$(echo "$line" | awk '{print $2}')
    read2=$(echo "$line" | awk '{print $3}')
    plasmid_name=$(echo "$line" | awk '{print $5}')
    output_subdir=$(echo "$line" | awk '{print $6}')
    type_of_run=$(echo "$line" | awk '{print $7}')
    reference=$(echo "$line" | awk '{print $8}')

    # Set reference to "No" if type_of_run is "n"
    if [[ "$type_of_run" == "n" ]]; then
        reference="No"
    fi

    # Construct the full path to the output subdirectory under the main output directory
    final_output_subdir="$output_dir/$output_subdir/$plasmid_name"

    # Construct the paths for input FASTQ files
    read1_path="$output_dir/$dir_name/$read1"
    read2_path="$output_dir/$dir_name/$read2"

    # Debug information
    echo "Values: read1=$read1_path, read2=$read2_path, plasmid_name=$plasmid_name, output_subdir=$output_subdir, type_of_run=$type_of_run, reference=$reference"

    start_time=$(date +"%d-%b-%Y %H:%M:%S")
    assembly_report=""

    # Convert type_of_run to lowercase for comparison
    type_of_run_lower=$(echo "$type_of_run" | tr '[:upper:]' '[:lower:]')

# For other types of runs ('y' and 'n'), use the existing code
if [[ "$type_of_run_lower" == "y" ]]; then
    if [[ -z "$reference" ]]; then
        echo "Error: Reference data file is required when a reference genome is available for the 'y' type of run."
        exit 1
    fi
    
    # Construct the path to the reference sequence
    reference_path="${reference_database}/${reference}"
    
    # Call run_plasmid_assembly_with_reference with the provided reference data file
    run_plasmid_assembly_with_reference "$reference_path" "$bcftools_PATH"
    #run_plasmid_assembly_with_reference "$reference_path"

    type_of_run_label="Reference"
    
    # Check for consensus sequence file in the specified path
    consensus_file="$final_output_subdir/${plasmid_name}_consensus.fa"
    if [[ -f "$consensus_file" ]]; then
        assembly_report="Consensus sequence file available"
        code_report="RP"  # Reference Success
    else
        assembly_report="Consensus sequence file not available"
        code_report="RF"  # Reference Unsuccess
    fi
    # Calculate the size of the consensus file
    size=$(python3 calculate_fasta.py "$consensus_file")

elif [[ "$type_of_run_lower" == "n" ]]; then
    # Call run_plasmid_assembly_without_reference
    run_plasmid_assembly_without_reference

    type_of_run_label="Denovo"
    
    # Check for the assembly.fasta file in the unicycler_output directory
    assembly_fasta="${unicycler_final}/assembly.fasta"
    if [[ -f "$assembly_fasta" ]]; then
        num_fasta_sequences=$(grep -c '^>' "$assembly_fasta")
        assembly_report="Unicycler: The assembly file contains $num_fasta_sequences sequence(s) in the file"
    else
        assembly_report="Unicycler: The assembly file not found"
    fi

    # Calculate the size of the assembly FASTA file
    size=$(python3 calculate_fasta.py "$assembly_fasta")

    
    # Check for the .html file in the plannotate_output directory
    plannotate_html="${final_output_subdir}/plannotate_output/*.gbk"
    if ls $plannotate_html 1> /dev/null 2>&1; then
        assembly_report="${assembly_report}, The annotation file is available in the plannotate output folder"
        code_report="DP"  # Denovo Success
    else
        assembly_report="${assembly_report}, The annotation file is not available"
        code_report="DF"  # Denovo Unsuccess
    fi

else
    echo "Error: Unknown value in Type_of_run column of the text file."
    exit 1
fi

# Get the end time
end_time=$(date +"%d-%b-%Y %H:%M:%S")

# Define the header
# header="Plasmid Name\tType of Run Label\tStart Time\tEnd Time\tAssembly Report"

# Append the run information to the report file 
# echo -e "$header" >> "$output_dir/Plasmid-Assemble-report.txt"
echo -e "$plasmid_name\t$type_of_run_label\t$start_time\t$end_time\t$size\t$assembly_report\t$code_report" >> "$output_dir/Plasmid-Assemble-report.txt"

done

# Check if the Reference directory exists
if [[ ! -d "$reference_database" ]]; then
    echo "Error: Reference directory not found."
    exit 1
fi

end_time_pro=$(date +"%d-%b-%Y %H:%M:%S")


# Function to generate the bar graph using Python script
generate_summary_plot() {
    python3 generate_summary_plot.py "$output_dir/Plasmid-Assemble-report.txt" "$output_dir" "$reference_success_count" "$reference_unsuccess_count" "$denovo_success_count" "$denovo_unsuccess_count"
}

# Read input values
output_dir="$output_dir"  # Set output directory path
total_running_data="$dir_name"
start_time="$start_time_pro"
end_time="$end_time_pro"

# Initialize counters
denovo_success_count=0
denovo_unsuccess_count=0
reference_success_count=0
reference_unsuccess_count=0

# Read and process the report file, skipping the header line
while IFS=$'\t' read -r line; do
    # Skip header
    if [ "$line" == "Plasmid Name	Type of Run Label	Start Time	End Time	Assembly Report	Code Report" ]; then
        continue
    fi
    
    # Extract the relevant columns from the line
    type_of_run_label=$(echo "$line" | awk '{print $2}')
    last_column=$(echo "$line" | awk '{print $NF}')
    
    # Print debug output
    echo "Type of Run Label: $type_of_run_label, Last Column: $last_column"

    # Update counts based on the extracted values
    if [ "$type_of_run_label" == "Reference" ]; then
        if [ "$last_column" == "RP" ]; then
            ((reference_success_count++))
        elif [ "$last_column" == "RF" ]; then
            ((reference_unsuccess_count++))
        fi
    elif [ "$type_of_run_label" == "Denovo" ]; then
        if [ "$last_column" == "DP" ]; then
            ((denovo_success_count++))
        elif [ "$last_column" == "DF" ]; then
            ((denovo_unsuccess_count++))
        fi
    fi
done < "$output_dir/Plasmid-Assemble-report.txt"

# Print counts outside of the loop
echo "Reference Success Count: $reference_success_count"
echo "Reference Unsuccess Count: $reference_unsuccess_count"
echo "Denovo Success Count: $denovo_success_count"
echo "Denovo Unsuccess Count: $denovo_unsuccess_count"
# Call the function to generate the summary plot
generate_summary_plot

# Define the path to the bar graph image
bar_graph_path="success_bar_graph.png"

# Define label meanings for colors
label_meanings=("Reference Success" "Reference Unsuccess" "Denovo Success" "Denovo Unsuccess" )
color_codes=("#ADF802" "#FFCE44" "#1589FF" "#797979")

# Generate the HTML report
echo "<html>
<head>
  <title>Plasmid Assembly Report</title>
  <style>
    .color-box {
      display: inline-block;
      width: 10px;
      height: 10px;
      margin-right: 5px;
    }
    .reference-success { background-color: ${color_codes[0]}; }
    .reference-unsuccess { background-color: ${color_codes[1]}; }
    .denovo-success { background-color: ${color_codes[2]}; }
    .denovo-unsuccess { background-color: ${color_codes[3]}; }
  </style>
</head>
<body>
  <h1>Plasmid Assembly Report</h1>
  <p>Total Running Data: $total_running_data</p>
  <p>Number of De Novo Success: $denovo_success_count</p>
  <p>Number of De Novo Unsuccess: $denovo_unsuccess_count</p>
  <p>Number of Reference Success: $reference_success_count</p>
  <p>Number of Reference Unsuccess: $reference_unsuccess_count</p>
  <p>Start Running Time: $start_time</p>
  <p>Finish Running Time: $end_time</p>
  <h2>Summary of Running</h2>
  <p>Label Meanings:</p>
  <ul>
    <li><span class='color-box reference-success'></span>${label_meanings[0]}</li>
    <li><span class='color-box reference-unsuccess'></span>${label_meanings[1]}</li>
    <li><span class='color-box denovo-success'></span>${label_meanings[2]}</li>
    <li><span class='color-box denovo-unsuccess'></span>${label_meanings[3]}</li>
  </ul>
  <img src='$bar_graph_path' alt='Summary Plot'>
</body>
</html>" > "$output_dir/report.html"

# Clean up the database
deleted_count=$(find "$reference_database" \( -name "*.fasta.fai" -o -name "*.bt2" \) -type f -delete -print | wc -l)
echo "Deleted $deleted_count files in the Reference database."
