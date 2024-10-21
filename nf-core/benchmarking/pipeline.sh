#!/bin/bash

# Sample files locations
declare -A forward_reads
declare -A reverse_reads

forward_reads[SRR6357070]="../data/SRR6357070_1.fastq.gz"
reverse_reads[SRR6357070]="../data/SRR6357070_2.fastq.gz"
forward_reads[SRR6357071]="../data/SRR6357071_1.fastq.gz"
reverse_reads[SRR6357071]="../data/SRR6357071_2.fastq.gz"
forward_reads[SRR6357072]="../data/SRR6357072_1.fastq.gz"
reverse_reads[SRR6357072]="../data/SRR6357072_2.fastq.gz"
forward_reads[SRR6357073]="../data/SRR6357073_1.fastq.gz"
reverse_reads[SRR6357073]="../data/SRR6357073_2.fastq.gz"
forward_reads[SRR6357074]="../data/SRR6357074_1.fastq.gz"
reverse_reads[SRR6357074]="../data/SRR6357074_2.fastq.gz"
forward_reads[SRR6357075]="../data/SRR6357075_1.fastq.gz"
reverse_reads[SRR6357075]="../data/SRR6357075_2.fastq.gz"
forward_reads[SRR6357076]="../data/SRR6357076_1.fastq.gz"
reverse_reads[SRR6357076]="../data/SRR6357076_2.fastq.gz"
forward_reads[SRR6357077]="../data/SRR6357077_1.fastq.gz"
reverse_reads[SRR6357077]="../data/SRR6357077_2.fastq.gz"

# Set pipeline parameters
max_memory='16GB'
max_cpus=4

# Activate a conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate comp_work_benchmark

# Run FastQC on each sample
echo "Running FastQC"
mkdir -p fastqc_results
for id in "${!forward_reads[@]}"; do
  fastq_1=${forward_reads[$id]}
  fastq_2=${reverse_reads[$id]}
  reads="$fastq_1 $fastq_2"
  fastqc $reads --threads $max_cpus -o fastqc_results -q
done

# Deactivate FastQC environment
conda deactivate

# Activate trimming environment
echo "Running Trim Galore"
mkdir -p trimmed_results
conda activate comp_work_benchmark_trim
for id in "${!forward_reads[@]}"; do
  prefix="../data/${id}"
  trim_galore --paired ${prefix}_1.fastq.gz ${prefix}_2.fastq.gz --fastqc -o trimmed_results --suppress_warn
done

# Deactivate trimming environment
conda deactivate

# Activate benchmarking environment
conda activate comp_work_benchmark

# Hisat2 build
echo "Creating HISAT2 index"
mkdir -p hisat_index
reference="../data/genome.fa"
name="hisat_index/benchmark_test_hisat"
hisat2-build ${reference} ${name} --quiet 

# Hisat2 align
echo "Running HISAT2 alignment"
hisat2_index="hisat_index/benchmark_test_hisat"
output_sam="output_hisat2_aligned_sam_file.sam"
fastq1_list=$(ls trimmed_results/*_1.fq.gz | tr '\n' ',')
fastq2_list=$(ls trimmed_results/*_2.fq.gz | tr '\n' ',')

hisat2 --fast -x ${hisat2_index} -1 ${fastq1_list} -2 ${fastq2_list} -S ${output_sam} \
  --rg-id "benchmark_test" \
  --rg "SM:benchmark_test" \
  --rg "LB:benchmark_test_lib" \
  --rg "PL:unknown" \
  --rg "PU:benchmark_test_unit" \
  --threads ${max_cpus}

# Samtools sort
echo "Sorting SAM file with Samtools"
mkdir -p sorted_results
input_sam="output_hisat2_aligned_sam_file.sam"
output_sorted_sam="sorted_results/sorted.sam"
samtools sort ${input_sam} -o ${output_sorted_sam}

# Picard mark duplicates
# Pull Picard Docker image
echo "Pulling Picard Docker image"
docker pull docker.io/broadinstitute/picard

echo "Marking duplicates with Picard"
mkdir -p picard_results
docker run --rm -v $(pwd):/data docker.io/broadinstitute/picard java -Xmx12g -jar /usr/picard/picard.jar MarkDuplicates \
  -I /data/${output_sorted_sam} \
  -O /data/picard_results/marked_duplicates.bam \
  -M /data/picard_results/metrics.txt

# Monitor time and resource usage for a few processes (run once)
echo "Monitoring processes..."
ps -eo pid,comm,%cpu,%mem --sort=-%cpu | head -n 5

