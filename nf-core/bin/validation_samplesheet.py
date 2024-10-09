import os
import csv

def validate_samplesheet(samplesheet_file):
    valid_strandedness = ["unstranded", "forward", "reverse", "auto"]
    valid = True

    # Open the samplesheet
    with open(samplesheet_file, 'r') as file:
        reader = csv.DictReader(file)

        # Check if header matches the expected format
        expected_columns = ["sample", "fastq_1", "fastq_2", "strandedness"]
        if reader.fieldnames != expected_columns:
            print(f"ERROR: Header does not match expected format: {expected_columns}")
            valid = False

        # Validate each row
        for row in reader:
            sample = row["sample"]
            fastq_1 = row["fastq_1"]
            fastq_2 = row["fastq_2"]
            strandedness = row["strandedness"]

            # Check if strandedness is valid
            if strandedness not in valid_strandedness:
                print(f"ERROR: Invalid strandedness value '{strandedness}' for sample {sample}")
                valid = False

            # Check if both fastq_1 and fastq_2 files exist (paired-end check)
            if not fastq_1 or not fastq_2:
                print(f"ERROR: Paired-end data missing for sample {sample}. Both fastq_1 and fastq_2 are required.")
                valid = False
            else:
                # Check if fastq_1 file exists
                if not os.path.exists(fastq_1):
                    print(f"ERROR: File '{fastq_1}' for sample {sample} does not exist")
                    valid = False

                # Check if fastq_2 file exists
                if not os.path.exists(fastq_2):
                    print(f"ERROR: File '{fastq_2}' for sample {sample} does not exist")
                    valid = False

    if valid:
        print("Samplesheet validation passed!")
    else:
        print("Samplesheet validation failed!")

# Run the validation
validate_samplesheet("C:\Users\tabat\MasterBioinformatik\Semester3\ComputationalWorkflows\Project\comp_work_project\nf-core\assets\samplesheet.csv")
