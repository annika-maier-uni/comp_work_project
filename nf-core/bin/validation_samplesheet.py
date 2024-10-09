import os
import csv

def validate_samplesheet(samplesheet_file):
    valid = True

    # Open the samplesheet
    with open(samplesheet_file, 'r') as file:
        reader = csv.DictReader(file)

        # Check if header matches the expected format
        expected_columns = ["sample", "fastq_1", "fastq_2", "strandedness"]
        if reader.fieldnames != expected_columns:
            print(f"ERROR: Header does not match expected format: {expected_columns}")
            valid = False

        # Validate each row in the samplesheet (only check fastq files)
        for row in reader:
            sample = row["sample"]
            fastq_1 = row["fastq_1"]
            fastq_2 = row["fastq_2"]
            strandedness = row["strandedness"]

            # Check if both fastq_1 and fastq_2 files exist (paired-end check)
            if not fastq_1 or not fastq_2:
                print(f"ERROR: Paired-end data missing for sample {sample}. Both fastq_1 and fastq_2 are required.")
                valid = False
            else:
                # Check if fastq_1 ends with .fastq.gz and file exists
                if not fastq_1.endswith(".fastq.gz"):
                    print(f"ERROR: File '{fastq_1}' for sample {sample} does not end with '.fastq.gz'.")
                    valid = False
                elif not os.path.exists(fastq_1):
                    print(f"ERROR: File '{fastq_1}' for sample {sample} does not exist.")
                    valid = False

                # Check if fastq_2 ends with .fastq.gz and file exists
                if not fastq_2.endswith(".fastq.gz"):
                    print(f"ERROR: File '{fastq_2}' for sample {sample} does not end with '.fastq.gz'.")
                    valid = False
                elif not os.path.exists(fastq_2):
                    print(f"ERROR: File '{fastq_2}' for sample {sample} does not exist.")
                    valid = False

    # Final output after validation
    if valid:
        print("Samplesheet validation passed!")
    else:
        print("Samplesheet validation failed!")

# Run the validation
if __name__ == "__main__":
    samplesheet_path = "path\to\comp_work_project\nf-core\assets\samplesheet.csv"
    validate_samplesheet(samplesheet_path)
