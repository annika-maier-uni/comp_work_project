import os
import csv
import argparse

def validate_samplesheet(samplesheet_file):
    valid = True

    # Check if the samplesheet file exists
    if not os.path.exists(samplesheet_file):
        print(f"ERROR: The file '{samplesheet_file}' does not exist.")
        return False

    # Open the samplesheet
    with open(samplesheet_file, 'r') as file:
        reader = csv.DictReader(file)

        # Check if header matches the expected format
        expected_columns = ["sample", "fastq_1", "fastq_2", "strandedness"]
        if reader.fieldnames != expected_columns:
            print(f"ERROR: Header does not match expected format. Expected columns are: {expected_columns}")
            print(f"Found columns are: {reader.fieldnames}")
            valid = False

        # Validate each row in the samplesheet (only check fastq files)
        for row in reader:
            sample = row["sample"]
            fastq_1 = row["fastq_1"]
            fastq_2 = row["fastq_2"]
            strandedness = row["strandedness"]

            # Check if both fastq_1 and fastq_2 files are provided (paired-end check)
            if not fastq_1 or not fastq_2:
                print(f"ERROR: Paired-end data missing for sample '{sample}'. Both fastq_1 and fastq_2 are required.")
                valid = False
            else:
                # Check if fastq_1 ends with .fastq.gz and file exists
                if not fastq_1.endswith(".fastq.gz"):
                    print(f"ERROR: File '{fastq_1}' for sample '{sample}' does not end with '.fastq.gz'.")
                    valid = False
                elif not os.path.exists(fastq_1):
                    print(f"ERROR: File '{fastq_1}' for sample '{sample}' does not exist.")
                    valid = False

                # Check if fastq_2 ends with .fastq.gz and file exists
                if not fastq_2.endswith(".fastq.gz"):
                    print(f"ERROR: File '{fastq_2}' for sample '{sample}' does not end with '.fastq.gz'.")
                    valid = False
                elif not os.path.exists(fastq_2):
                    print(f"ERROR: File '{fastq_2}' for sample '{sample}' does not exist.")
                    valid = False

    # Final output after validation
    if valid:
        print("Samplesheet validation passed!")
    else:
        print("Samplesheet validation failed!")

def main():
    # Set up argument parsing
    parser = argparse.ArgumentParser(description="Validate a samplesheet CSV file.")
    parser.add_argument('--samplesheet', required=True, help="Path to the samplesheet CSV file")
    args = parser.parse_args()

    # Run the validation
    validate_samplesheet(args.samplesheet)

if __name__ == "__main__":
    main()
