"""
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simple RNA Sequencing Pipeline
# Author: Weronika Jaśkowiak, Maike Nägele, Tabea Attig, Annika Maier
# Date: 11.10.2024
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""

"""
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Python Script: validate_samplesheet
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Description:
# This script validates a samplesheet CSV file used for processing paired-end FASTQ files.
# It checks for the correct format of the samplesheet, ensuring that all required fields 
# are present and that the specified FASTQ files exist and are correctly named.
#
# Input:
# - str samplesheet_file
#   - Path to the samplesheet CSV file to be validated.
#
# Output:
# - validation.txt
#   - A text file indicating the result of the validation process. It contains either
#     "Samplesheet validation passed!" or "Samplesheet validation failed!" based on the 
#     checks performed on the input samplesheet.
#
# Validation Checks Performed:
# - Ensure the samplesheet header matches the expected format: ["sample", "fastq_1", "fastq_2", "strandedness"].
# - Verify that no values are missing in any row of the samplesheet.
# - Check that both FASTQ files exist and are named correctly, with extensions '.fastq.gz'.
# - Confirm that the paths to the FASTQ files are correctly constructed based on the script's execution context.
#
# Error Messages:
# - Prints error messages to the console for any validation failures encountered during the 
#   validation process, detailing the specific issues found for each sample.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""


import os  # Module to interact with the operating system
import csv  # Module for reading and writing CSV files
import argparse  # Module for parsing command-line arguments


def validate_samplesheet(samplesheet_file):
    valid = True  # Variable to track overall validation status

    # Open the samplesheet CSV file
    with open(samplesheet_file, 'r') as file:
        reader = csv.DictReader(file)  # Create a CSV reader that maps the data into dictionaries

        # Define the expected header columns for the samplesheet
        expected_columns = ["sample", "fastq_1", "fastq_2", "strandedness"]

        # Check if the header matches the expected format
        if reader.fieldnames != expected_columns:
            print(f"ERROR: Header does not match expected format: {expected_columns}")
            valid = False  # Set valid to False if the header is incorrect

        # Validate each row in the samplesheet
        for row in reader:
            sample = row["sample"]  # Sample name
            fastq_1 = row["fastq_1"]  # Path to the first FASTQ file
            fastq_2 = row["fastq_2"]  # Path to the second FASTQ file
            strandedness = row["strandedness"]  # Strandedness information

            # Check for missing values in the row
            if not all([sample, fastq_1, fastq_2, strandedness]):
                print(f"ERROR: Missing values in the row for sample '{sample}'. All columns must be filled.")
                valid = False  # Set valid to False if any value is missing
                continue  # Skip further checks for this row

            # Construct the correct file paths (assuming the script is run from 'nf-core/bin')
            fastq_1_path = os.path.join("../../../", fastq_1)  # Adjust path for fastq_1
            fastq_2_path = os.path.join("../../../", fastq_2)  # Adjust path for fastq_2

            # Check if both fastq_1 and fastq_2 files exist (paired-end check)
            if not fastq_1 or not fastq_2:
                print(f"ERROR: Paired-end data missing for sample {sample}. Both fastq_1 and fastq_2 are required.")
                valid = False
            else:
                # Check if fastq_1 ends with .fastq.gz and if the file exists
                if not fastq_1.endswith(".fastq.gz"):
                    print(f"ERROR: File '{fastq_1}' for sample '{sample}' does not end with '.fastq.gz'.")
                    valid = False
                elif not os.path.exists(fastq_1_path):
                    print(
                        f"ERROR: File '{os.path.basename(fastq_1)}' for sample '{sample}' does not exist at '{fastq_1_path}'.")
                    valid = False

                # Check if fastq_2 ends with .fastq.gz and if the file exists
                if not fastq_2.endswith(".fastq.gz"):
                    print(f"ERROR: File '{fastq_2}' for sample '{sample}' does not end with '.fastq.gz'.")
                    valid = False
                elif not os.path.exists(fastq_2_path):
                    print(
                        f"ERROR: File '{os.path.basename(fastq_2)}' for sample '{sample}' does not exist at '{fastq_2_path}'.")
                    valid = False

    file_name = "validation.txt"  # Output file name for validation results

    # Final output after validation
    if valid:
        with open(file_name, 'w') as file:
            # Write success message to the validation file
            file.write("Samplesheet validation passed!")
    else:
        with open(file_name, 'w') as file:
            # Write failure message to the validation file
            file.write("Samplesheet validation failed!")


# Set up command-line argument parsing
parser = argparse.ArgumentParser(description="Validate samplesheet CSV file.")
parser.add_argument("-file", required=True,
                    help="Path to the samplesheet CSV file")  # Add argument for the samplesheet file
args = parser.parse_args()  # Parse the command-line arguments

# Call the validation function with the provided samplesheet file path
validate_samplesheet(args.file)