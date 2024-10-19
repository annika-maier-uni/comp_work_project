/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simple RNA Sequencing Pipeline
# Author: Weronika Jaśkowiak, Maike Nägele, Tabea Attig, Annika Maier
# Date: 11.10.2024
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Nextflow Process: SAMPLESHEET VALIDATION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Description:
# The `SAMPLESHEET_VALIDATION` process is responsible for checking the validity of
# a samplesheet using a Python validation script. The validation script verifies the
# format and checks the existence of the specified FASTQ files. If the validation
# fails, the pipeline will stop, and an error message will be returned.
#
#
# Input:
# - path(python): Path to the Python validation script.
# - path(samplesheet): Path to the samplesheet file that needs to be validated.
#
# Output:
# - path "validation.txt": A text file that contains the results of the validation.
#                          If the validation passes, it contains "Samplesheet validation passed!".
#                          If it fails, it contains "Samplesheet validation failed!".
#
# Execution:
# - The process runs the Python script, captures the output in "validation.txt", and
#   checks whether the validation passed or failed based on the contents of that file.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


process SAMPLESHEET_VALIDATION {

    debug true

    input:
    path python
    path samplesheet

    output:
    path samplesheet

    script:
    """
    python $python -file $samplesheet

    if grep -q "failed" "validation.txt"; then
        echo "Samplesheet validation failed!"
        exit 1
    else
        echo "Samplesheet validation passed!"
    fi
    """
}
