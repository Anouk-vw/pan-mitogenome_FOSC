# -*- coding: utf-8 -*-
"""
Fasta file combiner
Speaks for itself
10/01/23 Kyran Wissink
"""

import sys
import os
import glob

def combine_fasta(directory: str) -> str:
    """
    Combine the fasta files in the given directory.

    Parameters
    ----------
    directory : str
        The directory containing the fasta files.

    Returns
    -------
    str
        The path to the combined fasta file.
    """

    # Add forward slash if not present in the config file
    if not directory.endswith("/"):
        directory += "/"

    # Set the output file name
    output_file = directory + 'combined.fa'


    # Skip this step if it was already done
    if os.path.isfile(output_file):
        print(f"'{output_file}' already exists, skipping file combination.")
        return output_file
    else:
        print(f"Combining fasta files in {directory} to {output_file}")


    # Check if the directory exists
    if not os.path.isdir(directory):
        print(f"Directory '{directory}' does not exist.")
        exit()

    # Find all .fasta files in the directory
    fasta_files = glob.glob(os.path.join(directory, '*.fa*'))
    if not fasta_files:
        print(f"No .fasta files found in directory '{directory}'.")
        exit()



    # Open the output file in append mode
    with open(output_file, 'a') as out_file:
        # Iterate through all fasta files in the directory
        for filename in fasta_files:
            file_name = os.path.splitext(os.path.basename(filename))[0]
            # Open each fasta file and write its contents to the output file
            with open(filename, 'r') as fasta_file:
                for line in fasta_file:
                    # If the line is a header line, modify it to include the original file name
                    if line.startswith('>'):
                        out_file.write(f'>{file_name}_{line[1:]}')
                    else:
                        out_file.write(line)



    print(f"Combined {len(fasta_files)} .fasta files into '{output_file}'.")


if __name__ == "__main__":
    if os.path.isdir(sys.argv[1]):
        combine_fasta(sys.argv[1])

