# -*- coding: utf-8 -*-
"""
Created on Sat Jan 14 13:57:30 2023

@author: Kyran Wissink

Python script for initialising parameters if not supplied directly

"""

import os
import glob
import subprocess
from typing import Tuple
import yaml
from yaml.loader import SafeLoader


def load_config_file(file_path: str) -> dict:
    """
    Load a YAML configuration file using the SafeLoader from PyYAML.

    Args:
        file_path (str): Path to the configuration file.

    Returns:
        dict: The parsed configuration data.
    """
    with open(file_path) as f:
        config = yaml.load(f, Loader=SafeLoader)

    return config


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

    return output_file


def missing_config_params(config: dict) -> list:
    """
    Return a list of missing parameters from the configuration.

    Args:
        config (dict): The configuration data.

    Returns:
        list: A list of missing parameter names.
    """
    missing_params = []
    for param in config["pggb"]:
        if not config["pggb"][param]:
            missing_params.append(param)    

    return missing_params


def update_config_file(config: dict, update_list: list) -> dict:
    """
    Update the configuration file with automatically generated parameters if not supplied.

    Parameters
    ----------
    config : dict
        The configuration file of the workflow in YAML format.
    update_list : list
        The list of parameters to be updated.

    Returns
    -------
    dict
        The updated configuration file.
    """
    if "percent_identity" in update_list or "haplotypes" in update_list:
        print("Calculating percentage identity...")
        lowest_percent_identity, haplotypes = mash_triangle(config)

        if "haplotypes" in update_list:
            config["pggb"]["haplotypes"] = haplotypes
        if "percent_identity" in update_list:
            config["pggb"]["percent_identity"] = lowest_percent_identity

    if "poa_params" in update_list:
        poa_params = get_poa_params(config["pggb"]["percent_identity"])
        config["pggb"]["poa_params"] = poa_params

    return config


def mash_triangle(config: dict) -> Tuple[float, int]:
    """
    Calculate the lowest percent identity of genomes in the input data using the 'mash triangle' command.

    Parameters
    ----------
    config : dict
        The configuration file of the workflow in YAML format.

    Returns
    -------
    Tuple[float, int]
        The lowest percent identity of genomes, for use in PGGB data and the number of haplotypes found in the datafile.
    """
    # Run mash triangle to get percent identities
    mash_matrix = subprocess.run(["mash", "triangle", config["sample"]], stdout=subprocess.PIPE)

    # Decode and split Unix shell output into a list
    mash_list = mash_matrix.stdout.decode("utf-8").split()
    haplotypes = int(mash_list.pop(0))  # Number of genomes is first in the list

    # Get the max divergence from the mash_list and convert to lowest_pct_identity
    mash_list = [float(x) for x in mash_list if x.replace(".", "").isdigit()]
    max_divergence = max(mash_list)
    lowest_percent_identity = 100 - max_divergence * 100
    print(lowest_percent_identity)

    return lowest_percent_identity, haplotypes


def get_poa_params(percent_identity: float) -> str:
    """
    Get the poa parameters based on the percent identity.
    This is based on a loose recommendation by the PGGB
    documentation.

    Parameters
    ----------
    percent_identity : float
        The lowest percent identity of genomes.

    Returns
    -------
    str
        The poa parameters.
    """
    if percent_identity > 99:
        poa_params = "asm5"
    elif percent_identity > 90:
        poa_params = "asm10"
    else:
        poa_params = "asm20"

    return poa_params


def write_config_file(file_path: str, config: dict):
    """
    Write the configuration data to a YAML file.

    Args:
        file_path (str): Path to the configuration file.
        config (dict): The configuration data.
    """
    with open(file_path, 'w') as file:
        yaml.dump(config, file)


def main():
    """
    Main function that loads a configuration file, updates it with missing
    parameters, and writes the updated configuration data back to the file.

    """
    config = load_config_file('config.yaml')

    if os.path.isdir(config["sample"]):
        config["sample"] = combine_fasta(config["sample"])

    missing_params = missing_config_params(config)
    if missing_params:
        config = update_config_file(config, missing_params)

    write_config_file('config.yaml', config)


if __name__ == "__main__":
    main()
