# ATG_fastq_demultiplexer.py 
# Version: 20240413
# Author: Amir Taheri-Ghahfarokhi
# Contact: 
#        Email: Amir.Taheri.Ghahfarokhi@Gmail.com
#        Linkedin: https://www.linkedin.com/in/ghahfarokhi/
# 
# Acronyms:
# bol = begining of line
# eol = end of line

import pandas as pd
import os
import sys
import argparse
import subprocess
import shutil

def check_barcodes_info_file(barcodes_info_file):
    # Check if the input file exists
    if not os.path.isfile(barcodes_info_file):
        sys.exit(f"ERROR: Input file '{barcodes_info_file}' not found.")

    df = pd.read_csv(barcodes_info_file, delimiter="\t")
    
    # Lower case column names 
    df.columns = map(str.lower, df.columns)

    # Upper case barcodes 
    df['barcode_bol'] = df['barcode_bol'].str.upper()
    df['barcode_eol'] = df['barcode_eol'].str.upper()

    # Check the required columns
    required_columns = ["sample_name", "barcode_bol", "barcode_eol", "demuxed_name"]
    input_columns = df.columns.tolist()
    missing_columns = [col for col in required_columns if col not in input_columns]

    if missing_columns:
        sys.exit(f"ERROR: Input file is missing the following required columns: {', '.join(missing_columns)}")
    else:
        return df

def create_bcfiles_for_fastx_splitter(df, intermediate_dir):
    df['sample_bol'] = df['sample_name'] + '_' + df['barcode_bol']
    df['sample_bol_eol'] = df['sample_bol'] + '_' + df['barcode_eol']

    os.makedirs(intermediate_dir, exist_ok=True)

    # Group by sample_name and count unique barcode_bol values
    sample_bol_counts = df.groupby("sample_name")["barcode_bol"].nunique()

    for sample_name, bol_count in sample_bol_counts.items():
        sample_df = df[df["sample_name"] == sample_name]
        unique_bol_barcodes = sample_df["barcode_bol"].unique()

        # Write the unique barcode_bol values to a txt file
        with open(os.path.join(intermediate_dir, sample_name + "_bol_barcodes.txt"), "w") as bol_file:
            for barcode_bol in unique_bol_barcodes:
                bol_file.write(f"{sample_name}_{barcode_bol}\t{barcode_bol}\n")

    # Group by sample_bol and count unique barcode_eol values
    sample_bol_eol_counts = df.groupby("sample_bol")["barcode_eol"].nunique()

    for sample_bol, eol_count in sample_bol_eol_counts.items():
        sample_bol_df = df[df["sample_bol"] == sample_bol]
        unique_eol_barcodes = sample_bol_df["barcode_eol"].unique()

        # Write the unique barcode_eol values to a txt file
        
        with open(os.path.join(intermediate_dir, sample_bol + "_eol_barcodes.txt"), "w") as bol_eol_file:
            for barcode_eol in unique_eol_barcodes:
                demuxed_name = sample_bol_df[sample_bol_df["barcode_eol"] == barcode_eol]["demuxed_name"].iloc[0]
                bol_eol_file.write(f"{demuxed_name}\t{barcode_eol}\n")

    return df    

def split_fastq_files(df, fastq_dir, intermediate_dir, output_dir):
    if not os.path.exists(output_dir):
            os.makedirs(output_dir)

    log_file = open(f"{output_dir}fastx_splitter.log", 'a')

    bol_df = df.drop_duplicates(['sample_name'])
    bol_df.sort_values(by=['sample_name'])
    
    for index, row in bol_df.iterrows():
        sample_name = row['sample_name']
        bcfile = os.path.join(intermediate_dir, sample_name + "_bol_barcodes.txt")
        print(f"Splitting bol barcodes for {sample_name} ...")
        fastq_file = os.path.join(fastq_dir, sample_name + ".fastq")
        cat = subprocess.Popen(("cat", fastq_file), stdout=subprocess.PIPE)
        try:
            splitter = subprocess.check_output(("fastx_barcode_splitter.pl", "--bcfile", bcfile, "--prefix", intermediate_dir, "--suffix", ".fastq", "--bol"), 
                                               stdin=cat.stdout, universal_newlines=True)
            cat.wait()
        except subprocess.CalledProcessError as exc:
            log_file.write(f"{sample_name} : failed!")
            print(f"{sample_name} : failed!")
        else:
            log_file.write("\n{}\n".format(splitter))
            print("\n{}\n".format(splitter))
        
        eol_df = df.loc[df['sample_name'] == sample_name, ]
        eol_df = eol_df.drop_duplicates(['sample_bol'])
        eol_df.sort_values(by=['barcode_eol'])

        for index, row in eol_df.iterrows():
            sample_bol = row['sample_bol']
            bcfile = os.path.join(intermediate_dir, sample_bol + "_eol_barcodes.txt")
            print(f"Splitting eol barcodes for {sample_bol} ...")
            
            cat = subprocess.Popen(("cat", os.path.join(intermediate_dir, sample_bol + ".fastq")), stdout=subprocess.PIPE)

            try:
                splitter = subprocess.check_output(("fastx_barcode_splitter.pl", "--bcfile", bcfile, "--prefix", output_dir, "--suffix", ".fastq", "--eol"), 
                                                stdin=cat.stdout, universal_newlines=True)
                cat.wait()
            except subprocess.CalledProcessError as exc:
                log_file.write(f"{sample_bol} : failed!")
                print(f"{sample_bol} : failed!")
            else:
                log_file.write("\n{}\n".format(splitter))
                print("\n{}\n".format(splitter))

def print_help():
    help = """
    This script uses barcode splitter from fastx toolkit to demultiplex fastq files.
    Author: Amir.Taheri.Ghahfarokhi@Gmail.com
    Version: 20240413

    Requirements:
    - fastx toolkit (fastx_barcode_splitter.pl must be in the PATH)

    usage:
    atg_fastq_demultiplexer.py [-h]
    OR
    atg_fastq_demultiplexer.py [--fastq-dir PATH_TO_FASTQ_DIR] [--barcodes-file PATH_TO_BARCODES_FILE] [--out-dir PATH_TO_OUT_DIR]
    OR
    atg_fastq_demultiplexer.py [-i PATH_TO_FASTQ_DIR] [-f PATH_TO_BARCODES_FILE] [-o PATH_TO_OUT_DIR]

    Important note 1:
    Fastq files in the FASTQ_DIR are expected to have ".fastq" suffix and their names match the "sample_name" in BARCODES_FILE.

    Important note 2:
    BARCODES_FILE is expected to be tab separated and have at least these four columns:
    sample_name barcode_bol barcode_eol demuxed_name

    Important Note 3: 
    sample_name and demuxed_name must be alphanumeric (i.e., avoid &, $, @, -, %, * and empty space in the file names)

    Acronyms: 
    bol = begining of line
    eol = end of line

    Reference for fastx_toolkit
    http://hannonlab.cshl.edu/fastx_toolkit/commandline.html
    """
    print(help)

def parse_arguments():
    if len(sys.argv) == 1 or '-h' in sys.argv or '--help' in sys.argv:
        print_help()
        sys.exit(1)
    
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument('-i', "--fastq-dir")
    parser.add_argument('-f', "--barcodes-file")
    parser.add_argument('-o', "--out-dir")
    args = parser.parse_args()
    
    if args.fastq_dir is None:
        print("Please provide a fastq directory using [--fastq-dir PATH_TO_FASTQ_DIR] or [-i PATH_TO_FASTQ_DIR].")
        sys.exit(1)
    
    if args.barcodes_file is None:
        print("Please provide a barcode_info file using [--barcodes-file PATH_TO_BARCODES_FILE] or [-f PATH_TO_BARCODES_FILE].")
        sys.exit(1)
    
    if args.out_dir is None:
        print("Please provide a path to a directory for saving the outputs using [--out-dir PATH_TO_DIR] or [-o PATH_TO_OUT_DIR].")
        sys.exit(1)

    return args

def main():

    args = parse_arguments()

    # Ensure that directories have abs path
    intermediate_dir = os.path.abspath(os.path.join(os.getcwd(), 'atg_fastq_demultiplexer_intermediate_dir'))
    output_dir =  os.path.abspath(args.out_dir)
    fastq_dir = os.path.abspath(args.fastq_dir)

    # Ensure that directories ends with slash
    intermediate_dir = os.path.join(intermediate_dir, '')
    output_dir =  os.path.join(output_dir, '')
    fastq_dir = os.path.join(fastq_dir, '')

    barcodes_file = args.barcodes_file

    df = check_barcodes_info_file(barcodes_file)

    df = create_bcfiles_for_fastx_splitter(df, intermediate_dir)

    split_fastq_files(df, fastq_dir, intermediate_dir, output_dir)
    
    shutil.rmtree(intermediate_dir)

if __name__ == "__main__":
    main()