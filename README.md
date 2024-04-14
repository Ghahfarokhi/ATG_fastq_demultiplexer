## ATG_fastq_demultiplexer

This script uses `FASTX Toolkit` with default parameters to split a fastq file into multiple fastq files based on the sequences at the begining of line (bol) and end of line (eol). Merge the paired-end fastq files before using this script (alternatively use [ATG_fastq_merger](https://github.com/Ghahfarokhi/ATG_fastq_merger), which is a python script that uses `FLASH` for merging illumina NGS R1 and R2). 


## Requirements
* `fastx_toolkit` package: can be installed via [conda](https://anaconda.org/bioconda/fastx_toolkit), **note**: `fastx_barcode_splitter.pl` is expected to be in the path after the installation of fastx toolkit.
* `python3`
* Tested on a linux machine, no reason to not to work on Mac and Windows.

## Usage 

`python atg_fastq_demultiplexer.py [--fastq-dir PATH_TO_FASTQ_DIR] [--barcodes-file PATH_TO_BARCODES_FILE] [--out-dir PATH_TO_OUT_DIR]`

or

`python atg_fastq_demultiplexer.py [-i PATH_TO_FASTQ_DIR] [-f PATH_TO_BARCODES_FILE] [-o PATH_TO_OUT_DIR]`

Prepare a tab delimited `barcodes_info.tsv` file with the following four column names: `sample_name`, `barcode_bol`, `barcode_eol`, and `demuxed_name`. For example:

| sample_name | barcode_bol | barcode_eol | demuxed_name |
| ----------- | ----------- | ----------- | ------------ |
| mysample_1  | CAGTT | AACTG | mysample_1_CAGTT_AACTG_TRT_x_Rep_1 |
| mysample_1  | CAGTT | TCTCT | mysample_1_CAGTT_TCTCT_TRT_x_Rep_2 |
| mysample_1  | CAGTT | GATCA | mysample_1_CAGTT_GATCA_TRT_x_Rep_3 |
| mysample_1  | CAGTT | TTCGG | mysample_1_CAGTT_TTCGG_TRT_y_Rep_1 |
| mysample_1  | AGAGA | AACTG | mysample_1_AGAGA_AACTG_TRT_y_Rep_2 |
| mysample_1  | AGAGA | TCTCT | mysample_1_AGAGA_TCTCT_TRT_y_Rep_3 |
| mysample_1  | AGAGA | GATCA | mysample_1_AGAGA_GATCA_TRT_z_Rep_1 |
| mysample_1  | AGAGA | TTCGG | mysample_1_AGAGA_TTCGG_TRT_z_Rep_2 |
| mysample_1  | TGATC | AACTG | mysample_1_TGATC_AACTG_TRT_z_Rep_3 |

bol = begining of line (i.e., read), eol = end of line (i.e., read). 

**Important Note 1**: `sample_name` must be alphanumeric (*i.e.*, avoid &, $, @, -, %, * and empty space in the file names).

**Important Note 2**: Fastq files in the `FASTQ_DIR` are expected to have `.fastq` suffix and their names match the "sample_name" in BARCODES_FILE.


## Outputs
* `out_dir/demux_name.fastq` contains the barcode splitted fastq file(s). Note that barcodes from the begining and of reads are **not** trimmed. 
* `out_dir/fastx_splitter.log` contains the fastx statistics output (e.g., number of reads per splitted files and number of unmatched reads). 

### References
* **Fastx toolkit**: http://hannonlab.cshl.edu/fastx_toolkit/commandline.html.
 
### Bugs
Please report errors/bugs to: Amir.Taheri.Ghahfarokhi@gmail.com
