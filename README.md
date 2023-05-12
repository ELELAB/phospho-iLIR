# phospho_iLIR

## Overview

`phospho_iLIR` is a Python package to run the phospho_iLIR pipeline and analyze the results.

The phospho_iLIR pipeline has been developed to predict changes in secondary structure propensity that may be induced by phosphorylation in the core of putative LIR motifs or in flanking regions.

The pipeline takes in input a list of UniProt IDs and for each ID:

* gets the corresponding protein sequence from UniProt;

* predicts LIR motifs in the sequence using `iLIR`[^kalvari2014](if a custom LIRs list has not been provided with the `--lirs` option);
* predicts possible phosphorylation sites within those motifs using `NetPhos` [^blom1999]. The motifs that can also be extended C-term or N-term by a variable number of residues to include flanking regions;
* constructs all possible LIR sequence variants with a subset of the phosphorylation sites mutated to phosphomimetic residues (to mimic a phosphorylation event);
* predicts the secondary structure propensities of the wild-type LIR sequence and of all the variants to assess possible local changes in conformational propensities upon phosphorylation. This is done using `Spider3` [^heffernan2017] and `PSIPRED` [^mcguffin2000].

## Background

The LIR motif is a short linear motif playing a crucial role in autophagy by mediating interactions between the ATG8 protein family and autophagy receptors and adaptors [^sora2020]. Phosphorylation in the regions flanking the motif has been identifed as a fine-tuning mechanism for the binding affinity and specificity of several LIR-containing proteins to the ATG8 family members.

## Requirements

The user must have these programs installed before running the phospho-iLIR scripts:

* `iLIR`
* `NetPhos` v3.1
* `Spider3`
* `PSIPRED`

The following Python requirements must also be met:

* `python` v3.8 or higher

All required Python dependencies will be installed together with phospho_iLIR, if not already present.

## Installation

To install phospho_iLIR, just download this folder, unzip it and run the following command from inside the folder:

`python setup.py install`

## Usage

### phospho-iLIR

This is the executable responsible for running the phospho-iLIR pipeline.

#### Command line

`phospho_iLIR [-h] -i IDSFILE -c CONFIGFILE [-d WORKDIR] [--lirs LIRS_FILE] [-n NPROC]`

#### Options

| Option               | Meaning                                                      |
| -------------------- | ------------------------------------------------------------ |
| `-h`, `--help`       | Show the help message and exit.                              |
| `-i`, `--idsfile`    | File containing the list of UniProt IDs.                     |
| `-c`, `--configfile` | Configuration file.                                          |
| `-d`, `--workdir`    | Working directory. The default is the current working directory. |
| `--lirs`             | Custom list of LIRs to investigate, the iLIR step on wt sequence will be skipped.      |
| `-l`, `--logfile`    | Log file. The default is phospho-_LIR.log, The log messages will be printed both to the log file and the standard output. |
| `-n`, `--nproc`      | Number of processes to use. The default is one process.      |

#### Input files

##### UniProt IDs file

A file containig a newline-separated list of UniProt IDs.

##### Configuration file

A YAML file containing the script configuration (please look at the `config.yaml` file in the `phospho_iLIR/config` directory for an example of configuration file).

##### Custom list of LIRs

The user can submit a custom list of LIRs. The file needs to be comma separated and resembles the iLIR output (with some differences).
The header (first line) contains: the line index, the position of the starting residue, the position of the ending residue, the LIR
sequence (i.e., `,START,END,LIR sequence`).
Every other line contains the LIRs info (i.e., `0,515,520,LQFLET`).
The LIR sequence starts from the two residues N-term of the core plus the four core aminoacids.  

#### Outputs

Note: the generation of the outputs for Spider3 and/or PSIPRED depends on whether in the configuration file the option for running them was turned on/off.

Suppose we have a list containing only one UniProt ID named Q0000, and we find two LIRs at positions 10-13 and 20-23 and three phosphorylation sites at position S8, S9 and Y25. In the configuration file, we decided to include as "flanking regions" two residues at each side of the LIR, and to use glutamate as phosphomimetic residue for both serine and tyrosine. Running the pipeline will generate the following directory tree:

```
Q00000
----| Q00000.fasta
----| Q00000_SLIMfast_input.csv
----| Q00000_SLIMfast_input_core.csv
----| ilir
----| netphos
----| psipred
----| spider3
----| lir_8-15
----|----| lir_8_15.fasta
----|----| lir_8_15-ilir.csv
----|----| lir_8_15-phosphosites.csv
----|----| lir_8_15-psipred.csv
----|----| lir_8_15-psipred.html
----|----| lir_8_15-spider3.csv
----|----| lir_8_15-variants.csv
----|----| lir_8_15-variants.md
----|----| var_S8E
----|----|----| var_S8E.fasta
----|----|----| ilir
----|----|----| psipred
----|----|----| spider3
----|----| var_S9E
----|----|----| var_S9E.fasta
----|----|----| ilir
----|----|----| psipred
----|----|----| spider3
----|----| var_S8E_S9E
----|----|----| var_S9E.fasta
----|----|----| ilir
----|----|----| psipred
----|----|----| spider3
----| lir_18-25
----|----| lir_18_25-phosphosites.csv
----|----| lir_18_25.fasta
----|----| lir_18_25-ilir.csv
----|----| lir_18_25-psipred.csv
----|----| lir_18_25-psipred.html
----|----| lir_18_25-spider3.csv
----|----| lir_18_25-variants.csv
----|----| lir_18_25-variants.md
----|----| var_Y25E
----|----|----| var_Y25E.fasta
----|----|----| ilir
----|----|----| psipred
----|----|----| spider3
```

The directories `ilir`, `netphos`, `psipred` and `spider3` contain the results of iLIR, NetPhos, PSIPRED and Spider3 for the willd-type sequence and each LIR variant (the names of these directories can be changed in the configuration file). 

`.fasta` files are FASTA files containing the sequence of either the full-length protein corresponding to the UniProt ID (`Q00000.fasta`) or a wild-type LIR sequence (`lir_*.fasta`) or the sequence of a full-length protein variant containing phosphomimetic mutations in a LIR motif  (`var_*.fasta`).

`_SLIMfast_input` is the csv file needed for the SLIMfast run. It is generated starting from the `lir_*-variants.csv` files. It contains the whole LIR sequence.

`_SLIMfast_input_core` is the csv file needed for the SLIMfast run. It is generated starting from the `lir_*-variants.csv` files. It contains the core LIR sequence.

`lir_*-variants.csv` and `lir_*-variants.md` are files summarizing the information about all the phosphomimetic variants generated for a specific LIR, either as a dataframe in a CSV file or as Markdown file (useful as a report).

`lir-*-ilir.csv` files are CSV files containing a dataframe summarizing the iLIR results for all the variants of a specific LIR.

`lir_*-phosphosites.csv` are CSV files containing a dataframe of all the phosphosites found by NetPhos for a specific LIR.

`lir_*-spider3.csv` files are CSV files containing a dataframe summarizing the Spider3 results for all the variants of a specific LIR.

`lir_*-psipred.csv` and `lir_*-psipred.html` are files summarizing the PSIPRED results for all the variants of a specific LIR either as a CSV file containing a dataframe or as a HTML file where the protein sequence is color-coded according to the secondary structure propensity of each residue (useful as a report).

## References

[^blom1999]: Blom, Nikolaj, Steen Gammeltoft, and Søren Brunak. "Sequence and structure-based prediction of eukaryotic protein phosphorylation sites." *Journal of molecular biology* 294.5 (1999): 1351-1362.
[^heffernan2017]: Heffernan, Rhys, et al. "Capturing non-local interactions by long short-term memory bidirectional recurrent neural networks for improving prediction of protein secondary structure, backbone angles, contact numbers and solvent accessibility." *Bioinformatics* 33.18 (2017): 2842-2849.
[^kalvari2014]: Kalvari, Ioanna, et al. "iLIR: A web resource for prediction of Atg8-family interacting proteins." *Autophagy* 10.5 (2014): 913-925.
[^mcguffin200]: McGuffin, Liam J., Kevin Bryson, and David T. Jones. "The PSIPRED protein structure prediction server." *Bioinformatics* 16.4 (2000): 404-405.
[^sora2020]: Sora, Valentina, et al. "Structure and dynamics in the ATG8 family from experimental to computational techniques." *Frontiers in Cell and Developmental Biology* 8 (2020): 420.

