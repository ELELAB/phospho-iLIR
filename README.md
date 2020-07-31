# phospho-iLIR

## Overview

`phospho-iLIR` is a set of Python scripts to run the phospho-iLIR pipeline and analyze the results.

The phospho-iLIR pipeline takes in input a list of UniProt IDs and for each ID:

* gets the corresponding protein sequence;

* predicts LIR motifs in the sequence (using `iLIR`);
* predicts possible phosphorylation sites within those motifs (using `NetPhos`);
* constructs all possible LIR sequence variants with a subset of the phosphorylation sites mutated to phosphomimetic residues;
* predicts the secondary structure propensities of the wild-type LIR sequence and of all the variants to assess possible local changes in conformational propensities upon phosphorylation events (using `Spider3` and `PSIPRED`).

## Background

## Requirements

The user must have these programs installed before running the phospho-iLIR scripts:

* `iLIR`
* `NetPhos` v3.1
* `Spider3`
* `PSIPRED`

The following Python requirements must also be met:

* `python` v3.7 or higher

Finally, the following Python dependencies must be installed:

* `dask` and `dask distributed`
* `matplotlib`
* `pandas`

## Installation

The scripts require no installation.

## Usage

### phospho-ilir.py

This is the script responsible for running the phospho-iLIR pipeline.

#### Command line

`phospho-ilir.py [-h] -i IDSFILE -c CONFIGFILE [-d WORKDIR] [-n NPROC]`

#### Options

| Option               | Meaning                                                      |
| -------------------- | ------------------------------------------------------------ |
| `-h`, `--help`       | Show the help message and exit.                              |
| `-i`, `--idsfile`    | File containing the list of UniProt IDs.                     |
| `-c`, `--configfile` | Configuration file.                                          |
| `-d`, `--workdir`    | Working directory. Default is the current working directory. |
| `-n`, `--nproc`      | Number of processes to use. Default is one process.          |

#### Input files

##### UniProt IDs file

A file containig a newline-separated list of UniProt IDs.

##### Configuration file

A YAML file containing the script configuration (please see the `config.yaml` file in `examples` for an example of configuration file).