#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-

#    util.py
#
#    Utility functions used by the phospho-iLIR pipeline.
#
#    Copyright (C) 2022 Valentina Sora 
#                       <sora.valentina1@gmail.com>
#                       Matteo Tiberti 
#                       <matteo.tiberti@gmail.com> 
#                       Elena Papaleo
#                       <elenap@cancer.dk>
#
#    This program is free software: you can redistribute it and/or
#    modify it under the terms of the GNU General Public License as
#    published by the Free Software Foundation, either version 3 of
#    the License, or (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public
#    License along with this program. 
#    If not, see <http://www.gnu.org/licenses/>.



# Standard library
from io import StringIO
import itertools
import os
import re
import subprocess
import urllib
import urllib.request
# Third-party packages
from bs4 import BeautifulSoup as bs
import matplotlib.cm as cm
import matplotlib.colors as mplcolors
import pandas as pd
import requests as rq
from six import text_type
import yaml
# phospho-iLIR
from .dask_patches import reset_worker_logger
from .defaults import (
    ILIR_CSV_COLS,
    ILIR_SERVERS,
    PSI_OUT_SUFFIX,
    SP3_OUT_SUFFIX,
    UNIPROT_FASTA_URL,
    )



############################# RUN COMMANDS ############################



def run_ilir(server,
             fasta,
             out_html,
             out_csv,
             wd,
             log_prefix = "",
             log_file = None):
    """Run iLIR. The code has been adapted from
    the iLIR standalone tool developed by Matteo
    Tiberti <matteo.tiberti@gmail.com>.
    """
    
    # Reset the distributed.worker logger
    logger = reset_worker_logger(log_file)

    # Make sure that the specified directory exists.
    # If not, create it.
    os.makedirs(wd, exist_ok = True)

    # If the output files have already been produced
    if os.path.exists(out_html) and os.path.exists(out_csv):

        # Inform the user that iLIR will not be run again
        warnstr = \
            f"{log_prefix}:The iLIR output files {out_html} " \
            f"and {out_csv} have been found. iLIR will " \
            f"not be run again."
        logger.warning(warnstr)

        # Return the output CSV file
        return out_csv

    # Get the text from the FASTA file
    fasta = open(fasta, "r").read()

    # Try to send the request to iLIR
    response = \
        rq.post(ILIR_SERVERS[server]["url"],
                data = {"input" : fasta, "hidden" : 1})

    # If something went wrong
    if not response.ok:

        # Warn the user
        errstr = \
            f"{log_prefix}:Something went wrong when trying to " \
            f"get the iLIR results from " \
            f"{ILIR_SERVERS[server]['url']}. Status code: " \
            f"{response.status_code}."
        logger.error(errstr)

        # Return an empty string
        return ""

    # Inform the user that everything went fine
    infostr = \
        f"{log_prefix}:The iLIR results have been successfully " \
        f"retrieved from {ILIR_SERVERS[server]['url']}."
    logger.info(infostr)

    # Open the output HTML file
    with open(out_html, "w") as oh:

        # Save the results
        oh.write(response.text)

        # Inform the user that the HTML file has been written
        infostr = \
            f"{log_prefix}:The raw iLIR results have ben written " \
            f"to {out_html}."
        logger.info(infostr)

    # Try to read in the results as a data frame
    try:
        
        soup = bs(response.text, "html5lib")
        df = \
            pd.read_html(\
                StringIO(text_type(\
                    soup.find_all("table")[\
                        ILIR_SERVERS[server]["table_number"]])),
                    flavor = "html5lib",
                    header = 0)[0]

    # If something went wrong
    except Exception as e:

        # Inform the user
        errstr = \
            f"{log_prefix}:Could not parse the results from iLIR. " \
            f"An error may have occurred in the iLIR server " \
            f"{server}: {e}"
        logger.error(errstr)

        # Return an empty string
        return ""

    # Inform the user that the output file from iLIR has been
    # successfully parsed
    infostr = \
        f"{log_prefix}:The results from iLIR (server: {server}) " \
        f"have been successfully parsed."
    logger.info(infostr)

    # Drop all NA values
    df = df.dropna(how = "all")

    # Reset the index
    df = df.reset_index(drop = True)

    # If the data came from the 'cyprus' server
    if server == "cyprus":

        # Format the column containing the PSSM score
        df_pssm = \
            df[ILIR_CSV_COLS["pssm_score"][1]].str[:-1].str.split(\
                "(", expand = True)

        # Convert the data type to an integer
        df[ILIR_CSV_COLS["pssm_score"][1]] = df_pssm[0].astype(int)

        # Fill the e-value column
        df[ILIR_CSV_COLS["evalue"][1]] = df_pssm[1]

        # Fill the 'In PDB' column
        df[ILIR_CSV_COLS["in_pdb"][1]] = \
            ["No" if pd.isnull(x) else "Yes" \
             for x in df[ILIR_CSV_COLS["in_pdb"][1]]]

        # Replace NA values with empty strings
        df = df.fillna(value = "")

    # If the data came from the 'warwick' server
    elif server == "warwick":

        # For each column
        for col, (name_in_warwick, name) in ILIR_CSV_COLS.items():

            # If the column has data
            if name_in_warwick is not None:

                # Rename the column
                df = df.rename({name_in_warwick : name},
                               axis = 1)

                # If the column is the one defining the start of
                # the motif, the end of the motif, or the PSSN
                # score
                if col in ("start", "end", "pssm_score"):

                    # Convert the data type to integer
                    df[name] = df[name].astype(int)

            # Otherwise
            else:

                # Fill it with NA
                df[name] = "NA"
    
    # Sort the columns
    df = df[[item[1] for item in ILIR_CSV_COLS.values()]]

    # Save the data frame to the output CSV file
    df.to_csv(out_csv)

    # Inform the user
    infostr = \
        f"{log_prefix}:iLIR results successfully written to " \
        f"{out_csv}."
    logger.info(infostr)

    # Return the output CSV file
    return out_csv

def custom_list(log_prefix = "",
                log_file = None,
                lirs_list = ""):

    # Reset the distributed.worker logger
    logger = reset_worker_logger(log_file)

    # Inform the user that a custom list of LIRs has been provided
    infostr = \
        f"{log_prefix}:A custom LIRs list ({lirs_list}) " \
        f"has been provided by the user. The iLIR run will be skipped."
    logger.info(infostr)

def run_netphos(executable,
                fasta,
                wd,
                out_dat,
                log_prefix = "",
                log_file = None):
    """Run NetPhos 3.1.
    """

    # Reset the distributed.worker logger
    logger = reset_worker_logger(log_file)

    # Make sure that the specified directory exists.
    # If not, create it. 
    os.makedirs(wd, exist_ok = True)

    # If the output file has already been produced
    if os.path.exists(out_dat):

        # Inform the user that NetPhos will not be run again
        warnstr = \
            f"{log_prefix}:The NetPhos output file {out_dat} " \
            f"has been found. NetPhos will not be run again."
        logger.warning(warnstr)

        # Return the output file
        return out_dat
    
    # Set the command
    args = [executable, fasta]

    # Inform the user that the process has started
    infostr = \
        f"Running NetPhos on {fasta}. The output will be " \
        f"written to {out_dat}."
    logger.info(infostr)
    
    # Start the process
    p = subprocess.Popen(args,
                         stdout = open(out_dat, "w"),
                         cwd = wd)
    
    # Wait for the process to complete
    p.wait()

    # If everything went fine
    if p.returncode == 0:

        # Inform the user
        infostr = \
            f"NetPhos ran successfully on {fasta}."
        logger.info(infostr)

    # Otherwise
    else:

        # Warn the user
        errstr = \
            f"Something went wrong when running NetPhos " \
            f"on {fasta}."
        logger.error(infostr)

        # Return an empty string
        return ""

    # Return the output file
    return out_dat


def run_spider3(executable,
                fasta,
                out_prefix,
                wd,
                stdout = subprocess.DEVNULL,
                log_prefix = "",
                log_file = None):
    """Run Spider3.
    """
    
    # Reset the distributed.worker logger
    logger = reset_worker_logger(log_file)

    # Make sure that the specified directory exists.
    # If not, create it.  
    os.makedirs(wd, exist_ok = True)

    # Get the output file
    out_file = os.path.join(wd, out_prefix + SP3_OUT_SUFFIX)

    # If the output file has already been produced
    if os.path.exists(out_file):

        # Inform the user that Spider3 will not be run again
        warnstr = \
            f"{log_prefix}:The Spider3 output file {out_file} " \
            f"has been found. Spider3 will not be run again."
        logger.warning(warnstr)

        # Return the output file
        return out_file
    
    # Set the command
    args = [executable, out_prefix, fasta]
    
    # Start the process
    p = subprocess.Popen(args,
                         cwd = wd,
                         stdout = stdout)
    
    # Wait for the process to complete
    p.wait()

    # If everything went fine
    if p.returncode == 0:

        # Inform the user
        infostr = \
            f"{log_prefix}:Spider3 ran successfully on {fasta}."
        logger.info(infostr)

    # Otherwise
    else:

        # Warn the user
        errstr = \
            f"{log_prefix}:Something went wrong when running " \
            f"Spider3 on {fasta}."
        logger.error(infostr)

        # Return an empty string
        return ""

    # Return the output file
    return out_file


def run_psipred(executable,
                fasta,
                out_prefix,
                wd,
                stdout = subprocess.DEVNULL,
                stderr = subprocess.DEVNULL,
                log_prefix = "",
                log_file = None):
    """Run PSIPRED.
    """

    # Reset the distributed.worker logger
    logger = reset_worker_logger(log_file)
    
    # Make sure that the specified directory exists.
    # If not, create it.
    os.makedirs(wd, exist_ok = True)
    
    # Set the command
    args = [executable, fasta]

    # Get the output file
    out_file = os.path.join(wd, out_prefix + PSI_OUT_SUFFIX)

    # If the output file has already been produced
    if os.path.exists(out_file):

        # Inform the user that PSIPRED will not be run again
        warnstr = \
            f"{log_prefix}:The PSIPRED output file {out_file} " \
            f"has been found. PSIPRED will not be run again."
        logger.warning(warnstr)

        # Return the output file
        return out_file

    # Inform the user that the process is about to be run
    infostr = \
        f"{log_prefix}:Running PSIPRED on {fasta}. The results " \
        f"will be written to {out_file}."
    logger.info(infostr)

    # Start the process
    p = subprocess.Popen(args,
                         cwd = wd,
                         stdout = stdout,
                         stderr = stderr)
    
    # Wait for the process to complete
    p.wait()

    # If everything went fine
    if p.returncode == 0:

        # Inform the user
        infostr = \
            f"{log_prefix}:PSIPRED ran successfully on {fasta}."
        logger.info(infostr)

    # Otherwise
    else:

        # Warn the user
        errstr = \
            f"{log_prefix}:Something went wrong when running " \
            f"PSIPRED on {fasta}."
        logger.error(infostr)

        # Return an empty string
        return ""

    # Return the output file
    return out_file


def run_process_netphos_output(interpreter,
                               script,
                               np_out,
                               out_csv,
                               wd,
                               stdout = subprocess.DEVNULL,
                               log_prefix = "",
                               log_file = None):
    """Process the output obtained from NetPhos 3.1.
    """

    # Reset the distributed.worker logger
    logger = reset_worker_logger(log_file)

    # Make sure that the specified directory exists.
    # If not, create it.
    os.makedirs(wd, exist_ok = True)

    # If the output file has already been produced
    if os.path.exists(out_csv):

        # Inform the user that the script will not be run again
        warnstr = \
            f"{log_prefix}:The {script} output file {out_csv} " \
            f"has been found. {script} will not be run again."
        logger.warning(warnstr)

        # Return the output file
        return out_csv
    
    # Set the command
    args = [interpreter, script, "-f", np_out, "-o", out_csv]
    
    # Start the process
    p = subprocess.Popen(args,
                         cwd = wd,
                         stdout = stdout)
    
    # Wait for the process to complete
    p.wait()

    # If everything went fine
    if p.returncode == 0:

        # Inform the user
        infostr = \
            f"{log_prefix}:{script} ran successfully on {np_out}."
        logger.info(infostr)

    # Otherwise
    else:

        # Warn the user
        errstr = \
            f"{log_prefix}:Something went wrong when running " \
            f"{script} on {np_out}."
        logger.error(infostr)

        # Return an empty string
        return ""

    # Return the output file
    return out_csv



########################## READ/PROCESS DATA ##########################



def get_uniprot_ids(uniprot_ids_file):
    """Read a list of newline-separated UniProt IDs from a file.
    """

    # Parse the file and return the list
    with open(uniprot_ids_file, "r") as f:
        return [l.strip("\n") for l in f if not re.match(r"^\s*$", l)]


def get_sequence_from_fasta(fasta_path,
                            log_prefix = "",
                            log_file = None):
    """Get a protein sequence from a FASTA file.
    """
    
    # Reset the distributed.worker logger
    logger = reset_worker_logger(log_file)

    # Initialize the sequence to an empty string
    sequence = ""
    
    # Open the FASTA file
    with open(fasta_path, "r") as f:
        
        # For each line
        for l in f:
            
            # If the line is empty or the header
            if re.match(r"^\s*$", l) or l.startswith(">"):
                continue
            
            # If there are multiple sequences in the file,
            # only the first one will be returned
            sequence = l.rstrip("\n")
            break

    # Inform the user that the sequence retrieval was
    # successful
    infostr = \
        f"Sequence successfully retrieved from {fasta_path}."
    logger.info(infostr)

    # Return the sequence
    return sequence


def get_lirs_ilir(ilir_res,
                  full_seq,
                  log_prefix = "",
                  log_file = None,
                  source = ""):
    """Parse the results from iLIR or custom csv and return the list of LIRs
    found.
    """

    # NB: cannot return a generator if used with
    # dask Client.submit() or similar
    
    # Reset the distributed.worker logger
    logger = reset_worker_logger(log_file)

    # Try to read the data frame containing the iLIR or the custom LIRs results
    try:

        il_df = pd.read_csv(ilir_res,
                            sep = ",",
                            index_col = 0)

    # If something went wrong
    except Exception as e:

        # Warn the user
        errstr = \
            f"{log_prefix}:Could not load the {source} results " \
            f"from {ilir_res}."
        logger.error(errstr)

        # Return an empty list
        return []

    # Inform the user that the iLIR results or the custom LIRs have been
    # successfully loaded
    infostr = \
        f"The {source} results have been successfully loaded " \
        f"from {ilir_res}."
    logger.info(infostr)
    
    # Create an empty list to store the LIRs found
    lirs = []
    
    # Iterate over the rows of the data frame (each
    # row is a LIR)
    for numrow, row in il_df.iterrows():
        
        # Get the starting and ending point of the LIR
        # sequence. We need to add 1 to "START"
        # because Python indexes start from 0, but
        # sequence numbering starts from 1, and iLIR
        # already adds two extra residues to the core LIR
        # on the left side
        start, end = row["START"]+1, row["END"]
        
        # Get the LIR sequence from the provided
        # complete sequence
        lir_seq = full_seq[start:end]
        ilir_seq = full_seq[start-2:end]
        
        # If it does not correspond to the sequence found
        # by iLIR or in the custom list (it assuments iLIR was run with the
        # complete FASTA sequence)
        if not ilir_seq == row["LIR sequence"]:
            
            # Warn the user
            errstr = \
                f"{log_prefix}:The LIR found at position " \
                f"{start-2}-{end} in the sequence does not " \
                f"correspond to the one found in the {source} results. " \
                f"Please check both the complete protein " \
                f"sequence and the {source} output CSV file " \
                f"for inconsistencies."
            logger.error(errstr)
            
            # Return an empty list
            return []
        
        # Append the LIR sequence and its starting
        # and ending points to the list of LIRs
        lirs.append((lir_seq, start, end))

    # Inform the user that the LIRs have been successfully
    # parsed
    infostr = \
        f"{log_prefix}:The LIR motifs have been successfully " \
        f"retrieved from the {source} results."
    logger.info(infostr)
    
    # Return the list of LIRs
    return lirs

def get_phosphosites_netphos(netphos_res,
                             log_prefix = "",
                             log_file = None):
    """Parse the results from NetPhos 3.1 and
    return a set of sequence positions predicted
    to be phosphorylation sites.
    """

    # NB: cannot return a generator if used with
    # dask Client.submit() or similar
    
    # Reset the distributed.worker logger
    logger = reset_worker_logger(log_file)
    
    # Try to read the NetPhos results as a data frame
    try:

        np_df = pd.read_csv(netphos_res,
                            sep = ",")

    # If something went wrong
    except Exception as e:

        # Warn the user
        errstr = \
            f"{log_prefix}:Could not parse the NetPhos " \
            f"results from {netphos_res}. Exception: {e}"
        logger.error(errstr)

        # Return an empty set
        return set()
    
    # Get the phosphorylation sites, use sets because
    # lookup is faster and we do not need them in
    # order here
    return set(np_df["resnum"])


def get_extended_lir(lir,
                     full_seq,
                     l_context,
                     r_context,
                     log_prefix = "",
                     log_file = None):
    """Include a variable length residue context
    into the original LIR sequence.
    """
    
    # WARNING: full_seq is assumed to start from 1
    # since it should come from UniProt

    # Reset the distributed.worker logger
    logger = reset_worker_logger(log_file)
    
    # Get the LIR sequence and its starting and
    # ending point
    seq, start, end = lir
    
    # Compute the starting and ending point of the LIR
    # extended sequence, given a number of context
    # residues on both sides
    ext_start = start - l_context
    ext_end = end + r_context

    # Get the starting and ending points in UniProt
    # numbering
    ext_start_uniprot = ext_start
    ext_end_uniprot = ext_end
    
    # If the number of residues on the left goes
    # beyond the beginning of the full sequence
    if ext_start < 0:
        
        # The starting point will be the
        # beginning of the full sequence
        ext_start = 0

        # The starting point in UniProt numbering
        # will be 1
        ext_start_uniprot = 1
    
    # If the number of residues on the right goes
    # beyond the end of the full sequence
    if ext_end > len(full_seq)-1:
        
        # The ending point will be the
        # end of the full sequence
        ext_end = len(full_seq)-1

        # Update the ending point in UniProt
        # numbering as well
        ext_end_uniprot = ext_end
    
    # Get the extended LIR sequence
    ext_seq = full_seq[ext_start:ext_end]
    
    # Log information about the extended LIR sequence
    # (to check that it was built correctly)
    logger.info(\
        f"{log_prefix}:The original LIR sequence is " \
        f"{start}-{seq}-{end-1}.")
    logger.info(\
        f"{log_prefix}:The extended LIR sequence is " \
        f"{ext_start_uniprot}-{ext_seq}-{ext_end_uniprot}.")
    
    # Set the extended LIR name
    ext_name = f"lir_{ext_start_uniprot}_{ext_end_uniprot}"
    
    # Return the extended LIR
    return (ext_name, ext_seq, ext_start, ext_end)


def get_lir_phosphosites(lir,
                         p_sites,
                         log_prefix = "",
                         log_file = None):
    """Get the phosphorylation sites found in a LIR,
    given the LIR and a set of possible phosphosites.
    """

    # NB: cannot return a generator if used with
    # dask Client.submit() or similar

    # Reset the distributed.worker logger
    logger = reset_worker_logger(log_file)
    
    # Get the LIR name, sequence, starting and ending point
    name, seq, start, end = lir
    
    # Get the sequence range of the LIR
    seq_range = range(start, end)

    # If no phosphorylation sites were found in the LIR
    if not p_sites:

        # Inform the user
        infostr = \
            f"{log_prefix}:No phosphorylation sites found in " \
            f"the extended LIR sequence {start+1}-{seq}-{end}."
        logger.info(infostr)

        # Return an empty list of phosphorylation sites
        # for the LIR
        return []
    
    # Otherwise
    else:

        # Get the LIR phosphorylation sites looking up
        # the set of phosphorylation sites provided
        lir_p_sites = \
            [(i, ps, res) for i, (ps, res) in \
             enumerate(zip(seq_range, seq)) \
             if ps+1 in set(p_sites)]
        
        # Log information about the phosphorylation sites
        infostr = \
            f"{log_prefix}:{len(lir_p_sites)} phosphorylation " \
            f"sites found in {start+1}-{seq}-{end}."
        logger.info(infostr)
        
        # Real absolute position of a phosphosites is shifted by 1, 
        # since Python indexing starts from 0 but sequence numbering
        # starts from 1
        logstr = ", ".join(\
            [f"{ap+1} ({rt})" for rp, ap, rt in lir_p_sites])
        logger.info(\
            f"{log_prefix}:Phosphorylation sites found at " \
            f"positions: {logstr}.")
    
    # Return the list of phosphorylation sites
    return lir_p_sites


def get_variants(up_id,
                 lir,
                 full_seq,
                 lir_p_sites,
                 pres2pmim,
                 log_prefix = "",
                 log_file = None):
    """Get all combinations of the phosphomimetic variants of
    the protein given a LIR possible phosphorylation sites.
    """

    # NB: cannot return a generator if used with
    # dask Client.submit() or similar

    # Reset the distributed.worker logger
    logger = reset_worker_logger(log_file)

    # Get the LIR name, sequence, starting and ending point
    name, seq, start, end = lir

    # If no phosphosites were found
    if not lir_p_sites:

        # Inform the user
        infostr = \
            f"{log_prefix}:Since no phosphorylation sites were " \
            f"found in {start+1}-{seq}-{end}, no variants " \
            f"will be generated for this LIR."
        logger.info(infostr)

        # Return an empty list
        return []
    
    # Get the full protein sequence before and after the LIR
    before_lir, after_lir = full_seq[:start], full_seq[end:]
    
    # Get the relative positions (= positions in the LIR), the
    # absolute positions (= positions in the full sequence)
    # and residue types of the phosphorylation sites
    rel_pos, abs_pos, res_types = zip(*lir_p_sites)
    
    # Add 1 to the positions since residue numbering starts
    # from 1 but Python indexing starts from 0
    res_num = [pos+1 for pos in abs_pos]
    
    # Create a mapping of the relative positions to
    # the residue numbers
    rel2num = dict(zip(rel_pos, res_num))
    
    # Convert the relative positions into a set
    # (faster lookup compared to a list)
    rel_pos = set(rel_pos)
    
    # Generate a list of possible options for each position
    # of the LIR (phosphorylation sites will have two options,
    # one being the wild-type residue and one being the
    # phosphomimetic residue chosen for that residue type;
    # all other residues will have only the wild-type residue
    # as possible option) 
    options = \
        [(res, pres2pmim[res]) if i in rel_pos else (res,) \
         for i, res in enumerate(seq)]
    
    # Create an empty list to store the variants
    variants = []
    
    # For each variant (generated by a Cartesian product
    # over all the possible options for each LIR position)
    for var_num, var in enumerate(itertools.product(*options)):
        
        # Convert the variant sequence from a list to a string
        var_seq = "".join(var)
        
        # The full variant sequence will be the portion of the
        # sequence before the LIR plus the variant LIR sequence
        # plus the portion of the sequence after the LIR
        var_full_seq = before_lir + var_seq + after_lir
        
        # For the first variant (wild-type sequence)
        if var_num == 0:

            # Add it to the list and continue
            variants.append(\
                (up_id, seq, full_seq, start, end, [], []))
            continue
        
        # Create an empty lists to store the list of mutations
        # (e.g. ["S3E", "T4E"]) and positions (e.g. ["3", "4"])
        # of all the phosphorylation sites mutated in the variant
        mutations = []
        positions = []
        
        # For each residue
        for i, var_res in enumerate(var_seq):
            
            # If the variant residue is different from the
            # corresponding one in the wild-type sequence, it
            # is a phosphorylation site that has been mutated
            if var_res != seq[i]:
                
                # Update the list of mutations and positions
                mutations.append(f"{seq[i]}{rel2num[i]}{var_res}")
                positions.append(rel2num[i])

        # Generate the variant name
        var_name = f"var_{'_'.join(mutations)}"
        
        # Update the list of variants
        variants.append((var_name, var_seq, var_full_seq, start,
                         end, mutations, positions))
    
    # Pretty-print out the variants for debug purposes
    max_l_name = max([len(n) for n in list(zip(*variants))[1]])
    for var_n, var, var_full, start, end, muts, pos in variants[1:]:
        r_just = max_l_name+5 - len(var_n)
        logger.info(\
            f"{log_prefix}:Variant {var_n}: " \
            f"{start+1:>{r_just}}-{var}-{end}")
    
    # Return the list of variants, including the wild-type
    return variants


def aggregate_ss_results(up_id,
                         source,
                         groupby,
                         **kwargs):
    """Aggregate secondary structure prediction results.
    """

    # Reset the distributed.worker logger
    logger = reset_worker_logger(log_file = kwargs["log_file"])

    # Get the log prefix
    log_prefix = kwargs["log_prefix"]
    
    # If the results come from PSIPRED
    if source == "PSIPRED":

        # PSIPRED .ss2 file columns names
        cols = ["Seq", "SS", "Coil", "Helix", "Strand"]

        # Symbol used for comments
        comment = "#"

        # Field separator
        sep = r"\s+"
    
    # If the results come from Spider3
    elif source == "Spider3":

        # spider3 .i1 file columns names
        cols = ["SS", "SS8", "ASA", "Phi", "Psi", "Theta",
                "Tau", "HSE_alpha_up", "HSE_alpha_down", "CN13"]
        
        # Symbol used for comments
        comment = "#"

        # Field separator
        sep = r"\s+"
    
    # Create an empty dictionary to store the raw data frames
    raw_dfs = {}

    # Get all keyword arguments representing results
    # for the variants
    ss_res = \
        {k : v for k, v in kwargs.items() \
         if k not in ("up_id", "source", "groupby",
                      "log_prefix", "log_file")}

    # Get start and end positions of the LIR in
    # the full sequence
    start, end = ss_res[up_id]["start"], ss_res[up_id]["end"]
    
    # For each variant name and corresponding result
    for var_name, res in ss_res.items():

        # Get the output file containing the results
        results = res["results"]

        # If the file does not exist
        if not os.path.exists(results):

            # Warn the user and ignore the file
            errstr = \
                f"{log_prefix}:Could not find the file " \
                f"{results} for the variant {var_name}. " \
                f"Therefore, the variant will be excluded " \
                f"from the aggregation of the {source} results."
            logger.error(errstr)

            # Continue to the next variant
            continue
        
        # Read the results as a data frame
        raw_df = pd.read_csv(res["results"],
                             sep = sep,
                             names = cols,
                             comment = comment)

        # The data frame name will be the variant name
        raw_df.name = var_name

        # If the results are that of the wild-type
        # full sequence
        if var_name == up_id:
            
            # Get only the rows corresponding to
            # the LIR results
            raw_df = raw_df[start-1:end]

        # Store the data frame in the dictionary
        raw_dfs[var_name] = raw_df
    
    # If the results should be grouped by variant
    if groupby == "variant":
        
        # The final data frames will be the raw data frames
        dfs = raw_dfs
    
    # If the results should be grouped by secondary structure
    elif groupby == "secstruc":
        
        # Create a new dictionary of empty data frames
        dfs = {col : pd.DataFrame() for col in cols}
        
        # For each raw data frame 
        for var_name, raw_df in raw_dfs.items():
            
            # For each column in the data frame
            for col, col_data in raw_df.iteritems():
                
                # Append that column to the data frame
                # collecting all columns of the same type
                dfs[col][var_name] = col_data
                
                # Rename the data frame with the name of the column
                dfs[col].name = col

    # Inform the user that the aggregation was successful
    infostr = \
        f"{log_prefix}:The secondary structure predictions " \
        f"from {source} for the variants have been successfully " \
        f"aggregated (grouped by: {groupby})."
    logger.info(infostr)
        
    # Return the data frames
    return dfs



############################# WRITE FILES #############################



def get_and_write_fasta(uniprot_id,
                        fasta_path,
                        log_prefix = "",
                        log_file = None):
    """Write a FASTA file with the protein sequence corresponding
    to a given UniProt ID.
    """

    # Reset the distributed.worker logger
    logger = reset_worker_logger(log_file)

    # Get the path to the FASTA file and the name of the file
    path, fasta_file = os.path.split(fasta_path)
    
    # Make sure that the path up to the directory containing
    # the file exists
    os.makedirs(path, exist_ok = True)
    
    # Open the file
    with open(fasta_path, "w") as o:
        
        # Get the data
        response = rq.get(UNIPROT_FASTA_URL.format(uniprot_id))

        # If something went wrong
        if not response.ok:

            # Warn the user
            errstr = \
                f"{log_prefix}:Something went wrong when trying to " \
                f"get the FASTA sequence of the protein with " \
                f"UniProt ID {uniprot_id} from " \
                f"{UNIPROT_FASTA_URL.format(uniprot_id)}. Status " \
                f"code: {response.status_code}."
            logger.error(errstr)

            # Raise en exception
            raise Exception(errstr)

        # Inform the user that everything went fine
        infostr = \
            f"{log_prefix}:The FASTA sequence of the protein with " \
            f"UniProt ID {uniprot_id} has been successfully " \
            f"retrieved from {UNIPROT_FASTA_URL.format(uniprot_id)}."
        logger.info(infostr)
        
        # Read the data
        data = response.text.split("\n")
        
        # Write the data to the file 
        o.write(data[0] + "\n" + "".join(data[1:]))

    # Return the path to the FASTA file
    return fasta_path


def write_fasta(sequence,
                fasta_path,
                log_prefix = "",
                log_file = None):
    """Write a FASTA file with a given protein sequence.
    """

    # Reset the distributed.worker logger
    logger = reset_worker_logger(log_file)
    
    # Get the path to the FASTA file and the name of the file
    path, fasta_file = os.path.split(fasta_path)
    
    # Make sure that the path up to the directory containing
    # the file exists
    os.makedirs(path, exist_ok = True)
    
    # Open the FASTA file
    with open(fasta_path, "w") as o:
        
        # Get the file name without the extension
        name = fasta_file.rstrip(".fasta")
        
        # Write data to the file, using the file name as a header
        o.write(f">{name}\n{sequence}")

    # Inform the user that the file has been successfully written
    infostr = \
        f"{log_prefix}:The FASTA file {fasta_path} has been " \
        f"successfully written."
    logger.info(infostr)

    # Return the path to the FASTA file
    return fasta_path


def write_lir_phosphosites_csv(lir_p_sites,
                               out_csv,
                               log_prefix = "",
                               log_file = None):
    """Write the phosphorylation sites found in a LIR 
    to a CSV file.
    """
    
    # Reset the distributed.worker logger
    logger = reset_worker_logger(log_file)

    # If no phosphosites were passed
    if not lir_p_sites:

        # Inform the user that no file will be written
        infostr = \
            f"{log_prefix}:Since no phosphorylation sites " \
            f"were found for the current LIR, no CSV file " \
            f"containing them will be written."
        logger.info(infostr)

        # Return
        return

    # Get the relative positions (= positions in the LIR), the
    # absolute positions (= positions in the full sequence)
    # and residue types of the phosphorylation sites
    rel_pos, abs_pos, res_types = zip(*lir_p_sites)
    
    # Add 1 to the positions since residue numbering starts
    # from 1 but Python indexing starts from 0
    res_num = tuple([pos + 1 for pos in abs_pos])
    
    # Create a data frame with the residue numbers and
    # residue types of the phosphorylation sites
    df = pd.DataFrame({"resnum" : res_num, "restype" : res_types})

    # Inform the user that the data frame has been successfully
    # created
    infostr = \
        f"{log_prefix}:The data frame containing the " \
        f"phosphorylation sites has been successfully created."
    logger.info(infostr)
    
    # Write the data frame to the output CSV file
    df.to_csv(out_csv,
              sep = ",",
              index = False)

    # Inform the user that the data frame has been successfully
    # written
    infostr = \
        f"{log_prefix}:The data frame containing the " \
        f"phosphorylation sites has been successfully written " \
        f"to {out_csv}."
    logger.info(infostr)


def write_variants_csv(variants,
                       out_csv,
                       log_prefix = "",
                       log_file = None):
    """Write a CSV file containing a dataframe with
    information about the LIR variants.
    """

    # Reset the distributed.worker logger
    logger = reset_worker_logger(log_file)

    # If no variants were passed
    if not variants:

        # Inform the user that no file will be written
        infostr = \
            f"{log_prefix}:Since no variants were generated, " \
            f"no CSV file for the variants will be written."
        logger.info(infostr)

        # Return
        return

    # Function to join list elements into a string
    list2str = lambda x: ",".join(map(str, x))
    
    # Set the columns name
    cols = \
        ["name", "sequence", "fullsequence", "start",
         "end", "mutations", "positions"]
    
    # Generate the data frame
    df = pd.DataFrame(data = variants,
                      columns = cols)
    
    # Drop the fullsequence column since we are only interested
    # in the LIR portion of the variant sequence
    df = df.drop(["fullsequence"],
                 axis = 1)
    
    # The UniProt residue numbering starts from 1 but Python
    # indexing starts from 0 (does not affect the end
    # index because in Python indexing the end of the
    # interval is not included)
    df["start"] = df["start"] + 1
    
    # Convert the lists into strings
    df["mutations"] = df["mutations"].apply(list2str)
    df["positions"] = df["positions"].apply(list2str)
    
    # Save the data frame to the output CSV file 
    df.to_csv(path_or_buf = out_csv,
              index = False)

    # Inform the user that the data frame has been successfully
    # written
    infostr = \
        f"{log_prefix}:The variants have been successfully " \
        f"written to {out_csv}."
    logger.info(infostr)


def write_variants_markdown(variants,
                            out_md,
                            log_prefix = "",
                            log_file = None):
    """Write a Markdown file containing a table with
    information about the LIR variants (where positions
    with mutated residues are shown in bold).
    """

    # Reset the distributed.worker logger
    logger = reset_worker_logger(log_file)

    # If no variants were passed
    if not variants:

        # Inform the user that no file will be written
        infostr = \
            f"{log_prefix}:Since no variants were generated, " \
            f"no MarkDown file for the variants will be written."
        logger.info(infostr)

        # Return
        return

    # Open the output file
    with open(out_md, "w") as o:

        # Write the header of the table
        o.write("| Name | Sequence | Start | End |\n")
        o.write("|---|---|---|---|\n")
        
        # For each variant
        for var_n, var, var_full, start, end, muts, pos in variants:
            
            # Make mutated positions bold and remember to add 1 to the
            # starting position because they are Python indexes
            # and they need to be converted to residue numbers
            var = \
                "".join([f"**{c}**" if start+1+i in set(pos) else c \
                         for i, c in enumerate(var)])
            
            # Strip consecutive "****" if consecutive phosphorylation
            # sites (does not render correctly otherwise)
            var = var.replace("****", "") 
            
            # Write a table entry for the current variant
            o.write(f"| {var_n} | {var} | {start+1} | {end} |\n")

    # Inform the user that the data frame has been successfully
    # written
    infostr = \
        f"{log_prefix}:The variants have been successfully " \
        f"written to {out_md}."
    logger.info(infostr)


def write_ilir_csv(lir_start,
                   lir_end,
                   out_csv,
                   **kwargs):   
    """Write a CSV file with iLIR data for the different 
    variants of a LIR.
    """
    
    # Reset the distributed.worker logger
    logger = reset_worker_logger(log_file = kwargs["log_file"])

    # Get the log prefix
    log_prefix = kwargs["log_prefix"]

    # Representation of missing values
    NA_REP = "NA"
    
    # CSV file field separator
    SEP = ","

    # Create an empty dictionary to store the processed
    # iLIR results.
    df_dict = {}

    # Get all keyword arguments representing results
    # for the variants
    ilir_res = \
        {k : v for k, v in kwargs.items() \
         if k not in ("lir_start", "lir_end", 
         "out_csv", "log_prefix", "log_file")}
    
    # For each (variant name, result) pair in the dictionary
    # collecting iLIR results for all variants
    for var_name, res in ilir_res.items():

        # If the file does not exist
        if not os.path.exists(res):

            # Warn the user and ignore the file
            errstr = \
                f"{log_prefix}:Could not find the file " \
                f"{res} for the variant {var_name}. " \
                f"Therefore, the variant will be excluded " \
                f"from the aggregation of the iLIR results."
            logger.error(errstr)

            # Continue to the next variant
            continue
        
        # Read the data frame containing the results
        il_df = pd.read_csv(res,
                            sep = ",",
                            index_col = 0)
        
        # Iterate over the rows of the dataframe (each
        # row is a LIR)
        for num_row, row in il_df.iterrows():
            
            # Get the starting and ending point of the LIR
            # sequence. We need to add 1 to "START"
            # because Python indexes start from 0, but
            # sequence numbering starts from 1, and iLIR
            # already adds two extra residues to the core LIR
            # on the left side
            start, end = row["START"]+1, row["END"]
            
            # If the current LIR starts at the same position
            # at the LIR of interest and also ends at the same
            # position (it IS the LIR of interest)
            if start == lir_start and end == lir_end:
                
                # Add te LIR data to the dictionary
                df_dict[var_name] = row
                
                # Stop parsing the current iLIR result
                # since we have already found the LIR
                break
        
        # If the LIR of interest was not found in the iLIR
        # results for the current variant
        if not var_name in df_dict.keys():
            
            # Create a row filled with NA values
            empty_row = [NA_REP]*len(row)
            
            # Add a NA-only series to the dictionary
            df_dict[var_name] = pd.Series(empty_row)
    
    # Create a new data frame built from the dictionary
    # (will have the same columns as the iLIR output, but
    # rows will be named after the variant corresponding
    # to each iLIR result)
    df = pd.DataFrame.from_dict(df_dict,
                                orient = "index")

    # Inform the user that the aggregation was successful
    infostr = \
        f"{log_prefix}:The iLIR results for the variants " \
        f"have been successfully aggregated."
    logger.info(infostr)
    
    # Save the data frame to the output CSV file
    df.to_csv(out_csv,
              sep = SEP,
              na_rep = NA_REP)

    # Inform the user that the aggregated results were
    # successfully written to the output file
    infostr = \
        f"{log_prefix}:The aggregated iLIR results for the " \
        f"variants have been successfully written to {out_csv}."
    logger.info(infostr)


def write_ss_csv(ss_dfs,
                 wt_seq,
                 start,
                 end,
                 out_csv,
                 source,
                 log_prefix = "",
                 log_file = None):
    """Write a CSV file where rows represent the different variants
    and columns represent the secondary structure predictions for
    each residue of the LIR sequence.
    """

    # Reset the distributed.worker logger
    logger = reset_worker_logger(log_file)

    # Wild type LIR sequence, start and end points are needed to set
    # the column names (residue type and number) and to select only
    # the portion of the results corresponding to the LIR sequence
    
    # The column names will be residue names in the form {type}{number}
    # add 1 to the starting index since Python indexing starts from
    # 0 but residue numbering starts from 1
    columns = [f"{res}{start+1+i}" for i, res in enumerate(wt_seq)]
    
    # The row names will be the variant names, while data will be the
    # secondary structure predictions for each position;
    # take only data corresponding to the LIR portion of the sequence
    index, data = zip(*[(n, df["SS"].tolist()[start:end]) \
                        for n, df in ss_dfs.items()])
    
    # Create the data frame
    df = pd.DataFrame(data = data,
                      index = index,
                      columns = columns)
    
    # Write the output CSV file
    df.to_csv(out_csv,
              sep = ",")

    # Inform the user that the aggregated results were
    # successfully written to the output file
    infostr = \
        f"{log_prefix}:The aggregated {source} results for the " \
        f"variants have been successfully written to {out_csv}."
    logger.info(infostr)


def write_psipred_html(psipred_dfs,
                       start,
                       end,
                       out_html,
                       cmaps,
                       chunk_size,
                       log_prefix = "",
                       log_file = None):
    """Write an HTML file with the sequences of the variants
    color-coded according to their propensity to be in
    different secondary structures.
    """

    # Reset the distributed.worker logger
    logger = reset_worker_logger(log_file)

    # Get the sequence starting and ending points
    seq_start = start
    seq_end = end
    
    # Adjust the chunk size if the sequence is shorter
    if (seq_end - seq_start) < chunk_size:
        chunk_size = seq_end - seq_start
    
    # Column where the sequence is stored
    seq_name = "Seq"
    
    # Colums where secondary structure propensities are stored
    ss_names = ["Coil", "Helix", "Strand"]
    
    # Width of the space dedicated to the secondary
    # structure name (names shorter than the longest
    # one will be padded with white spaces)
    name_width = max([len(name) for name in ss_names])+2
    
    # Title of the HTML file
    title = "Secondary structure propensities"
    
    # Colors will be normalized between 0 and 1 (the range of
    # secondary structure propensities)
    norm = mplcolors.Normalize(vmin = 0.0,
                               vmax = 1.0)
    
    # Open the HTML output file
    with open(out_html, "w") as o:
        
        # Write the header
        o.write("<!DOCTYPE html>\n")
        o.write("<html>\n")
        o.write(f"<head>\n<title>{title}</title>\n</head>")
        o.write("<body>\n<article>\n<header>\n")
        
        # Define the style (use a monospace font)
        style_str = "<p style=\"font-family:'Courier'\">{:s}</p>\n"
        
        # Define the string for padding
        pad_str = "&nbsp;"
        
        # For each data frame
        for var_name, df in psipred_dfs.items():
            
            # Select only the data frame slice
            # corresponding to the LIR
            df = df.iloc[start:end]
            
            # Write the variant name as header
            o.write(f"<h1>{var_name}</h1>\n")
            
            # Create a list where the lenght of each chunk
            # of text will be stored
            len_chunks = []
            
            # Create an empty list where each chunk
            # of text will be stored
            chunks = [[] for ss_name in ss_names]
            
            # Create an iterator over the columns
            # containing the propensities mapped
            # to the corresponding color maps
            iter_df = zip(df[ss_names].iteritems(), cmaps)
            
            # Iterate over the columns
            for i, ((ss_name, column), cmap) in enumerate(iter_df):

                # Create a scalar mappable from the normalized
                # color map
                smap = cm.ScalarMappable(norm = norm,
                                         cmap = cmap)
                
                # Convert each RGBA color of the color map to the
                # corresponding HEX code
                colors = \
                    [mplcolors.to_hex(smap.to_rgba(v)) for v in column]
                
                # Map each residue of the sequence to the corresponding
                # color
                res_colors = zip(df[seq_name], colors)
                
                # Iterate over the residues and colors
                for j, (res, color) in enumerate(res_colors):
                    
                    # Start a new chunk of text every 'chunk_size'
                    # characters
                    if j % chunk_size == 0:
                        chunks[i].append("")
                        len_chunks.append(0)
                    
                    # Add the string representing the colored
                    # residue to the chunk
                    chunks[i][-1] += \
                        f'<span style="color:{color}">{res}</span>'
                    
                    # Update the counter for the length of the chunk
                    len_chunks[-1] += 1
            
            # For each chunk
            for chunk, len_chunk in zip(zip(*chunks), len_chunks):
                
                # For each secondary structure sub-chunk
                for ss_chunk, ss_name in zip(chunk, ss_names):

                    # The ending point of the chunk is the
                    # starting point plus the chunk size
                    end = start + len_chunk-1
                    
                    # Calculate the padding
                    padding = pad_str*(name_width-len(ss_name))
                    
                    # Write the formatted and stylized chunk
                    # add 1 to both start and end indexes to convert
                    # between Python indexes and residue numbering
                    chunk_str = \
                        f"{ss_name}{padding}{start+1} - " \
                        f"{ss_chunk} - {end+1}"
                    o.write(style_str.format(chunk_str))
                
                # Separate each chunk with an extra new line
                o.write("\n")
                
                # Set the starting point of the next chunk
                start = end + 1
            
            # Reset the starting and ending points of the sequence
            start = seq_start
            end = seq_end
        
        # Close the remaining HTML tags
        o.write("</header>\n</article>\n</body>")  
        o.write("</html>")

    # Inform the user that the aggregated results were
    # successfully written to the output file
    infostr = \
        f"{log_prefix}:The aggregated PSIPRED results for the " \
        f"variants have been successfully written to {out_html}."
    logger.info(infostr)