#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-

#    phospho_iLIR.py
#
#    Run the phospho-iLIR pipeline to generate all possible
#    phosphomimetic variants of a LIR motif and predict changes
#    in the secondary structure possibly due to phosphorylation
#    events.
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
import argparse
import functools
import logging
import os
from pkg_resources import resource_filename, Requirement
import sys
# Third-party packages
import dask
import distributed
import matplotlib.pyplot as plt
import yaml
# phospho-iLIR
from . import util
from .defaults import (
    LOG_FILE_DEFAULT,
    PSI_OUT_SUFFIX,
    SP3_OUT_SUFFIX,
    )



def run(logger):



    ######################### ARGUMENT PARSER #########################


    
    # Create the parser
    parser = argparse.ArgumentParser()

    # Add the arguments
    i_helpstr = "File containing the list of UniProt IDs."
    parser.add_argument("-i", "--idsfile",
                        type = str,
                        required = True,
                        help = i_helpstr)

    c_helpstr = "Configuration file."
    parser.add_argument("-c", "--configfile",
                        type = str,
                        required = True,
                        help = c_helpstr)

    d_helpstr = \
        "Working directory. The default is the current working " \
        "directory."
    parser.add_argument("-d", "--workdir",
                        type = str,
                        default = os.getcwd(),
                        help = d_helpstr)

    l_helpstr = \
        f"Log file. The default is: {LOG_FILE_DEFAULT}. The log " \
        f"messages will be printed both to the log file and " \
        f"the standard output."
    parser.add_argument("-l", "--logfile",
                        type = str,
                        default = LOG_FILE_DEFAULT,
                        help = l_helpstr)

    n_helpstr = \
        "Number of processes to use. The default is one process."
    parser.add_argument("-n", "--nproc",
                        type = int,
                        default = 1,
                        help = n_helpstr)

    # Parse the arguments
    args = parser.parse_args()



    ###################### GENERAL CONFIGURATION ######################
 

    
    # Get the UniProt IDs file
    IDS_FILE = args.idsfile
    
    # Try to parse the configuration
    try:

        CONFIG = yaml.full_load(open(args.configfile, "r"))

    # If something went wrong
    except Exception as e:

        # Warn the user
        errstr = \
            f"Could not parse the configuration file " \
            f"{args.configfile}. Exception: {e}"
        logger.error(errstr)

        # Raise an exception
        raise Exception(errstr)
    
    # Get the top-level working directory
    top_wd = args.workdir
    
    # If only a directory name was passed, it will be a directory
    # created inside the current working directory
    WD = \
        os.path.abspath(top_wd) if os.path.basename(top_wd) != top_wd \
        else os.path.join(os.getcwd(), top_wd)

    # Log file
    LOG_FILE = args.logfile
    
    # Number of processes to be used when running
    NPROC = args.nproc
    
    # Configuration - LIRs
    LIR_CONFIG = CONFIG["lirs"]
    LIR_P_SITES_CSV = LIR_CONFIG["out_p_sites_suffix"] + ".csv"
    L_CONTEXT = LIR_CONFIG["l_context"]
    R_CONTEXT = LIR_CONFIG["r_context"]
    
    # Configuration - variants
    VAR_CSV = LIR_CONFIG["variants"]["out_suffix"] + ".csv"
    VAR_MD = LIR_CONFIG["variants"]["out_suffix"] + ".md"
    PRES2PMIM = LIR_CONFIG["variants"]["substitutions"]
    
    # Configuration - iLIR
    IL_CONFIG = CONFIG["ilir"]
    IL_SERVER = IL_CONFIG["server"]
    IL_DIR = IL_CONFIG["dir_name"]
    IL_HTML = IL_CONFIG["out_suffix"] + ".html"
    IL_CSV = IL_CONFIG["out_suffix"] + ".csv"
    
    # Configuration - NetPhos
    NP_CONFIG = CONFIG["netphos"]
    NP_EXEC = NP_CONFIG["executable"]
    NP_SCRIPT = resource_filename(\
                    Requirement("phospho_iLIR"),
                    "phospho_iLIR/process_netphos_output.py")
    NP_DIR = NP_CONFIG["dir_name"]
    NP_DAT = NP_CONFIG["raw_out_suffix"] + ".dat"
    NP_CSV = NP_CONFIG["proc_out_suffix"] + ".csv"
    
    # Configuration - Spider3
    SP3_CONFIG = CONFIG["spider3"]
    SP3_RUN = SP3_CONFIG["run"]
    SP3_EXEC = SP3_CONFIG["executable"]
    SP3_DIR = SP3_CONFIG["dir_name"]
    SP3_CSV = SP3_CONFIG["aggregation"]["out_suffix"] + ".csv"
    
    # Configuration - PsiPred
    PSI_CONFIG = CONFIG["psipred"]
    PSI_RUN = PSI_CONFIG["run"]
    PSI_EXEC = PSI_CONFIG["executable"]
    PSI_DIR = PSI_CONFIG["dir_name"]
    PSI_CSV = PSI_CONFIG["aggregation"]["out_suffix"] + ".csv"
    PSI_HTML = PSI_CONFIG["aggregation"]["out_suffix"] + ".html"
    CMAPS = PSI_CONFIG["aggregation"]["cmaps"]
    PSI_HTML_CMAPS = \
        [plt.get_cmap(cmap.strip()) for cmap in CMAPS.split(",")]
    PSI_HTML_CHUNK_SIZE = PSI_CONFIG["aggregation"]["chunk_size"]
    
    # Generate callables from the functions with the executables set
    part_sp3 = functools.partial(util.run_spider3,
                                 executable = SP3_EXEC)
    part_psi = functools.partial(util.run_psipred,
                                 executable = PSI_EXEC)
    part_np = functools.partial(util.run_netphos,
                                executable = NP_EXEC)
    
    # Use the same Python interpreter in use for the processing script
    part_np_script = functools.partial(util.run_process_netphos_output,
                                       interpreter = sys.executable,
                                       script = NP_SCRIPT)



    ######################### RUN THE PIPELINE ########################



    # Change the scheduler if you want to use threads or a single core
    with dask.config.set(scheduler = "processes"):
        
        # Create a local cluster with the desired number of workers
        cluster = distributed.LocalCluster(n_workers = NPROC,
                                           silence_logs = "INFO",
                                           processes = True,
                                           threads_per_worker = 1)
        
        # Set a client to submit the jobs to
        client = distributed.Client(cluster)
        
        # Try to get the UniProt IDs
        try:
            
            up_ids = util.get_uniprot_ids(IDS_FILE)

        # If something went wrong
        except Exception as e:

            # Warn the user
            errstr = \
                f"Could not get the UniProt IDs from {IDS_FILE}. " \
                f"Exception: {e}."
            logging.error(errstr)

            # Raise an exception
            raise Exception(errstr)

        # Create a list for orphan futures that need to be collected
        # before exiting
        futures = []

        
        #------------------------ UniProt IDs ------------------------#

        
        # For each UniProt ID
        for up_id in up_ids:
            
            # Set a path for the directory corresponding
            # to the current UniProt ID
            up_id_dir = os.path.join(WD, up_id)

            # Create the directory
            os.makedirs(up_id_dir, exist_ok = True)

            # Set the logging options
            log_opts = \
                {"log_prefix" : up_id,
                 "log_file" : LOG_FILE}
            
            # Write the FASTA file corresponding to the
            # UniProt sequence
            fasta_path = os.path.join(up_id_dir, up_id + ".fasta")
            fasta = \
                client.submit(util.get_and_write_fasta,
                              uniprot_id = up_id,
                              fasta_path = fasta_path,
                              **log_opts)

            # Get that sequence from the FASTA file
            full_seq = \
                client.submit(util.get_sequence_from_fasta,
                              fasta_path = fasta,
                              **log_opts)


            #----------------------- Run iLIR ------------------------#


            # Set the path to the directory where iLIR will be run
            il_dir = os.path.join(up_id_dir, IL_DIR)

            # Set the path to the output HTML file that iLIR will write
            il_html = os.path.join(il_dir, up_id + IL_HTML)

            # Set the path to the output CSV file that iLIR will write
            il_csv = os.path.join(il_dir, up_id + IL_CSV)

            # Launch iLIR
            il_csv = client.submit(util.run_ilir,
                                   server = IL_SERVER,
                                   fasta = fasta,
                                   out_html = il_html,
                                   out_csv = il_csv,
                                   wd = il_dir,
                                   **log_opts)


            #---------------------- Run NetPhos ----------------------#


            # Set the path to the directory where NetPhos will be run
            np_dir = os.path.join(up_id_dir, NP_DIR)

            # Set the path to the output .dat file that NetPhos
            # will write
            np_dat_path = os.path.join(np_dir, up_id + NP_DAT)

            # Launch NetPhos
            np_dat = client.submit(part_np,
                                   fasta = fasta,
                                   wd = np_dir,
                                   out_dat = np_dat_path,
                                   **log_opts)
            
            # Set the path to the output CSV file that the processing
            # script will write
            np_csv_path = os.path.join(np_dir, up_id + NP_CSV)

            # Process the NetPhos .dat file
            np_csv = client.submit(part_np_script,
                                   np_out = np_dat,
                                   out_csv = np_csv_path,
                                   wd = np_dir,
                                   **log_opts)


            #--------------- Get LIRs and phosphosites ---------------#
            

            # Get the LIRs
            lirs = client.submit(util.get_lirs_ilir,
                                 full_seq = full_seq,
                                 ilir_res = il_csv,
                                 **log_opts)
            
            # Get the phosphorylation sites
            p_sites = client.submit(util.get_phosphosites_netphos,
                                    netphos_res = np_csv,
                                    **log_opts)


            #---------------------- Run Spider3 ----------------------#
            

            # If Spider3 needs to be run
            if SP3_RUN:
                
                # Set the path to the directory where Spider3 will
                # be run
                sp3_dir = os.path.join(up_id_dir, SP3_DIR)

                # Launch Spider3
                futures.append(client.submit(part_sp3,
                                             fasta = fasta,
                                             out_prefix = up_id,
                                             wd = sp3_dir,
                                             **log_opts))
            
            # If PSIPRED needs to be run
            if PSI_RUN:
                
                # Set the path to the directory where PSIPRED will
                # be run
                psi_dir = os.path.join(up_id_dir, PSI_DIR)

                # Launch PSIPRED
                futures.append(client.submit(part_psi,
                                             fasta = fasta,
                                             out_prefix = up_id,
                                             wd = psi_dir,
                                             **log_opts))


            #------------------------- LIRs --------------------------#
            
            
            # For each LIR (use result() here since it is a very
            # fast calculation and it does not block the workers
            # for a long time)
            for lir in lirs.result():

                # Get the raw LIR attributes
                raw_lir_seq, raw_lir_start, raw_lir_end = lir
                
                # Get the extended LIR sequence (use result() here
                # since it is a very fast calculation and it does
                # not block the workers for a long time)
                ext_lir = \
                    client.submit(util.get_extended_lir,
                                  lir = lir,
                                  full_seq = full_seq,
                                  l_context = L_CONTEXT,
                                  r_context = R_CONTEXT).result()

                # Get the LIR name, sequence, starting and ending point
                lir_name, lir_seq, lir_start, lir_end = ext_lir

                # Set a path for the LIR directory
                lir_dir = os.path.join(up_id_dir, lir_name)

                # Create the directory
                os.makedirs(lir_dir, exist_ok = True)

                # Set the options for logging
                lir_log_opts = \
                    {"log_prefix" : f"{up_id}:{lir_name}",
                     "log_file" : LOG_FILE}

                # Set the path for the FASTA file which will contain
                # the LIR sequence
                lir_fasta_path = \
                    os.path.join(lir_dir, lir_name + ".fasta")

                # Generate a FASTA file with the LIR sequence
                lir_fasta = client.submit(util.write_fasta,
                                          sequence = lir_seq,
                                          fasta_path = lir_fasta_path,
                                          **lir_log_opts)
                
                # Get the LIR phosphorylation sites
                lir_p_sites = client.submit(util.get_lir_phosphosites,
                                            lir = ext_lir,
                                            p_sites = p_sites,
                                            **lir_log_opts)

                # Set the path to the output CSV file that will contain
                # the phosphorylation sites found in the LIR
                lir_p_sites_csv = \
                    os.path.join(lir_dir, lir_name + LIR_P_SITES_CSV)

                # Write out the file
                futures.append(\
                    client.submit(util.write_lir_phosphosites_csv,
                                  lir_p_sites = lir_p_sites,
                                  out_csv = lir_p_sites_csv,
                                  **lir_log_opts))
                
                # Get the LIR phosphomimetic variants
                variants = client.submit(util.get_variants,
                                         up_id = up_id,
                                         lir = ext_lir,
                                         full_seq = full_seq,
                                         lir_p_sites = lir_p_sites,
                                         pres2pmim = PRES2PMIM,
                                         **lir_log_opts)


                #--------------------- Variants ----------------------#


                # Set the path to the output CSV file that will contain
                # the variants
                var_csv = os.path.join(lir_dir, lir_name + VAR_CSV)

                # Write out the file
                futures.append(\
                    client.submit(util.write_variants_csv,
                                  variants = variants,
                                  out_csv = var_csv,
                                  **lir_log_opts))
                
                # Set the path to the output MarkDown file that will
                # contain the variants
                var_md = os.path.join(lir_dir, lir_name + VAR_MD)

                # Write out the file
                futures.append(\
                    client.submit(util.write_variants_markdown,
                                  variants = variants,
                                  out_md = var_md,
                                  **lir_log_opts))
                
                # Create dictionaries to store the results for
                # all variants (output files)

                # Add a dictionary for iLIR's results 
                var_il_res = {up_id : il_csv}

                # If Spider3 needs to be run
                if SP3_RUN:

                    # Add a dictionary for its results
                    var_sp3_res = \
                        {up_id : \
                            {"results" : \
                                os.path.join(sp3_dir, up_id + ".i1"),
                             "start" : raw_lir_start,
                             "end" : raw_lir_end}}

                # If PSIPRED needs to be run
                if PSI_RUN:

                    # Add a dictionary for its results
                    var_psi_res = \
                        {up_id : \
                            {"results" : \
                                os.path.join(psi_dir, up_id + ".ss2"),
                             "start" : raw_lir_start,
                             "end" : raw_lir_end}}

                # For each LIR variant (apart from the wild-type)
                # (use result() here since it is a very fast
                # calculation and it does not block the workers
                # for a long time)
                for variant in variants.result()[1:]:
                    
                    # Get the variant name, partial and full sequence,
                    # starting and ending point, mutations present and
                    # positions of such mutations in the sequence
                    var_name, var_seq, var_full_seq, var_start, \
                        var_end, var_muts, var_pos = variant
                    
                    # Set a path for the variant directory
                    var_dir = os.path.join(lir_dir, var_name)

                    # Create the directory
                    os.makedirs(var_dir, exist_ok = True)

                    # Set the options for logging
                    var_log_opts = \
                        {"log_prefix" : \
                            f"{up_id}:{lir_name}:{var_name}",
                         "log_file" : \
                            LOG_FILE}
                    
                    # Set the path to a FASTA file with the variant
                    # sequence within the context of the full UniProt
                    # sequence
                    var_fasta_path = \
                        os.path.join(var_dir, var_name + ".fasta")

                    # Write the file
                    var_fasta = \
                        client.submit(util.write_fasta,
                                      sequence = var_full_seq,
                                      fasta_path = var_fasta_path,
                                      **var_log_opts)                 
                    
                    # Set the path to the directory that will contain
                    # the results from iLIR for the current variant
                    var_il_dir = os.path.join(var_dir, IL_DIR)

                    # Set the path to the the output HTML file that
                    # will contain the results from iLIR for the
                    # current variant
                    var_il_html = \
                        os.path.join(var_il_dir, var_name + IL_HTML)

                    # Set the path to the the output CSV file that
                    # will contain the results from iLIR for the
                    # current variant
                    var_il_csv = \
                        os.path.join(var_il_dir, var_name + IL_CSV)

                    # Launch iLIR
                    var_il_res[var_name] = \
                        client.submit(util.run_ilir,
                                      server = IL_SERVER,
                                      fasta = var_fasta,
                                      out_html = var_il_html,
                                      out_csv = var_il_csv,
                                      wd = var_il_dir,
                                      **var_log_opts)

                    # If Spider3 needs to be run
                    if SP3_RUN:

                        # Set the path to the directory that will
                        # contain the results from Spider3 for
                        # the current variant
                        var_sp3_dir = os.path.join(var_dir, SP3_DIR)

                        # Launch Spider3
                        var_sp3_res[var_name] = {}
                        var_sp3_res[var_name]["results"] = \
                            client.submit(part_sp3,
                                          fasta = var_fasta,
                                          out_prefix = var_name,
                                          wd = var_sp3_dir,
                                          **var_log_opts)

                    # If PSIPRED needs to be run
                    if PSI_RUN:

                        # Set the path to the directory that will
                        # contain the results from PSIPRED for
                        # the current variant
                        var_psi_dir = os.path.join(var_dir, PSI_DIR)

                        # Launch PSIPRED
                        var_psi_res[var_name] = {}
                        var_psi_res[var_name]["results"] = \
                            client.submit(part_psi,
                                          fasta = var_fasta,
                                          out_prefix = var_name,
                                          wd = var_psi_dir,
                                          **var_log_opts)

           
                #------------------- Aggregation ---------------------#
       

                # If Spider3 was run
                if SP3_RUN:

                    # Aggregate the results
                    sp3ss_dfs = client.submit(\
                                    util.aggregate_ss_results,
                                    up_id = up_id,
                                    source = "Spider3",
                                    groupby = "variant",
                                    **{**var_sp3_res, **lir_log_opts})
                    
                    # Set the path to the output CSV file that will
                    # contain the Spider3 results for all the variants
                    out_sp3_csv = \
                        os.path.join(lir_dir, lir_name + SP3_CSV)

                    # Write the file
                    futures.append(client.submit(\
                                   util.write_ss_csv,
                                   ss_dfs = sp3ss_dfs,
                                   wt_seq = lir_seq,
                                   start = lir_start,
                                   end = lir_end,
                                   out_csv = out_sp3_csv,
                                   source = "Spider3",
                                   **lir_log_opts))

                # If PSIPRED was run
                if PSI_RUN:

                    # Aggregate the results
                    psiss_dfs = client.submit(\
                                    util.aggregate_ss_results,
                                    up_id = up_id,
                                    source = "PSIPRED",
                                    groupby = "variant",
                                    **{**var_psi_res, **lir_log_opts})
                    
                    # Set the path to the output CSV file that will
                    # contain the PSIPRED results for all the variants
                    out_psi_csv = \
                        os.path.join(lir_dir, lir_name + PSI_CSV)

                    # Write the file
                    futures.append(\
                        client.submit(\
                            util.write_ss_csv,
                            ss_dfs = psiss_dfs,
                            wt_seq = lir_seq,
                            start = lir_start,
                            end = lir_end,
                            out_csv = out_psi_csv,
                            source = "PSIPRED",
                            **lir_log_opts))
                    
                    # Set the path to the output HTML file that will
                    # contain the PSIPRED results for all the variants
                    out_psi_html = \
                        os.path.join(lir_dir, lir_name + PSI_HTML)

                    # Write the file
                    futures.append(\
                        client.submit(\
                            util.write_psipred_html,
                            psipred_dfs = psiss_dfs,
                            start = lir_start,
                            end = lir_end,
                            out_html = out_psi_html,
                            cmaps = PSI_HTML_CMAPS,
                            chunk_size = PSI_HTML_CHUNK_SIZE,
                            **lir_log_opts))

                
                # Set the path to the output CSV file that will contain
                # the iLIR results for all the variants
                out_il_csv = os.path.join(lir_dir, lir_name + IL_CSV)
                
                # Write the file
                futures.append(\
                    client.submit(\
                        util.write_ilir_csv,
                        lir_start = raw_lir_start,
                        lir_end = raw_lir_end,
                        out_csv = out_il_csv,
                        **{**var_il_res, **lir_log_opts}))


        # Gather all orphan futures still running
        client.gather(futures)



def main():

    # Configure the logging
    logging.basicConfig(level = logging.INFO)

    # Get the module logger
    logger = logging.getLogger(__name__)

    # Try running the pipeline
    try:
        
        run(logger)

    # If something went wrong
    except Exception as e:

        # Warn the user
        errstr = \
            f"Could not run phospho_iLIR. Exception: {e}"
        logging.error(errstr)

        # Exit
        sys.exit(errstr)
