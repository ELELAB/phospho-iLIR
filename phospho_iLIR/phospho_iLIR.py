#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-

# standard library
import argparse
import functools
import logging
import os
import sys
# dask
import dask
import distributed
from distributed import fire_and_forget
# others
import matplotlib.pyplot as plt
import yaml

from . import util



def main():


    ######################### ARGUMENT PARSER #########################

    # create the parser
    parser = argparse.ArgumentParser()

    # add arguments
    i_helpstr = "File containing the list of UniProt IDs."
    parser.add_argument("-i", "--idsfile", \
                        type = str, \
                        required = True, \
                        help = i_helpstr)

    c_helpstr = "Configuration file."
    parser.add_argument("-c", "--configfile", \
                        type = str, \
                        required = True, \
                        help = c_helpstr)

    d_helpstr = \
        "Working directory. Default is the current working directory."
    parser.add_argument("-d", "--workdir", \
                        type = str, \
                        default = os.getcwd(), \
                        help = d_helpstr)

    n_helpstr = "Number of processes to use. Default is one process."
    parser.add_argument("-n", "--nproc", \
                        type = int, \
                        default = 1, \
                        help = n_helpstr)

    # parse the arguments
    args = parser.parse_args()


    ###################### LOGGING CONFIGURATION ######################
    
    logging.basicConfig(level = logging.INFO)


    ###################### GENERAL CONFIGURATION ######################
 
    # UniProt IDs file
    IDSFILE = args.idsfile
    # load and parse the configuration
    CONFIG = yaml.full_load(open(args.configfile, "r"))
    # top-level working directory
    _wd = args.workdir
    # if only a directory name was passed, it will be a directory
    # created inside the current working directory
    WD = os.path.abspath(_wd) if os.path.basename(_wd) != _wd \
         else os.path.join(os.getcwd(), _wd)
    # number of processes
    NPROC = args.nproc
    
    # LIRs
    LIRCONFIG = CONFIG["lirs"]
    LIRPSITESCSV = LIRCONFIG["outpsitessuffix"] + ".csv"
    LCONTEXT = LIRCONFIG["lcontext"]
    RCONTEXT = LIRCONFIG["rcontext"]
    # variants
    VARCSV = LIRCONFIG["variants"]["outsuffix"] + ".csv"
    VARMD = LIRCONFIG["variants"]["outsuffix"] + ".md"
    PRES2PMIM = LIRCONFIG["variants"]["substitutions"]
    # iLIR
    ILCONFIG = CONFIG["ilir"]
    ILEXEC = ILCONFIG["executable"]
    ILDIR = ILCONFIG["dirname"]
    ILHTML = ILCONFIG["outsuffix"] + ".html"
    ILCSV = ILCONFIG["outsuffix"] + ".csv"
    # NetPhos
    NPCONFIG = CONFIG["netphos"]
    NPEXEC = NPCONFIG["executable"]
    NPSCRIPT = NPCONFIG["procscript"]
    NPDIR = NPCONFIG["dirname"]
    NPDAT = NPCONFIG["rawoutsuffix"] + ".dat"
    NPCSV = NPCONFIG["procoutsuffix"] + ".csv"
    # Spider3
    SP3CONFIG = CONFIG["spider3"]
    SP3RUN = SP3CONFIG["run"]
    SP3EXEC = SP3CONFIG["executable"]
    SP3DIR = SP3CONFIG["dirname"]
    SP3CSV = SP3CONFIG["aggregation"]["outsuffix"] + ".csv"
    # PsiPred
    PSICONFIG = CONFIG["psipred"]
    PSIRUN = PSICONFIG["run"]
    PSIEXEC = PSICONFIG["executable"]
    PSIDIR = PSICONFIG["dirname"]
    PSICSV = PSICONFIG["aggregation"]["outsuffix"] + ".csv"
    PSIHTML = PSICONFIG["aggregation"]["outsuffix"] + ".html"
    CMAPS = PSICONFIG["aggregation"]["cmaps"]
    PSIHTMLCMAPS = \
        [plt.get_cmap(cmap.strip()) for cmap in CMAPS.split(",")]
    PSIHTMLCHUNKSIZE = PSICONFIG["aggregation"]["chunksize"]
    
    # generate callables from the functions with the executables set
    partsp3 = functools.partial(util.run_spider3, \
                                executable = SP3EXEC)
    partil = functools.partial(util.run_ilir, \
                               executable = ILEXEC)
    partpsi = functools.partial(util.run_psipred, \
                                executable = PSIEXEC)
    partnp = functools.partial(util.run_netphos, \
                               executable = NPEXEC)
    # use the same Python interpreter in use for the processing script
    partnpscript = functools.partial(util.run_process_netphos_output, \
                                     interpreter = sys.executable, \
                                     script = NPSCRIPT)


    ######################### RUN THE PIPELINE ########################


    # change scheduler if you want to use threads or a single core
    with dask.config.set(scheduler = "processes"):
        # create a local cluster with the desired number of workers
        cluster = distributed.LocalCluster(n_workers = NPROC, \
                                           silence_logs = "INFO", \
                                           processes = True, \
                                           threads_per_worker = 1)
        
        # set a client to submit the jobs to
        client = distributed.Client(cluster)
        
        # get the UniProt IDs 
        upids = client.submit(util.get_uniprotids, IDSFILE)

        
        #------------------------ UniProt IDs ------------------------#

        
        for upid in upids.result():
            
            # create a path for the directory corresponding
            # to the current UniProt ID
            upiddir = os.path.join(WD, upid)
            
            # write the FASTA file corresponding to the
            # UniProt sequence
            fasta = os.path.join(upiddir, upid + ".fasta")
            fastaproc = client.submit(util.get_and_write_fasta, \
                                      uniprotid = upid, \
                                      fastapath = fasta).result()

            # get that sequence from the FASTA file
            fullseq = client.submit(util.get_sequence_from_fasta, \
                                    fastapath = fasta)

            # run iLIR
            ildir = os.path.join(upiddir, ILDIR)
            ilhtml = os.path.join(ildir, upid + ILHTML)
            ilcsv = os.path.join(ildir, upid + ILCSV)
            ilproc = client.submit(partil, \
                                   fasta = fasta, \
                                   outhtml = ilhtml, \
                                   outcsv = ilcsv, \
                                   wd = ildir).result()
            
            # run NetPhos
            npdir = os.path.join(upiddir, NPDIR)
            npdat = os.path.join(npdir, upid + NPDAT)
            npproc = client.submit(partnp, \
                                   fasta = fasta, \
                                   wd = npdir, \
                                   outdat = npdat).result()
            
            # process NetPhos output
            npcsv = os.path.join(npdir, upid + NPCSV)
            npscriptproc = client.submit(partnpscript, \
                                         npout = npdat, \
                                         outcsv = npcsv, \
                                         wd = npdir).result()
            
            # get the LIRs
            lirs = client.submit(util.get_lirs_ilir, \
                                 fullseq = fullseq, \
                                 ilirres = ilcsv)
            
            # get the phosphorylation sites
            psites = client.submit(util.get_phosphosites_netphos, \
                                   netphosres = npcsv)
            
            if SP3RUN:
                # launch Spider3 and forget about it
                sp3dir = os.path.join(upiddir, SP3DIR)
                fire_and_forget(client.submit(partsp3, \
                                              fasta = fasta, \
                                              outprefix = upid, \
                                              wd = sp3dir))
            
            if PSIRUN:
                # launch spider3 and forget about it
                psidir = os.path.join(upiddir, PSIDIR)
                fire_and_forget(client.submit(partpsi, \
                                              fasta = fasta, \
                                              wd = psidir))


            #------------------------- LIRs --------------------------#
            
            
            for lir in lirs.result():

                # get the raw LIR attributes
                rawlirseq, rawlirstart, rawlirend = lir
                
                # get the extended LIR sequence
                extlir = client.submit(util.get_extended_lir, \
                                       lir = lir, \
                                       fullseq = fullseq, \
                                       lcontext = LCONTEXT, \
                                       rcontext = RCONTEXT).result()
                
                # get the LIR phosphorylation sites
                lirpsites = client.submit(util.get_lir_phosphosites, \
                                          lir = extlir, \
                                          psites = psites).result()
                
                # if no phosphorylation sites were found in the LIR,
                # go on to the next LIR
                if not lirpsites:
                    continue

                # get the LIR name, sequence, starting and ending point
                lirname, lirseq, lirstart, lirend = extlir

                # create a path for the LIR directory
                lirdir = os.path.join(upiddir, lirname)
                
                # generate a FASTA file with the LIR sequence
                lirfasta = os.path.join(lirdir, lirname + ".fasta")
                lirfastaproc = client.submit(util.write_fasta, \
                                             sequence = lirseq, \
                                             fastapath = lirfasta)
                
                # get the LIR phosphomimetic variants
                variants = client.submit(util.get_variants, \
                                         lir = extlir, \
                                         fullseq = fullseq, \
                                         lirpsites = lirpsites, \
                                         pres2pmim = PRES2PMIM).result()
                
                # write the LIR phosphorylation sites to a CSV file
                lirpsitescsv = os.path.join(lirdir, lirname + LIRPSITESCSV)
                fire_and_forget(\
                    client.submit(util.write_lir_phosphosites_csv, \
                                  lirpsites = lirpsites, \
                                  outcsv = lirpsitescsv))

                # write a CSV file with all the variants
                varcsv = os.path.join(lirdir, lirname + VARCSV)
                fire_and_forget(\
                    client.submit(util.write_variants_csv, \
                                  variants = variants, \
                                  outcsv = varcsv))
                
                # write a Markdown file with all the variants
                varmd = os.path.join(lirdir, lirname + VARMD)
                fire_and_forget(\
                    client.submit(util.write_variants_markdown, \
                                  variants = variants, \
                                  outmd = varmd))
                
                # create empty lists to store the results for
                # all variants (output files and futures)
                varilres = {}
                varsp3res = {}
                varpsires = {}
                varilfutures = []
                varsp3futures = []
                varpsifutures = []

                
                #--------------------- Variants ----------------------#
                

                for variant in variants:
                    
                    # get the variant name, sequences, starting
                    # and ending point, mutations present and
                    # positions of such mutations in the sequence
                    varname, varseq, varfullseq, varstart, \
                        varend, varmuts, varpos = variant
                    
                    # create a path for the variant directory
                    vardir = os.path.join(lirdir, varname)
                    
                    # generate a FASTA file with the variant sequence
                    # within the context of the full UniProt sequence
                    varfasta = os.path.join(vardir, varname + ".fasta")
                    varfastaproc = \
                        client.submit(util.write_fasta, \
                                      sequence = varfullseq, \
                                      fastapath = varfasta).result()                   
                    
                    # run iLIR
                    varildir = os.path.join(vardir, ILDIR)
                    varilhtml = os.path.join(varildir, varname + ILHTML)
                    varilcsv = os.path.join(varildir, varname + ILCSV)
                    varilres[varname] = varilcsv
                    varilfutures.append(\
                        client.submit(partil, \
                                      fasta = varfasta, \
                                      outhtml = varilhtml, \
                                      outcsv = varilcsv, \
                                      wd = varildir))

                    if SP3RUN:                  
                        # run Spider3
                        varsp3dir = os.path.join(vardir, SP3DIR)
                        varsp3res[varname] = \
                            os.path.join(varsp3dir, varname + ".i1")
                        varsp3futures.append(\
                            client.submit(partsp3, \
                                          fasta = varfasta, \
                                          outprefix = varname, \
                                          wd = varsp3dir))

                    if PSIRUN:
                        # run PSIPRED
                        varpsidir = os.path.join(vardir, PSIDIR)
                        varpsires[varname] = \
                            os.path.join(varpsidir, varname + ".ss2")
                        varpsifutures.append(\
                            client.submit(partpsi, \
                                          fasta = varfasta, \
                                          wd = varpsidir))

           
                #------------------- Aggregation ---------------------#
       

                if SP3RUN:
                    # gather Spider3 results for all variants
                    client.gather(varsp3futures)
                    sp3ssdfs = client.submit(\
                                    util.aggregate_ss_results, \
                                    ssres = varsp3res, \
                                    source = "spider3", \
                                    groupby = "variant")
                    
                    # write a summary CSV file of the Spider3 results
                    outsp3csv = os.path.join(lirdir, lirname + SP3CSV)
                    fire_and_forget(client.submit(\
                                    util.write_ss_csv, \
                                    ssdfs = sp3ssdfs, \
                                    wtseq = lirseq, \
                                    start = lirstart, \
                                    end = lirend, \
                                    outcsv = outsp3csv))

                if PSIRUN:
                    # gather PSIPRED results for all variants
                    client.gather(varpsifutures)
                    psissdfs = client.submit(\
                                    util.aggregate_ss_results, \
                                    ssres = varpsires, \
                                    source = "psipred", \
                                    groupby = "variant")
                    
                    # write a summary CSV file of the PSIPRED results
                    outpsicsv = os.path.join(lirdir, lirname + PSICSV)
                    fire_and_forget(client.submit(\
                                    util.write_ss_csv, \
                                    ssdfs = psissdfs, \
                                    wtseq = lirseq, \
                                    start = lirstart, \
                                    end = lirend, \
                                    outcsv = outpsicsv))
                    
                    # write a summary HTML file of the PSIPRED results
                    outpsihtml = os.path.join(lirdir, lirname + PSIHTML)
                    fire_and_forget(client.submit(\
                                    util.write_psipred_html, \
                                    psipreddfs = psissdfs, \
                                    start = lirstart, \
                                    end = lirend, \
                                    outhtml = outpsihtml, \
                                    cmaps = PSIHTMLCMAPS, \
                                    chunksize = PSIHTMLCHUNKSIZE))

                # gather iLIR results for all variants
                client.gather(varilfutures)
                # write the aggregated iLIR results for all variants
                outilcsv = os.path.join(lirdir, lirname + ILCSV)
                fire_and_forget(client.submit(\
                                util.write_ilir_csv, \
                                ilirres = varilres, \
                                lirstart = rawlirstart, \
                                lirend = rawlirend, \
                                outcsv = outilcsv))

if __name__ == "__main__":
    main()

