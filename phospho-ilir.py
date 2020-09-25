#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-

# standard library
import argparse
import functools
import itertools
import logging
import os
import re
import subprocess
import sys
import urllib
import urllib.request
# dask
import dask
import distributed
from distributed import fire_and_forget
import distributed.utils
# others
import matplotlib.cm as cm
import matplotlib.colors as mplcolors
import matplotlib.pyplot as plt
import pandas as pd
import yaml

# to address a bug that resets the distributed.worker
# logger to WARNING level when a task is launched on
# the worker, no matter what the configuration was
def reset_worker_logger():
    """Utility function to reset a Dask logger handlers
    and level to desired values.
    """

    # new level
    NEWLEVEL = logging.INFO
    # get the logger
    logger = logging.getLogger("distributed.worker")
    # define the handlers to keep
    htokeep = [h for h in logger.handlers if type(h).__name__ == \
               distributed.utils.DequeHandler.__name__]
    # remove all the handlers
    for h in logger.handlers:
        logger.removeHandler(h)
    # add the handlers to keep
    for h in htokeep:
        # set the new level
        h.setLevel(NEWLEVEL)
        # add the handler to the logger
        logger.addHandler(h)
    # reset the logger level to the new level
    logger.setLevel(NEWLEVEL)
    # return the new logger
    return logger


############################# RUN COMMANDS ############################


# giving subprocess.DEVNULL to stdout 
# raises OSError(9, "Bad file descriptor")
def run_psipred(executable, \
                fasta, \
                wd, \
                stdout = None):
    """Run PSIPRED.
    """

    # reset the distributed.worker logger
    logger = reset_worker_logger()
    
    # make sure that the specified directory exists.
    # If not, create it.
    os.makedirs(wd, exist_ok = True)
    # set the command
    args = [executable, fasta]
    # start the process
    p = subprocess.Popen(args, cwd = wd, stdout = stdout)
    # wait for the process to complete
    return p.wait()


def run_ilir(executable, \
             fasta, \
             outhtml, \
             outcsv, \
             wd, \
             stdout = subprocess.DEVNULL):
    """Run iLIR.
    """
    
    # reset the distributed.worker logger
    logger = reset_worker_logger()

    # make sure that the specified directory exists.
    # If not, create it.
    os.makedirs(wd, exist_ok = True)
    # set the command
    args = [executable, fasta, "-l", outhtml, "-o", outcsv]
    # start the process
    p = subprocess.Popen(args, cwd = wd, stdout = stdout)
    # wait for the process to complete
    return p.wait()


def run_netphos(executable, \
                fasta, \
                wd, \
                outdat):
    """Run NetPhos 3.1.
    """

    # reset the distributed.worker logger
    logger = reset_worker_logger()

    # make sure that the specified directory exists.
    # If not, create it. 
    os.makedirs(wd, exist_ok = True)
    # set the command
    args = [executable, fasta]
    # start the process
    p = subprocess.Popen(args, stdout = open(outdat, "w"), cwd = wd)
    # wait for the process to complete
    return p.wait()


def run_spider3(executable, \
                fasta, \
                outprefix, \
                wd, \
                stdout = subprocess.DEVNULL):
    """Run Spider3.
    """
    
    # reset the distributed.worker logger
    logger = reset_worker_logger()

    # make sure that the specified directory exists.
    # If not, create it.  
    os.makedirs(wd, exist_ok = True)
    # set the command
    args = [executable, outprefix, fasta]
    # start the process
    p = subprocess.Popen(args, cwd = wd, stdout = stdout)
    # wait for the process to complete
    return p.wait()


def run_process_netphos_output(interpreter, \
                               script, \
                               npout, \
                               outcsv, \
                               wd, \
                               stdout = subprocess.DEVNULL):
    """Process the output obtained from NetPhos 3.1.
    """

    # reset the distributed.worker logger
    logger = reset_worker_logger()

    # make sure that the specified directory exists.
    # If not, create it.
    os.makedirs(wd, exist_ok = True)
    # set the command
    args = [interpreter, script, "-f", npout, "-o", outcsv]
    # start the process
    p = subprocess.Popen(args, cwd = wd, stdout = stdout)
    # wait for the process to complete
    return p.wait()


########################## READ/PROCESS DATA ##########################


def get_uniprotids(uniprotidsfile):
    """Read a list of newline-separated UniProt IDs from a file.
    """
    
    # reset the distributed.worker logger
    logger = reset_worker_logger()

    # parse the file and return the list
    with open(uniprotidsfile, "r") as f:
        return [l.strip("\n") for l in f if not re.match(r"^\s*$", l)]


def get_sequence_from_fasta(fastapath):
    """Get a protein sequence from a FASTA file.
    """
    
    # reset the distributed.worker logger
    logger = reset_worker_logger()
    
    # parse the FASTA file
    with open(fastapath, "r") as f:
        for l in f:
            # ignore empty lines and the header line
            if re.match(r"^\s*$", l) or l.startswith(">"):
                continue
            # if there are multiple sequences in the file,
            # only the first one will be returned
            return l.rstrip("\n")


def get_lirs_ilir(ilirres, \
                  fullseq):
    """Parse the results from iLIR and return the list of LIRs
    found.
    """

    # NB: cannot return a generator if used with
    # dask Client.submit() or similar
    
    # reset the distributed.worker logger
    logger = reset_worker_logger()

    # read the dataframe containing the iLIR results
    ildf = pd.read_csv(ilirres, sep = ",", index_col = 0)
    # create an empty list to store the LIRs found
    lirs = []
    # iterate over the rows of the dataframe (each
    # row is a LIR)
    for numrow, row in ildf.iterrows():
        # get the starting and ending point of the LIR
        # sequence. We need to subtract 1 to "START"
        # because Python indexes start from 0, but
        # sequence numbering starts from 1
        start, end = row["START"]-1, row["END"]
        # get the LIR sequence from the provided
        # complete sequence
        lirseq = fullseq[start:end]
        # make sure it corresponds to the one found
        # by iLIR (it assuments iLIR was run with the
        # complete FASTA sequence)
        if not lirseq == row["LIR sequence"]:
            # raise an error
            errstr = f"The LIR found in the sequence provided at " \
                     f"position {start}-{end} does not correspond " \
                     f"to the one found in the iLIR CSV file. " \
                     f"Please check both the complete protein " \
                     f"sequence and the CSV file for inconsistencies."
            raise ValueError(errstr)
        # append the LIR sequence and its starting
        # and ending points to the list of LIRs
        lirs.append((lirseq, start, end))
    # return a list of LIRs
    return lirs


def get_phosphosites_netphos(netphosres):
    """Parse the results from NetPhos 3.1 and
    return a set of sequence positions predicted
    to be phosphorylation sites.
    """

    # NB: cannot return a generator if used with
    # dask Client.submit() or similar
    
    # reset the distributed.worker logger
    logger = reset_worker_logger()
    
    # read the CSV file as a dataframe
    npdf = pd.read_csv(netphosres, sep = ",")
    # get the phosphorylation sites, use sets because
    # lookup is faster and we do not need them in
    # order here
    return set(npdf["resnum"])


def get_extended_lir(lir, \
                     fullseq, \
                     lcontext, \
                     rcontext):
    
    """Include a variable length residue context
    into the original LIR sequence.
    """
    
    # WARNING: fullseq is assumed to start from 1
    # since it should come from UniProt

    # reset the distributed.worker logger
    logger = reset_worker_logger()
    
    # get the LIR sequence and its starting and
    # ending point
    seq, start, end = lir
    # compute the starting and ending point of the LIR
    # extended sequence, given a number of context
    # residues on both sides
    extstart = start - lcontext
    extend = end + rcontext
    # if the number of residues on the left goes
    # beyond the beginning of the full sequence
    if extstart < 0:
        # the starting point will be the
        # beginning of the full sequence
        extstart = 0
    # if the number of residues on the right goes
    # beyond the end of the full sequence
    if extend > len(fullseq)-1:
        # the ending point will be the
        # end of the full sequence
        extend = len(fullseq)-1
    # get the extended LIR sequence
    extseq = fullseq[extstart:extend]
    # log information about the extended LIR sequence
    # (to check that it was built correctly)
    logger.info(f"Original LIR sequence is " \
                f"{start}-{seq}-{end-1}")
    logger.info(f"Extended LIR sequence is " \
                f"{extstart}-{extseq}-{extend-1}")
    # set the extended LIR name
    extname = f"lir_{extstart}_{extend-1}"
    # return the extended LIR
    return (extname, extseq, extstart, extend)


def get_lir_phosphosites(lir, psites):
    """Get the phosphorylation sites found in a LIR,
    given the LIR and a set of possible phosphosites.
    """

    # NB: cannot return a generator if used with
    # dask Client.submit() or similar

    # reset the distributed.worker logger
    logger = reset_worker_logger()
    
    # get the LIR name, sequence, starting and ending point
    name, seq, start, end = lir
    # get the sequence range of the LIR
    seqrange = range(start, end)
    # get the LIR phosphorylation sites looking up
    # the set of phosphorylation sites provided
    lirpsites = [(i,ps,res) for i, (ps,res) in \
                 enumerate(zip(seqrange, seq)) \
                 if ps+1 in set(psites)]
    # if no phosphorylation sites were found in the LIR
    if not lirpsites:
        logger.info("No phosphorylation sites found in " \
                    "the extendend LIR sequence.")
    else:
        # log information about the phosphorylation sites
        logger.info(f"Found {len(lirpsites)} phosphorylation " \
                    f"sites in {seq}.")
        # real absolute position of a phosphosites is shifted by 1, 
        # since Python indexing starts from 0 but sequence numbering
        # starts from 1
        logstr = ", ".join(\
            [f"{ap+1} ({rt})" for rp, ap, rt in lirpsites])
        logger.info(f"Phosphorylation sites at positions: {logstr}.")
    # return the list of phosphorylation sites
    return lirpsites


def get_variants(lir, fullseq, lirpsites, pres2pmim):
    """Get all combinations of the phosphomimetic variants of
    the protein given a LIR possible phosphorylation sites.
    """

    # NB: cannot return a generator if used with
    # dask Client.submit() or similar

    # reset the distributed.worker logger
    logger = reset_worker_logger()
    
    # get the LIR name, sequence, starting and ending point
    name, seq, start, end = lir
    # get the full protein sequence before and after the LIR
    beforelir, afterlir = fullseq[:start], fullseq[end:]
    # get relative positions (= positions in the LIR),
    # absolute positions (= positions in the full sequence)
    # and residue types of the phosphorylation sites
    relpos, abspos, restypes = zip(*lirpsites)
    # add 1 to the positions since residue numbering starts
    # from 1 but Python indexing starts from 0
    resnum = [pos+1 for pos in abspos]
    # create a mapping of the relative positions to
    # the absolute positions
    rel2num = dict(zip(relpos, resnum))
    # convert the relative positions into a set
    # (faster lookup compared to a list)
    relpos = set(relpos)
    # generate a list of possible options for each position
    # of the LIR (phosphorylation sites will have two options,
    # one being the wild-type residue and one being the
    # phosphomimetic residue chosen for that residue type;
    # all other residues will have only the wild-type residue
    # as possible option) 
    options = [(res, pres2pmim[res]) if i in relpos else (res,) \
                for i, res in enumerate(seq)]
    # create an empty list to store the variants
    variants = []
    # for each variant (generated by a Cartesian product
    # over all the possible options for each LIR position)
    for varnum, var in enumerate(itertools.product(*options)):
        # convert the variant sequence from a list to a string
        varseq = "".join(var)
        # the full variant sequence will be the portion of the
        # sequence before the LIR plus the variant LIR sequence
        # plus the portion of the sequence after the LIR
        varfullseq = beforelir + varseq + afterlir
        # the variant generated first is the wild-type
        # sequence, so skip it and go the next one
        if varnum == 0:
            continue
        # create empty lists to store the list of mutations
        # (e.g. ["S3E", "T4E"]) and positions (e.g. ["3", "4"])
        # of all the phosphorylation sites mutated in the variant
        mutations = []
        positions = []
        # for each residue
        for i, varres in enumerate(varseq):
            # if the variant residue is different from the
            # corresponding one in the wild-type sequence, it
            # is a phosphorylation site that has been mutated
            if varres != seq[i]:
                # update the list of mutations and positions
                mutations.append(f"{seq[i]}{rel2num[i]}{varres}")
                positions.append(rel2num[i])
        # generate the variant name
        varname = f"var_{'_'.join(mutations)}"
        # update the list of variants
        variants.append((varname, varseq, varfullseq, start, \
                         end, mutations, positions))
    # pretty-print out the variants for debug purposes
    maxlname = max([len(n) for n in list(zip(*variants))[1]])
    for varn, var, varfull, start, end, muts, pos in variants:
        rjust = maxlname+5 - len(varn)
        logger.info(f"Variant {varn}: {start+1:>{rjust}}-{var}-{end}")
    # return the list of variants
    return variants


def aggregate_ss_results(ssres, source, groupby):
    """Aggregate secondary structure prediction results."""

    # reset the distributed.worker logger
    logger = reset_worker_logger()

    if source == "psipred":
        # PSIPRED .ss2 file columns names
        cols = ["Seq", "SS", "Coil", "Helix", "Strand"]
        comment = "#"
        sep = r"\s+"
    elif source == "spider3":
        # spider3 .i1 file columns names
        cols = ["SS", "SS8", "ASA", "Phi", "Psi", "Theta", \
                "Tau", "HSE_alpha_up", "HSE_alpha_down", "CN13"]
        comment = "#"
        sep = r"\s+"
    # create an empty dictionary to store the raw dataframes
    rawdfs = {}
    # for each variant name and corresponding result
    for varname, res in ssres.items():
        # read the results as a DataFrame
        rawdf = pd.read_csv(res, \
                            sep = sep, \
                            names = cols, \
                            comment = comment)
        # the dataframe name will be the variant name
        rawdf.name = varname
        # store the dataframe in the dictionary
        rawdfs[varname] = rawdf
    # if the results should be grouped by variant
    if groupby == "variant":
        # return the raw dataframes
        return rawdfs
    # if the results should be grouped by secondary structure
    elif groupby == "secstruc":
        # create a new dictionary of empty dataframes
        dfs = {col : pd.DataFrame() for col in cols}
        # for each raw dataframe 
        for varname, rawdf in rawdfs.items():
            # for each column in the dataframe
            for col, coldata in rawdf.iteritems():
                # append that column to the dataframe
                # collecting all columns of the same type
                dfs[col][varname] = coldata
                # rename the dataframe with the name of the column
                dfs[col].name = col
        # return the dataframes
        return dfs


############################# WRITE FILES #############################


def get_and_write_fasta(uniprotid, fastapath):
    """Write a FASTA file with the protein sequence corresponding
    to a given UniProt ID.
    """

    # reset the distributed.worker logger
    logger = reset_worker_logger()

    # get the path to the FASTA file and the name of the file
    path, fastafile = os.path.split(fastapath)
    # make sure that the path up to the directory containing
    # the file exists
    os.makedirs(path, exist_ok = True)
    # URL where to retrieve the FASTA data
    url = "https://www.uniprot.org/uniprot/{:s}.fasta"
    # open the file
    with open(fastapath, "w") as o:
        # open the URL
        response = urllib.request.urlopen(url.format(uniprotid))
        # read and decode the data
        data = response.read().decode("utf-8").split("\n")
        # write the data to the file 
        o.write(data[0] + "\n" + "".join(data[1:]))


def write_fasta(sequence, fastapath):
    """Write a FASTA file with a given protein sequence.
    """

    # reset the distributed.worker logger
    logger = reset_worker_logger() 
    
    # get the path to the FASTA file and the name of the file
    path, fastafile = os.path.split(fastapath)
    # make sure that the path up to the directory containing
    # the file exists
    os.makedirs(path, exist_ok = True)
    with open(fastapath, "w") as o:
        # get the file name without the extension
        name = fastafile.rstrip(".fasta")
        # write data to the file, using the file name as a header
        o.write(f">{name}\n{sequence}")


def write_lir_phosphosites_csv(lirpsites, outcsv):
    """Write the phosphorylation sites found in a LIR 
    to a CSV file.
    """
    
    # reset the distributed.worker logger
    logger = reset_worker_logger()

    # get relative positions (= positions in the LIR),
    # absolute positions (= positions in the full sequence)
    # and residue types of the phosphorylation sites
    relpos, abspos, restypes = zip(*lirpsites)
    # add 1 to the positions since residue numbering starts
    # from 1 but Python indexing starts from 0
    resnum = tuple([pos + 1 for pos in abspos])
    # create a dataframe with the residue numbers and
    # residue types of the phosphorylation sites
    df = pd.DataFrame({"resnum" : resnum, "restype" : restypes})
    # write the dataframe to a CSV file
    df.to_csv(outcsv, sep = ",", index = False)   


def write_variants_csv(variants, outcsv):
    """Write a CSV file containing a dataframe with
    information about the LIR variants.
    """

    # reset the distributed.worker logger
    logger = reset_worker_logger()

    # function to join list elements into a string
    list2str = lambda x: ",".join(map(str, x))
    # set the columns name
    cols = ["name", "sequence", "fullsequence", "start", \
            "end", "mutations", "positions"]
    # generate the dataframe
    df = pd.DataFrame(data = variants, columns = cols)
    # drop the fullsequence column since we are only interested
    # in the LIR portion of the variant sequence
    df = df.drop(["fullsequence"], axis = 1)
    # UniProt residue numbering starts from 1 but Python
    # indexing starts from 0 (does not affect the end
    # index because in Python indexing the end of the
    # interval is not included)
    df["start"] = df["start"] + 1
    # convert lists into strings
    df["mutations"] = df["mutations"].apply(list2str)
    df["positions"] = df["positions"].apply(list2str)
    # save the dataframe to a CSV file 
    df.to_csv(path_or_buf = outcsv, index = False)


def write_variants_markdown(variants, outmd):
    """Write a Markdown file containing a table with
    information about the LIR variants (where positions
    with mutated residues are shown in bold).
    
    TODO: write HTML file instead of Markdown.
    """

    # reset the distributed.worker logger
    logger = reset_worker_logger()

    with open(outmd, "w") as o:
        # write the header of the table
        o.write("| Name | Sequence | Start | End |\n")
        o.write("|---|---|---|---|\n")
        # for each variant
        for varn, var, varfull, start, end, muts, pos in variants:
            # make mutated positions bold
            # remeber to add 1 to start because they are Python indexes
            # and they need to be converted to residue numbers
            var = \
                "".join([f"**{c}**" if start+1+i in set(pos) else c \
                         for i, c in enumerate(var)])
            # strip consecutive "****" if consecutive
            # phosphorylation sites (does not render
            # correctly otherwise)
            var = var.replace("****", "") 
            # write a table entry for the current variant
            o.write(f"| {varn} | {var} | {start+1} | {end} |\n")


def write_ilir_csv(ilirres, \
                   lirstart, \
                   lirend, \
                   outcsv):
    
    """Write a CSV file with iLIR data for the different 
    variants of a LIR.
    """
    
    # reset the distributed.worker logger
    logger = reset_worker_logger()

    # create an empty dictionary to store the processed
    # iLIR results
    dfdict = {}
    # for each (variant name, result) pair in the dictionary
    # collecting iLIR results for all variants
    for varname, res in ilirres.items():
        # read the dataframe containing the results
        ildf = pd.read_csv(res, sep = ",", index_col = 0)
        # iterate over the rows of the dataframe (each
        # row is a LIR)
        for numrow, row in ildf.iterrows():
            # get the starting and ending point of the LIR
            # sequence. We need to subtract 1 to "START"
            # because indexes start from 0, but sequence
            # numbering starts from 1
            start, end = row["START"]-1, row["END"]
            # if the current LIR starts at the same position
            # at the LIR of interest and also ends at the same
            # position (it IS the LIR of interest)
            if start == lirstart and end == lirend:
                # add te LIR data to the dictionary
                dfdict[varname] = row
                # stop parsing the current iLIR result
                # since we have already found the LIR
                break  
    # create a new dataframe built from the dictionary
    # (will have the same columns as the iLIR output, but
    # rows will be named after the variant corresponding
    # to each iLIR result)
    df = pd.DataFrame.from_dict(dfdict, orient = "index")
    # save the dataframe to a CSV file
    df.to_csv(outcsv, sep = ",", na_rep = "NA")


def write_ss_csv(ssdfs, \
                 wtseq, \
                 start, \
                 end, \
                 outcsv):
    """Write a CSV file where rows represent the different variants
    and columns represent the secondary structure predictions for
    each residue of the LIR sequence.
    """

    # reset the distributed.worker logger
    logger = reset_worker_logger()

    # wild type LIR sequence, start and end points are needed to set
    # the column names (residue type and number) and to select only
    # the portion of the results corresponding to the LIR sequence
    
    # the column names will be residue names in the form {type}{number}
    # add 1 to the starting index since Python indexing starts from
    # 0 but residue numbering starts from 1
    columns = [f"{res}{start+1+i}" for i, res in enumerate(wtseq)]
    # the row names will be the variant names, while data will be the
    # secondary structure predictions for each position;
    # take only data corresponding to the LIR portion of the sequence
    index, data = zip(*[(n, df["SS"].tolist()[start:end]) \
                        for n, df in ssdfs.items()])
    # create the dataframe
    df = pd.DataFrame(data = data, index = index, columns = columns)
    # write the CSV file
    df.to_csv(outcsv, sep = ",")


def write_psipred_html(psipreddfs, \
                       start, \
                       end, \
                       outhtml, \
                       cmaps, \
                       chunksize):
    
    """Write an HTML file with the sequences of the variants
    color-coded according to their propensity to be in
    different secondary structures.
    """

    # reset the distributed.worker logger
    logger = reset_worker_logger()

    # get the sequence starting and ending points
    seqstart = start
    seqend = end
    # adjust the chunk size if the sequence is shorter
    if seqend-seqstart < chunksize:
        chunksize = seqend-seqstart
    # column where the sequence is stored
    seqname = "Seq"
    # colums where secondary structure propensities are stored
    ssnames = ["Coil", "Helix", "Strand"]
    # width of the space dedicated to the secondary
    # structure name (names shorter than the longest
    # one will be padded with white spaces)
    namewidth = max([len(name) for name in ssnames])+2
    # title of the HTML file
    title = "Secondary structure propensities"
    # colors will be normalized between 0 and 1 (the range of
    # secondary structure propensities)
    norm = mplcolors.Normalize(vmin = 0.0, vmax = 1.0)
    # open the HTML output file
    with open(outhtml, "w") as o:
        # write the header
        o.write("<!DOCTYPE html>\n")
        o.write("<html>\n")
        o.write(f"<head>\n<title>{title}</title>\n</head>")
        o.write("<body>\n<article>\n<header>\n")
        # define style (use a monospace font)
        stylestr = "<p style=\"font-family:'Courier'\">{:s}</p>\n"
        # define string for padding
        padstr = "&nbsp;"
        # for each dataframe
        for varname, df in psipreddfs.items():
            # select only the dataframe slice corresponding to the LIR
            df = df.iloc[start:end]
            # write the variant name as header
            o.write(f"<h1>{varname}</h1>\n")
            # create a list where the lenght of each chunk
            # of text will be stored
            lenchunks = []
            # create an empty list where each chunk
            # of text will be stored
            chunks = [[] for ssname in ssnames]
            # create an iterator over the columns
            # containing the propensities mapped
            # to the corresponding color maps
            iterdf = zip(df[ssnames].iteritems(), cmaps)
            # iterate over the columns
            for i, ((ssname, column), cmap) in enumerate(iterdf):
                # create a scalar mappable from the normalized
                # color map
                smap = cm.ScalarMappable(norm = norm, cmap = cmap)
                # convert each RGBA color of the color map to the
                # corresponding HEX code
                colors = \
                    [mplcolors.to_hex(smap.to_rgba(v)) for v in column]
                # map each residue of the sequence to the corresponding
                # color
                rescolors = zip(df[seqname], colors)
                # iterate over the residues and colors
                for j, (res, color) in enumerate(rescolors):
                    # start a new chunk of text every 'chunksize'
                    # characters
                    if j % chunksize == 0:
                        chunks[i].append("")
                        lenchunks.append(0)
                    # add the string representing the colored
                    # residue to the chunk
                    chunks[i][-1] += \
                        f'<span style="color:{color}">{res}</span>'
                    # update the counter for the length of the chunk
                    lenchunks[-1] += 1
            # for each chunk
            for chunk, lenchunk in zip(zip(*chunks), lenchunks):
                # for each secondary structure sub-chunk
                for sschunk, ssname in zip(chunk, ssnames):
                    # the ending point of the chunk is the
                    # starting point plus the chunk size
                    end = start + lenchunk-1
                    # calculate the padding
                    padding = padstr*(namewidth-len(ssname))
                    # write the formatted and stylized chunk
                    # add 1 to both start and end indexes to convert
                    # between Python indexes and residue numbering
                    chunkstr = \
                        f"{ssname}{padding}{start+1} - {sschunk} - {end+1}"
                    o.write(stylestr.format(chunkstr))
                # separate each chunk with an extra new line
                o.write("\n")
                # set the starting point of the next chunk
                start = end + 1
            # reset the starting and ending points of the sequence
            start = seqstart
            end = seqend
        # close the remaining tags
        o.write("</header>\n</article>\n</body>")  
        o.write("</html>")



if __name__ == "__main__":


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
    WD = args.workdir
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
    partsp3 = functools.partial(run_spider3, \
                                executable = SP3EXEC)
    partil = functools.partial(run_ilir, \
                               executable = ILEXEC)
    partpsi = functools.partial(run_psipred, \
                                executable = PSIEXEC)
    partnp = functools.partial(run_netphos, \
                               executable = NPEXEC)
    # use the same Python interpreter in use for the processing script
    partnpscript = functools.partial(run_process_netphos_output, \
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
        upids = client.submit(get_uniprotids, IDSFILE)

        
        #------------------------ UniProt IDs ------------------------#

        
        for upid in upids.result():
            
            # create a path for the directory corresponding
            # to the current UniProt ID
            upiddir = os.path.join(WD, upid)
            
            # write the FASTA file corresponding to the
            # UniProt sequence
            fasta = os.path.join(upiddir, upid + ".fasta")
            fastaproc = client.submit(get_and_write_fasta, \
                                      uniprotid = upid, \
                                      fastapath = fasta).result()

            # get that sequence from the FASTA file
            fullseq = client.submit(get_sequence_from_fasta, \
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
            lirs = client.submit(get_lirs_ilir, \
                                 fullseq = fullseq, \
                                 ilirres = ilcsv)
            
            # get the phosphorylation sites
            psites = client.submit(get_phosphosites_netphos, \
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
                extlir = client.submit(get_extended_lir, \
                                       lir = lir, \
                                       fullseq = fullseq, \
                                       lcontext = LCONTEXT, \
                                       rcontext = RCONTEXT).result()
                
                # get the LIR phosphorylation sites
                lirpsites = client.submit(get_lir_phosphosites, \
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
                lirfastaproc = client.submit(write_fasta, \
                                             sequence = lirseq, \
                                             fastapath = lirfasta)
                
                # get the LIR phosphomimetic variants
                variants = client.submit(get_variants, \
                                         lir = extlir, \
                                         fullseq = fullseq, \
                                         lirpsites = lirpsites, \
                                         pres2pmim = PRES2PMIM).result()
                
                # write the LIR phosphorylation sites to a CSV file
                lirpsitescsv = os.path.join(lirdir, lirname + LIRPSITESCSV)
                fire_and_forget(client.submit(write_lir_phosphosites_csv, \
                                              lirpsites = lirpsites, \
                                              outcsv = lirpsitescsv))

                # write a CSV file with all the variants
                varcsv = os.path.join(lirdir, lirname + VARCSV)
                fire_and_forget(client.submit(write_variants_csv, \
                                              variants = variants, \
                                              outcsv = varcsv))
                
                # write a Markdown file with all the variants
                varmd = os.path.join(lirdir, lirname + VARMD)
                fire_and_forget(client.submit(write_variants_markdown, \
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
                        client.submit(write_fasta, \
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
                                    aggregate_ss_results, \
                                    ssres = varsp3res, \
                                    source = "spider3", \
                                    groupby = "variant")
                    
                    # write a summary CSV file of the Spider3 results
                    outsp3csv = os.path.join(lirdir, lirname + SP3CSV)
                    fire_and_forget(client.submit(\
                                    write_ss_csv, \
                                    ssdfs = sp3ssdfs, \
                                    wtseq = lirseq, \
                                    start = lirstart, \
                                    end = lirend, \
                                    outcsv = outsp3csv))

                if PSIRUN:
                    # gather PSIPRED results for all variants
                    client.gather(varpsifutures)
                    psissdfs = client.submit(\
                                    aggregate_ss_results, \
                                    ssres = varpsires, \
                                    source = "psipred", \
                                    groupby = "variant")
                    
                    # write a summary CSV file of the PSIPRED results
                    outpsicsv = os.path.join(lirdir, lirname + PSICSV)
                    fire_and_forget(client.submit(\
                                    write_ss_csv, \
                                    ssdfs = psissdfs, \
                                    wtseq = lirseq, \
                                    start = lirstart, \
                                    end = lirend, \
                                    outcsv = outpsicsv))
                    
                    # write a summary HTML file of the PSIPRED results
                    outpsihtml = os.path.join(lirdir, lirname + PSIHTML)
                    fire_and_forget(client.submit(\
                                    write_psipred_html, \
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
                                write_ilir_csv, \
                                ilirres = varilres, \
                                lirstart = rawlirstart, \
                                lirend = rawlirend, \
                                outcsv = outilcsv))

