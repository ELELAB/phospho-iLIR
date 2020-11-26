#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
#
#    process_netphos_output.py
#
#    Write a .csv file with data parsed from a NetPhos 3.1 output file.
#    The script can also perform some filtering on the data parsed
#    before writing the .csv file.
#
#    Copyright (C) 2020 Valentina Sora <sora.valentina1@gmail.com>
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

import argparse
import pandas as pd
import re

def get_netphos_data(infile):
    """Get data from a NetPhos 3.1 output file.
    Parameters
    ----------
    infile : `str`
        NetPhos 3.1 output file.
    
    Returns
    -------
    `pandas.DataFrame`
        Dataframe with data parsed from the file.
    """
    
    with open(infile, "r") as f:
        # function to identify non-empty strings
        noempty = lambda x: x != ""
        # function to identify lines containing useful data
        isdata = \
            lambda x: x.startswith("#") \
                      and not x.startswith("#\n") \
                      and not x.startswith("#  prediction results") \
                      and not x.startswith("# Sequence") \
                      and not x.startswith("# ---")
        # create an empty list to store parsed data
        data = []
        # for each line
        for l in f:
            # ignore empty lines
            if re.match(r"^\s*$", l):
                continue
            # ignore lines not containing useful data
            if not isdata(l):
                continue
            # get data
            hashtag, sequence, resnum, restype, \
            context, score, kinase, answer = \
                filter(noempty, l.rstrip("\n").split(" "))
            # convert the residue number into an integer
            resnum = int(resnum)
            # convert the score into a float
            score = float(score)
            # append data parsed from the current line
            data.append({"sequence" : sequence, \
                         "resnum" : resnum, \
                         "restype" : restype, \
                         "context" : context, \
                         "score" : score, \
                         "kinase" : kinase, \
                         "answer" : answer})
        # create a pandas dataframe and return it
        return pd.DataFrame(data)


def filter_by_highest_score(df):
    """Keep only the highest-score prediction for
    each phosphporylation site.
    Parameters
    ----------
    df : `pandas.DataFrame`
        Input dataframe.
    Returns
    -------
    newdf : `pandas.DataFrame`
        Output dataframe.
    """

    # create a temporary column to uniquely identify each
    # phosphorylation site
    sequence_resnum = df["sequence"] + "_" + df["resnum"].astype(str)
    # create a copy of the original dataframe to modify
    newdf = df.copy(deep = True)
    # add the new column to the dataframe
    newdf["sequence_resnum"] = sequence_resnum
    # first, sort rows according to the prediction score,
    # highest scores first
    newdf = newdf.sort_values(by = "score", ascending = False)
    # drop duplicates rows according to the new column,
    # keeping only the first entry (that will be the one
    # with highest score, thanks to the sorting)
    newdf = newdf.drop_duplicates(subset = "sequence_resnum", \
                                  keep = "first")
    # re-sort the index so that rows will be sorted
    # as in the original dataframe
    newdf = newdf.sort_index()
    # drop the temporary column
    newdf = newdf.drop("sequence_resnum", axis = 1)
    # return the filtered dataframe
    return newdf


def filter_by_restypes(df, restypes, keepin):
    """Keep only the highest-score prediction for
    each phosphporylation site.
    Parameters
    ----------
    df : `pandas.DataFrame`
        Input dataframe.
    restypes : `list`
        List of residue types.
    keepin : `bool`
        Whether the residue type of the phosphorylation 
        site must kept if in the list or must be
        discarded if in the list.
    Returns
    -------
    newdf : `pandas.DataFrame`
        Output dataframe.
    """
    
    # residue types that can be phosphorylated
    RESTYPES = {"S", "T", "Y"}
    # which ones are included in the list?
    inters = RESTYPES.intersection(set(restypes))
    # keep/discard selected residue types
    tokeep = inters if keepin else RESTYPES - inters
    # filter the dataframe
    newdf = df[df["restype"].isin(tokeep)]
    # return the filtered dataframe
    return newdf


def filter_by_prediction(df, keeppositive):
    """Keep only positive/negative predictions.
    Parameters
    ----------
    df : `pandas.DataFrame`
        Input dataframe.
    keeppositive : `bool`
        Whether to keep only positive or negative
        predictions.
    Returns
    -------
    newdf : `pandas.DataFrame`
        Output dataframe.
    """

    # keep either the positive or the negative predictions
    tokeep = "YES" if keeppositive else "."
    # filter the dataframe
    newdf = df[df["answer"] == tokeep]
    # return the filtered dataframe
    return newdf


def filter_by_number(df, column, minval = None, maxval = None):
    """Filter by numerical columns (being higher/lower than
    certain thresholds).
    Parameters
    ----------
    df : `pandas.DataFrame`
        Input dataframe.
    column : `str`
        Name of the column to filter upon.
    minval : `int`, default: `None`
        Minimum value.
    maxval : `int`, default: `None`
        maximum value.
    Returns
    -------
    newdf : `pandas.DataFrame`
        Output dataframe.
    """

    # set minimum and maximum values to infinites if not specified
    minval = float("-inf") if minval is None else minval
    maxval = float("inf") if maxval is None else maxval
    # filter the dataframe
    newdf = df[((df[column] >= minval) & (df[column] <= maxval))]
    # return the filtered dataframe
    return newdf



if __name__ == "__main__":

    #---------------------- Set argument parser ----------------------#

    # create the parser
    parser = argparse.ArgumentParser()

    f_helpstr = "NetPhos 3.1 output file."
    parser.add_argument("-f", "--netphos-file", \
                        dest = "netphos_file", \
                        type = str, \
                        required = True, \
                        help = f_helpstr)

    o_helpstr = "Output dataframe."
    parser.add_argument("-o", "--output-df", \
                        dest = "output_df", \
                        type = str, \
                        required = True, \
                        help = o_helpstr)

    hs_helpstr = \
        "For each phosphorylation site, keep only the " \
        "highest-scoring prediction"
    parser.add_argument("-hs", "--only-highest-scores", \
                        dest = "only_highest_scores", \
                        action = "store_true", \
                        default = False, \
                        help = hs_helpstr)

    pp_helpstr = \
        "Keep only positive predictions."
    parser.add_argument("-pp", "--only-positive-predictions", \
                        dest = "only_positive_predictions", \
                        action = "store_true", \
                        default = False, \
                        help = pp_helpstr)

    rt_helpstr = \
        "Keep only phosphorylation sites that are of " \
        "specific residue types."
    parser.add_argument("-rt", "--restypes", \
                        dest = "restypes", \
                        type = str, \
                        nargs = "+", \
                        default = None, \
                        help = rt_helpstr)

    minr_helpstr = \
        "Keep only phosphorylation sites whose residue " \
        "number is higher than the minimum."
    parser.add_argument("-minr", "--minimum-resnum", \
                        dest = "minimum_resnum", 
                        type = int, \
                        default = None, \
                        help = minr_helpstr)

    maxr_helpstr = \
        "Keep only phosphorylation sites whose residue " \
        "number is lower than the maximum."
    parser.add_argument("-maxr", "--maximum-resnum", \
                        dest = "maximum_resnum", 
                        type = int, \
                        default = None, \
                        help = maxr_helpstr)

    mins_helpstr = \
        "Keep only predictions whose score " \
        "is higher than the minimum."
    parser.add_argument("-mins", "--minimum-score", \
                        dest = "minimum_score", 
                        type = float, \
                        default = None, \
                        help = mins_helpstr)

    maxs_helpstr = \
        "Keep only predictions whose score " \
        "is lower than the maximum."
    parser.add_argument("-maxs", "--maximum-score", \
                        dest = "maximum_score", 
                        type = float, \
                        default = None, \
                        help = maxs_helpstr)

    # parse the arguments
    args = parser.parse_args()

    #---------------------------- Get data ---------------------------#
    
    # get data from the NetPhos output file
    df = get_netphos_data(args.netphos_file)

    #-------------------------- Filter data --------------------------#

    # keep only the highest-scoring predictions
    if args.only_highest_scores:
        df = filter_by_highest_score(df)
    # keep only positive predictions
    if args.only_positive_predictions:
        df = filter_by_prediction(df = df, \
                                  keeppositive = True)
    # keep only phosphorylation sites being of certain
    # residue types
    if args.restypes:
        df = filter_by_restypes(df = df, \
                                restypes = args.restypes, \
                                keepin = True)
    # keep only phosphorylation sites having a residue
    # number lower/higher than certain thresholds
    if args.minimum_resnum or args.maximum_resnum:
        df = filter_by_number(df = df, \
                              column = "resnum", \
                              minval = args.minimum_resnum, \
                              maxval = args.maximum_resnum)
    # keep only predictions having a score lower/higher
    # than certain thresholds
    if args.minimum_score or args.maximum_score:
        df = filter_by_number(df = df, \
                              column = "score", \
                              minval = args.minimum_score, \
                              maxval = args.maximum_score)

    #--------------------------- Write data --------------------------#
    
    # write the dataframe as a .csv file
    df.to_csv(path_or_buf = args.output_df, \
              index = False, \
              float_format = "%.3f")