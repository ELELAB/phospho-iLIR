#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-

#    process_netphos_output.py
#
#    Write a .csv file with data parsed from a NetPhos 3.1 output file.
#    The script can also perform some filtering on the data parsed
#    before writing the .csv file.
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
import re
# Third-party packages
import pandas as pd



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
    
    # Open the NetPhos 3.1 output file
    with open(infile, "r") as f:
        
        # Function to identify non-empty strings
        noempty = lambda x: x != ""
        
        # Function to identify lines containing useful data
        isdata = \
            lambda x: x.startswith("#") \
                      and not x.startswith("#\n") \
                      and not x.startswith("#  prediction results") \
                      and not x.startswith("# Sequence") \
                      and not x.startswith("# ---")
        
        # Create an empty list to store parsed data
        data = []
        
        # For each line
        for l in f:
            
            # If the line is empty
            if re.match(r"^\s*$", l):

                # Ignore it and continue
                continue
            
            # If the line does not contain useful data
            if not isdata(l):

                # Ignore it and continue
                continue
            
            # Get data
            hashtag, sequence, resnum, restype, \
            context, score, kinase, answer = \
                filter(noempty, l.rstrip("\n").split(" "))
            
            # Convert the residue number into an integer
            resnum = int(resnum)
            
            # Convert the score into a float
            score = float(score)
            
            # Append the data parsed from the current line
            # to the list
            data.append({"sequence" : sequence,
                         "resnum" : resnum,
                         "restype" : restype,
                         "context" : context,
                         "score" : score,
                         "kinase" : kinase,
                         "answer" : answer})
        
        # Create a pandas data frame and return it
        return pd.DataFrame(data)


def filter_by_highest_score(df):
    """Keep only the highest-score prediction for
    each phosphporylation site.
    
    Parameters
    ----------
    df : `pandas.DataFrame`
        Input data frame.
    
    Returns
    -------
    new_df : `pandas.DataFrame`
        Output data frame.
    """

    # Create a temporary column to uniquely identify each
    # phosphorylation site
    sequence_resnum = df["sequence"] + "_" + df["resnum"].astype(str)
    
    # Create a copy of the original data frame to modify
    new_df = df.copy(deep = True)
    
    # Add the new column to the data frame
    new_df["sequence_resnum"] = sequence_resnum
    
    # First, sort rows according to the prediction score,
    # highest scores first
    new_df = new_df.sort_values(by = "score",
                                ascending = False)
    
    # Drop duplicates rows according to the new column,
    # keeping only the first entry (that will be the one
    # with highest score, thanks to the sorting)
    new_df = new_df.drop_duplicates(subset = "sequence_resnum",
                                    keep = "first")
    
    # Re-sort the index so that rows will be sorted
    # as in the original data frame
    new_df = new_df.sort_index()
    
    # Drop the temporary column
    new_df = new_df.drop("sequence_resnum",
                         axis = 1)
    
    # Return the filtered data frame
    return new_df


def filter_by_restypes(df,
                       res_types,
                       keep_in):
    """Keep only the highest-score prediction for
    each phosphporylation site.
    
    Parameters
    ----------
    df : `pandas.DataFrame`
        Input data frame.
    
    res_types : `list`
        List of residue types.
    
    keep_in : `bool`
        Whether the residue type of the phosphorylation 
        site must kept if in the list or must be
        discarded if in the list.
    
    Returns
    -------
    new_df : `pandas.DataFrame`
        Output data frame.
    """
    
    # Residue types that can be phosphorylated
    RES_TYPES = {"S", "T", "Y"}
    
    # Which ones are included in the list?
    inters = RES_TYPES.intersection(set(res_types))
    
    # Keep/discard selected residue types
    to_keep = inters if keep_in else RES_TYPES - inters
    
    # Filter the data frame
    new_df = df[df["restype"].isin(to_keep)]
    
    # Return the filtered dataframe
    return new_df


def filter_by_prediction(df,
                         keep_positive):
    """Keep only positive/negative predictions.
    
    Parameters
    ----------
    
    df : `pandas.DataFrame`
        Input data frame.
    
    keep_positive : `bool`
        Whether to keep only positive or negative
        predictions.
    
    Returns
    -------
    new_df : `pandas.DataFrame`
        Output data frame.
    """

    # Keep either the positive or the negative predictions
    to_keep = "YES" if keep_positive else "."
    
    # Filter the data frame
    new_df = df[df["answer"] == to_keep]
    
    # Return the filtered data frame
    return new_df


def filter_by_number(df,
                     column,
                     min_val = None,
                     max_val = None):
    """Filter by numerical columns (being higher/lower than
    certain thresholds).
    
    Parameters
    ----------
    
    df : `pandas.DataFrame`
        Input data frame.
    
    column : `str`
        Name of the column to filter upon.
    
    min_val : `int`, default: `None`
        Minimum value.
    
    max_val : `int`, default: `None`
        Maximum value.
    
    Returns
    -------
    new_df : `pandas.DataFrame`
        Output data_frame.
    """

    # Set minimum and maximum values to infinites if not specified
    min_val = float("-inf") if min_val is None else min_val
    max_val = float("inf") if max_val is None else max_val
    
    # Filter the data frame
    new_df = df[((df[column] >= min_val) & (df[column] <= max_val))]
    
    # Return the filtered data frame
    return new_df



if __name__ == "__main__":



    #-------------------- Set the argument parser --------------------#



    # Create the parser
    parser = argparse.ArgumentParser()


    # Add the arguments
    f_helpstr = "NetPhos 3.1 output file."
    parser.add_argument("-f", "--netphos-file",
                        dest = "netphos_file",
                        type = str,
                        required = True,
                        help = f_helpstr)

    o_helpstr = "Output dataframe."
    parser.add_argument("-o", "--output-df",
                        dest = "output_df",
                        type = str,
                        required = True,
                        help = o_helpstr)

    hs_helpstr = \
        "For each phosphorylation site, keep only the " \
        "highest-scoring prediction"
    parser.add_argument("-hs", "--only-highest-scores",
                        dest = "only_highest_scores",
                        action = "store_true",
                        default = False,
                        help = hs_helpstr)

    pp_helpstr = \
        "Keep only positive predictions."
    parser.add_argument("-pp", "--only-positive-predictions",
                        dest = "only_positive_predictions",
                        action = "store_true",
                        default = False,
                        help = pp_helpstr)

    rt_helpstr = \
        "Keep only phosphorylation sites that are of " \
        "specific residue types."
    parser.add_argument("-rt", "--restypes",
                        dest = "restypes",
                        type = str,
                        nargs = "+",
                        default = None,
                        help = rt_helpstr)

    minr_helpstr = \
        "Keep only phosphorylation sites whose residue " \
        "number is higher than the minimum."
    parser.add_argument("-minr", "--minimum-resnum",
                        dest = "minimum_resnum", 
                        type = int,
                        default = None,
                        help = minr_helpstr)

    maxr_helpstr = \
        "Keep only phosphorylation sites whose residue " \
        "number is lower than the maximum."
    parser.add_argument("-maxr", "--maximum-resnum",
                        dest = "maximum_resnum", 
                        type = int,
                        default = None,
                        help = maxr_helpstr)

    mins_helpstr = \
        "Keep only predictions whose score " \
        "is higher than the minimum."
    parser.add_argument("-mins", "--minimum-score",
                        dest = "minimum_score", 
                        type = float,
                        default = None,
                        help = mins_helpstr)

    maxs_helpstr = \
        "Keep only predictions whose score " \
        "is lower than the maximum."
    parser.add_argument("-maxs", "--maximum-score",
                        dest = "maximum_score", 
                        type = float,
                        default = None,
                        help = maxs_helpstr)

    # Parse the arguments
    args = parser.parse_args()



    #-------------------------- Get the data -------------------------#


    
    # Get data from the NetPhos output file
    df = get_netphos_data(args.netphos_file)



    #------------------------ Filter the data ------------------------#



    # Keep only the highest-scoring predictions
    if args.only_highest_scores:
        df = filter_by_highest_score(df)
    
    # Keep only positive predictions
    if args.only_positive_predictions:
        df = filter_by_prediction(df = df,
                                  keep_positive = True)
    
    # Keep only phosphorylation sites being of certain
    # residue types
    if args.restypes is not None:
        df = filter_by_restypes(df = df,
                                res_types = args.restypes,
                                keep_in = True)
    
    # Keep only phosphorylation sites having a residue
    # number lower/higher than certain thresholds
    if args.minimum_resnum is not None \
    or args.maximum_resnum is not None:
        df = filter_by_number(df = df,
                              column = "resnum",
                              min_val = args.minimum_resnum,
                              max_val = args.maximum_resnum)
    
    # Keep only predictions having a score lower/higher
    # than certain thresholds
    if args.minimum_score is not None \
    or args.maximum_score is not None:
        df = filter_by_number(df = df,
                              column = "score",
                              min_val = args.minimum_score,
                              max_val = args.maximum_score)



    #------------------------- Write the data ------------------------#
    


    # Write the dataframe as a .csv file
    df.to_csv(path_or_buf = args.output_df,
              index = False,
              float_format = "%.3f")