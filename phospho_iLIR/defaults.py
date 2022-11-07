#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-

#    defaults.py
#
#    Hard-coded values that should only be changed for development
#    purposes.
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



UNIPROT_FASTA_URL = "https://www.uniprot.org/uniprot/{:s}.fasta"


#------------------------------- iLIR --------------------------------#



# iLIR servers
ILIR_SERVERS = \
    {"cyprus" : \
        {"url" : "http://repeat.biol.ucy.ac.cy/cgi-bin/iLIR/iLIR_cgi",
         "description" : "University of Cyprus",
         "table_number" : 1},
     "warwick" : \
         {"url" : "https://ilir.warwick.ac.uk/lirpredict.php",
          "description" : "University of Warwick",
          "table_number" : 0}}

# Columns in the output CSV file created from iLIR results mapped
# to the names they have in the results when yielded by the 'warwick'
# server and they name they should have in the final CSV file
ILIR_CSV_COLS = \
    {"motif" : ("Motif", "MOTIF"),
     "start" : ("Start", "START"),
     "end" : ("End", "END"),
     "seq" : ("Pattern", "LIR sequence"),
     "evalue" : (None, "e-value"),
     "pssm_score" : ("PSSM Score", "PSSM score"),
     "similar_lirs" : (None, "Similar LIRs"),
     "in_pdb" : (None, "In PDB"),
     "anchor" : ("LIR in Anchor", "Anchor")}



#------------------------------ Spider3 ------------------------------#



SP3_OUT_SUFFIX = ".i1"



#------------------------------ PSIPRED ------------------------------#


PSI_OUT_SUFFIX = ".ss2"