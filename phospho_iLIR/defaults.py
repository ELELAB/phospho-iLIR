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

# Columns in the output CSV file created from iLIR results
ILIR_CSV_COLS = \
    {"motif" : "MOTIF",
     "start" : "START",
     "end" : "END",
     "seq" : "LIR sequence",
     "evalue" : "e-value",
     "pssm_score" : "PSSM score",
     "similar_lirs" : "Similar LIRs",
     "in_pdb" : "In PDB",
     "anchor" : "Anchor"}