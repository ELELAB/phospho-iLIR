#### phospho-iLIR run with a custom LIRs list ####

# In thid folder we ran phospho-iLIR submitting it a custom
# LIRs list (lirs.csv). 
# The protein to be investigated is NDP52 (Uniprot ID: Q13137)
# and it presents a non canonical LIR without the aromatic 
# residue at position HP1 (334-LVV-336) which the iLIR software cannot detect,
# thus it is a perfect example for a run with a custom LIRs list.

######### INPUT FILES #########

# config.yaml --> configuration file with the settings for the run
		  (8 residues N-term flanking region,
		   8 residues C-term flanking region,
		   spider3 for the secondary structure,
		   Glu as phosphomimetic for Thr, Tyr and Ser)
# lirs.csv --> custom LIRs list containing the NDP52 LIR
# list.txt --> txt file containing the Uniprot IDs (in this case,
               only Q13137)

######### USAGE ##########

phosphoiLIR -i list.txt -c config.yaml --lirs

######### DIRECTORY TREE ##########

.
├── Q13137
│   ├── Q13137.fasta
│   ├── Q13137_SLIMfast_input.csv
│   ├── Q13137_SLIMfast_input_core.csv
│   ├── lir_124_144
│   │   ├── lir_124_144-ilir.csv
│   │   ├── lir_124_144-phosphosites.csv
│   │   ├── lir_124_144-spider3.csv
│   │   ├── lir_124_144-variants.csv
│   │   ├── lir_124_144-variants.md
│   │   ├── lir_124_144.fasta
│   │   ├── var_T137E
│   │   │   ├── ilir
│   │   │   │   ├── var_T137E-ilir.csv
│   │   │   │   └── var_T137E-ilir.html
│   │   │   ├── spider3
│   │   │   │   ├── var_T137E.i0c
│   │   │   │   ├── var_T137E.i0r
│   │   │   │   ├── var_T137E.i0s
│   │   │   │   ├── var_T137E.i1
│   │   │   │   ├── var_T137E.i1c
│   │   │   │   ├── var_T137E.i1r
│   │   │   │   └── var_T137E.i1s
│   │   │   └── var_T137E.fasta
│   │   ├── var_T137E_T138E
│   │   │   ├── ilir
│   │   │   │   ├── var_T137E_T138E-ilir.csv
│   │   │   │   └── var_T137E_T138E-ilir.html
│   │   │   ├── spider3
│   │   │   │   ├── var_T137E_T138E.i0c
│   │   │   │   ├── var_T137E_T138E.i0r
│   │   │   │   ├── var_T137E_T138E.i0s
│   │   │   │   ├── var_T137E_T138E.i1
│   │   │   │   ├── var_T137E_T138E.i1c
│   │   │   │   ├── var_T137E_T138E.i1r
│   │   │   │   └── var_T137E_T138E.i1s
│   │   │   └── var_T137E_T138E.fasta
│   │   └── var_T138E
│   │       ├── ilir
│   │       │   ├── var_T138E-ilir.csv
│   │       │   └── var_T138E-ilir.html
│   │       ├── spider3
│   │       │   ├── var_T138E.i0c
│   │       │   ├── var_T138E.i0r
│   │       │   ├── var_T138E.i0s
│   │       │   ├── var_T138E.i1
│   │       │   ├── var_T138E.i1c
│   │       │   ├── var_T138E.i1r
│   │       │   └── var_T138E.i1s
│   │       └── var_T138E.fasta
│   ├── netphos
│   │   ├── Q13137-netphos.dat
│   │   └── Q13137-netphos_processed.csv
│   └── spider3
│       ├── Q13137.i0c
│       ├── Q13137.i0r
│       ├── Q13137.i0s
│       ├── Q13137.i1
│       ├── Q13137.i1c
│       ├── Q13137.i1r
│       └── Q13137.i1s
├── config.yaml
├── lirs.csv
├── list.txt
├── phospho_iLIR.log
└── readme.txt

13 directories, 53 files
