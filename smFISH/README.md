March 7th, 2025:
The Python program "smFISH_intermolecular_pairs.py" is platform independent, and in addition to designing smFISH probes it generates an output file containing intermolecular free energy changes for duplex formation between unlike probes ("*combined_output.csv"). Please read the program's printed message for more information.


# smFISH
Design smFISH probes (take into consideration target structure)
This is a work in progress, aiming to improve the design of efficient smFISH probes for RNA targets by considering the RNA secondary structure.
REQUIREMENTS:
1. Python 3.x.x with pandas, biopython == 1.83 (or earlier)
2. RNAstructure text interface version 6.5 (earlier versions do not work! and please see installation notes for PinMol for help on how to configure this)
3. input file = "ct" format (containing MFE AND suboptimal structures)

CONSTRAINTS implemented ($ indicates thhat a constraint matches the one in the Biosearch Technologies smFISH probe design tool:
1. $smFISH probe length = 20 nucleotides
2. $smFISH probe GC content is between 45 and 60%
3. $two or more nucleotides difference between regions targeted by smFISH probes
4. target region = full length (note: CDS is more commonly used to maximize signal when more than one mRNA variant is expressed) - check carefully if the probes match the mRNA variant(s) of interest
5. it takes into account the minimum fee energy (MFE) AND sub-optimal structures predicted using an energy minimization algorithm such as RNAstructure by Dave Mathews at the University of Rochester y calculating the hybridization efficiency of each probe in presence of 10% formamide and at 37C.
6. manual probe selection may be required if CDS only is targeted or not enough probes are selected (48 recommended) when using these constraints
