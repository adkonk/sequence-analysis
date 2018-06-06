# sequence-analysis
_ Programs to analyze protein sequences

Requires Python 3.6 and Biopython. Install Biopython using "pip install biopython" in terminal if necessary

# disorderalignment.py
_ Creates a graphic of aligned predicted disorder

Perform disorder prediction using SPOT-disorder on full dataset: http://sparks-lab.org/server/SPOT-disorder/

Save folder of predictions in the same folder as disorderalignment.py (requires unzipping twice)

Save .FASTA file of protein sequences as basename_protein_translation.py in the same folder as disorderalignment.py

Save .FA file of protein sequence MSA in the same folder as disorderalignment.py


Key: 

File naming method: programname_dataset(conservation setting)(identity/type setting)(window size)

Blue: Aligned ordered sequence

Red: Aligned disordered sequence

Dark: Poorly conserved

# distance.py
_ Calculates distance between residues of the same identity and compares the protein of interest to the full proteome

Must be in the same folder as the TAIR10 FASTA protein file

Takes about 90 seconds the first time it is run, but each time after should be much faster

