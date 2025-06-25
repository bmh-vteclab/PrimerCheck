# PrimerCheck
PrimerCheck is a python script that uses BLAST+ to detect if DNA sequences bind to given forward primers, reverse primers, and probes and whether they produce a valid PCR product (≤500 bp with proper primer orientation and probe positioning).

Files:
PrimerCheck.py: python script
forward_vtx.fasta: forward primers to test PrimerCheck
reverse_vtx.fasta: reverse primers to test PrimerCheck
probe_vtx.fasta: probes to test PrimerCheck
PrimerCheck_test_sequences.fasta: sequences to test PrimerCheck
Summary_expected_test_sequences: Expected results with testing PrimerCheck with vtx primers, probes and test sequences

Installation
1. Install Python
Download from python.org.
Make sure to check "Add Python to PATH" during installation.

2. Install Required Python Packages
biopython pandas click tqdm

3. Install NCBI BLAST+
Download and install from: https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
Add the BLAST+ folder to your system PATH.

Setup for Command Line Access (Optional but Recommended)
On Windows
Make sure the folder where your PrimerCheck.py file is located is in your system PATH.
To find the path, open the folder and copy the address bar (e.g. C:\Users\YourName\scripts), then:
Press ⊞ Win, search for "Environment Variables"
Under "System Variables", find Path, click Edit → New → paste the folder path.
Click OK and restart your terminal.

On Linux
If using pip install --user, ensure this is in your shell config (~/.bashrc, ~/.zshrc, etc.):
export PATH="$HOME/.local/bin:$PATH"
Run:
source ~/.bashrc  # or ~/.zshrc

Usage
python PrimerCheck.py -f input_sequences.fasta -d forward.fasta -r reverse.fasta -p probe.fasta
All input files must be in FASTA format. BLAST databases are built automatically.

Output Files
FinalResults.csv
All successful primer+probe+reverse matches forming valid PCR products (≤500 bp).

TopHits.csv
One best result per input sequence (longest product with highest identity).

Summary.txt
Summary report of:
Total sequences processed
Sequences forming valid products
Sequences binding only primers
Sequences binding one or no primers
Product size ranges
Oversized binding distances

Notes
PCR product must form with forward and reverse primers on opposite strands and the probe in between them (++-, +--, -++ or --+ orientation).
Primers must bind to sequence with a binding length of ≥17 bp and match ≥88% identity.
Output folders are overwritten each run — back up results if needed.
