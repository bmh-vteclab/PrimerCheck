import os
import platform
import subprocess
import click
import pandas as pd
from Bio import SeqIO
from tqdm import tqdm

# Minimum allowed primer length
MIN_PRIMER_LENGTH = 17

# Get files from users
@click.command(help="PrimerCheck: Check primer/probe matches against sequences using BLAST. Primer and probe sequences should be in 3 separate fasta files.\
One file for forward primers, one file for reverse primers and one file for probes. PrimerCheck was designed for Windows and Linux platforms and has not been \
tested on others.")
@click.option('-f', '--fasta', type=click.Path(exists=True), help='Path to the fasta file.')
@click.option('-d', '--forward', type=click.Path(exists=True), help='Path to the forward primers fasta file.')
@click.option('-r', '--reverse', type=click.Path(exists=True), help='Path to the reverse primers fasta file.')
@click.option('-p', '--probe', type=click.Path(exists=True), help='Path to the probe primers fasta file.')

def main(fasta, forward, reverse, probe):
    # Determine system platform (Windows or Linux) and confirm Blast is installed
    system = platform.system().lower()
    check_blast_installed(system)

    # Put sequences to be tested in a list
    sequences = list(SeqIO.parse(fasta, 'fasta'))

    # make a Blast database for the forward primers, reverse primers and probes. 
    make_blast_db(system, forward, reverse, probe, 'forward_db', 'reverse_db',  'probe_db')
    
    # Make result and summary lists to collect data for the result and summary output files
    results = []
    summary = {
        'total': len(sequences),
        'successful_ids': [],
        'primer_only_ids': [],
        'incomplete_ids': [],
        'failed_ids': [],
        'valid_sizes': [],
        'oversized': []
    }

    # list for top Blast hits for Top hits output file
    top_hits = []

    # Go through the sequences one at a time and perform Blast against each database (forward primer, reverse primer and probe).
    # Also makes a progress bar for the number of sequences.
    for seq_record in tqdm(sequences, desc="Processing sequences"):
        seq_id = seq_record.id
        temp_query = 'temp_query.fasta'
        with open(temp_query, 'w') as f:
            SeqIO.write(seq_record, f, 'fasta')

        # Blast against each database
        forward_hits = blast_query(system, temp_query, 'forward_db')
        reverse_hits = blast_query(system, temp_query, 'reverse_db')
        probe_hits = blast_query(system, temp_query, 'probe_db')
        
        # list of all matches to then get top hits
        matched_result = []

        # If there are matches for forward and reverse primers, find the length of the product and confirm they are on opposite strands
        if not forward_hits.empty and not reverse_hits.empty:
            matched = False
            for fwd in forward_hits.itertuples(index=False):
                for rev in reverse_hits.itertuples(index=False):
                    fwd_start, fwd_end = sorted([fwd.qstart, fwd.qend])
                    rev_start, rev_end = sorted([rev.qstart, rev.qend])
                    amp_start, amp_end = min(fwd_start, rev_start), max(fwd_end, rev_end)
                    product_size = amp_end - amp_start
                    opposite_strand = fwd.sstrand != rev.sstrand

                    if not opposite_strand:
                        continue

                    # If the prodcut size is 500bp or less, confirm that the probe is in between the forward and reverse primers.
                    # If the probe is in between the primers, record the data into output files.
                    if product_size <= 500:
                        for prb in probe_hits.itertuples(index=False):
                            probe_start, probe_end = sorted([prb.qstart, prb.qend])
                            in_between = amp_start < probe_start < probe_end < amp_end
                            if in_between:
                                row = [
                                    seq_id,
                                    fwd.sacc, fwd.length, fwd.qstart, fwd.qend, fwd.pident,
                                    prb.sacc, prb.length, prb.qstart, prb.qend, prb.pident,
                                    rev.sacc, rev.length, rev.qstart, rev.qend, rev.pident
                                ]
                                results.append(row)
                                matched_result.append(row)
                                
                                # Gather data on IDs that would make a pcr product and sizes of those products
                                summary['successful_ids'].append(seq_id)
                                summary['valid_sizes'].append(product_size)
                                matched = True

            # Gather data on IDs that would not make a pcr product and products that would be over 500 bp
            if not matched:
                summary['primer_only_ids'].append(seq_id)
                for fwd in forward_hits.itertuples(index=False):
                    for rev in reverse_hits.itertuples(index=False):
                        dist = abs(fwd.qstart - rev.qend)
                        if dist > 500:
                            summary['oversized'].append(dist)
        elif not forward_hits.empty or not reverse_hits.empty:
            summary['incomplete_ids'].append(seq_id)
        else:
            summary['failed_ids'].append(seq_id)

        # Get the best Blast hit for each sequence that would make a pcr product under 500bp
        if matched_result:
            best_match = max(matched_result, key=lambda x: (x[2]+x[7]+x[12], x[5]+x[10]+x[15]))
            top_hits.append(best_match)

        # Delete query temperary file
        os.remove(temp_query)

    # Prepare headers for result file
    columns = ['Strain', 'Forward', 'Fwd length', 'Fwd Start', 'Fwd End', 'Fwd %indent',
               'Probe',  'Probe length', 'Probe Start', 'Probe End', 'Probe %indent',
               'Reverse', 'Rev length', 'Rev Start', 'Rev End', 'Rev %indent']

    # Make the results file
    results_df = pd.DataFrame(results, columns=columns)
    results_df.to_csv('FinalResults.csv', index=False)

    pd.DataFrame(top_hits, columns=columns).to_csv('TopHits.csv', index=False)

    # Remove duplicate IDs before computing summary file
    successful_ids = sorted(set(summary['successful_ids']))
    primer_only_ids = sorted(set(summary['primer_only_ids']))
    incomplete_ids = sorted(set(summary['incomplete_ids'] + summary['failed_ids']))

    with open('Summary.txt', 'w') as f:
        f.write(f"Total sequences processed: {summary['total']}\n")
        f.write(f"Sequences with valid PCR products: {len(successful_ids)}\n")
        f.write(f"Sequences with primer binding but no product: {len(primer_only_ids)}\n")
        f.write(f"Sequences binding to one or no primers: {len(incomplete_ids)}\n")
        f.write(f"\nIDs with valid PCR products:\n" + '\n'.join(successful_ids) + '\n')
        f.write(f"\nIDs with primer binding but no product:\n" + '\n'.join(primer_only_ids) + '\n')
        f.write(f"\nIDs with one or no primer binding:\n" + '\n'.join(incomplete_ids) + '\n')
        if summary['valid_sizes']:
            f.write(f"\nPCR product size range: {min(summary['valid_sizes'])} - {max(summary['valid_sizes'])} bp\n")
        else:
            f.write("\nPCR product size range: None\n")
        if summary['oversized']:
            f.write(f"Oversized product distances (>500 bp): {min(summary['oversized'])} - {max(summary['oversized'])} bp\n")
        else:
            f.write("Oversized product distances (>500 bp): None\n")

    # Delete all temperary files
    cleanup_temp_files()
    click.echo("\nAll finished!")

# Confirm that Blast is installed 
def check_blast_installed(system):
    try:
        command = ['blastn', '-version'] if system != 'windows' else 'blastn -version'
        subprocess.run(command, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True, shell=(system=='windows'))
    except (subprocess.CalledProcessError, FileNotFoundError):
        click.echo("BLAST+ is not installed or not added to PATH.")
        raise SystemExit(1)

# Make Blast Databases for forward primers, reverse primers and probes
def make_blast_db(system, forward, reverse, probe, db_fwd, db_rev, db_probe):
    if system == 'windows':
        subprocess.check_output(f'makeblastdb -in "{forward}" -parse_seqids -dbtype nucl -out "{db_fwd}"', shell=True)
    else:
        subprocess.check_output(f'makeblastdb -in {forward} -parse_seqids -dbtype nucl -out {db_fwd}', shell=True)
    if system == 'windows':
        subprocess.check_output(f'makeblastdb -in "{reverse}" -parse_seqids -dbtype nucl -out "{db_rev}"', shell=True)
    else:
        subprocess.check_output(f'makeblastdb -in {reverse} -parse_seqids -dbtype nucl -out {db_rev}', shell=True)
    if system == 'windows':
        subprocess.check_output(f'makeblastdb -in "{probe}" -parse_seqids -dbtype nucl -out "{db_probe}"', shell=True)
    else:
        subprocess.check_output(f'makeblastdb -in {probe} -parse_seqids -dbtype nucl -out {db_probe}', shell=True)

# Perform Blast on forward database, reverse database and probe database
def blast_query(system, query_file, db_name):
    output_file = 'temp_result.csv'
    if system == 'windows':
        subprocess.check_output(f'blastn -task blastn -query "{query_file}" -db "{db_name}" -out "{output_file}" -outfmt "6 sacc qseqid qstart qend length pident sstrand" -word_size 7 -evalue 1000 -perc_identity 88', shell=True)
    else:
        subprocess.check_output(f'blastn -task blastn -query {query_file} -db {db_name} -out {output_file} -outfmt "6 sacc qseqid qstart qend length pident sstrand" -word_size 7 -evalue 1000 -perc_identity 88', shell=True)

    if os.path.exists(output_file):
        df = pd.read_csv(output_file, sep='\t', header=None, names=['sacc', 'qseqid', 'qstart', 'qend', 'length', 'pident', 'sstrand'])
        os.remove(output_file)
        return df[df.length >= MIN_PRIMER_LENGTH]
    else:
        return pd.DataFrame(columns=['sacc', 'qseqid', 'qstart', 'qend', 'length', 'pident', 'sstrand'])

# Delete temporary files
def cleanup_temp_files():
    for prefix in ['forward_db', 'reverse_db', 'probe_db']:
        for ext in ['nhr', 'nin', 'nsq', 'nos', 'njs', 'not', 'ntf', 'nto', 'ndb', 'nog']:
            file = f'{prefix}.{ext}'
            if os.path.exists(file):
                os.remove(file)

if __name__ == '__main__':
    main()
