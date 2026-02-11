#!/usr/bin/env python3
"""
Bootstrap analysis script for gene family evolutionary rates
Performs bootstrap resampling and MEGA analysis on multiple gene families
"""

import os
import sys
import argparse
import pickle
from pathlib import Path
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
import subprocess
import tempfile
import random
import re
import numpy as np
import pandas as pd
from datetime import datetime

def setup_logging():
    """Setup basic logging to stdout with timestamps"""
    def log(message):
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print(f"[{timestamp}] {message}")
        sys.stdout.flush()  # Ensure immediate output in HPC environment
    return log

def parse_meg_file(path: Path) -> pd.DataFrame:
    """Return a lower-triangle DataFrame of pairwise rates from MEGA .meg output."""
    lines = path.read_text().splitlines()
    # extract labels
    labels = []
    for L in lines:
        m = re.match(r"^\[\s*(\d+)\]\s*#\s*(.+)$", L)
        if m: 
            labels.append(m.group(2).strip())
        if re.match(r"^\[\s*\d+\]", L) and '#' not in L:
            break
    n = len(labels)
    M = np.zeros((n,n))
    # fill matrix
    for L in lines:
        m = re.match(r"^\[\s*(\d+)\]\s*(.*)$", L)
        if not m: 
            continue
        i = int(m.group(1))-1
        vals = re.findall(r"[-+]?\d*\.\d+(?:[eE][-+]?\d+)?", m.group(2))
        if len(vals)==i:
            for k,v in enumerate(vals):
                j = k
                M[i,j] = M[j,i] = float(v)
    np.fill_diagonal(M, 0)
    df = pd.DataFrame(M, index=labels, columns=labels)
    return df.where(np.tril(np.ones(df.shape), k=-1).astype(bool))

def build_cluster_list(families, data_dir, log):
    """Build the cluster list per family"""
    log("Building cluster list per family...")
    cluster_list_per_family = {}
    
    for family in families:
        # ensure the alignments directory exists
        cluster_alignments = f"{data_dir}/sequences_x_updated/{family}_selected_isoform/blastdb/cluster_alignments"
        os.makedirs(cluster_alignments, exist_ok=True)

        # grab every .fa basename in the cluster_fastas dir
        cluster_dir = f"{data_dir}/sequences_x_updated/{family}_selected_isoform/blastdb/cluster_fastas"
        
        if not os.path.exists(cluster_dir):
            log(f"WARNING: Cluster directory does not exist for {family}: {cluster_dir}")
            cluster_list_per_family[family] = []
            continue
            
        all_clusters = [
            os.path.splitext(fn)[0]
            for fn in os.listdir(cluster_dir)
            if fn.endswith(".fa")
        ]

        # filter out FASTAs with only one sequence
        filtered = []
        for name in all_clusters:
            path = os.path.join(cluster_dir, f"{name}.fa")
            try:
                with open(path) as f:
                    nseq = sum(1 for line in f if line.startswith(">"))
                if nseq > 1:
                    filtered.append(name)
            except FileNotFoundError:
                log(f"WARNING: File not found: {path}")
                continue

        # optional sanity-check for duplicate IDs
        for name in filtered:
            path = os.path.join(cluster_dir, f"{name}.fa")
            seen, dups = set(), set()
            with open(path) as f:
                for line in f:
                    if line.startswith(">"):
                        seqid = line[1:].split()[0]
                        if seqid in seen:
                            dups.add(seqid)
                        else:
                            seen.add(seqid)
            if dups:
                log(f"[{family}] {name}.fa has duplicate IDs: {', '.join(dups)}")

        # store the filtered list for later
        cluster_list_per_family[family] = filtered
        log(f"{family}: keeping {len(filtered)} clusters")
    
    return cluster_list_per_family

def run_bootstrap_analysis(families, cluster_list_per_family, data_dir, n_reps, megacc, ds_model, dn_model, fmt, log, output_dir):
    """Main bootstrap analysis function"""
    results = []
    total_clusters = sum(len(clusters) for clusters in cluster_list_per_family.values())
    processed_clusters = 0
    
    log(f"Starting bootstrap analysis with {n_reps} replicates")
    log(f"Total clusters to process: {total_clusters}")
    
    for fam, clusters in cluster_list_per_family.items():
        log(f"Processing family: {fam} ({len(clusters)} clusters)")
        
        for clust in clusters:
            processed_clusters += 1
            log(f"Processing cluster {processed_clusters}/{total_clusters}: {fam}/{clust}")
            
            aln_fp = (
                data_dir/'sequences_x_updated'/f"{fam}_selected_isoform"/
                'blastdb'/'cluster_alignments'/f"{clust}_NT.fa"
            )
            if not aln_fp.exists():
                log(f"WARNING: missing {aln_fp}")
                continue

            try:
                aln = AlignIO.read(aln_fp, fmt)
                L = aln.get_alignment_length()
                assert L % 3 == 0, f"{clust} not codon-aligned"
                n_codons = L // 3

                # stats will hold one mean per replicate per species-pair
                stats = {}  # (sp1,sp2) → {'ds': [], 'dn': []}

                for rep in range(n_reps):
                    if rep % 100 == 0:
                        log(f"  Replicate {rep+1}/{n_reps}")
                    
                    # 1) bootstrap codon blocks
                    idx = [random.randrange(n_codons) for _ in range(n_codons)]
                    boot_recs = []
                    for r in aln:
                        seq = ''.join(str(r.seq[i*3:(i+1)*3]) for i in idx)
                        boot_recs.append(SeqRecord(Seq(seq), id=r.id, description=''))
                    boot = MultipleSeqAlignment(boot_recs)

                    # 2) write temp FASTA & run MEGACC
                    with tempfile.NamedTemporaryFile(suffix='.fa', delete=False) as tf:
                        AlignIO.write(boot, tf.name, fmt)
                        ftmp = tf.name
                    tmp_ds = Path(tempfile.mktemp(suffix='.meg'))
                    tmp_dn = Path(tempfile.mktemp(suffix='.meg'))
                    
                    try:
                        subprocess.run(f"{megacc} -a {ds_model} -d {ftmp} -o {tmp_ds}",
                                     shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                        subprocess.run(f"{megacc} -a {dn_model} -d {ftmp} -o {tmp_dn}",
                                     shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

                        df_ds = parse_meg_file(tmp_ds)
                        df_dn = parse_meg_file(tmp_dn)

                        # 3) collect per-replicate values
                        rep_vals = {}  # (sp1,sp2) → {'ds':[], 'dn':[]}
                        for t1 in df_ds.index:
                            for t2 in df_ds.columns:
                                v_ds = df_ds.loc[t1, t2]
                                if pd.isna(v_ds): 
                                    continue
                                v_dn = df_dn.loc[t1, t2]
                                # extract species suffix
                                sp1 = t1.rsplit('_',1)[-1]
                                sp2 = t2.rsplit('_',1)[-1]
                                key = tuple(sorted((sp1, sp2)))
                                rep_vals.setdefault(key, {'ds':[], 'dn':[]})
                                rep_vals[key]['ds'].append(v_ds)
                                rep_vals[key]['dn'].append(v_dn)

                        # 4) compute means and append
                        for key, vv in rep_vals.items():
                            mean_ds = np.mean(vv['ds'])
                            mean_dn = np.mean(vv['dn'])
                            stats.setdefault(key, {'ds':[], 'dn':[]})
                            stats[key]['ds'].append(mean_ds)
                            stats[key]['dn'].append(mean_dn)

                    except subprocess.CalledProcessError as e:
                        log(f"ERROR: MEGACC failed for {fam}/{clust} replicate {rep}: {e}")
                        continue
                    finally:
                        # cleanup
                        Path(ftmp).unlink(missing_ok=True)
                        tmp_ds.unlink(missing_ok=True)
                        tmp_dn.unlink(missing_ok=True)

                # 5) dump per-cluster
                for (sp1, sp2), vv in stats.items():
                    results.append({
                        'family':   fam,
                        'cluster':  clust,
                        'species1': sp1,
                        'species2': sp2,
                        'dS_rates': vv['ds'],   # length == n_reps
                        'dN_rates': vv['dn']
                    })
                log(f"Completed {fam}/{clust}")
                
                # Save intermediate results every 10 clusters
                if processed_clusters % 10 == 0:
                    intermediate_file = output_dir / f"intermediate_results_{processed_clusters}.pkl"
                    with open(intermediate_file, 'wb') as f:
                        pickle.dump(results, f)
                    log(f"Saved intermediate results to {intermediate_file}")
                
            except Exception as e:
                log(f"ERROR processing {fam}/{clust}: {e}")
                continue
    
    return results

def main():
    parser = argparse.ArgumentParser(description='Bootstrap analysis for gene family evolutionary rates')
    parser.add_argument('--data-dir', required=True, help='Base data directory')
    parser.add_argument('--output-dir', required=True, help='Output directory for results')
    parser.add_argument('--n-reps', type=int, default=1000, help='Number of bootstrap replicates')
    parser.add_argument('--megacc', default='megacc', help='Path to MEGACC executable')
    parser.add_argument('--ds-model', required=True, help='Path to dS model file')
    parser.add_argument('--dn-model', required=True, help='Path to dN model file')
    parser.add_argument('--families', nargs='+', help='List of gene families to process')
    parser.add_argument('--families-file', help='File containing list of gene families (one per line)')
    
    args = parser.parse_args()
    
    # Setup logging
    log = setup_logging()
    
    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Determine families to process
    if args.families_file:
        with open(args.families_file, 'r') as f:
            families = [line.strip() for line in f if line.strip()]
    elif args.families:
        families = args.families
    else:
        # families list for X
        families = ['CSF2RA', 'SPANX', 'TBL1X', 'VCX', 'TMSB', 'MAGEB',
                   'TCEAL8', 'H2A', 'endogenous', 'SPACA5',
                   'SSX', 'GAGE', 'NUDT10', 'CENPVL',
                   'FLJ39060', 'XAGE1', 'FAM156', 'SPIN',
                   'ZXD', 'CXorf49', 'DMRTC1', 'FAM236', 'PABPC', 'RPL36A', 'ARMCX', 'NXF',
                   'TCP11X2', 'GPRASP', 'RAB40A', 'H2BW', 'CT47', 'RHOXF2', 'SMIM10', 'ETD',
                   'INTS6L', 'CT45A', 'CXorf51', 'EOLA', 'HSFX', 'TMEM185A', 'CSAG', 'PNMA',
                   'PWWP4', 'OPN1LW', 'TEX28', 'LAGE3', 'IKBKG', 'F8A1',
                   'collagen', 'LOC129475109', 'LOC115932372']
        # families list for Y
        # families = ['CDY1', 'glutamate', 'TSPY8' ,'DAZ1', 
        #              'BPY2', 'RBMY1B', 'MTRNR2', 'proline', 'VCY1B', 
        #              'HSFY1', 'keratin' ,'FRG1', 
        #              'centriole','FAM47A', 'zinc','isoenzyme', 'retrovirus','TATA-box']
    
    log(f"Processing {len(families)} gene families")
    log(f"Output directory: {output_dir}")
    log(f"Bootstrap replicates: {args.n_reps}")
    
    # Convert paths
    data_dir = Path(args.data_dir)
    ds_model = Path(args.ds_model)
    dn_model = Path(args.dn_model)
    
    # Verify model files exist
    if not ds_model.exists():
        log(f"ERROR: dS model file not found: {ds_model}")
        sys.exit(1)
    if not dn_model.exists():
        log(f"ERROR: dN model file not found: {dn_model}")
        sys.exit(1)
    
    # Build cluster list
    cluster_list_per_family = build_cluster_list(families, data_dir, log)
    
    # Save cluster list
    cluster_file = output_dir / "cluster_list_per_family.pkl"
    with open(cluster_file, 'wb') as f:
        pickle.dump(cluster_list_per_family, f)
    log(f"Saved cluster list to {cluster_file}")
    
    # Run bootstrap analysis
    results = run_bootstrap_analysis(
        families, cluster_list_per_family, data_dir, args.n_reps, 
        args.megacc, ds_model, dn_model, 'fasta', log, output_dir
    )
    
    # Save final results
    log("Saving final results...")
    df = pd.DataFrame(results)
    
    # Save as both pkl and CSV
    results_pkl = output_dir / "bootstrap_results.pkl"
    results_csv = output_dir / "bootstrap_results.csv"
    
    with open(results_pkl, 'wb') as f:
        pickle.dump(results, f)
    
    # For CSV, handle the list columns
    df_csv = df.copy()
    df_csv['dS_rates'] = df_csv['dS_rates'].apply(lambda x: ','.join(map(str, x)))
    df_csv['dN_rates'] = df_csv['dN_rates'].apply(lambda x: ','.join(map(str, x)))
    df_csv.to_csv(results_csv, index=False)
    
    log(f"Results saved to {results_pkl} and {results_csv}")
    log(f"Total results: {len(results)} records")
    log("Analysis complete!")

if __name__ == "__main__":
    main()
