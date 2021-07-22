#!/usr/bin/env python3

import argparse
from collections import defaultdict
from fastdist import fastdist
from functools import lru_cache
import pandas as pd
import numpy as np
from scipy.spatial import distance
from sklearn.decomposition import PCA
import sys
from time import time
import logging
import os
import json, gzip

###################
# PARSE ARGUMENTS #
###################

# Create the parser
parser = argparse.ArgumentParser(
    description='Group genes by co-abundance'
)

# Add the arguments
parser.add_argument(
    '--assembly-hdf',
    type=str,
    help='Details of gene assembly in HDF5 format'
)
# Add the arguments
parser.add_argument(
    '--abundance-hdf',
    type=str,
    help='Details of gene abundance in HDF5 format'
)
parser.add_argument(
    '--metric',
    type=str,
    default='cosine',
    help='Distance metric to be used for clustering'
)
parser.add_argument(
    '--threshold',
    type=float,
    default=0.15,
    help='Maximum distance used to join genes'
)
parser.add_argument(
    '--min-contig-size',
    type=int,
    default=3,
    help='Only cluster genes which are found on contigs with at least this number of genes on them'
)
parser.add_argument(
    '--min-contig-depth',
    type=int,
    default=5,
    help='Minimum depth of sequencing per contig'
)
parser.add_argument(
    '--min-specimens',
    type=int,
    default=2,
    help='Only cluster genes which are detected by alignment in at least this number of specimens'
)
parser.add_argument(
    '--max-n-cags',
    type=int,
    default=250000,
    help='Maximum number of CAGs -- halt if this point is reached'
)
parser.add_argument(
    '--output-folder',
    type=str,
    default='./',
    help='Folder for output files'
)
parser.add_argument(
    '--output-prefix',
    type=str,
    default="CAGs",
    help='Prefix for output files'
)

# Parse the arguments
args = parser.parse_args()

# Make sure that the input file exists
for fp in [args.abundance_hdf, args.assembly_hdf, args.output_folder]:
    msg = f"Input HDF {fp} does not exist"
    assert os.path.exists(fp), msg

# Set up logging
logFormatter = logging.Formatter(
    '%(asctime)s %(levelname)-8s [makeCAGs] %(message)s'
)
rootLogger = logging.getLogger()
rootLogger.setLevel(logging.INFO)

# Write to file
log_fp = os.path.join(
    args.output_folder,
    f"{args.output_prefix}.makeCAGs.log"
)
fileHandler = logging.FileHandler(log_fp)
fileHandler.setFormatter(logFormatter)
rootLogger.addHandler(fileHandler)

# Also write to STDOUT
consoleHandler = logging.StreamHandler()
consoleHandler.setFormatter(logFormatter)
rootLogger.addHandler(consoleHandler)

logging.info("Joining genes by co-abundance")
logging.info(f"Input Assembly HDF: {args.assembly_hdf}")
logging.info(f"Input Abundance HDF: {args.abundance_hdf}")
logging.info(f"Metric: {args.metric}")
logging.info(f"Threshold: {args.threshold}")
logging.info(f"Minimum Contig Size: {args.min_contig_size}")
logging.info(f"Minimum Contig Depth: {args.min_contig_depth}")
logging.info(f"Minimum Specimens per Gene: {args.min_specimens}")
logging.info(f"Output Folder: {args.output_folder}")
logging.info(f"Output Prefix: {args.output_prefix}")


######################################
# DEFINE FUNCTIONS USED TO MAKE CAGS #
######################################

# Get the list of specimens which have assembly information
@lru_cache(maxsize=1)
def specimen_list(prefix="/abund/allele/assembly/"):
    with pd.HDFStore(args.assembly_hdf, 'r') as store:
        return [
            p[len(prefix):]
            for p in store
            if p.startswith(prefix)
        ]
    
# Read the contig information for a single assembly
def read_contigs(specimen):
    # Read the table
    with pd.HDFStore(args.assembly_hdf, 'r') as store:
        df = pd.read_hdf(store, "/abund/allele/assembly/%s" % specimen)
        
    # Remove genes which don't match the gene catalog
    df = df.loc[df["catalog_gene"].isnull().apply(lambda v: v is False)]
    
    # Only return the contig and gene lists
    df = df.reindex(
        columns=["contig", "catalog_gene", "multi"]
    ).drop_duplicates(
    ).rename(
        columns=dict(multi='depth')
    )

    # Filter by depth
    df = df.query(
        f"depth >= {args.min_contig_depth}"
    )
    
    # Count the number of genes 
    contig_vc = df.contig.value_counts()
    
    # Assign the contig size and specimen name
    df = df.assign(
        n_genes = df.contig.apply(contig_vc.get)
    ).assign(
        specimen=specimen
    )
    
    # Filter by size
    df = df.query(
        f"n_genes >= {args.min_contig_size}"
    )

    logging.info(f"Read in {df.shape[0]:,} genes on {df.contig.unique().shape[0]:,} contigs from specimen '{specimen}'")
        
    return df.reset_index(drop=True)


# Get the abundance of each gene/allele in a single specimen
def read_abund(specimen, verbose=False):
    if verbose:
        logging.info(f"Reading abundances from {specimen}")

    # Read the table
    with pd.HDFStore(args.abundance_hdf, 'r') as store:
        df = pd.read_hdf(store, f"/abund/gene/long/{specimen}")
    
    # Only return the gene name and the depth of sequencing, formatted as a Series
    df = df.reindex(
        columns=["id", "depth"]
    ).set_index(
        "id"
    )['depth']
    
    # Divide by the total depth of sequencing across all genes
    df = df / df.sum()
    
    return df

# Automatically incrementing index
class IncrementingIndex:
    
    def __init__(self):
        
        self.ix = dict()
        self.ix_name = dict()
        
    def get(self, key):
        
        if self.ix.get(key) is None:
            self.ix_name[len(self.ix)] = key
            self.ix[key] = len(self.ix)
        return self.ix.get(key)
    
    def name(self, ix):
        return self.ix_name.get(ix)

class GeneAbund:
    
    def __init__(self, verbose=True):
        
        # Keep a numeric index of each gene
        self.gene_ix = IncrementingIndex()
        
        # Keep a numeric index of each specimen
        self.specimen_ix = IncrementingIndex()
        
        # Store the abundance of each gene in a dict of dicts
        # The first key is the gene, the second key is the specimen
        self.abund = defaultdict(dict)
        
        # Read in the abundances across all specimens
        for n, specimen in enumerate(specimen_list()):
            
            if verbose and n % 10 == 0 and n > 0:
                logging.info(f"Read in abundances from {n:,} specimens")
                
            # Iterate over every gene observed
            for gene, abund in read_abund(specimen).items():
                
                # Add the measurement to the dict of dicts
                self.abund[
                    self.gene_ix.get(gene)
                ][
                    self.specimen_ix.get(specimen)
                ] = abund
                
        if verbose:
            logging.info(f"Read in abundances from {(n + 1):,} specimens")


class IAC:
    
    def __init__(
        self,
        max_dim=len(specimen_list())
    ):

        # Keep track of the dense abundance vectors for the active set
        self.abunds = np.zeros((0, max_dim), dtype=float)
        
        # Keep track of the members of each group
        self.groups = list()
        
    def add(self, gene_ix_set, gene_abund):
        
        # If genes have already been added to this object
        if len(self.groups) > 0:
        
            # Calculate the co-abundance distance to all gene (groups) in the active set
            dists = 1 - fastdist.vector_to_matrix_distance(
                gene_abund,
                self.abunds,
                fastdist.__dict__[args.metric],
                args.metric
            )

            # Find the lowest distance
            min_dist = None
            closest_group = None
            for ix, d in enumerate(dists):
                if min_dist is None or d < min_dist:
                    closest_group = ix
                    min_dist = d

            # If the distance is not above the threshold
            if min_dist <= args.threshold:

                # Combine the new gene with the previous group
                self.combine(
                    gene_ix_set, 
                    gene_abund, 
                    closest_group
                )
                return
            
        # Otherwise

        # Add this gene to the active set
        self.abunds = np.concatenate((self.abunds, gene_abund[np.newaxis, :]))
        self.groups.append(gene_ix_set)

    def combine(
        self,
        gene_ix_set,
        gene_abund,
        closest_group
    ):
        
        # Compute the combined abundance as a weighted average
        comb_abund = (
            (self.abunds[closest_group] * len(self.groups[closest_group])) + 
            (gene_abund * len(gene_ix_set))
        ) / (
            len(self.groups[closest_group]) + len(gene_ix_set)
        )
        
        # Make the new set of genes
        new_gene_set = gene_ix_set | self.groups[closest_group]
        
        # Remove the closest group from the active set
        del self.groups[closest_group]
        self.abunds = np.delete(self.abunds, closest_group, 0)
        
        # Add the new combined set of genes back to the array
        self.add(
            new_gene_set,
            comb_abund
        )

    def items(self):
        """Yield the members (set) and abundances (array) of the gene groups."""

        # Iterate over each gene group
        for group_ix, group_set in enumerate(self.groups):

            # Yield the set of genes, and the vector of abundances
            yield group_set, self.abunds[group_ix]
        

def make_dense_abund(a):
    # Make a dense abundance vector
    abund = np.zeros((len(specimen_list()),))
    for i, v in a.items():
        abund[i] = v
        
    # If the mean abundance is 0
    if np.sum(abund) == 0:
        
        # Return None
        return
    
    # Return the abundance vector
    return abund


def report_progress(iac):
    
    # Calculate the total number of CAGs across all specimens
    
    # If the user provided a dict of IACs
    if isinstance(iac, dict):
        tot_cags = 0
        tot_genes = 0
        for specimen in iac.keys():
            tot_cags += iac[specimen].abunds.shape[0]
            tot_genes += sum(map(len, iac[specimen].groups))
            
    # Otherwise, it should be a single IAC
    else:
        
        tot_cags = iac.abunds.shape[0]
        tot_genes = sum(map(len, iac.groups))

    # Report the progress
    logging.info(f"Added {tot_genes:,} genes in total -- {tot_cags:,} current CAGs")

    # Reset the clock
    return time()


class SpecimenPCA:
    """Run PCA on a subset of genes to generate embeddings on specimens."""

    def __init__(self, gene_abund, max_genes=1000000, variance_cutoff=0.001):

        # Get the total list of genes which were detected
        gene_key_list = list(gene_abund.gene_ix.ix.keys())
        # Randomly sample up to `max_genes` of those genes
        gene_key_list = pd.Series(gene_key_list).sample(
            min(max_genes, len(gene_key_list))
        ).tolist()

        n_genes = len(gene_key_list)
        n_specimens = len(specimen_list())
        logging.info(f"Running PCA on a subset of {n_genes:,} genes across all {n_specimens:,} specimens")

        # Set up an empty array
        abund_array = np.zeros((len(specimen_list()), len(gene_key_list)))

        # Iterate over each of the genes which were selected
        for gene_ix, gene_key in enumerate(gene_key_list):

            # Get the abundance for this gene
            a = make_dense_abund(
                gene_abund.abund[
                    gene_abund.gene_ix.get(
                        gene_key
                    )
                ]
            )

            # Add the abundances for that gene to the array
            abund_array[:, gene_ix] = a

        # Set up PCA
        self.pca = PCA()
        
        # Fit the data
        self.pca.fit(abund_array.T)

        # Calculate the number of dimensions which capture 1-variance_cutoff variance
        cumulative_exp = np.cumsum(self.pca.explained_variance_ratio_)

        # Set the maximum number of dimensions to use
        self.max_dim = len(specimen_list())
        for i, v in enumerate(cumulative_exp):
            if v > (1 - variance_cutoff):
                self.max_dim = (i + 1)
                break

        logging.info(f"Using {self.max_dim:,} / {len(specimen_list()):,} PCA dimensions")


def co_assembly_boosted_clustering(
    contig_df,
    gene_abund,
    log_interval=300,
):
    
    # Start the clock
    start_time = time()

    # Start by running PCA on specimens using a random subset of genes
    specimen_pca = SpecimenPCA(gene_abund)

    # Now let's use the approach where we add genes which are found on the same contig first

    # Set up the IAC object for the overall assembly
    iac = IAC(max_dim=specimen_pca.max_dim)
    
    # Keep track of what genes have been added
    all_genes = set([])

    # Keep track of which contig we are considering at each round
    current_contig = None
    current_contig_iac = None

    # Iterate over every gene, starting with the longest contigs (due to previous sorting)
    for _, r in contig_df.iterrows():
    
        # If we've already added this gene
        if r.catalog_gene in all_genes:

            # Skip it
            continue
            
        # Mark that we've added this gene now
        all_genes.add(r.catalog_gene)

        # Get the index of this gene
        gene_ix = gene_abund.gene_ix.get(r.catalog_gene)

        # If this gene does not meet the --min-specimens prevalence threshold
        if len(gene_abund.abund[gene_ix]) < args.min_specimens:

            # Skip it
            continue

        # Make a dense abundance vector
        abund = make_dense_abund(
            gene_abund.abund[
                gene_ix
            ]
        )

        # Transform the abundance based on PCA weightings
        embedded_abund = specimen_pca.pca.transform(
            abund[np.newaxis, :]
        )[0, :specimen_pca.max_dim]

        # If this is a new contig
        if current_contig is None or current_contig != r["contig"]:

            # If we have grouped genes for the previous contig
            if current_contig_iac is not None:

                # Iterate over the previous set of groups
                for prev_genes, prev_abund in current_contig_iac.items():

                    # Add them to the global IAC
                    iac.add(prev_genes, prev_abund)

                    # If we have exceeded the args.max_n_cags threshold
                    if iac.abunds.shape[0] > args.max_n_cags:

                        # Bail out of the process
                        logging.info(f"Reached limit of {args.max_n_cags:,}, stopping.")
                        return

            # Set up an object for the genes in this new contig
            current_contig = r["contig"]
            current_contig_iac = IAC(max_dim=specimen_pca.max_dim)

        # Add this gene to the IAC for this contig
        current_contig_iac.add(
            set([gene_ix]),
            embedded_abund
        )

        # If more than log_interval seconds have passed
        if (time() - start_time) >= log_interval:

            # Report the progress
            start_time = report_progress(iac)

    # Iterate over the set of groups from the last contig
    for prev_genes, prev_abund in current_contig_iac.items():

        # Add them to the global IAC
        iac.add(prev_genes, prev_abund)

        # If we have exceeded the args.max_n_cags threshold
        if iac.abunds.shape[0] > args.max_n_cags:

            # Bail out of the process
            logging.info(f"Reached limit of {args.max_n_cags:,}, stopping.")
            return

    # Report the end of the entire process
    logging.info("FINISHED")
    report_progress(iac)

    # Before returning the results, add a dict with the size-ordered index of each CAG
    # Keys are the index position of CAGs in iac.groups, Values are the sorted CAG IDs
    iac.cag_id = pd.DataFrame([
        {
            "ix": ix,
            "size": len(l)
        }
        for ix, l in enumerate(iac.groups)
    ]).sort_values(
        by="size",
        ascending=False
    ).assign(
        cag_id = range(len(iac.groups))
    ).set_index(
        "ix"
    )[
        "cag_id"
    ]

    # Return the final IAC that has been constructed
    return iac


##############
# ENTRYPOINT #
##############

if __name__ == "__main__":

    # Get the list of specimens
    logging.info(f"There are {len(specimen_list()):,} specimens present.")

    # Read in assembly information for all specimens
    assembly_start = time() # Start the clock
    contig_df = pd.concat(
        read_contigs(specimen)
        for specimen in specimen_list()
    ).sort_values(
        by=["n_genes", "contig"], # Keep the contig labels sorted within each size rank
        ascending=False
    ).reset_index(
        drop=True
    )

    # Report the total number of genes
    logging.info(f"Read in {contig_df.catalog_gene.unique().shape[0]:,} unique genes on {contig_df.contig.unique().shape[0]:,} contigs.")

    # Record the total time for reading in assembly information
    assembly_elapsed = time() - assembly_start

    # Read in the gene abundances across all specimens
    gene_abund_start = time()
    gene_abund = GeneAbund()
    gene_abund_elapsed = time() - gene_abund_start

    # Start the clock for the IAC process
    iac_start = time()

    iac = co_assembly_boosted_clustering(
        contig_df,
        gene_abund,
    )
    iac_elapsed = time() - iac_start
    logging.info(f"DONE CLUSTERING - {round(iac_elapsed, 1):,} seconds")

    # Record the CAG assignments as CSV
    output_fp = os.path.join(
        args.output_folder, 
        f"{args.output_prefix}.csv.gz"
    )
    output_start = time()
    logging.info(f"Writing out CAG assignments to {output_fp}")

    # Compress the output with gzip
    with gzip.open(output_fp, 'wt') as handle:

        # Write out the header
        handle.write('gene,CAG\n')

        # Iterate over every CAG
        for groups_ix, cag_id in iac.cag_id.items():

            # If the gropu is a single unclustered gene
            if len(iac.groups[groups_ix]) == 1:

                # Skip it
                continue

            # Otherwise, iterate over every gene in the CAG
            for gene_id in iac.groups[groups_ix]:

                # Get the name from the numeric gene ID
                gene_name = gene_abund.gene_ix.name(gene_id)

                # Write out the pairing of the gene with the CAG in long format
                handle.write(f"{gene_name},{cag_id}\n")

    output_elapsed = time() - output_start
    logging.info(f"Writing out results took {round(output_elapsed, 2)} seconds")

    # Record the total time elapsed
    with open(
        os.path.join(
            args.output_folder, 
            f"{args.output_prefix}.time.json"
        ), 
        'w'
    ) as handle:
        handle.write(json.dumps(dict(
            assembly=assembly_elapsed,
            abundance=gene_abund_elapsed,
            clustering=iac_elapsed,
            output=output_elapsed,
        )))

    logging.info("DONE SAVING RESULTS")
