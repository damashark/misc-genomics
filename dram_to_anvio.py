#!/usr/bin/env python

# Julian Damashek
# Hamilton College
# jdamashe@hamilton.edu, juliandamashek@gmail.com
#
# Parse genes.faa file produced by DRAM into a functions-txt file to import into anvi'o.
# This assumes you have already run DRAM on a genome.
# Additionally, you need to have made an anvi'o contigs database of the same genome, and exported the genes using anvi-export-gene-calls. 
# 
# NOTE: when testing this, I found out .faa file of genes exported from anvi'o had a wrap after a certain number of characters,
#		so the amino acid sequence of long genes were wrapped over >1 lines. The genes.faa file from DRAM did not do this, and only
#		had one line of amino acids per gene. 
#		
#		This was an issue when trying to combine dataframes from both programs.
#
#		It's probably pretty easy to fix this in the script, but that's for another day. I got exams to grade today :)
#		So for now, run anvi-script-reformat-fasta on the exported genes .faa file from anvi'o and check to make sure it's not wrapped.
#		(Or any other tool to make each sequence in a fasta file only 1 line of sequence.)
#
# Inputs needed: genes.faa file outputted by DRAM, GENOME.db_genes.faa file exported by anvi-export-gene-calls.
#				Both need to be 1 description and 1 sequence line per sequence.
#				Script assumes you are in the DRAM output directory, e.g., 'genes.faa' already exists in the working directory.
#
#
#
#
#
# Usage: dram_to_anvio.py [GENOME.db_genes.faa file]
#
# This creates an output file 'dram_annotations_for_anvio.txt' that can be uploaded to an anvi'o database with anvi-import-functions
#

import sys
import pandas as pd
import re

# save separate lists of anvi'o gene names and their amino acid sequences
anvio_gene_ids = []
anvio_gene_seqs = []

with open(str(sys.argv[1]),'r') as anvio_genes:
	for line in anvio_genes: 
		if line.startswith('>'): 
			anvio_gene_ids.append(line.strip('\r\n').split('>')[1])
			anvio_gene_seqs.append(next(anvio_genes).strip('\r\n'))#[:35])

	# combine into single table
	anvio_gene_table = pd.DataFrame(zip(anvio_gene_ids, anvio_gene_seqs), columns=['gene_callers_id', 'AA_sequence'])
	
#sanity check of anvi'o table. Should have gene_caller_ids numbers and amino acid sequences
#anvio_gene_table.to_csv('anvio_only_table.txt', sep='\t', index=False)

# Now read in and process DRAM outputs.
with open('genes.faa', 'r') as dram_file:
		
	gene_count = 0
		
	# lists to store column data for output spreadsheet
	dram_gene_seqs = []
	dram_gene_annot = []
	dram_db = []
	dram_accession = []
	dram_E = []
	
	for line in dram_file:
		if line.startswith('>'):
		
			# save AA sequence, removing * from end if it's there (anvi'o doesn't have them)
			raw_seq = str(next(dram_file).strip('\r\n'))#[:35])
			if raw_seq.endswith('*'):
				dram_gene_seqs.append(re.sub(r'\*$','',raw_seq))
			else:
				dram_gene_seqs.append(raw_seq)
				
			# save blanks for unannotated genes (rank=E from DRAM). Should this just be skipped? Probably. Test later :)
			if line.strip('\r\n').split('rank: ')[1] == 'E':
				dram_gene_annot.append('')
				dram_db.append('')
				dram_accession.append('')
				dram_E.append('')
				
			# save annotations, database, accessions for annotated genes
			else: 
				#annotation
				dram_gene_annot.append(str(line.strip('\r\n').split('; ',1)[1])+" ("+str(line.strip('\r\n').split('rank: ')[1].split(';')[0])+")")
				#database - commented line would add the database from with DRAM, e.g., 'pfam.' 
				#dram_db.append(line.strip('\r\n').split('db=')[1].split(')')[0])
				dram_db.append(str('DRAM'))
				
				#accession number
				if ' [' in line.strip('\r\n'):
					dram_accession.append(line.strip('\r\n').split(' [')[1].split(']')[0])
				else: 
					dram_accession.append('')
				
				#confidence score
				#dram_E.append(line.strip('\r\n').split('rank: ')[1].split(';')[0])
				dram_E.append('')
		gene_count += 1
			
	# combine all DRAM data into one table
	dram_gene_table = pd.DataFrame(zip(dram_gene_seqs, dram_db, dram_accession, dram_gene_annot, dram_E), columns=['AA_sequence', 'source', 'accession', 'function', 'e_value'])
	#sanity check of anvi'o table. Should have gene_caller_ids numbers and amino acid sequences
	#dram_gene_table.to_csv('dram_only_table.txt', sep='\t', index=False)
	
	# join with anvi'o gene names to make final dataframe, save as tsv file
	output_table = pd.merge(anvio_gene_table, dram_gene_table, on='AA_sequence', how='left').drop(columns='AA_sequence')
	output_table.to_csv('dram_annotations_for_anvio.txt', sep='\t', index=False)

print("\nDone! Saved DRAM annotations of %i genes to \'dram_annotations_for_anvio.txt\'\nLoad 'em into anvi'o: https://merenlab.org/2016/06/18/importing-functions/" % (gene_count))		
print("\n~*~ Goodbye. ~*~\n")