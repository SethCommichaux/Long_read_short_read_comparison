# INPUT = the reference cgmlst or wgmlst alleles in fasta format as well as the blast results 
# (-outfmt '6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore gaps')
# where the reference alleles were queried against the assembly.
# OUTPUT: the blast distance between the reference cgmlst or wgmlst and the queried assembly 
# we define the blast distance as the total number of mismatches and indels between the blast alignment of the referenece
# alleles and the assembly. if there was not blast hit to a reference allele in the assembly, the whole length of the missing
# gene is counted towards the blast distance

import sys
from Bio import SeqIO

ref_fasta = sys.argv[1]
blast = sys.argv[2]

refs = {str(i.id):len(str(i.seq).replace('-','')) for i in SeqIO.parse(ref_fasta,'fasta')}
results = {}


for i in open(blast):
	tmp = i.strip().split('\t')
	allele = tmp[0]
	aln_len = float(tmp[3])
	mismatch = float(tmp[4])
	gaps = float(tmp[12])
	if allele not in results:
		if aln_len >= refs[allele]:
			results[allele] = gaps+mismatch
		else:
			results[allele] = gaps+mismatch+(refs[allele]-aln_len)

for k,v in refs.items():
	if k in results: refs[k]=results[k]
	else: refs[k] = refs[k]

print int(round(sum(refs.values())))
