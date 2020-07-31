# input is contigs in fasta format. outputs assembly size, the number of contigs, the N50 and the length of the longest contig 

import sys
from Bio import SeqIO

def calc_n50(contig_lens,assembly_size):
	n = 0
	for i in sorted(contig_lens,reverse=True):
		n += i
		if n >= assembly_size/2.:
			return i

contig_lens = [len(i.seq) for i in SeqIO.parse(sys.argv[1],'fasta')]

assembly_size = sum(contig_lens)
n50 = calc_n50(contig_lens,assembly_size)
longest_contig = max(contig_lens)
number_contigs = len(contig_lens)

print assembly_size,number_contigs,n50,longest_contig
