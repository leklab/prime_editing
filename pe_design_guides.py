import pysam
import argparse


#padding = 250
padding = 30
#padding = 1


pbs_len = 13
rtt_len = 13

pam_len = 20

def reverse_complement(seq):
	rc_seq = ""

	for i in range(0,len(seq)):
		if (seq[i]=="A"):
			rc_seq="T"+rc_seq
		elif (seq[i]=="T"):
			rc_seq="A"+rc_seq
		elif (seq[i]=="G"):
			rc_seq="C"+rc_seq
		elif (seq[i]=="C"):
			rc_seq="G"+rc_seq
		else:
			rc_seq=seq[i]+rc_seq

	return rc_seq


def get_sequence(chrom,pos):

	f1 = pysam.FastaFile('/home/ml2529/Prime_Editing/Homo_sapiens_assembly19.fasta')

	variant_region = '%s:%d-%d' % (chrom,pos-padding, pos+padding)

	#print(variant_region)
	seq = f1.fetch(region=variant_region)

	#print(seq)	

	return seq

def nearest_pam(seq, pos, edit):


	mod_seq = seq[0:padding] + edit + seq[padding+1:len(seq)] 
	rc_seq = reverse_complement(mod_seq)


	pos_min_pos = padding

	for i in range(0, len(seq)-1):
		if (seq[i]=="G") and (seq[i+1]=="G"):

			relative_pos = i-1-padding

			if (relative_pos >= -4 and relative_pos < pos_min_pos):
				pos_min_pos = relative_pos

			'''
			relative_pos = i-1-padding
			relative_pos = i-1-padding
			nick_pos = relative_pos - 4

			pam_end = i-2
			pam_seq = mod_seq[pam_end-pam_len+1:pam_end+1]
			
			pbs_end = (i-1)-4
			rtt_start = (i-1)-3

			pbs_seq = reverse_complement(mod_seq[pbs_end-pbs_len+1:pbs_end+1])
			rtt_seq = reverse_complement(mod_seq[rtt_start:rtt_start+rtt_len])

			print("Seq: %s Pos: %d Nick_pos: %d PAM_seq: %s RTT_seq: %s PBS_seq: %s" % (rc_seq[i-3:i+4],relative_pos, nick_pos, pam_seq, rtt_seq, pbs_seq))
			'''
			

			#print("Seq: %s Pos: %d" % (seq[i-3:i+4],relative_pos))



	#print(rc_seq)
	#print(mod_seq)

	neg_min_pos = padding

	for i in range(0, len(rc_seq)-1):
		if (rc_seq[i]=="G") and (rc_seq[i+1]=="G"):

			relative_pos = i-1-padding

			if (relative_pos >= -4 and relative_pos < neg_min_pos):
				neg_min_pos = relative_pos

			'''
			relative_pos = i-1-padding
			nick_pos = relative_pos - 4

			pam_end = i-2
			pam_seq = rc_seq[pam_end-pam_len+1:pam_end+1]
			
			pbs_end = (i-1)-4
			rtt_start = (i-1)-3

			pbs_seq = reverse_complement(rc_seq[pbs_end-pbs_len+1:pbs_end+1])
			rtt_seq = reverse_complement(rc_seq[rtt_start:rtt_start+rtt_len])

			print("Seq: %s Pos: %d Nick_pos: %d PAM_seq: %s RTT_seq: %s PBS_seq: %s" % (rc_seq[i-3:i+4],relative_pos, nick_pos, pam_seq, rtt_seq, pbs_seq))
			'''

	#print("Min pos: %d Min pos: %d" % (pos_min_pos,neg_min_pos))

	if(pos_min_pos < neg_min_pos):
		print_oligos(mod_seq,pos_min_pos,edit)
	elif(neg_min_pos < padding):
		print_oligos(rc_seq,neg_min_pos,edit)


	#print_oligos(seq,neg_min_pos,edit)



def print_oligos(seq,relative_pos,edit):

	#mod_seq = seq[0:padding] + edit + seq[padding+1:len(seq)] 
	#rc_seq = reverse_complement(mod_seq)

	i = relative_pos + 1 + padding	
	nick_pos = relative_pos - 4

	pam_end = i-2
	pam_seq = seq[pam_end-pam_len+1:pam_end+1]
	
	pbs_end = (i-1)-4
	rtt_start = (i-1)-3

	pbs_seq = reverse_complement(seq[pbs_end-pbs_len+1:pbs_end+1])
	rtt_seq = reverse_complement(seq[rtt_start:rtt_start+rtt_len])

	print("Nick_pos: %d PAM_seq: %s RTT_seq: %s PBS_seq: %s" % (nick_pos, pam_seq, rtt_seq, pbs_seq))




def main(args):
	pos = int(args.pos)
	seq = get_sequence(args.chr,pos)

	nearest_pam(seq,pos,args.edit)




if __name__ == '__main__':
	parser = argparse.ArgumentParser()

	#parser.add_argument('--vcf', '--input', '-i', help='Input VCF file (.vcf) or gzipped VCF file (.vcf.gz)', required=True)
	#parser.add_argument('--out', '-o', help='Output tsv file', default="")

	parser.add_argument('--pos', '-p', help='Variant Position', default="", required=True)
	parser.add_argument('--chr', '-c', help='Variant Position', default="", required=True)
	parser.add_argument('--edit', '-e', help='Variant Position', default="*")


	parser.add_argument('--pbs_len', help='Variant Position', default=13)
	parser.add_argument('--rtt_len', help='Variant Position', default=13)
	parser.add_argument('--fasta', help='Variant Position', default="")


	args = parser.parse_args()

	main(args)

