import sys

input=sys.argv[1]
info=sys.argv[2]
output=sys.argv[3]

from Bio import SeqIO

circular_contigs=[]
with open(info) as f:
	for line in f:
		values = line.split("\t")
		if values[3] =="+":
			circular_contigs.append(values[0])		
print(circular_contigs)	

records=[]
for record in SeqIO.parse(input, "fasta"):
	#print("%s %i" % (record.id, len(record)))
	#print(record.seq)
	#rotate if circular contig
	if record.id in circular_contigs:
		middle=round(len(record)/2)
		record.seq=record.seq[middle:]+record.seq[:middle]
		records.append(record)
	else:
		records.append(record)

#output new fasta
SeqIO.write(records, output, "fasta")

