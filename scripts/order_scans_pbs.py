import sys

out = open(sys.argv[3],'w')
out.write('CHROM'+'\t'+'SCAFF'+'\t'+'BIN_START'+'\t'+'BIN_END'+'\t'+'pbs'+'\n')

chrom_name = sys.argv[1].split('list.')[1]
chrom_name = chrom_name.split('.txt')[0]

for line in open(sys.argv[1],'r'):
	chrom = line.split()[0]
	matches = []
	order = []
	with open(sys.argv[2],'r') as pbs:
		next(pbs)
		for p in pbs:
			scaff = p.split()[0]
			start = p.split()[1]
			if str(scaff) == str(chrom):
				matches.append(p)
				order.append(int(start))
	sort_order = sorted(order)
	for o in sort_order:
		for m in matches:
			if int(m.split()[1]) == int(o):
				scaff = m.split()[0]
				start = m.split()[1]
				end = m.split()[2]
				pbs = m.split()[14]
				out.write(str(chrom_name)+'\t'+str(scaff)+'\t'+str(start)+'\t'+str(end)+'\t'+str(pbs)+'\n')
