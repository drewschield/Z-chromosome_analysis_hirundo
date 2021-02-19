import sys

out = open(sys.argv[3],'w')
out.write('scaffold'+'\t'+'start'+'\t'+'end'+'\t'+'mid'+'\t'+'sites'+'\t'+'lnL'+'\t'+'topo1'+'\t'+'topo2'+'\t'+'topo3'+'\n')

dat = []

with open(sys.argv[1],'r') as weights:
	next(weights)
	for line in weights:
		datum = line.strip()
		dat.append(datum)

for line in open(sys.argv[2],'r'):
	chrom = line.split()[0]
	start = int(line.split()[1])
	end = int(line.split()[2])
	
	for d in dat:
		if d.split()[3] != 'nan':
			if str(d.split()[0]) == str(chrom) and int(d.split()[3]) >= start and int(d.split()[3]) <= end:
				out.write(str(d)+'\n')
