ped=$1
bootstraps=$2
for K in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16; do
	admixture --cv --B$bootstraps $ped $K | tee log${K}.out
done
