python parseDNAFile.py --memory=globalMem --device=0 -m 10 > 0_global_10.out &
python parseDNAFile.py --memory=globalMem --device=1 -m 10 > 1_global_10.out &
python parseDNAFile.py --memory=constantMem --device=2 -m 10 > 0_constant_10.out &
python parseDNAFile.py --memory=constantMem --device=3 -m 10 > 1_constant_10.out &
