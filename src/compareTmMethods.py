#!/usr/bin/env python3

import numpy as np
import random
import sys
import subprocess
from Bio.SeqUtils import MeltingTemp as mt


if (__name__ == "__main__"):
    num_samples = 300
    Tm_diffs = []
    Tm_RNAs = []
    Tm_MELTINGs = []
    for _ in range(num_samples):
        seq = ''.join([random.choice("ACGT") for _ in range(50)])

        input_str = "%s\n" % seq
        cmd = ["/usr/local/src/ViennaRNA-2.4.3/src/bin/RNAplex", "-p", "-l", "50", "-N", "0.05", "-M", "0", "-K", "0", "-U", "0.1", "-Q", "0.00000003"]
        output = subprocess.check_output(cmd, universal_newlines=True, input=input_str)
        Tm_RNA = float(output.strip().split("\n")[-1].split()[-1])
        Tm_RNAs.append(Tm_RNA)
        
        Tm_MELTING = mt.Tm_NN(seq=seq, dnac1=30, nn_table=mt.RNA_NN2, Na=50, Tris=100, Mg=0, saltcorr=2) # 2, 1, 3, 7, 6, 4, 5 are closest to Vienna RNA, in that order
        Tm_MELTINGs.append(Tm_MELTING)

        Tm_diff = Tm_RNA - Tm_MELTING
        Tm_diffs.append(Tm_diff)

    diff_mean = np.mean(Tm_diffs)
    diff_stddev = np.std(Tm_diffs)
    print("Num samples = %d\tMean Tm diff = %4.2f\tStddev Tm diff = %4.2f" % (num_samples, diff_mean, diff_stddev))
    print("Mean Tm ViennaRNA = %4.2f\tMean Tm MELTING = %4.2f" % (np.mean(Tm_RNAs), np.mean(Tm_MELTINGs)))

    sys.exit(0)
