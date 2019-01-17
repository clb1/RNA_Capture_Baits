#!/usr/bin/env python3

import sys

if (__name__ == "__main__"):
    op = open(sys.argv[1], 'w')

    initial_start = None
    prev_chrom = None
    prev_start = None
    prev_stop = None
    prev_isoform_group = None
    strands = set()
    for line in sys.stdin:
        fields = line.strip().split("\t")
        chrom = fields[0]
        start = int(fields[1])
        stop = int(fields[2])
        isoform_group = fields[3]
        if (prev_isoform_group==None and prev_start==None):
            prev_chrom = chrom
            initial_start = start
            prev_start = start
            prev_stop = stop
            prev_isoform_group = isoform_group
            strands.add(fields[5])
        elif (start - prev_start == 1 and isoform_group == prev_isoform_group):
            prev_chrom = chrom
            prev_start = start
            prev_stop = stop
            strands.add(fields[5])
        else:
            assert(len(strands)==1)
            strand = list(strands)[0]
            segment_len = prev_stop - initial_start
            op.write("%s\t%d\t%d\t%s\t%d\t%s\n" % (prev_chrom, initial_start, prev_stop, prev_isoform_group, segment_len, strand))
            strands.clear()
            strands.add(fields[5])
            prev_chrom = chrom
            initial_start = start
            prev_start = start
            prev_stop = stop
            prev_isoform_group = isoform_group
    
    assert(len(strands)==1)
    strand = list(strands)[0]
    segment_len = prev_stop - initial_start
    op.write("%s\t%d\t%d\t%s\t%d\t%s\n" % (prev_chrom, initial_start, prev_stop, prev_isoform_group, segment_len, strand))

    op.close()
    sys.exit(0)
    
