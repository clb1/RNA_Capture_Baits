class Bait(object):
    def __init__(self, src_segment, chromosome, start, stop, strand, seq, Tm, frac_hairpin, frac_homodimer, opt_position_deviation):
        self.label = None
        self.src_segment = src_segment
        self.chromosome = chromosome
        self.start = start
        self.stop = stop

        self.target_isoform_strand = self.src_segment.getStrand()
        self.bait_strand = strand
        assert (self.target_isoform_strand != self.bait_strand), "Isoform and bait strands must be opposite"

        self.target_isoforms = set() # TODO: delete
        self.num_offtarget_hybs = None
        self.sum_duplex_penalty = None
        self.sum_frac_identity = None

        self.opt_position_deviation = opt_position_deviation
        
        self.Tm = Tm
        #self.opt_Tm_dev = opt_Tm_dev # absolute deviation from optimal Tm
        self.num_lc = sum([seq[i].islower() for i in range(len(seq))])
        self.frac_hairpin = frac_hairpin
        self.frac_homodimer = frac_homodimer
        self.comparitor_tuple = None

        self.bait_sequence = seq    # Reverse complement of the subsequence of the target RNA


    def __str__(self):
        return "%s:%d-%d %s" % (self.chromosome, self.start, self.stop, self.bait_strand)

    def __hash__(self):
        return id(self)
    
    def __eq__(self, other_bait):
        return self.comparitor_tuple == other_bait.getComparitorTuple()

    def __ne__(self, other_bait):
        return self.comparitor_tuple != other_bait.getComparitorTuple()

    def __lt__(self, other_bait):
        return self.comparitor_tuple < other_bait.getComparitorTuple()

    def __le__(self, other_bait):
        return self.comparitor_tuple <= other_bait.getComparitorTuple()

    def __gt__(self, other_bait):
        return self.comparitor_tuple > other_bait.getComparitorTuple()

    def __ge__(self, other_bait):
        return self.comparitor_tuple >= other_bait.getComparitorTuple()


    def setComparitorTuple(self, Ta):
        if (self.sum_frac_identity is None):
            self.comparitor_tuple = (self.num_lc, self.frac_hairpin, self.frac_homodimer, abs(Ta-self.Tm), self.opt_position_deviation)
        else:
            self.comparitor_tuple = (self.sum_frac_identity, self.num_lc, self.frac_hairpin, self.frac_homodimer, abs(Ta-self.Tm), self.opt_position_deviation)


    def getComparitorTuple(self):
        return self.comparitor_tuple

            
    def getStart(self):
        return self.start

    
    def getChromStartStopStrand(self):
        return (self.chromosome, self.start, self.stop, self.bait_strand)


    def overlapsWith(self, other_bait):
        other_chrom, other_start, other_stop, other_strand = other_bait.getChromStartStopStrand()
        return (self.chromosome==other_chrom and self.bait_strand==other_strand and not(other_stop < self.start or other_start > self.stop))


    def setLabel(self, label):
        assert (self.label is None), "Bait label already set"
        self.label = label


    def getLabel(self):
        assert (self.label is not None), "Bait label not set"
        return self.label


    def getSegmentLabel(self):
        return self.src_segment.getLabel()


    def getSequence(self):
        assert (self.bait_sequence is not None), "Bait sequence not set"
        return self.bait_sequence


    def getTm(self):
        return self.Tm
    

    def getAsBED(self):
        assert (self.label is not None), "Bait label not set"
        return (self.chromosome, self.start-1, self.stop, self.label, len(self.target_isoforms), self.bait_strand)


    def setNumOfftargetData(self, num_hybs, sum_duplex_penalty, sum_frac_identity):
        self.num_offtarget_hybs = num_hybs
        self.sum_duplex_penalty = sum_duplex_penalty
        self.sum_frac_identity = sum_frac_identity
        

    def getTargetIsoformIDs(self):
        self.src_segment.getTargetIsoformIDs()

