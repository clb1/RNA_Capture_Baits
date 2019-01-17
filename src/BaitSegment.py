
class BaitSegment(object):
    def __init__(self, chromosome, start, stop, strand, isoform_exons, bait_length):
        self.label = None
        self.chromosome = chromosome
        self.start = start  # 1-based
        self.stop = stop    # 1-based
        self.strand = strand # Strand of the isoforms for this segment
        self.comparitor_tuple = (self.chromosome, self.start, self.stop)
        
        self.bait_length = bait_length
        
        # List of (isoform_ID, exon_num for targeted exon in isoform_ID, number of exons for isoform_ID)
        self.isoform_exon_num = {}
        self.isoform_num_exons = {}
        self.isoform_is_a_target = {}
        self.target_locus = None

        # Record whether this segment targets the last exon in multi-exon isoforms
        all_is_last_exon = []
        n, d = 0, 0
        sum_exon_num = 0
        target_loci = []
        for isoform_ID, exon_num, num_isoform_exons, is_a_target in isoform_exons:
            self.isoform_exon_num[isoform_ID] = exon_num
            self.isoform_num_exons[isoform_ID] = num_isoform_exons
            self.isoform_is_a_target[isoform_ID] = is_a_target
            if (is_a_target):
                CG_locus = isoform_ID.split('.')[0]
                target_loci.append(CG_locus)
                d += 1
                all_is_last_exon.append(exon_num == num_isoform_exons and exon_num > 1)
                if (exon_num == num_isoform_exons and exon_num > 1):
                    n += 1
                sum_exon_num += exon_num

        self.target_locus = frozenset(target_loci)
        self.frac_last_exon = float(n)/float(d)
        self.avg_exon_num = round(sum_exon_num/float(d), 2)
        self.is_last_exon = all(all_is_last_exon)
                                         
        self.ranked_candidate_baits = None
        self.num_pass = None
        self.num_fail_Tm = None
        self.num_fail_hairpin = None
        self.num_fail_homodimer = None


    def __len__(self):
        return self.stop - self.start + 1

    def __str__(self):
        return "%s:%d-%d %s" % (self.chromosome, self.start, self.stop, self.strand)

    def __hash__(self):
        return id(self)

    def __eq__(self, other_segment):
        return self.comparitor_tuple == other_segment.getComparitorTuple()

    def __ne__(self, other_segment):
        return self.comparitor_tuple != other_segment.getComparitorTuple()

    def __lt__(self, other_segment):
        return self.comparitor_tuple < other_segment.getComparitorTuple()

    def __le__(self, other_segment):
        return self.comparitor_tuple <= other_segment.getComparitorTuple()

    def __gt__(self, other_segment):
        return self.comparitor_tuple > other_segment.getComparitorTuple()

    def __ge__(self, other_segment):
        return self.comparitor_tuple >= other_segment.getComparitorTuple()


    def getComparitorTuple(self):
        return self.comparitor_tuple


    def setLabel(self, label):
        assert (self.label is None), "Segment label already set"
        self.label = label


    def getLabel(self):
        assert (self.label is not None), "Segment label not set"
        return self.label


    def getAsBED(self):
        assert (self.label is not None), "Segment label not set"
        return (self.chromosome, self.start-1, self.stop, self.label, len(self.isoform_is_a_target), self.strand)


    def getTargetIsoformIDs(self):
        return set(filter(lambda i: self.isoform_is_a_target[i], self.isoform_is_a_target.keys()))


    def getTargetIsoformIDsAndExons(self):
        ret_list = []
        for isoform_ID in filter(lambda i: self.isoform_is_a_target[i], self.isoform_is_a_target.keys()):
            ret_list.append( (isoform_ID, self.isoform_exon_num[isoform_ID]) )
        return ret_list


    def targetsSimpleIsoform(self, min_num_exons):
        does_target_simple = False
        for isoform_ID in filter(lambda i: self.isoform_is_a_target[i], self.isoform_is_a_target.keys()):
            if (self.isoform_num_exons[isoform_ID] <= min_num_exons):
                does_target_simple = True
                break
        return does_target_simple


    def targetsOnlySimpleIsoforms(self, max_num_exons):
        does_target_simple = True
        for isoform_ID in filter(lambda i: self.isoform_is_a_target[i], self.isoform_is_a_target.keys()):
            if (self.isoform_num_exons[isoform_ID] > max_num_exons):
                does_target_simple = False
                break
        return does_target_simple


    def targetsUntargetedExon(self, exons_w_segments_per_isoform):
        does_target_untargeted = False
        for isoform_ID in filter(lambda i: self.isoform_is_a_target[i], self.isoform_is_a_target.keys()):
            if (self.isoform_exon_num[isoform_ID] not in exons_w_segments_per_isoform[isoform_ID]):
                does_target_untargeted = True
                break
        return does_target_untargeted


    def targetsLastExonOf(self, isoform_ID):
        return self.isoform_exon_num[isoform_ID] == self.isoform_num_exons[isoform_ID]


    def targetsExonToAvoid(self, exons_want_to_avoid):
        return len(self.isoform_exon_num.items() & exons_want_to_avoid) > 0
    
        
    def getChromosome(self):
        return self.chromosome

    
    def getStrand(self):
        return self.strand


    def getStartStop(self):
        return (self.start, self.stop)


    def getChromStartStop(self):
        return (self.chromosome, self.start, self.stop)


    def getChromStartStopStrand(self):
        return (self.chromosome, self.start, self.stop, self.strand)


    def resetStartStop(self, new_start, new_stop):
        self.start = new_start
        self.stop = new_stop
        self.comparitor_tuple = (self.chromosome, self.start, self.stop)
        
        
    def getTargetLocus(self):
        return self.target_locus


    def fromSameLocus(self, locus_set):
        return not self.target_locus.isdisjoint(locus_set)
    

    def overlapsWith(self, other_segment):
        other_chrom, other_start, other_stop, other_strand = other_segment.getChromStartStopStrand()
        return (self.chromosome==other_chrom and self.strand==other_strand and not(other_stop < self.start or other_start > self.stop))


    def fracCoverageBy(self, other_segment):
        '''Assumes that this segment and other_segment have been predetermined to overlap'''
        other_chrom, other_start, other_stop, other_strand = other_segment.getChromStartStopStrand()
        assert (self.chrom == other_chrom and self.strand == other_strand)

        if (self.start <= other_start <= self.stop):
            frac_coverage = float(self.stop - other_start) / float(self.stop - self.start)
        elif (self.start <= other_stop <= self.stop):
            frac_coverage = float(other_stop - self.start) / float(self.stop - self.start)
        else:
            print("ERROR: segments %s and %s do not overlap" % (str(self), str(other_segment)), file=sys.stderr)
            sys.exit(1)

        return frac_coverage


    def noncoverageBy(self, other_segment):
        '''Assumes that this segment and other_segment have been predetermined to overlap'''
        other_chrom, other_start, other_stop, other_strand = other_segment.getChromStartStopStrand()
        assert (self.chromosome == other_chrom and self.strand == other_strand)

        if (self.start <= other_start <= self.stop):
            nonoverlap_len = other_start - self.start
        elif (self.start <= other_stop <= self.stop):
            nonoverlap_len = self.stop - other_stop
        else:
            print("ERROR: segments %s and %s do not overlap" % (str(self), str(other_segment)), file=sys.stderr)
            sys.exit(1)

        return nonoverlap_len
    

    def getMidpoint(self):
        return int((self.stop - self.start)/2)


    def getExonNumberForEachIsoform(self, consider_only_targets):
        ret_list = []
        if (consider_only_targets):
            for isoform_ID, exon_num in filter(lambda x: x[0] in self.isoform_is_a_target, self.isoform_exon_num.items()):
                is_last_exon = exon_num == self.isoform_num_exons[isoform_ID]
                ret_list.append( (isoform_ID, exon_num, is_last_exon) )
        else:
            for isoform_ID, exon_num in self.isoform_exon_num.items():
                is_last_exon = exon_num == self.isoform_num_exons[isoform_ID]
                ret_list.append( (isoform_ID, exon_num, is_last_exon) )

        return ret_list
    

    def getBaitLength(self):
        return self.bait_length

    
    def getBaitAtIndex(self, index):
        return self.ranked_candidate_baits[index]

    
    def isLastExonSegment(self):
        return self.is_last_exon


    def getFracLastExonTargets(self):
        return self.frac_last_exon


    def getAvgExonNum(self):
        return self.avg_exon_num

    
    def getCandidateBaits(self):
        assert (self.ranked_candidate_baits is not None), "Must first evaluate and set the candidate baits"
        return self.ranked_candidate_baits


    def hasCandidateBaits(self):
        assert (self.ranked_candidate_baits is not None), "Must first evaluate and set the candidate baits"
        return len(self.ranked_candidate_baits) > 0


    def getNumCandidateBaits(self):
        assert (self.ranked_candidate_baits is not None), "Must first evaluate and set the candidate baits"
        return len(self.ranked_candidate_baits)


    def getAllPossibleBaits(self, genome_ref):
        baits = {}
        segment_seq = genome_ref[self.chromosome][self.start-1:self.stop]
        for offset in range(self.stop - self.start - self.bait_length + 2):
            bait_seq = segment_seq[offset:offset+self.bait_length]
            assert (len(bait_seq) == self.bait_length)
            if (self.strand == '+'): # Bait strand is on '-' if segment's isoform is on '+'
                bait_seq = bait_seq[::-1].translate(DNA_complement_table)
                bait_strand_name = "minus"
            else:
                bait_strand_name = "plus"
            label = "%s-%d-%d-%s" % (self.chromosome, self.start+offset, self.start+offset+self.bait_length-1, bait_strand_name)  # positions in label are 1-based
            baits[label] = bait_seq
        return baits


    def getTopRankedBait(self):
        self.ranked_candidate_baits.sort()
        return self.ranked_candidate_baits[0]
    

    def getAdditionalBait(self, not_olap_baits):
        secondary_bait = None
        for bait in self.ranked_candidate_baits:
            if (not any(map(lambda b: b.overlapsWith(bait), not_olap_baits))):
                secondary_bait = bait
                break
        return secondary_bait


    def getLowestBaitTm(self):
        return min(map(methodcaller("getTm"), self.ranked_candidate_baits))


    def getOptBaitPosition(self):
        covers_a_first_exon, covers_a_last_exon = (False, False)
        for isoform_ID, is_a_target in self.isoform_is_a_target.items():
            if (is_a_target):
                if (self.isoform_exon_num[isoform_ID] == 1):
                    covers_a_first_exon = True
                if (self.isoform_exon_num[isoform_ID] == self.isoform_num_exons[isoform_ID]):
                    covers_a_last_exon = True

        if (covers_a_first_exon == covers_a_last_exon):
            opt_position = self.start + int((self.stop-self.start)/2)
        elif (covers_a_first_exon):
            opt_position = self.stop if (self.strand == '+') else self.start # Want baits furthest from RNA 5'
        else: # covers_a_last_exon == True
            opt_position = self.start if (self.strand == '+') else self.stop # Want baits furthest from RNA 3'

        return opt_position

    
    def setCandidateBaits(self, bait_data, num_pass, num_fail_Tm, num_fail_hairpin, num_fail_homodimer):
        self.num_pass = num_pass
        self.num_fail_Tm = num_fail_Tm
        self.num_fail_hairpin = num_fail_hairpin
        self.num_fail_homodimer = num_fail_homodimer

        opt_bait_pos = self.getOptBaitPosition()
        self.ranked_candidate_baits = []
        for chromosome, start, stop, bait_strand, bait_seq, Tm, frac_hairpin, frac_homodimer in bait_data:
            opt_position_deviation = abs(opt_bait_pos - int(start + (stop-start)/2))
            new_bait = Bait(self, chromosome, start, stop, bait_strand, bait_seq, Tm, frac_hairpin, frac_homodimer, opt_position_deviation)
            self.ranked_candidate_baits.append(new_bait)
        

    def rankCandidateBaits(self, Ta):
        for bait in self.ranked_candidate_baits:
            bait.setComparitorTuple(Ta)
        self.ranked_candidate_baits.sort()
        

    def trimCandidateBaits(self, Ta, max_num=10):
        had_lc = False

        # Trim off baits that are very similar
        if (len(self.ranked_candidate_baits) > 1):
            self.rankCandidateBaits(Ta)
            to_keep = [self.ranked_candidate_baits[0]]
            to_keep_starts = [self.ranked_candidate_baits[0].getStart()]
            for bait in self.ranked_candidate_baits[1:]:
                this_bait_start = bait.getStart()
                if all(map(lambda another_start: abs(this_bait_start-another_start) >= 10, to_keep_starts)):
                    to_keep.append(bait)
                    to_keep_starts.append(this_bait_start)

                if (len(to_keep) == max_num):
                    break

            self.ranked_candidate_baits = to_keep
            

    def getOrderLabeledBaitSequences(self):
        bait_seqs = []
        for counter, bait in enumerate(self.ranked_candidate_baits):
            bait_seqs.append( ("%s_%d" % (self.label,counter), bait.getSequence()) )
        return bait_seqs
