#!/usr/bin/env python3

import sys
import copy
import gzip
import numpy as np
import os
import tempfile
import subprocess
import multiprocessing
from multiprocessing.pool import Pool
from subprocess import check_output, CalledProcessError, Popen
from itertools import chain, combinations
from operator import itemgetter, methodcaller
from collections import defaultdict, Counter
import pickle

# Uses BED-style indexing.
# Start index is 0-based and stop index is the 1-based index of the last nucleotide position that you want (or, 0-based position just after the last wanted nucleotide).
from pyfaidx import Fasta

from RNAThermodynamics import RNAThermodynamics
import Bait, BaitSegment

DNA_complement_table = str.maketrans("ACGTNacgtn","TGCANtgcan")


def readTargetIsoformIDs(target_isoforms_file):
    target_isoform_IDs = set()

    with open(target_isoforms_file, 'r') as ip:
        for line in ip:
            isoform_ID = line.strip()
            assert (isoform_ID.startswith("CG"))
            target_isoform_IDs.add(isoform_ID)

    return target_isoform_IDs


def getBaitSegments(bait_segments_bed_gz, target_isoform_IDs, bait_length):
    all_segments = set()

    isoform_map_segments = defaultdict(list)
    #exons_want_to_avoid = set()
    with gzip.open(bait_segments_bed_gz, 'rb') as ip:
        for line in ip:
            chromosome, start, stop, annot_isoform_IDs, segment_len, isoforms_strand = line.decode().strip().split("\t")
            assert (int(segment_len) >= bait_length), "Found segment shorter than bait length (%s < %d)" % (segment_len, bait_length)
            isoform_exons = []
            min_num_exons = 1e10
            for annot_isoform_ID in annot_isoform_IDs.split(','):
                isoform_ID, exon_info = annot_isoform_ID.rsplit('_',1)
                exon_num, num_exons_for_isoform = map(int, exon_info.split('/'))
                min_num_exons = min(min_num_exons, num_exons_for_isoform)
                isoform_exons.append( (isoform_ID, exon_num, num_exons_for_isoform, isoform_ID in target_isoform_IDs) )

            assert (any(map(lambda i: i in target_isoform_IDs, map(itemgetter(0), isoform_exons)))), "Shouldn't have segments without any target isoforms"
            new_segment = BaitSegment(chromosome, int(start)+1, int(stop), isoforms_strand, isoform_exons, bait_length)
            all_segments.add(new_segment)
            for isoform_ID in map(itemgetter(0), isoform_exons):
                isoform_map_segments[isoform_ID].append(new_segment)

    # For isoforms with a low number of segments, try splitting segments so that
    # they can have more than one segment, enabling more baits for the isoform.
    # Don't split segments that target only last exons.
    segments_to_try_splitting = set()
    for isoform_ID, segments in filter(lambda t: len(t[1])<20 and t[0] in target_isoform_IDs, isoform_map_segments.items()):
        segments_to_try_splitting.update( filter(lambda s: s.getFracLastExonTargets() == 0.0 or s.targetsOnlySimpleIsoforms(7), segments) )
            
    all_replaced_segments = set()
    all_replacement_segments = list()
    for segment_to_split in segments_to_try_splitting:
        max_splits = 2 if (segment_to_split.getFracLastExonTargets() > 0.0) else 5
        replacement_segments = splitSegment(segment_to_split, bait_length, max_splits)
        if (replacement_segments is not None):
            all_replaced_segments.add(segment_to_split)
            all_replacement_segments.extend(replacement_segments)

    print("INFO: split %d segments into %d segments" % (len(all_replaced_segments), len(all_replacement_segments)), file=sys.stderr, flush=True)

    all_segments -= all_replaced_segments
    all_segments.update( all_replacement_segments )

    # Confirm that all targets have at least one segment
    no_segments, num_nonlast_segments_per_isoform = checkIsoformSegmentCoverage(target_isoform_IDs, all_segments)
    if (len(no_segments) > 0):
        print("WARNING: After reading and splitting segments, %d isoforms have no segments" % len(no_segments), file=sys.stderr)

    return all_segments


def splitSegment(segment_to_split, bait_length, max_splits):
    replacement_segments = None

    segment_len = len(segment_to_split)
    for bait_len_multiplier in list(map(lambda x: float(x)/100.0, range(300, 1001, 50))):
        num_subsegments = int(segment_len/(bait_len_multiplier*bait_length)) # Calc number of segments made when segment length is 10x->3x the bait length
        if (num_subsegments >= 2):
            num_subsegments = max_splits if (num_subsegments > max_splits) else num_subsegments
            break

    if (num_subsegments > 1):
        replacement_segments = []
        start, stop = segment_to_split.getStartStop()
        for new_start, new_stop in map(lambda x: (min(x),max(x)), np.array_split(range(start,stop+1), num_subsegments)):
            new_subsegment = copy.deepcopy(segment_to_split)
            new_subsegment.resetStartStop(int(new_start), int(new_stop))
            replacement_segments.append(new_subsegment)

    return replacement_segments


def setPerIsoformData(all_segments, target_isoform_IDs):
    isoform_next_exons = {}
    num_segments_per_isoform = {}
    max_num_baits_per_isoform = {}

    for isoform_ID in target_isoform_IDs:
        isoform_next_exons[isoform_ID] = []
        num_segments_per_isoform[isoform_ID] = 0
        
    for segment in all_segments:
        for isoform_ID, exon_num, is_last_exon in segment.getExonNumberForEachIsoform(consider_only_targets=True):
            isoform_next_exons[isoform_ID].append(exon_num)
            num_segments_per_isoform[isoform_ID] += 1

    for isoform_ID in isoform_next_exons.keys():
        if (len(isoform_next_exons[isoform_ID]) == 0):
            print("WARNING: zero segments for %s" % isoform_ID, file=sys.stderr, flush=True)
        else:
            isoform_next_exons[isoform_ID].sort(reverse=True)

    for isoform_ID in target_isoform_IDs:
        max_num_baits_per_isoform[isoform_ID] = 0 if (len(isoform_next_exons[isoform_ID]) == 0) else max(30, max(isoform_next_exons[isoform_ID]))

    return (isoform_next_exons, num_segments_per_isoform, max_num_baits_per_isoform)


def discardLastExonSegments(segments, num_nonlast_segments_per_isoform, min_num_segments_per_isoform):
    '''Preserves ordering of segments'''
    remaining_segments = []
    for segment in segments:
        if (all(map(lambda isoform_ID: num_nonlast_segments_per_isoform[isoform_ID] >= min_num_segments_per_isoform, segment.getTargetIsoformIDs()))):
            # Then it is okay to discard this segment if it is a last exon segment
            if (not segment.isLastExonSegment()):
                remaining_segments.append(segment)
        else:
            remaining_segments.append(segment)
    return remaining_segments


def selectSegmentlevelSolution(all_segments, target_isoform_IDs, num_baits, bait_length, antitargets, genome_repeats_fasta, genome_ref, transcriptome_ref, tempdir, thermo_analysis):
    working_segments = set()
    selected_segments = []
    exons_w_segments_per_isoform = defaultdict(set)
    needed_segment_count_per_isoform = 1
    num_working_segments_for_isoform = {}
    curr_segment_count_per_isoform = defaultdict(int)
    last_exon_has_segment = defaultdict(lambda: False)

    unused_ordered_segments, segments_with_no_baits, segment_pass_fail_counts = evalAndPreorderSegments(all_segments, bait_length, genome_ref, transcriptome_ref,
                                                                                                        antitargets, genome_repeats_fasta, tempdir, thermo_analysis)

    no_segments, num_nonlast_segments_per_isoform = checkIsoformSegmentCoverage(target_isoform_IDs, unused_ordered_segments)
    if (len(no_segments) > 0):
        print("WARNING: After evaluating and preordering segments, %d isoforms have no segments" % len(no_segments), file=sys.stderr)
    unused_ordered_segments = discardLastExonSegments(unused_ordered_segments, num_nonlast_segments_per_isoform, min_num_segments_per_isoform=10)
    isoform_next_exons, num_segments_per_isoform, max_num_baits_per_isoform = setPerIsoformData(unused_ordered_segments, target_isoform_IDs)
    
    segment_locus_order = {}
    locus_segment_counter = defaultdict(int)
    for segment in unused_ordered_segments:
        segment_locus = segment.getTargetLocus()
        segment_locus_order[segment] = locus_segment_counter[segment_locus]
        locus_segment_counter[segment_locus] += 1

    for isoform_ID in target_isoform_IDs:
        num_working_segments_for_isoform[isoform_ID] = 0

    num_new_baits_used = 1e10
    prev_num_baits_used = 0
    iteration_number = 0
    max_frac_last_exon = 0.0
    while (len(selected_segments) < num_baits and num_new_baits_used > 0):
        prev_num_baits_used = len(selected_segments)
        
        # Get next exon needed for each isoform
        isoform_next_exon = {}
        for isoform_ID, next_exons in isoform_next_exons.items():
            if (len(next_exons) > 0):
                isoform_next_exon[isoform_ID] = next_exons.pop() # Next exons are in reverse order
            else:
                isoform_next_exon[isoform_ID] = 1e10
                
        # Get new segments that target the next exon that needs to be targeted for each isoform
        # First get the segments that overlap an exon that is needed to be targeted for an isoform.
        # unused_ordered_segments are preordered by their 5' proximity, mainly for the purpose of selecting the more 5' segment in last exon
        new_working_segments = set()
        retired_segments = set()
        for segment in filter(lambda s: s.getFracLastExonTargets() <= max_frac_last_exon, unused_ordered_segments): #  or s.targetsSimpleIsoform(5)
            if (any(map(lambda i: curr_segment_count_per_isoform[i] < max_num_baits_per_isoform[i], segment.getTargetIsoformIDs()))):
                # If any segment's isoforms are below their max segment limit...
                for isoform_ID, exon_num, is_last_exon in segment.getExonNumberForEachIsoform(consider_only_targets=True):
                    if (exon_num <= isoform_next_exon[isoform_ID] or num_working_segments_for_isoform[isoform_ID]==0):
                        new_working_segments.add(segment)
                        for isoform_ID in segment.getTargetIsoformIDs():
                            num_working_segments_for_isoform[isoform_ID] += 1
                        break
            else:
                retired_segments.add(segment)
                
        unused_ordered_segments = list(filter(lambda s: s not in new_working_segments and s not in retired_segments, unused_ordered_segments))
        working_segments.update(new_working_segments)

        iteration_number += 1
        print("Iteration %d with %d working segments and %d unused segments" % (iteration_number, len(working_segments), len(unused_ordered_segments)), file=sys.stderr)

        all_isoforms_got_new_segment = False
        chrom_of_best = None
        ranked_segments = []
        while (not all_isoforms_got_new_segment and len(selected_segments) < num_baits):
            print("%d segments used, needed number segments per isoform = %d" % (len(selected_segments), needed_segment_count_per_isoform), file=sys.stderr, flush=True)

            # For efficiency, only evaluate segments from affected chromosome. If affected_chrom is None, then eval all
            # If the exon that a segment targets has already been targeted by another segment, for all segment isoforms, then skip the segment
            retired_segments = set()
            for segment in filter(lambda s: chrom_of_best is None or (chrom_of_best==s.getChromosome() and s.fromSameLocus(best_segment_locus)), working_segments):
                if (any(map(lambda i: curr_segment_count_per_isoform[i] < max_num_baits_per_isoform[i], segment.getTargetIsoformIDs()))):
                    num_unfulfilled_isoforms = 0 # Number of target isoforms whose segment count deficient would be reduced by including this segment
                    # Implements check to determine if this segment overlaps its isoforms in their last exon and whether all their last exons have
                    # been targeted. If so remove this segment from further use.
                    for isoform_ID in segment.getTargetIsoformIDs():
                        if (num_segments_per_isoform[isoform_ID] < needed_segment_count_per_isoform):
                            if (curr_segment_count_per_isoform[isoform_ID] < num_segments_per_isoform[isoform_ID]):
                                num_unfulfilled_isoforms += 1
                        elif (curr_segment_count_per_isoform[isoform_ID] < needed_segment_count_per_isoform):
                            num_unfulfilled_isoforms += 1

                    overlaps_with_a_used_segment = any(map(lambda s:chrom_of_best==s.getChromosome() and s.overlapsWith(segment), selected_segments))
                    olap_sorting_criterion = 0 if (overlaps_with_a_used_segment) else 1
                    lowest_bait_Tm = segment.getLowestBaitTm()
                    avg_exon_num = segment.getAvgExonNum()
                    frac_last_exon = segment.getFracLastExonTargets()
                    perc_last_exon = int(round(100 * frac_last_exon))
                    ranked_segments.append( (segment, olap_sorting_criterion, -avg_exon_num, 100-perc_last_exon,
                                             num_unfulfilled_isoforms, -segment_locus_order[segment], -lowest_bait_Tm) )
                else:
                    retired_segments.add(segment)

            working_segments -= retired_segments
            
            # Select one best segment, remove it from further use, and increment the segment count for all of its target isoforms
            # Ideally want a segment that doesn't overlap and targets most unfulfilled isoforms. Disregard segments that don't fulfill counts for any isoforms
            ranked_segments.sort(key=itemgetter(1,2,3,4,5,6), reverse=True)
            filt_ranked_segments = list(filter(lambda s: s[4]>0, ranked_segments))

            # Find the highest ranking segment whose target isoforms haven't all met their segment limit
            best_segment = None
            while (len(filt_ranked_segments) > 0 and best_segment is None):
                if (any(map(lambda i: curr_segment_count_per_isoform[i] < max_num_baits_per_isoform[i], filt_ranked_segments[0][0].getTargetIsoformIDs()))):
                    best_segment = filt_ranked_segments[0][0]
                else:
                    segment_to_retire = filt_ranked_segments.pop(0)
                    working_segments.remove(segment_to_retire)
                    
            if (best_segment is not None):
                chrom_of_best = best_segment.getChromosome()
                best_segment_locus = best_segment.getTargetLocus()

                selected_segments.append(best_segment)
                working_segments.remove(best_segment)

                for isoform_ID, exon_num in best_segment.getTargetIsoformIDsAndExons():
                    if (best_segment.targetsLastExonOf(isoform_ID)):
                        last_exon_has_segment[isoform_ID] = True
                    exons_w_segments_per_isoform[isoform_ID].add(exon_num)
                    curr_segment_count_per_isoform[isoform_ID] += 1
                    num_working_segments_for_isoform[isoform_ID] -= 1

                ranked_segments = list(filter(lambda t:chrom_of_best != t[0].getChromosome() or not t[0].fromSameLocus(best_segment_locus), ranked_segments))
            else:
                all_isoforms_got_new_segment = True

        needed_segment_count_per_isoform += 1
        max_frac_last_exon += 0.2
        num_new_baits_used = len(selected_segments) - prev_num_baits_used

    print("Finished with %d baits" % len(selected_segments), file=sys.stderr, flush=True)
    if (len(selected_segments) != num_baits):
        print("WARNING: Only collected %d segments, needed %s" % (len(selected_segments), num_baits), file=sys.stderr, flush=True)

    selected_segments.sort(key=methodcaller("getChromStartStop"))

    print("Trimming candidate baits", file=sys.stderr)
    Ta = thermo_analysis.getTa()
    for counter, segment in enumerate(selected_segments, 1):
        segment.setLabel("S%d" % counter)

        num_before = segment.getNumCandidateBaits()
        segment.trimCandidateBaits(Ta)
        num_after = segment.getNumCandidateBaits()
        print("Segment %s\t%d -> %d baits" % (str(segment), num_before, num_after), file=sys.stderr)

    return selected_segments


def evalAndPreorderSegments(all_segments, bait_length, genome_ref, transcriptome_ref, antitargets, genome_repeats_fasta, tempdir, thermo_analysis):
    plus_strand_segments = []
    minus_strand_segments = []

    antitargets_fasta = tempfile.NamedTemporaryFile(mode='wt', suffix=".fa", dir=tempdir, delete=True)
    with open(antitargets, 'r') as ip:
        for line in ip:
            isoform_ID = line.strip()
            if (isoform_ID != ''):
                antitargets_fasta.write(">%s\n%s\n" % (isoform_ID, str(transcriptome_ref[isoform_ID])))
    antitargets_fasta.flush()
    makeblastdb_output = check_output(['makeblastdb', '-dbtype', 'nucl', '-in', antitargets_fasta.name], stderr=subprocess.STDOUT)

    # Evaluate the segments for their bait potential in parallel 
    segment_and_bait_data, segment_pass_fail_counts = evaluateSegmentsForBaitsInParallel(all_segments, genome_ref, antitargets_fasta.name, genome_repeats_fasta, tempdir, thermo_analysis)
    segments_with_no_baits = set()
    for segment, bait_data in segment_and_bait_data.items():
        if (len(bait_data) > 0):
            segment.setCandidateBaits(bait_data, *segment_pass_fail_counts[segment])
            if (segment.getStrand() == '+'):
                plus_strand_segments.append(segment)
            else:
                minus_strand_segments.append(segment)
        else:
            segments_with_no_baits.add(segment)

    antitargets_fasta.close()
    os.remove("%s.nhr" % antitargets_fasta.name)
    os.remove("%s.nin" % antitargets_fasta.name)
    os.remove("%s.nsq" % antitargets_fasta.name)

    superior_plus_strand_segments = removeInferiorSegments(plus_strand_segments, bait_length)
    superior_minus_strand_segments = removeInferiorSegments(minus_strand_segments, bait_length)
    num_removed = len(plus_strand_segments) + len(minus_strand_segments) - len(superior_plus_strand_segments) - len(superior_minus_strand_segments)
    print("INFO: removed %d inferior segments" % num_removed, file=sys.stderr)
    
    superior_plus_strand_segments.sort(key=methodcaller("getChromStartStop"))
    superior_minus_strand_segments.sort(key=methodcaller("getChromStartStop"), reverse=True)
    unused_ordered_segments = superior_plus_strand_segments + superior_minus_strand_segments

    return unused_ordered_segments, segments_with_no_baits, segment_pass_fail_counts


def evaluateSegmentsForBaitsInParallel(segments, genome_ref, antitargets_fasta_name, genome_repeats_fasta, tempdir, thermo_analysis):
    num_workers = min(19, len(segments))
    segment_and_bait_data = {}
    segment_pass_fail_counts = {}
    segments_to_eval = list(segments)
    
    # Do parallel processing
    jobs_queue = multiprocessing.Queue()
    results_queue = multiprocessing.Queue()
    seg_index = 0
    workers = []
    for _ in range(num_workers):
        proc = multiprocessing.Process(target=evaluateSegmentsForBaitsInParallel_Worker, args=(jobs_queue, results_queue, antitargets_fasta_name,
                                                                                               genome_repeats_fasta, tempdir, thermo_analysis))
        proc.start()
        workers.append(proc)
        segment_baits = segments_to_eval[seg_index].getAllPossibleBaits(genome_ref)
        jobs_queue.put( (seg_index, segment_baits) )
        seg_index += 1

    while (len(segment_and_bait_data) < len(segments_to_eval)):
        results_seg_index, bait_data, eval_str, pass_fail_counts = results_queue.get()
        eval_str = "Segment %s, %s" % (str(segments_to_eval[results_seg_index]), eval_str)

        if (seg_index < len(segments_to_eval)):
            segment_baits = segments_to_eval[seg_index].getAllPossibleBaits(genome_ref)
            jobs_queue.put( (seg_index, segment_baits) )
            seg_index += 1

        print(eval_str, file=sys.stderr, flush=True)
        segment = segments_to_eval[results_seg_index]
        segment_and_bait_data[segment] = bait_data
        segment_pass_fail_counts[segment] = pass_fail_counts

    for proc in workers:
        proc.terminate()
    
    return (segment_and_bait_data, segment_pass_fail_counts)


def evaluateSegmentsForBaitsInParallel_Worker(jobs_queue, results_queue, antitargets_fasta_name, genome_repeats_fasta, tempdir, thermo_analysis):
    while True:
        seg_index, baits = jobs_queue.get(block=True)

        # Filter baits that have out-of-range Tm or stable intramolecular structure or self-self interaction
        # Record thermo values for others
        num_baits_start = len(baits)
        thermo_vals, num_pass, num_fail_Tm, num_fail_hairpin, num_fail_homodimer = evalBaitSequenceThermodynamics(baits, thermo_analysis)
        num_baits_after_thermo = len(baits)

        # Filter baits that could hybridize to antitargets
        if (len(baits) > 0):
            scanBaitsAgainstAntitargets(baits, antitargets_fasta_name, tempdir, max_identity=20, only_masked=False)
            num_baits_after_1st_antitarget_filt = len(baits)
        else:
            num_baits_after_1st_antitarget_filt = 0
            
        if (len(baits) > 0):
            scanBaitsAgainstAntitargets(baits, genome_repeats_fasta, tempdir, max_identity=20, only_masked=True)
        
        eval_str = "thermo eval reduced number of baits %d -> %d\t(%d Tm, %d hairpin, %d homodimer) ->\t%d after RNA antitargets ->\t%d after genome repeats antitargets" % \
                   (num_baits_start, num_baits_after_thermo, num_fail_Tm, num_fail_hairpin, num_fail_homodimer, num_baits_after_1st_antitarget_filt, len(baits))

        good_baits = []
        for bait_label, bait_seq in baits.items():
            chromosome, start, stop, strand_name = bait_label.split('-')
            bait_strand = '+' if (strand_name == "plus") else '-'
            Tm, frac_hairpin, frac_homodimer = thermo_vals[bait_label]
            good_baits.append( (chromosome, int(start), int(stop), bait_strand, bait_seq, Tm, frac_hairpin, frac_homodimer) )

        results_queue.put( (seg_index, good_baits, eval_str, (num_pass, num_fail_Tm, num_fail_hairpin, num_fail_homodimer)) )

    
def removeInferiorSegments(segments, bait_length):
    '''Remove segments that overlap a segment that targets a superset of isoforms.
    There should be no overlapping with a segment that has identical target isoforms.'''
    segments_to_remove = set()
    for s1, s2 in combinations(segments,2):
        if (s1.overlapsWith(s2)):
            s1_target_isoforms = s1.getTargetIsoformIDs()
            s2_target_isoforms = s2.getTargetIsoformIDs()
                
            if (s1_target_isoforms < s2_target_isoforms):
                if (s1.noncoverageBy(s2) < bait_length):
                    segments_to_remove.add(s1)
            elif (s2_target_isoforms < s1_target_isoforms):
                if (s2.noncoverageBy(s1) < bait_length):
                    segments_to_remove.add(s2)
                
    return list(filter(lambda s: s not in segments_to_remove, segments))


def evalBaitSequenceThermodynamics(baits, thermo_analysis):
    thermo_vals = {}
        
    min_Tm, opt_Tm, max_Tm, Ta = thermo_analysis.getTemperatures()

    num_pass, num_fail_Tm, num_fail_hairpin, num_fail_homodimer = (0,0,0,0)
    baits_to_remove = set()
    for bait_label, bait_seq in baits.items():
        bait_seq = bait_seq.upper()
        Tm = thermo_analysis.calcTm(bait_seq)
        #print("-----\n%s\t%3.1f" % (bait_label, Tm), file=sys.stderr, flush=True)
        if (min_Tm <= Tm <= max_Tm):
            frac_hairpin = thermo_analysis.calcFracSelfStructure(bait_seq)
            #print(bait_seq, file=sys.stderr, flush=True)
            if (frac_hairpin < 0.0005):
                frac_homodimer = thermo_analysis.calcFracHomodimer(bait_seq)
                if (frac_homodimer < 0.75):
                    frac_hairpin = round(frac_hairpin,5)
                    frac_homodimer = round(frac_homodimer,4)
                    thermo_vals[bait_label] = (int(round(Tm,0)), frac_hairpin, frac_homodimer)
                    num_pass += 1
                    #print("\tfrac_hairpin = %6.5f\tfrac_homodimer = %5.4f" % (frac_hairpin, frac_homodimer), file=sys.stderr, flush=True)
                else:
                    #print("\tfrac_homodimer = %6.5f\tFAIL HOMODIMER" % frac_homodimer, file=sys.stderr, flush=True)
                    baits_to_remove.add(bait_label)
                    num_fail_homodimer += 1
            else:
                #print("\tfrac_hairpin = %6.5f\tFAIL HAIRPIN" % frac_hairpin, file=sys.stderr, flush=True)
                baits_to_remove.add(bait_label)
                num_fail_hairpin += 1
        else:
            baits_to_remove.add(bait_label)
            num_fail_Tm += 1

    for bait_label in baits_to_remove:
        del baits[bait_label]
            
    return (thermo_vals, num_pass, num_fail_Tm, num_fail_hairpin, num_fail_homodimer)
    

def scanBaitsAgainstAntitargets(baits, antitargets_fasta_name, tempdir, max_identity, only_masked):
    baits_fasta = tempfile.NamedTemporaryFile(mode='w+t', suffix=".fa", dir=tempdir, delete=True)
    num_seqs = 0
    if (only_masked):
        for bait_label, bait_seq in baits.items():
            if (sum([bait_seq[i].islower() for i in range(len(bait_seq))]) > 0):
                num_seqs += 1
                baits_fasta.write(">%s\n%s\n" % (bait_label, bait_seq))
    else:
        num_seqs += 1
        for bait_label, bait_seq in baits.items():
            baits_fasta.write(">%s\n%s\n" % (bait_label, bait_seq))
    baits_fasta.flush()
            
    if (num_seqs > 0):
        blastn_cmd = "blastn -task blastn-short -query %s -db %s -out - -max_target_seqs 1000 -num_threads 2 -strand minus -dust no \
                     -word_size 13 -evalue 1000 -qcov_hsp_perc 20 -perc_identity 80 -gapopen 2 -gapextend 2 -reward 1 -penalty -1 \
                     -outfmt '6 qseqid qlen length qend pident nident mismatch gaps qseq sseq'" % (baits_fasta.name, antitargets_fasta_name)
    
        blastn_process = Popen(blastn_cmd, shell=True, stdout=subprocess.PIPE)
        line = blastn_process.stdout.readline().decode()

        ids_of_matching_baits = set()
        while (line != '' or blastn_process.poll() is None):
            if (line != ''):
                qseqid, qlen, length, qend, pident, nident, mismatch, gaps, qseq, sseq = line.strip().split('\t')
                nident = int(nident)
                if (int(nident)-int(mismatch)-int(gaps) >= max_identity):
                    ids_of_matching_baits.add(qseqid)

            line = blastn_process.stdout.readline().decode()

        rc = blastn_process.poll()
        assert (rc == 0), "blastn returned with exit code %d" % rc

        for bait_label in ids_of_matching_baits:
            del baits[bait_label]

    baits_fasta.close()


def countBaitsOfftargets(baits_fasta_name, transcriptome_fasta, target_isoform_IDs, thermo_analysis):
    blastn_cmd = "blastn -task blastn-short -query %s -db %s -out - -max_target_seqs 50000 -num_threads 1 -strand minus -dust no \
                  -word_size 9 -evalue 1000 -qcov_hsp_perc 20 -perc_identity 60 -gapopen 2 -gapextend 2 -reward 1 -penalty -1 \
                  -outfmt '6 qseqid sseqid qlen length qend pident nident mismatch gaps qseq sseq'" % (baits_fasta_name, transcriptome_fasta)
    
    blastn_process = Popen(blastn_cmd, shell=True, stdout=subprocess.PIPE)
    line = blastn_process.stdout.readline().decode()

    all_hyb_frac_duplexed = defaultdict(list)
    while (line != '' or blastn_process.poll() is None):
        if (line != ''):
            qseqid, sseqid, qlen, length, qend, pident, nident, mismatch, gaps, qseq, sseq = line.strip().split('\t')
            nident = int(nident)
            if (nident >= 18 and sseqid not in target_isoform_IDs):
                frac_duplexed = thermo_analysis.calcFracHeteroDimer(qseq, sseq[::-1].translate(DNA_complement_table))
                frac_ident = float(nident)/float(qlen)
                if (frac_duplexed >= 0.01):
                    all_hyb_frac_duplexed[qseqid].append( (frac_duplexed,frac_ident) )

        line = blastn_process.stdout.readline().decode()

    rc = blastn_process.poll()
    assert (rc == 0), "blastn returned with exit code %d" % rc

    return all_hyb_frac_duplexed
        

def countBaitCrossHybs_DEPRACATED(segmentlevel_solution, tempdir):
    baits_fasta = tempfile.NamedTemporaryFile(mode='w+t', suffix=".fa", dir=tempdir, delete=True)
    for segment in segmentlevel_solution:
        for bait_label, bait_seq in segment.getOrderLabeledBaitSequences():
            baits_fasta.write(">%s\n%s\n" % (bait_label, bait_seq))
    baits_fasta.flush()
    makeblastdb_output = check_output(['makeblastdb', '-dbtype', 'nucl', '-in', baits_fasta.name], stderr=subprocess.STDOUT)

    blastn_cmd = "blastn -task blastn-short -query %s -db %s -out - -max_target_seqs 50000 -num_threads 2 -strand minus -dust no \
                  -word_size 9 -evalue 1000 -qcov_hsp_perc 20 -perc_identity 60 -gapopen 2 -gapextend 2 -reward 1 -penalty -1 \
                  -outfmt '6 qseqid sseqid qlen length qend pident nident mismatch gaps qseq sseq'" % (baits_fasta.name, baits_fasta.name)

    blastn_process = Popen(blastn_cmd, shell=True, stdout=subprocess.PIPE)
    line = blastn_process.stdout.readline().decode()

    crosshybs_penalties = {}
    while (line != '' or blastn_process.poll() is None):
        if (line != ''):
            qseqid, sseqid, qlen, length, qend, pident, nident, mismatch, gaps, qseq, sseq = line.strip().split('\t')
            nident = int(nident)
            if (nident >= 10):
                qseqid_segment, qseqid_baitnum = qseqid.split('_')
                sseqid_segment, sseqid_baitnum = sseqid.split('_')
                if (qseqid_segment != sseqid_segment):
                    crosshybs_penalties[((qseqid_segment, int(qseqid_baitnum)), (sseqid_segment, int(sseqid_baitnum)))] = nident # (qseqid,sseqid)
                    crosshybs_penalties[((sseqid_segment, int(sseqid_baitnum)), (qseqid_segment, int(qseqid_baitnum)))] = nident # (sseqid,qseqid)

        line = blastn_process.stdout.readline().decode()

    rc = blastn_process.poll()
    assert (rc == 0), "blastn returned with exit code %d" % rc

    baits_fasta.close()
    os.remove("%s.nhr" % baits_fasta.name)
    os.remove("%s.nin" % baits_fasta.name)
    os.remove("%s.nsq" % baits_fasta.name)

    return crosshybs_penalties


def selectBestBaitlevelSolution(segmentlevel_solution, num_needed_baits, transcriptome_fasta, target_isoform_IDs, tempdir, thermo_analysis):
    '''Segments in segmentlevel_solution should be in same order that they were selected, so that the best segments are first. This is because
    secondary baits, if needed, are chosen in a greedy fashion from it. Also, baits in each segment must be in sorted in order of "goodness"'''
    baitlevel_solution = defaultdict(list)
    arg_tuples = []
    baits_fastas = []
    for indices in np.array_split(range(len(segmentlevel_solution)), 19):
        baits_fasta = tempfile.NamedTemporaryFile(mode='w+t', suffix=".fa", dir=tempdir, delete=True)
        baits_fastas.append(baits_fasta)
        for seg_index in indices:
            segment = segmentlevel_solution[seg_index]
            segment_baits = segment.getCandidateBaits()
            for bait_index, bait in enumerate(segment_baits):
                bait_seq = bait.getSequence()
                bait_seq_label = "%s.%s" % (seg_index, bait_index)
                baits_fasta.write(">%s\n%s\n" % (bait_seq_label, bait_seq))
        baits_fasta.flush()

        arg_tuples.append( (baits_fasta.name, transcriptome_fasta, target_isoform_IDs, thermo_analysis) )

    all_hyb_frac_duplexed = {}
    with Pool(processes=len(arg_tuples)) as pool:
        counting_results = pool.starmap_async(countBaitsOfftargets, arg_tuples)
        counting_results.wait()
        for D in counting_results.get():
            all_hyb_frac_duplexed.update(D)
        
    for baits_fasta in baits_fastas:
        baits_fasta.close()

    for seg_index, segment in enumerate(segmentlevel_solution):
        segment_baits = segment.getCandidateBaits()
        for bait_index, bait in enumerate(segment_baits):
            bait_seq_label = "%s.%s" % (seg_index, bait_index)
            frac_duplexed_vals = list(map(itemgetter(0), all_hyb_frac_duplexed[bait_seq_label]))
            frac_identity_vals = list(map(itemgetter(1), all_hyb_frac_duplexed[bait_seq_label]))
            bait.setNumOfftargetData(len(all_hyb_frac_duplexed[bait_seq_label]), sum(frac_duplexed_vals), sum(frac_identity_vals))

    overlapping_segments = {}
    for segment in segmentlevel_solution:
        overlapping_segments[segment] = [segment]
    for s1,s2 in combinations(segmentlevel_solution,2):
        if (s1.overlapsWith(s2)):
            overlapping_segments[s1].append(s2)
            overlapping_segments[s2].append(s1)

    # Get baits from segments until the target number of baits is reached
    # Avoid new Baits that overlap an already selected Bait
    Ta = thermo_analysis.getTa()
    baitlevel_solution_size = 0
    for segment in segmentlevel_solution:
        segment.rankCandidateBaits(Ta)
        existing_baits = list(chain.from_iterable(map(lambda s: baitlevel_solution[s], overlapping_segments[segment])))
        bait = segment.getAdditionalBait(existing_baits)
        if (bait is not None):
            src_segment_label = bait.getSegmentLabel()
            bait_label = "%s-1" % src_segment_label
            bait.setLabel(bait_label)
            baitlevel_solution[segment].append(bait)
            baitlevel_solution_size += 1

    # Get additional baits from existing segments if deficient in num needed
    i = 0
    while (baitlevel_solution_size < num_needed_baits):
        segment = segmentlevel_solution[i]
        if (segment.getFracLastExonTargets() < 0.75):
            existing_baits = list(chain.from_iterable(map(lambda s: baitlevel_solution[s], overlapping_segments[segment])))
            bait = segment.getAdditionalBait(existing_baits)
            if (bait is not None):
                src_segment_label = bait.getSegmentLabel()
                new_num_segment_baits = len(baitlevel_solution[segment]) + 1
                bait_label = "%s-%d" % (src_segment_label, new_num_segment_baits)
                bait.setLabel(bait_label)
                baitlevel_solution[segment].append(bait)
                print("Added additional bait for Segment %s, now has %d" % (str(segment), new_num_segment_baits), file=sys.stderr)
                baitlevel_solution_size += 1
        i = 0 if(i==len(segmentlevel_solution)-1) else i+1
            
    print("Finished bait-level solution with %d baits" % baitlevel_solution_size, file=sys.stderr)

    return list(chain.from_iterable(baitlevel_solution.values()))


def checkBaitOverlaps(baits):
    for bait1, bait2 in combinations(baits, 2):
        if (bait1.overlapsWith(bait2)):
            print("Bait %s overlaps with Bait %s\n" % (str(bait1), str(bait2)), file=sys.stderr)


# Confirm that all targets have at least one segment
def checkIsoformSegmentCoverage(target_isoform_IDs, segments):
    num_nonlast_segments_per_isoform = defaultdict(int)
    segments_per_isoform = defaultdict(int)
    for segment in segments:
        is_last_exon = segment.isLastExonSegment()
        for isoform_ID in segment.getTargetIsoformIDs():
            segments_per_isoform[isoform_ID] += 1
            if (not is_last_exon):
                num_nonlast_segments_per_isoform[isoform_ID] += 1
    no_segments = target_isoform_IDs - segments_per_isoform.keys()
    return (no_segments, num_nonlast_segments_per_isoform)


def histogramSegmentsPerTargetIsoform(segments, target_isoform_IDs):
    op = open("baits_per_target_hist.tsv", 'w')
    baits_per_isoform = defaultdict(int)
    for segment in segments:
        for isoform_ID in segment.getTargetIsoformIDs():
            baits_per_isoform[isoform_ID] += 1

    # Confirm that all targets have at least one bait
    no_baits = target_isoform_IDs - baits_per_isoform.keys()
    if (len(no_baits) > 0):
        print("WARNING: %d isoforms have no baits" % len(no_baits), file=sys.stderr)
        for isoform_ID in no_baits:
            print("\t%s" % isoform_ID, file=sys.stderr)
            
    C = Counter(baits_per_isoform.values())
    op.write("NumberOfBaits\tIsoformCount\n")
    num_bait_and_count = C.most_common()
    num_bait_and_count.sort(key=itemgetter(0))
    for num_baits, num_isoforms_w_num_baits in num_bait_and_count:
        op.write("%d\t%d\n" % (num_baits, num_isoforms_w_num_baits))

    baits_per_isoform_counts = list(baits_per_isoform.items())
    baits_per_isoform_counts.sort(key=itemgetter(1))
    op.write("\nBaitCount\tIsoformID\n")
    for isoformID, bait_count in baits_per_isoform_counts:
        op.write("%d\t%s\n" % (bait_count, isoformID))

    op.close()
    

def writeSegmentsBed(segments, output_segments_bed):
    op = open(output_segments_bed, 'w')
    op.write("track type=bed name=BaitSegments description=BaitSegments visibility=full\n")
    for segment in segments:
        bed_tuple = segment.getAsBED()
        op.write("%s\t%d\t%d\t%s\t%d\t%s\n" % bed_tuple)
    op.close()


def writeBaitsBed(baits, output_baits_bed):
    op = open(output_baits_bed, 'w')
    op.write("track type=bed name=Baits description=Baits visibility=full\n")
    for bait in baits:
        bed_tuple = bait.getAsBED()
        op.write("%s\t%d\t%d\t%s\t%d\t%s\n" % bed_tuple)
    op.close()


def writeBaitFasta(baits, output_baits_fasta):
    op = open(output_baits_fasta, 'w')
    for bait in baits:
        label = bait.getLabel()
        seq = bait.getSequence()
        op.write(">%s\n%s\n" % (label, seq))
    op.close()


def setThermodynamics():
    rxn_temp_C = 65
    monovalent_salt_conc_mM = 100
    divalent_salt_conc_mM = 0
    tris_conc_mM = 100
    input_bait_conc_nM = 30

    min_Tm = 76
    opt_Tm = 85
    max_Tm = 100
    Ta = 65

    thermo_analysis = RNAThermodynamics(rxn_temp_C, monovalent_salt_conc_mM, divalent_salt_conc_mM, tris_conc_mM, input_bait_conc_nM)
    thermo_analysis.setTemperatures(min_Tm, opt_Tm, max_Tm, Ta)

    return thermo_analysis
    

if (__name__ == "__main__"):
    num_baits, bait_length, transcriptome_fasta, cgdb_bed, genome_fasta, genome_repeats_fasta, tempdir_root, antitargets, \
        bait_segments_bed_gz, target_isoforms_file, output_segments_bed, output_baits_bed, output_baits_fasta = sys.argv[1:]
    num_baits = int(num_baits)
    bait_length = int(bait_length)
    
    tempdir = "%s/buildBaitSet_%d" % (tempdir_root, os.getpid())
    os.mkdir(tempdir)
    print("Temporary directory is %s" % tempdir, file=sys.stderr, flush=True)

    genome_ref = Fasta(genome_fasta, as_raw=True, sequence_always_upper=False)
    transcriptome_ref = Fasta(transcriptome_fasta, as_raw=True, sequence_always_upper=False)
    thermo_analysis = setThermodynamics()

    # TODO: Confirm 0/1 bed intervals are being correctly read, used internally, and written.

    target_isoform_IDs = readTargetIsoformIDs(target_isoforms_file)
    all_segments = getBaitSegments(bait_segments_bed_gz, target_isoform_IDs, bait_length)
    
    segmentlevel_solution = selectSegmentlevelSolution(all_segments, target_isoform_IDs, num_baits, bait_length, antitargets,
                                                       genome_repeats_fasta, genome_ref, transcriptome_ref, tempdir, thermo_analysis)

    #pickle.dump(segmentlevel_solution, open("segmentlevel_soln.pkl", "wb"))
    #segmentlevel_solution = pickle.load(open("segmentlevel_soln.pkl", "rb"))
    histogramSegmentsPerTargetIsoform(segmentlevel_solution, target_isoform_IDs)
    writeSegmentsBed(segmentlevel_solution, output_segments_bed)
    
    baitlevel_solution = selectBestBaitlevelSolution(segmentlevel_solution, num_baits, transcriptome_fasta, target_isoform_IDs, tempdir, thermo_analysis)
    #pickle.dump(baitlevel_solution, open("baitlevel_soln.pkl", "wb"))
    #baitlevel_solution = pickle.load(open("bait_soln.pkl", "rb"))

    writeBaitsBed(baitlevel_solution, output_baits_bed)
    writeBaitFasta(baitlevel_solution, output_baits_fasta)
    
    sys.exit(0)

