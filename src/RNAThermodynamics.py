#!/usr/bin/env python3

from math import exp
import random
import re
import sys
import subprocess
import RNA
import pdb
from Bio.SeqUtils import MeltingTemp as mt


# Interesting paper, not really related though since it is about mRNA transcript concentration in cells:
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4325848/

# For thermodynamics of baits, see https://academic.oup.com/nar/article/31/17/4989/1192547

    # See http://biophysics.idtdna.com/HelpMelt.html
    #     https://blog-biosyn.com/2013/05/10/thermodynamics-of-dna-hybridization/
    #
    # From https://www.tbi.univie.ac.at/RNA/tutorial/#sec6,  "Finding potential binding sites with RNAduplex"
    # 
    # For bait self-structure:
    #
    # echo "CUACGGCGCGGCGCCCUUGGCGA" | /usr/local/src/ViennaRNA-2.4.3/src/bin/RNAfold --noPS -T 66
    # CUACGGCGCGGCGCCCUUGGCGA
    # ....(((.....)))........ ( -1.82)
    #
    #
    #
    # For self-other  bait vs bait evaluations:
    #      Prescreen with RNAduplex, then use RNAup on those that appear to have stable interactions
    # For self-self  bait vs bait evaluations:
    #      Use RNAup, since more realistic and number of baits is relatively small
    #
    #
    # /usr/local/src/ViennaRNA-2.4.3/src/bin/RNAduplex -T 66 < duplex.seqs
    # .(((((((((((((((((((((....(((.&.)))..))))))))))))))))))))).  14,43  :  10,37  (-24.46)
    #
    # /usr/local/src/ViennaRNA-2.4.3/src/bin/RNAup -b -w 50 -T 66 --no_header --no_output_file < duplex.seqs
    # (((((((((((((((((((((&)))))))))))))))))))))  15,35  :  16,36  (-9.33 = -22.63 + 1.02 + 12.28)

# These might be the two that I need for hairping and homo/heterodimer
# /usr/local/src/ViennaRNA-2.4.3/tests/python/test-RNA-mfe_eval.py
# seq1          = "CGCAGGGAUACCCGCG"
#    def test_mfe(self):
#        print "test_mfe"
#        fc= RNA.fold_compound(seq1)
#        (ss,mfe) = fc.mfe()
#        print ss, "[ %6.2f ]" % mfe
#        self.assertEqual(ss,struct1)

# seq1Dimer     = "CGCAGGGA&ACCCGCG"
#    def test_mfe_Dimer(self):
#        print "test_mfe_Dimer"
#        fc=RNA.fold_compound(seq1Dimer)
#        (ss,mfe) = fc.mfe_dimer()
#        print ss, "[ %6.2f ]" % mfe
#        self.assertEqual(ss,struct1Dimer)

# Below are from /usr/local/src/ViennaRNA-2.4.3/tests/python/test-RNA.py
#    def test_model_details_structure(self):
#        print "test_model_details_parameter_structure"

#        # check model details structure
#        md = RNA.md(); # default values
#        self.assertEqual(int(md.dangles), 2)
#        self.assertEqual(md.temperature, 37.0)

#        RNA.cvar.dangles     = 0
#        RNA.cvar.temperature = 40.1
#        md = RNA.md("global") # global values
#        self.assertEqual(int(md.dangles), 0)
#        self.assertEqual(int(RNA.cvar.dangles), 0)
#        self.assertEqual(md.temperature, 40.1)

#        #reset globals to default
#        RNA.cvar.dangles = 2
#        RNA.cvar.temperature= 37.0

#        # check parameter structures
#        params = RNA.param()
#        self.assertEqual(params.get_temperature(),37.0)
#        params = RNA.param(md)
#        self.assertEqual(params.get_temperature(),40.1)
#
#        pf_params = RNA.exp_param()
#        self.assertEqual(pf_params.get_temperature(),37.0)
#        pf_params = RNA.exp_param(md)
#        self.assertEqual(pf_params.get_temperature(),40.1)
#        md = None

#   def test_duplexfold(self):
#        print "testing duplexfold()"
#        duplex = RNA.duplexfold(seq1, seq2)
#        self.assertEqual(duplex.structure, ".(((.....(((.&)))))).")

#    def test_foldASequence(self):
#        print "test_foldASequence"
#        # new better interface
#        (struct, mfe) = RNA.fold(seq1)
#        self.assertEqual(struct,struct1)
#        # check energy
#        self.assertTrue(abs(RNA.energy_of_struct(seq1, struct1) - mfe) < 0.0001)



class RNAThermodynamics(object):
    def __init__(self, rxn_temp_C, monovalent_salt_conc_mM, divalent_salt_conc_mM, tris_conc_mM, input_bait_conc_nM):
        self.rxn_temp_C = rxn_temp_C
        self.rxn_temp_K = 273.15 + rxn_temp_C
        self.monovalent_salt_conc_mM = monovalent_salt_conc_mM
        self.divalent_salt_conc_mM = divalent_salt_conc_mM
        self.input_bait_conc_nM = input_bait_conc_nM
        self.input_bait_conc_M = input_bait_conc_nM * 0.000000001
        self.tris_conc_mM = tris_conc_mM
        
        RNA.cvar.temperature = rxn_temp_C
        md = RNA.md("global")
        md.noLP = 1

        self.RT = 1.98717 * self.rxn_temp_K  # 1.98717 is the Gas Constant in cal/(mol*K)

        # dG units are expected to be cal/mol and NOT kcal/mol
        #self.calc_frac_duplexed = lambda dG: (designParams.input_primer_conc_M * exp(-dG/RT))/(1 + designParams.input_primer_conc_M * exp(-dG/RT))

        #self.reRNAfold_ensemble_dG = re.compile(r"^\S+\s+\[\s*(.+)\]$")
        #self.reRNAup = re.compile(r"^\S+\s+\S+\s+.\s+\S+\s+\(\s*(.+)\s+\=\s+$")
        #self.reRNAduplex = re.compile(r"^\S+\s+\S+\s+.\s+\S+\s+\(\s*(.+)\)$")

        # Temperature constraints
        self.min_Tm = None
        self.opt_Tm = None
        self.max_Tm = None
        self.Ta = None


    def setTemperatures(self, min_Tm, opt_Tm, max_Tm, Ta):
        self.min_Tm = min_Tm
        self.opt_Tm = opt_Tm
        self.max_Tm = max_Tm
        self.Ta = Ta


    def getTemperatures(self):
        return (self.min_Tm, self.opt_Tm, self.max_Tm, self.Ta)
        

    def getTa(self):
        return self.Ta
    

    def calcTm(self, seq):
        #input_str = "%s\n" % seq
        #cmd = ["/usr/local/src/ViennaRNA-2.4.3/src/bin/RNAplex", "-p", "-l", "50", "-N", str(self.monovalent_salt_conc_M),
        #       "-M", str(self.divalent_salt_conc_M), "-K", "0", "-U", str(self.tris_conc_M), "-Q", "%10.9f" % self.input_bait_conc_M]
        #output = subprocess.check_output(cmd, universal_newlines=True, input=input_str)
        #Tm = float(output.strip().split("\n")[-1].split()[-2])
        Tm = mt.Tm_NN(seq=seq, dnac1=self.input_bait_conc_nM, nn_table=mt.RNA_NN2, Na=self.monovalent_salt_conc_mM,
                      Tris=self.tris_conc_mM, Mg=self.divalent_salt_conc_mM, saltcorr=2)
        return Tm

        
    def calcFracSelfStructure(self, seq):
        #input_str = "%s\n" % seq
        #cmd = ["/usr/local/src/ViennaRNA-2.4.3/src/bin/RNAfold", "-p", "--noPS", "-T",  str(self.rxn_temp_C)]
        #output = subprocess.check_output(cmd, universal_newlines=True, input=input_str)
        #mo = self.reRNAfold_ensemble_dG.match(output.strip().strip().split("\n")[-3])
        #dG_kcal = float(mo.group(1))

        #fc= RNA.fold_compound(seq)
        #dG_kcal = float(fc.pf()[1])
        # OR
        (struct, dG_kcal) = RNA.fold(seq)
        
        dG_cal = 1000.0 * dG_kcal
        K = exp(-dG_cal/self.RT)
        frac_unstructured = 1.0/(1.0 + self.input_bait_conc_M * K)
        frac_structured = 1.0 - frac_unstructured
        
        return frac_structured


    def calcFracHomodimer(self, seq):
        '''Assumes seq1 and seq2 have both been prescreened to have negligible self structure
        so that RNAduplex (which only considers intermolecular interactions and not the energy required
        for seq1/seq2 to melt and intramolecular interactions'''

        duplex = RNA.duplexfold(seq, seq)
        dG_kcal = float(duplex.energy)
        # OR
        #fc=RNA.fold_compound(dimer)
        #dG_kcal = float(fc.mfe_dimer()[-1])

        dG_cal = 1000.0 * dG_kcal
        K = exp(-dG_cal/self.RT)
        frac_unduplexed = 1.0/(1.0 + self.input_bait_conc_M * K)
        frac_duplexed = 1.0 - frac_unduplexed
        #print("dG = %s, frac_duplexed = %5.3f" % (dG, frac_duplexed), file=sys.stderr)

        return frac_duplexed


    def calcFracHeteroDimer(self, seq1, seq2):
        '''Assumes seq1 and seq2 have both been prescreened to have negligible self structure
        so that RNAduplex (which only considers intermolecular interactions and not the energy required
        for seq1/seq2 to melt and intramolecular interactions'''

        duplex = RNA.duplexfold(seq1, seq2)
        dG_kcal = float(duplex.energy)

        dG_cal = 1000.0 * dG_kcal
        
        K = exp(-dG_cal/self.RT)
        frac_unduplexed = 1.0/(1.0 + self.input_bait_conc_M * K)
        frac_duplexed = 1.0 - frac_unduplexed
        #print("dG = %s, frac_duplexed = %5.3f" % (dG, frac_duplexed), file=sys.stderr)

        return frac_duplexed
    

    def shuffleSeq(self, seq):
        l = list(seq)
        random.shuffle(l)
        return ''.join(l)

        
if (__name__ == "__main__"):

    seqs = ["CUACGGCGCGGCGCCCUUGGCGA", "CCUGAAAUACCAAAGCCAAUGAUGAGUGCAAUUGCUCCAGCCAGAAUAAUGAUGAUGCUA",
            "GAUUUGCCAGGGCGCGCAAUUGCACUCAUCAUUGGCAUGAUGAGUGCAAUUAUGAUGCUA"]

    rxn_temp_C = 95
    monovalent_salt_conc = 50
    divalent_salt_conc = 1.5
    input_bait_conc_nM = 20
    thermo_analysis = RNAThermodynamics(rxn_temp_C, monovalent_salt_conc, divalent_salt_conc, input_bait_conc_nM)
    #thermo_analysis.getFracIntramolcularStructured(seqs[0])
    #new_seq1 = thermo_analysis.shuffleSeq(seqs[1])
    #thermo_analysis.getFracSelfDimerized(new_seq1, seqs[2])
    pdb.set_trace()
    Tm = thermo_analysis.calcTm(seqs[2])

    sys.exit(0)
