#! /usr/bin/env python
###############################################################################
# parses pileup base string and returns the counts for all possible alleles 
# for each position
# reads input (mpileup output) from sys.stdin
###############################################################################

import os
import sys

class parseString(object):
    
    def __init__(self, ref, string):
        self.ref = ref.upper()
        self.string = string.upper()
        self.types = {'A':0,'G':0,'C':0,'T':0,'-':[],'*':0,'+':[],'X':[]}
        self.process()
        
    def process(self):
        # remove end of read character
        self.string = self.string.replace('$','')
        while self.string != '':
            if self.string[0] == '^':
                # skip two characters when encountering '^' as it indicates
                # a read start mark and the read mapping quality
                self.string = self.string[2:]
            elif self.string[0] == '*':
                self.types['*'] += 1
                # skip to next character
                self.string = self.string[1:]
            
            elif self.string[0] in ['.',',']:
                if (len(self.string)== 1) or (self.string[1] not in ['+','-']):
                    # a reference base
                    self.types[self.ref] += 1
                    self.string = self.string[1:]
                elif self.string[1] == '+': 
                    insersionLength = int(self.string[2])
                    insersionSeq = self.string[3:3+ insersionLength]
                    self.types['+'].append(insersionSeq)
                    self.string = self.string[3+insersionLength:]
                elif self.string[1] == '-':
                    deletionLength = int(self.string[2])
                    deletionSeq = self.string[3:3+deletionLength]
                    self.types['-'].append(deletionSeq)
                    self.string = self.string[3+deletionLength:]
                    
            elif self.types.has_key(self.string[0]) and\
                 ((len(self.string)==1) or (self.string[1] not in ['-','+'])):
                # one of the four bases
                self.types[self.string[0]] += 1
                self.string = self.string[1:]
            else:
                # unrecognized character
                # or a read that reports a substitition followed by an insersion/deletion
                self.types['X'].append(self.string[0])
                self.string = self.string[1:]
        return
    def __repr__(self):
        types = self.types
        return '\t'.join(map(str,[types['A'], types['C'], types['G'],types['T'],\
                                  types['*']]) +\
                         map(','.join, [types['-'],types['+'],types['X']]))
        

def main():
    print >>sys.stdout, "chrom\tpos\tref\tcov\tA\tC\tG\tT\t*\t-\t+\tX"
    for line in sys.stdin:
        toks = line.strip('\n').split('\t')
        ref = toks[2].upper()
        cov = toks[3]
        print >>sys.stdout, '\t'.join([toks[0], toks[1],ref, cov]) + '\t' + \
            parseString(ref, toks[4]).__repr__()

if __name__ == '__main__':
    main()
