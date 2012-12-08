#!/usr/bin/env python

import os
import re
import sys
import bz2
import gzip
import math

import argparse
import itertools
import logging as lg
import subprocess as sp

from peekorator import Peekorator

lg.basicConfig(level=lg.DEBUG)

nucleotides =  ['A', 'C', 'G', 'T']

class SampleParseError(Exception):
    pass

class SampleStats(object):

    def __init__(self, format_def, format_string, sample_name):
        self.format_def = format_def
        self.format_string = format_string
        self.sample = sample_name
        self.fields = self.format_def.split(':')
        self.values = self.format_string.split(':')
        self._nf = []
        self.info = {}
        
        if format_string == '.':
            pass
            #no data
        else:
            self.interpret()


    def interpret(self):
        try:
            assert(len(self.fields) == len(self.values))
        except AssertionError:
            lg.critical("parse problem '%s' - '%s'" % (self.fields, self.values))
            raise SampleParseError()
        
        for i, f in enumerate(self.fields):
            if f == 'DP':
                self.dp = int(self.values[i])
            elif f == 'NF':
                self._nf = map(int, self.values[i].split('|'))
            else:
                try:
                    self.info[f] = self.values[i]
                    self.__dict__[f] = self.values[i]
                except:
                    print 'error'
                    print 'key', f
                    print 'val', self.values[i]
                    raise
                
    def get_format_string(self):
        newval = []
        for i, fieldname in enumerate(self.fields):
            if fieldname == 'DP':
                newval.append(str(sum(self.nf)))
            elif fieldname == 'NF':
                newval.append("|".join(map(str, self.nf)))
            else:
                newval.append(self.values[i])
                
        return ":".join(newval)

    def set_nf(self, val):
        self._nf = val
        self.db = sum(val)

    def get_nf(self):
        return self._nf

    nf = property(get_nf, set_nf)
        
    def __str__(self):
        return "".join(['A', 'C', 'G', 'T'][i]
                       for i,x
                       in enumerate(self.nf) if x)

    def alleles(self):
        return [['A', 'C', 'G', 'T'][i]
                for i,x
                in enumerate(self.nf) if x]

    def get_major_allele(self, d=0.2):
        """
        return the major allele - but only if there is a fractional
        difference of d between the major & secondary allele        
        """
        alfs = [(x,i) for i,x in enumerate(self.allele_freqs())]
        alfs.sort()
        alfs.reverse()
        lg.warning('%s' % self.seq_logo())
        diff = alfs[0][0] - alfs[1][0]

        if diff < d:
            return None
        else:
            return ['A', 'C', 'G', 'T'][alfs[0][1]]

    def allele_freqs(self):
        if self.dp > 0:
            rv = [(float(x) / self.dp) for  x in self.nf]
        else:
            rv = [0.0, 0.0, 0.0, 0.0]
        assert(max(rv) <= 1)
        return rv

    def allele_freqs_sort(self):
        af = self.allele_freqs()
        af.sort()
        return af
    

    def get_absolute_distance(self, other):
        if str(self) != str(other):
            return 1
        else:
            return 0
        
    def get_euclidian_distance(self, other):
        taf= self.allele_freqs()
        oaf = other.allele_freqs()
        squadiff = 0
        for i in range(4):
            squadiff += (taf[i] - oaf[i]) ** 2

        return math.sqrt(squadiff)

    def get_max_allelefreq_distance(self, other):
        """
        Use the maximal allelefreq difference between self & other as
        a distance measure
        """
        return max([abs(x-y) for x,y in zip(self.allele_freqs(), other.allele_freqs())])
        
    def no_alleles(self):
        return len([x for x in self.nf if x > 0])
        
    def seq_logo(self, logo_len=10, dottify='#'):
        """
        Try to make a ascii seq logo
        """
        if self.dp == 0:
            return '-' + (' ' * (logo_len-1))
        
        gt = [{True : 1, False: 0}[x>0] for x in self.nf]
        ntleft = logo_len - sum(gt)
        gtplus = [int(round( ntleft * (float(max(x,1)-1) / (self.dp-sum(gt)))))
                  for x in self.nf ]
        gt = [sum(x) for x in zip(gt, gtplus)]
        
        if sum(gt) < logo_len:
            mx = gt.index(max(gt))
            gt[mx] += (logo_len - sum(gt))

        #assert(len(gt) == 4)
        #assert(sum(gt) == logo_len)
        
        if sum(gt) > logo_len:
            mx = gt.index(max(gt))
            gt[mx] -= (sum(gt) - logo_len)
            
        logo = ('A' * gt[0]) + \
               ('C' * gt[1]) + \
               ('G' * gt[2]) + \
               ('T' * gt[3])
        logo = logo.replace(dottify, '.')
        return logo

class RefSampleStats(SampleStats):

    def __init__(self, ref):
        format_str = {
            'a' : '50:50|00|00|00',
            'c' : '50:00|50|00|00',
            'g' : '50:00|00|50|00',
            't' : '50:00|00|00|50',
            'n' : '00:00|00|00|00'}.get(ref.lower(), '00:00|00|00|00')

        super(RefSampleStats, self).__init__("DP:NF", format_str, "Reference")

        
class Locus(object):
    def __init__(self, vcf_line, sample_names):
        
        self.vcf_line = vcf_line
        self.sample_names = sample_names

        ls = self.vcf_line.split("\t")
        self.ls = ls
        self.seq = ls[0]
        self.pos = int(ls[1])
        self.snpid = ls[2]
        self.ref = ls[3]
        self.refStats = RefSampleStats(self.ref)
        self.vars = ls[4].split(',')
        self.score = float(ls[5])
        self.filter = ls[6]

        #parse info fields
        self.info = {}
        for i in ls[7].split(';'):
            if not '=' in i:
                k = i
                v = True
            else:
                k,v = i.split('=', 1)

            self.info[k] = v
            if k == 'NS':
                self.ns = int(v)
            elif k == 'DP':
                pass
                #self._dp = int(v)
            elif k == 'AN':
                self.an = int(v)

        #parse format fields
        self.format_def = ls[8]
        assert(len(ls[9:]) == len(self.sample_names))
        self.samples = []
        for i, format_string in enumerate(ls[9:]):
            self.samples.append(SampleStats(self.format_def, format_string,
                                           sample_name = self.sample_names[i]))


    def set_dp(self, val):
        pass

    def get_dp(self):
        return sum([sum(x.nf) for x in self.samples])

    dp=property(get_dp, set_dp)
        

    def simple(self):
        """
        Return True if this locus has only one allele in all samples,
        and if that allele corresponds to the reference genome
        """
        for i in self.samples:
            if str(i) != self.ref: return False
        return True

    def has_sample_variation(self, measure='absolute', maxval = 0, internal_only=True):
        dm = self.dist_matrix(measure = measure)
        #print 'x' * 80
        #print self
        #print dm

        #get the submatrix that displays just the internal vaiation
        #the leftmost column & bottom row show variation outer
        im = [x[:-1] for x in dm[:-1]]

        if measure=='absolute':
            dist = sum([sum(x) for x in im])
            #print 'ab', dist

        elif measure=='allelefreq':
            #print self.format_dist_matrix(range(len(im)), im)
            dist = max([max(x) for x in im])

        return dist > maxval

        
    def format_dist_matrix(self, samples, matrix):
        """
        Nicely format a dist matrix
        """
        rv = ""
        rv += "%20s" % ''
        rv += " | ".join(["%10s" % x for x in samples])
        rv += "\n"
        for i, rw in enumerate(matrix):
            rv += "%20s" % samples[i]
            for j, cl  in enumerate(rw):
                if j > 0:
                    rv += " - "
                if i <= j:
                    rv += "%10.4f" % cl
                else:
                    rv += "%10s" % '*'
            rv += "\n"
        return rv

    def dist_matrix(self, measure='absolute'):
        
        nosamples = len(self.samples)
        dm = [ [ 0 for x in range(nosamples+1)]
               for x in range(nosamples+1) ]

        for i in range(nosamples+1):
            if i < nosamples: 
                s1 = self.samples[i]
            else: 
                s1 = self.refStats
        
            for j in range(nosamples + 1):
                if j < nosamples:
                    s2 = self.samples[j]
                else: 
                    s2 = self.refStats
                if i == j:
                    dm[i][j] = 0
                if i > j:
                    continue
                if measure == 'absolute':
                    dm[i][j] = s1.get_absolute_distance(s2)
                if measure == 'euclidian':
                    dm[i][j] = s1.get_euclidian_distance(s2)
                if measure == 'allelefreq':
                    dm[i][j] = s1.get_max_allelefreq_distance(s2)

        return dm

    def build_vcf_line(self):
        ivals = []
        for k in self.info.keys():
            if k == 'DP':
                ivals.append(self.dp)
            elif k == 'AN':
                ivals.append(self.an)
            else:
                ivals.append(self.info[k])
        
        infostr = ';'.join(
            [  '%s=%s' % (k,v) for (k,v) in zip(self.info.keys(), ivals)])
        
        groupstr = "\t".join([x.get_format_string() for x in self.samples])
        
        return "\t".join(
            map(str,
                [ self.seq,
                  self.pos,
                  self.snpid,
                  self.ref,
                  ','.join(self.vars),
                  self.score,
                  self.filter,
                  infostr,
                  self.format_def,
                  groupstr
                  ]))
                  
    def min_depth(self):
        return min([x.dp for x in self.samples])

    def max_no_alleles(self):
        return max([x.no_alleles() for x in self.samples])

    def min_no_alleles(self):
        return min([x.no_alleles() for x in self.samples])
        
    def nice_str(self, t='genotype', sep="\t"):

        if t == 'genotype':
            return sep.join(
                map(str, [
                    self.seq, self.pos, self.ref, self.dp]
                    + self.samples
                    ))
        elif t == 'nuccount':
            return sep.join(
                map(str, [
                    self.seq, self.pos, self.ref, self.dp] +
                    [sep.join(["%5d" % d for d in x.nf])for x in self.samples]
                    ))
        elif t == 'allfreq':
            return sep.join(
                map(str, [
                    self.seq, self.pos, self.ref, self.dp] +
                    [sep.join(["%5d" % int(d*100) for d in x.allele_freqs_sort()])
                     for x in self.samples]
                    ))

        elif t[:8] == 'logoplus':
            logo_len = int(t[8:])
            return sep.join(
                map(str,
                    [ self.seq, self.pos, self.ref, self.dp] +
                    [x.seq_logo(logo_len=logo_len, dottify=self.ref)
                     for x in self.samples] +
                    [ 'ACGT' ] +
                    ["|".join(["%2d" % d for d in x.nf])for x in self.samples]
            ))

        elif t[:4] == 'logo':
            logo_len = int(t[4:])
            return sep.join(
                map(str,
                    [ self.seq, self.pos, self.ref, self.dp] +
                    [x.seq_logo(logo_len=logo_len, dottify=self.ref)
                     for x in self.samples]
            ))

    def __str__(self):
        return "%-20s %6d %2s - %s" % (
            self.seq, self.pos, self.ref,
            " ".join(["%4s" % x for x in self.samples])
            )


class PSVCF(object):
    
    def __init__(self, filename, region=None):
        
        self.filename = filename
        self.region = False
        self.region_str = region
        self.region_has_pos = False
        self.been_in_region = False


        self.lg = lg.getLogger('PSVCF')

        if region:
            self.region = True
            _x = self.region_str.split(':')
            self.region_seq = _x[0]
            if len(_x) > 1:
                self.region_has_pos = True
                _y = _x[1].split('-')
                self.region_start = int(_y[0])
                self.region_end = int(_y[1])
            
        self.file_mode = 'file'
        if filename == '-':
            self.file_mode = 'stdin'
            self.F = Peekorator(sys.stdin)
        elif filename[-3:] == '.gz' and os.path.exists(filename + '.tbi') \
                and self.region:
            self.lg.debug('tabix mode')
            self.file_mode = 'tabix'
            #TABIX mode!
            cl = 'tabix -h %s %s' % (filename, self.region_str)
            self.lg.debug("running tabis: %s" % cl)
            self.TABIX_PROCESS = sp.Popen(cl.split(), stdout=sp.PIPE)
            self.F = Peekorator(self.TABIX_PROCESS.stdout)
        elif filename[-3:] == '.gz':
            self.lg.debug('Opening as bz2')
            self.F = Peekorator(gzip.open(self.filename))
        elif filename[-4:] == '.bz2':
            self.lg.debug('Opening as bz2')
            self.F = Peekorator(bz2.BZ2File(self.filename))
        else:
            self.lg.debug('normal file mode')
            self.F = Peekorator(open(self.filename))

        self.meta_header_lines = []
        self.header_line = ""

        #read header
        while True:
            line = self.F.peek
            if not line:
                raise StopIteration
            elif not line.strip():
                self.F.next()
                pass
            elif line[:2] == '##':         
                self.F.next()
                self.meta_header_lines.append(line.strip())
            elif line[:6] == '#CHROM':
                #header:
                self.interpret_header(line)
                break
            elif line[0] != '#':
                #bummer - no header??                
                self.create_dummy_header()
                break
            
        self.lg.info('Opening %s' % self.filename)
        
    def __iter__(self):
        return self

    def create_dummy_header(self):        
        self.header_line = "\t".join(
            ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'sample'])
        ls = self.header_line.split("\t")
        self.sample_names = ls[9:]


    def interpret_header(self, line): 
        self.header_line = line.strip()
        ls = self.header_line.split("\t")
        assert(ls[:9] == ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'])
        self.sample_names = ls[9:]

    def simple_names(self):
        """
        Return simplified sample names
        """
        def _simple_name(s):
            s = os.path.basename(s)
            s = s.replace('.vcf', '')
            s = s.replace('.bam', '')
            s = s.replace('.bz2', '')
            return s
        return [_simple_name(x) for x in self.sample_names]

    def next(self):
        while True:
            line = self.F.next()
            if not line:
                #EOF
                raise StopIteration
            line = line.strip()
            if not line:
                #empty line
                continue            
            if line[:2] == '##':
                self.meta_header_lines.append(line.strip())
                continue
            if line[:6] == '#CHROM':
                #header:
                self.interpret_header(line)
                continue

            try:
                loc = Locus(line, sample_names = self.sample_names)
            except SampleParseError:
                lg.critical("could not parse line")
                lg.critical(line)
                raise SampleParseError()

            if not self.region:
                #No region is specified - return regardless
                return loc

            #figure out if this locus falls within the specified region
            if loc.seq != self.region_seq:
                if self.been_in_region:
                    self.fin()

            if not self.region_has_pos:
                self.been_in_region = True
                return loc
            
            if loc.pos < self.region_start:
                continue
            if loc.pos > self.region_end:
                raise StopIteration

            self.been_in_region = True
            return loc

    def fin():
        """
        Finish iterations
        """
        if self.file_mode == 'file':
            self.F.close()
        elif self.file_mode == 'stdin':
            pass
        elif self.file_mode == 'tabix':
            self.F.close()
        raise StopIteration


    def add_meta(self, k, v):
        """
        Set a meta key/value pair - will be added to the header
        """
        self.meta_header_lines.append("##%s=%s" % (k,v))
        
    def build_header(self):
        """
        Return a new header
        """
        rv = self.meta_header_lines
        rv.append(self.header_line)
        return "\n".join(rv) + "\n"
 
