#! /usr/bin/env python

import os
import sys
import math
import csv
import collections
import docopt

import peakzilla_qnorm_mapq_patched as pz


__doc__ = '''
Usage: join_peaks.py [options] PEAKS CHIP INPUT [ (PEAKS CHIP INPUT) ... ]

This script finds peaks in common between multiple ChIP experiments determined
by peakzilla. For each ChIP experiment, input a PEAKS file as otuput by
peakzilla, and 2 BED files (CHIP and INPUT) as input to peakzilla.

This will output a table with 3 columns identifying the peaks (Chromosome,
Start, End, Name,'NPeaks','Spread','ChipSE','EnrichSE'). NPeaks signifies the 
number of peaks that were called among all the ChIP experiments, Spread is the
difference between the biggest and smallest ChIP peak, ChipSE and EnrichSE are
the standard error on the mean among the ChIP and Enrich values for the peaks.
For each experinent "X", information about the peaks are output: 'XPZName','XPZScore',
'XPZChip','XPZInput','XPZEnrich','XPZFDR','XChip','XInput','XEnrich','XMapq'.
All 'PZ' columns are the original output from peakzilla and the remaining
columns are re-calculated in this script (also output regardless of the presence
of a peak).

Options:
    --max-distance=DIST     maximum summit distance to join peaks [default: 10]
'''

args = docopt.docopt(__doc__)

#np.set_printoptions(precision=1,suppress=True)

def stddev(l):
    mean = sum(l)/float(len(l))
    variance = sum((x-mean)**2 for x in l)/(len(l)-1)
    return math.sqrt(variance)

def std_err(l):
    return stddev(l)/math.sqrt(len(l))

class Peak(object):
    def dist(self,other):
        if self.chrom == other.chrom:
            return abs(self.center-other.center)
        else:
            return -1
    def compute_fold_enrichment(self):
        self.computed_fold_enrichment = float(self.computed_chip
                                             )/self.computed_control

class SlavePeak(Peak):
    def __init__(self,set_name,center):
        self.name = 'Slave'
        self.set_name = set_name
        self.center = center

class PZPeak(Peak):
    def __init__(self,set_name,chrom,start,end,name,summit,score,chip,control,
                 fold_enrichment,distribution_score,fdr):
        self.set_name = set_name
        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)
        self.name = name
        self.center = int(summit)
        self.score = float(score)
        self.chip = float(chip)
        self.control = float(control)
        self.fold_enrichment = float(fold_enrichment)
        self.distribution_score = float(distribution_score)
        self.fdr = float(fdr)
    def width(self):
        return self.end-self.start+1

class JoinedPeak(Peak):
    WIDTH = 0
    HEADER = ['#Chromosome','Start','End','Name','NPeaks','Spread','ChipSE','EnrichSE']
    HEADER_TYPES = set()
    def __init__(self,pzpeak):
        self.chrom = pzpeak.chrom
        self.peaks = {}
        self.center = self.add(pzpeak) #pzpeak.center
    def can_add(self,pzpeak):
        return not pzpeak.set_name in self.peaks
    def add(self,pzpeak):
        self.HEADER_TYPES.add(pzpeak.set_name)
        self.peaks[pzpeak.set_name] = pzpeak
        return sum(p.center for p in self.peaks.values())/len(self.peaks)
    def name(self):
        return '%s_%d' % (self.chrom,self.center)
    @classmethod
    def header(cls):
        s = '\t'.join(cls.HEADER) + '\t'
        #'#Chromosome\tPosition\tNPeaks\tSpread\t'
        for htype in cls.HEADER_TYPES:
            s += '\t'.join(
                htype + '_' + x for x in [
                    'PZName','PZScore','PZChip','PZInput','PZEnrich','PZFDR','Chip','Input','Enrich','Mapq']
            ) + '\t'
        return s
    def __str__(self):
        s = ''
        called_peaks = 0
        peak_signals = []
        peak_enrichs = []
        for set_name,peak in self.peaks.items():
            if hasattr(peak,'score'):
                s += peak.name + '\t' + '\t'.join('%.2f' % x for x in
                               [peak.score,peak.chip,peak.control,peak.fold_enrichment,peak.fdr]) + '\t'
                called_peaks += 1
                #s += '%.1f\t%.1f\t%.1f\t%.1f\t' % (
                    #peak.score,peak.chip,peak.control,peak.fold_enrichment)
            else:
                s += 'NA\tNA\tNA\tNA\tNA\tNA\t'
            if hasattr(peak,'pzpeak'):
                s += '\t'.join('%.2f' % x for x in [
                    peak.pzpeak.nrom_signal,peak.pzpeak.norm_background,peak.pzpeak.fold_enrichment,peak.pzpeak.mapq_score
                ]) + '\t'
                peak_signals.append(peak.pzpeak.nrom_signal)
                peak_enrichs.append(peak.pzpeak.fold_enrichment)
            else:
                s += 'NA\tNA\tNA\tNA\tNA\t'

                #peak.computed_chip,peak.computed_control,peak.computed_fold_enrichment
            #s += '%.1f\t%.1f\t%.1f\t' % (
                #peak.computed_chip,peak.computed_control,peak.computed_fold_enrichment)
            #s += '\t'.join([str(x) for x in
                            #[peak.score,peak.chip,peak.fold_enrichment]])

        try:
            if len(peak_signals):
                s = '\t'.join([self.chrom,str(self.center-self.WIDTH/2),str(self.center+self.WIDTH/2),
                               self.chrom+'_'+str(self.center),str(called_peaks)]) +\
                    '\t%.2f\t%.2f\t%.2f\t' % (
                        max(peak_signals)/(min(peak_signals) + sys.float_info.epsilon),
                        std_err(peak_signals), std_err(peak_enrichs),
                    ) + s
            else:
                s = '\t'.join([self.chrom,str(self.center),
                               self.chrom+'_'+str(self.center),str(called_peaks)]) +\
                    '\tNA\tNA\tNA\t' + s

        except:
            print max(peak_signals),min(peak_signals)
            raise

        return s

class PeakScorer(pz.PeakContainer):
    def __init__(self, ip_tags, control_tags, peak_size, plus_model, minus_model):
        self.ip_tags = ip_tags
        self.control_tags = control_tags
        self.peak_size = peak_size
        self.peak_shift = (peak_size - 1) / 2
        self.score_threshold = 10
        self.plus_model = plus_model
        self.minus_model = minus_model
        self.peaks = collections.defaultdict(list)
        self.peak_count = 0
        self.plus_window = collections.deque([])
        self.minus_window = collections.deque([])
        self.position = 0

    def fill_scores(self,chrom,libtype,scoretype):
        plus_tags = collections.deque(getattr(self,'%s_tags' % libtype).get_tags(chrom, '+'))
        plus_mapq = collections.deque(getattr(self,'%s_tags' % libtype).get_mapq(chrom, '+'))
        minus_tags = collections.deque(getattr(self,'%s_tags' % libtype).get_tags(chrom, '-'))
        minus_mapq = collections.deque(getattr(self,'%s_tags' % libtype).get_mapq(chrom, '-'))
        self.plus_window = collections.deque([])
        self.minus_window = collections.deque([])
        self.plus_mapq = collections.deque([])
        self.minus_mapq = collections.deque([])

        for peak in self.peaks[chrom]:
            # fill windows
            while plus_tags and plus_tags[0] <= (peak.position + self.peak_shift):
                self.plus_window.append(plus_tags.popleft())
                self.plus_mapq.append(plus_mapq.popleft())
            while minus_tags and minus_tags[0] <= (peak.position + self.peak_shift):
                self.minus_window.append(minus_tags.popleft())
                self.minus_mapq.append(minus_mapq.popleft())
            # get rid of old tags not fitting in the window any more
            while self.plus_window and self.plus_window[0] < (peak.position - self.peak_shift):
                self.plus_window.popleft()
                self.plus_mapq.popleft()
            while self.minus_window and self.minus_window[0] < (peak.position - self.peak_shift):
                self.minus_window.popleft()
                self.minus_mapq.popleft()
            # calculate normalized background level
            # add position to region if over threshold
            self.position = peak.position
            if libtype == 'ip':
                peak.mapq_score = float(sum(self.plus_mapq) + sum(self.minus_mapq)
                          )/max(1,(len(self.plus_mapq) + len(self.minus_mapq)))
                #if peak.name == 'Peak_12869': 
                    #print zip(self.plus_window,self.plus_mapq)
                    #print zip(self.minus_window,self.minus_mapq)
                    #print sum(self.plus_mapq) , sum(self.minus_mapq), len(self.plus_mapq) , len(self.minus_mapq)
                    #print peak.mapq_score
            setattr(peak,scoretype,self.calculate_score())

    def score_peaks(self,peak_dict):
        for chrom,peaks in peak_dict.items():
            for jp in peaks:
                jp.pzpeak = pz.Peak()
                jp.pzpeak.size = self.peak_size
                jp.pzpeak.shift = self.peak_shift
                jp.pzpeak.position = jp.center
                jp.pzpeak.name = jp.name
                self.peaks[chrom].append(jp.pzpeak)
                self.peak_count += 1
        for chrom,peaks in self.peaks.items():
            self.peaks[chrom] = sorted(self.peaks[chrom], 
                                       lambda a,b: cmp(a.position,b.position))
            self.fill_scores(chrom,'ip','score')
            self.fill_scores(chrom,'control','background')
            self.determine_fold_enrichment(chrom)
            self.determine_signal_over_background(chrom)

class FileSet(object):
    def __init__(self,peakfile,chipfile,controlfile):
        self.peakfile = peakfile
        self.chip_file = chipfile
        self.chip_tags = pz.TagContainer(store_mapq=True)
        self.chip_tags(chipfile,True)
        self.control_file = controlfile
        self.control_tags = pz.TagContainer(store_mapq=True)
        self.control_tags(controlfile,True)
 
        #print self.chip_tags, self.control_tags
    def get_file(self,type):
        return getattr(self, '%s_file' % type)
    def get_tagcount(self,type):
        return getattr(self, '%s_tags' % type)



maxdist = int(args['--max-distance'])
peaksets = {}
filesets = {}
for peakfile,chipfile,controlfile in zip(args['PEAKS'],args['CHIP'],args['INPUT']):
    set_name = os.path.basename(peakfile).split('.')[0]
    peaksets[set_name] = collections.defaultdict(list)
    filesets[set_name] = FileSet(peakfile,chipfile,controlfile)
    r = csv.reader(open(peakfile),delimiter='\t')
    r.next() # header
    '''
    #XXX: limit peaks
    maxpeaks = 20
    peakcounter = 0
    for row in r:
        if float(row[5]) >= 100 and float(row[8]) >= 10:
            peakcounter += 1
            if peakcounter > maxpeaks:
                break
            peaksets[set_name][row[0]].append(PZPeak(set_name,*row))
            '''
    for row in r:
        peaksets[set_name][row[0]].append(PZPeak(set_name,*row))

    JoinedPeak.WIDTH += peaksets[set_name].itervalues().next()[0].width()

JoinedPeak.WIDTH /= len(peaksets)

# find closest peak to each peak in the new set
#  make new peaks when there's no qualifying one
npeaks = 0
joined_peaks = collections.defaultdict(list)
for set_name,peakset in peaksets.items():
    for chrom,peaks in peakset.items():
        for peak in peaks:
            closest = None
            for jp in joined_peaks[chrom]:
                dist = jp.dist(peak)
                if dist >= 0 and dist <= maxdist:
                    if closest is None or closest.dist(peak) > dist:
                        closest = jp
            if closest is None or not closest.can_add(peak):
                npeaks += 1
                joined_peaks[chrom].append(JoinedPeak(peak))
            else:
                closest.add(peak)


plus_model,minus_model = pz.generate_ideal_model(JoinedPeak.WIDTH)

for set_name,fileset in filesets.items():
    scorer = PeakScorer(fileset.chip_tags,fileset.control_tags,
                        JoinedPeak.WIDTH,plus_model,minus_model)
    peaks_to_score = collections.defaultdict(list)
    for chrom,peaks in joined_peaks.items():
        for jp in peaks:
            if set_name not in jp.peaks:
                jp.peaks[set_name] = SlavePeak(set_name,jp.center)
            peaks_to_score[chrom].append(jp.peaks[set_name])
    scorer.score_peaks(peaks_to_score)

print JoinedPeak.header()
for chrom,peaks in joined_peaks.items():
    for peak in peaks:
        print peak

#plus_model,minus_model = pz.generate_ideal_model(JoinedPeak.WIDTH)
#def get_coverage(fileset,type,jp,pseudocount=0):
    #score = 0
    #start = max(0,jp.center-JoinedPeak.WIDTH/2)
    #for aln in fileset.get_file(type).fetch(
        #reference = jp.chrom, start = start,
        #end = jp.center+JoinedPeak.WIDTH/2):
        #if aln.is_reverse:
            #score += minus_model[aln.pos-start]
        #else:
            #score += plus_model[aln.pos-start]
    #return (score+pseudocount)*10.**6/fileset.get_tagcount(type)

    #return 10.**6*fileset.get_file(type).count(
        #reference = jp.chrom,
        #start = max(0,jp.center-JoinedPeak.WIDTH/2),
        #end = jp.center+JoinedPeak.WIDTH/2)/fileset.get_tagcount(type)

        #start = jp.center,
        #end = jp.center+1)

#matrix = np.zeros((npeaks,len(peaksets)*2))
#i = 0
#for chrom,peaks in joined_peaks.items():
    #for jp in peaks:
        #for j,set_name in enumerate(peaksets.keys()):
            #control_coverage = get_coverage(filesets[set_name],'control',jp,pseudocount=1)
            #chip_coverage = get_coverage(filesets[set_name],'chip',jp)
            #matrix[i][j] = float(chip_coverage)
            #matrix[i][j+len(peaksets)] = float(control_coverage)
        #i += 1

#quantile_normalize.quantile_norm(matrix)

#i = 0
#for chrom,peaks in joined_peaks.items():
    #for jp in peaks:
        #for j,set_name in enumerate(peaksets.keys()):
            #if set_name not in jp.peaks:
                #jp.peaks[set_name] = SlavePeak(
                    #set_name,matrix[i][j],matrix[i][j + len(peaksets)])
            #else:
                #jp.peaks[set_name].computed_chip = matrix[i][j]
                #jp.peaks[set_name].computed_control = matrix[i][j+len(peaksets)]
                #jp.peaks[set_name].compute_fold_enrichment()
        #print jp
        #i += 1


'''

i = 0
for chrom,peaks in joined_peaks.items():
    for jp in peaks:
        for j,set_name in enumerate(filesets.keys()):
            matrix[i][j] = float(jp.peaks[set_name].computed_chip)
            matrix[i][j+len(peaksets)] = float(jp.peaks[set_name].computed_control)
        i += 1
        '''
