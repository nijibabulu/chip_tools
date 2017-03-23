#! /usr/bin/env python

import collections
import itertools
import csv
import docopt
import IntervalTree

__doc__ = '''
Usage:  associate_genes.py [options] PEAKFILE BEDFILE BLASTFILE [DESEQ ...]

Options:
    --max-distance=DIST         maximum distance to associate a gene with [default: 5000]
    --effective-width=WIDTH     width to consider when deciding if exonic [default: 20]
    --better-descrs=FILE        specify a 2-column tab-separated file with manual annotations
    --start-only                only consider the 5' end of genes when considering intergenic peaks
'''

args = docopt.docopt(__doc__)
max_distance = int(args['--max-distance'])
eff_width = int(args['--effective-width'])

peakfile = args['PEAKFILE']
bedfile = args['BEDFILE']

manual_descrs = {}
if args['--better-descrs'] is not None:
    manual_descrs = dict([l.strip().split('\t') for l in open(args['--better-descrs'])])

class Gene(object):
    def __init__(self,chrom,start,end,name,score,strand,cds_start,cds_end, 
                 item_rgb,block_count,block_sizes,block_starts,*args):
        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)
        self.name = name
        self.score = int(score)
        self.strand = strand
        self.cds_start = int(cds_start)
        self.cds_end = int(cds_end)
        self.item_rgb = int(item_rgb)
        self.block_count = int(block_count)
        self.block_sizes = [int(s) for s in block_sizes.split(',') if len(s)]
        self.block_starts = [int(s) for s in block_starts.split(',') if len(s)]
        self.blast = ''

        self.cds_exons = []

        for start,size in zip(self.block_starts,self.block_sizes):
            start += self.start
            end = start+size
            sc = sorted([start,end,self.cds_start,self.cds_end])
            if (sc[0] == start and sc[1] == end) or (
                sc[0] == self.cds_start and sc[1] == self.cds_end):
                continue
            self.cds_exons.append((max(self.cds_start,start),min(self.cds_end,end)))

        if strand == '+':
            self.fivep = self.start
            self.cds_fivep = self.cds_start
        else:
            self.fivep = self.end
            self.cds_fivep = self.cds_end


    # TODO: should allow for multiple overlapping genes
    def dist(self,coord,max_dist):
        if coord >= self.start and coord <= self.end:
            return 0
        if args['--start-only']:
            return min(abs(coord-self.cds_fivep),abs(coord-self.fivep))
        else:
            return min(abs(coord-self.start),abs(coord-self.end))

        #if self.strand == '+':
            #if (len(self.cds_exons) > 1 and coord >= self.cds_exons[1][0]) or (
                #len(self.cds_exons) == 1 and coord >= self.cds_exons[0][1]):
                #return max_dist+1
        #else:
            #if (len(self.cds_exons) > 1 and coord <= self.cds_exons[-2][1]) or (
                #len(self.cds_exons) == 1 and coord <= self.cds_exons[0][0]):
                #return max_dist+1
        #return min(abs(coord-self.cds_fivep),abs(coord-self.fivep))

    def __str__(self):
        blast = self.blast
        if '|' in self.blast:
            blast = blast[blast.rindex('|')+1:]
        return '%s\t%s' % (self.name,blast)

class Peak(object):
    def __init__(self,chrom,start,end,name,*args):
        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)
        self.center = self.start+(self.end-self.start)/2
        self.name = name
        
        self.rest = list(args)

        self.closest_gene = ''
        self.multiple = False
        self.exonic = False

bf = open(args['BLASTFILE'])
bf.next()
blast = {l.split('\t')[1]: l.split('\t')[4] for l in bf}

pf = open(args['PEAKFILE'])
new_header = pf.next().strip() + '\tMultiple\tExonic\tGene\tBlast Hit'
for deseq_group in args['DESEQ']:
    new_header += '\t' + '\t'.join([deseq_group+'_'+type for type in
                             ['log2FoldChange','padj']])
print new_header
peaks = collections.defaultdict(list)
for line in pf:
    l = line.split()
    peaks[l[0]].append(Peak(*l))

gf = open(args['BEDFILE'])
genes = collections.defaultdict(list)
intervals = collections.defaultdict(list)
for line in gf:
    l = line.split()
    if len(l) < 10:
        continue
    gene = Gene(*l)
    genes[l[0]].append(gene)

    if gene.name in blast:
        gene.blast = blast[gene.name]

    if gene.name in manual_descrs:
        gene.blast = manual_descrs[gene.name]

    intervals[l[0]] += [IntervalTree.Interval(s,e) for s,e in gene.cds_exons]
it = {c: IntervalTree.IntervalTree(i) for c,i in intervals.items()}

DeseqResult = collections.namedtuple('DeseqResult',[
    'name','baseMean','log2FoldChange','lfcSE','stat','pvalue','padj'])

deseqs = collections.defaultdict(dict)
for deseq_tab in args['DESEQ']:
    f = open(deseq_tab)
    r = csv.reader(f,delimiter='\t')
    r.next()
    for row in r:
        deseqs[deseq_tab][row[0]] = DeseqResult(*row)

nogene = 'NA\tNA'
for chrom,p in peaks.items():
    for peak in p:
        track = False
        closest = nogene
        multiple = False
        for gene in genes[chrom]:
            dist = gene.dist(peak.center,max_distance)

            #print gene.name,dist
            if dist < max_distance:
                if closest != nogene:
                    multiple = True
                if closest == nogene or closest.dist(peak.center,max_distance) > dist:
                    closest = gene

        exonic = chrom in it and len(it[chrom].search(peak.center-eff_width/2,
                                                      peak.center+eff_width/2)) > 0

        new_row = [peak.chrom,peak.start,peak.end,peak.name] + peak.rest + \
                [int(multiple),int(exonic),closest]
        for deseq_group in args['DESEQ']:
            if closest != nogene and closest.name in deseqs[deseq_group]:
                new_row += [deseqs[deseq_group][closest.name].log2FoldChange,
                            deseqs[deseq_group][closest.name].padj]
            else:
                new_row += ['','']

        print '\t'.join(str(x) for x in new_row)
