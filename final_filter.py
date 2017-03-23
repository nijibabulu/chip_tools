#! /usr/bin/env python

start_columns = [
    'Chromosome','Start','End','Name','NPeaks','Spread','ChipSE','EnrichSE']
recurring_columns = [
    'PZName','PZScore','PZChip','PZInput','PZEnrich','PZFDR','Chip','Input','Enrich','Mapq']
end_columns = [
    'Multiple','Exonic','Gene','Blast Hit']

final_start_columns = [ 'Chromosome','Start','End','Name','NPeaks','Spread','ChipSE','EnrichSE' ]
final_recurring_columns = [ 'PZScore','PZIPDist','Chip','Enrich','Mapq', ]
final_end_columns = [ 'Multiple','Gene','Blast Hit' ]
final_columns = final_start_columns  + final_recurring_columns + final_end_columns
import sys
import collections

if len(sys.argv) < 2:
    raise SystemExit, 'usage: final_filter.py TSV'

f = open(sys.argv[1])
header_line = f.next().strip()
header = collections.deque(header_line.split('\t'))
fields = []
pools = []

if len(set(start_columns + recurring_columns + end_columns)) < \
   len(start_columns)+len(recurring_columns)+len(end_columns):
    raise ValueError, 'repeated column'

def parse_columns(header,columns):
    for col,head_col in zip(columns,header):
        if col != head_col.strip('#'):
            raise ValueError, 'Non standard column: %s, should be %s' %(head_col,col)
        fields.append(header.popleft())

def parse_next_pool(header):
    for col,head_col in zip(recurring_columns,header):
        if not '_' in head_col:
            return False
        pool,colname = head_col.split('_')
        if pool not in pools:
            pools.append(pool)
        if colname != col:
            if col == recurring_columns[0]:
                return False
            else:
                print colname,col
                raise ValueError, 'Non standard recurring column: %s, should be %s_%s' % (head_col,pool,col)
        fields.append(header.popleft())
    return True

while parse_columns(header,start_columns):
    pass
while parse_next_pool(header):
    pass
while parse_columns(header,end_columns):
    pass

fields +=  header # extra columns
final_columns += header

#print '\t'.join(fields)
#sys.stdout.write('\t'.join(final_start_columns) + '\t')
#for pool in pools:
    #sys.stdout.write('\t'.join('%s_%s' % (pool,field) for field in final_recurring_columns) + '\t')
#sys.stdout.write('\t'.join(final_end_columns) + '\t')
#sys.stdout.write('\t'.join(header) + '\t') #extra columns
#print

score_min = 80
enrich_min = 10

def check_score(d):
    for pool in pools:
        score = d['%s_PZScore' % pool]
        enrich = d['%s_PZEnrich' % pool]
        if score != 'NA' and float(score) >= score_min and float(enrich) >= enrich_min:
            return True
    return False

def check_nonstrict_score(d):
    for k,v in d.items():
        if 'PZScore' in k and v != 'NA' and float(v) >= score_min:
            return True
        if 'PZEnrich' in k and v != 'NA' and float(v) >= enrich_min:
            return True
    return False

def strict_score(d):
    for k,v in d.items():
        if 'PZScore' in k and v != 'NA' and float(v) < score_min:
            return False
        if 'PZEnrich' in k and v != 'NA' and float(v) < enrich_min:
            return False
    return True

print header_line
for line in f:
    items = line.strip().strip('#').split('\t')
    d = collections.OrderedDict(zip(fields,items))
    if any([int(d['Exonic']) == 1, int(d['NPeaks']) < 2, not check_score(d), 
            d['Gene'] == 'NA' ]):
        continue
    print line,
    #for f,i in d.items():
        #for ff in final_columns:
            #if ff == f.strip('#').split('_')[-1]:
                #sys.stdout.write(i+'\t')
                #break
    #print



#for col,head_col in zip(end_columns,header):
  #if col != head_col:
    #raise ValueError, 'Non standard column: %s, should be %s' %(head_col,col)
  #header.popleft()


