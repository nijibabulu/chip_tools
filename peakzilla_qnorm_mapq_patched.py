#!/usr/bin/env python

# Copyright (c) Jonas Steinmann, 2010-2013
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License version 2 as 
# published by the Free Software Foundation.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see http://www.gnu.org/licenses/.

import sys
import csv
import os
import math
import itertools
from operator import add
from time import strftime, localtime
from collections import deque
from array import array
from optparse import OptionParser
from copy import copy
from math import exp, sqrt, pi

def main():
	# option parser
	usage = 'python peakzilla.py [OPTIONS] chip.bed control.bed > results.tsv'
	parser = OptionParser(usage=usage)
	
	parser.add_option("-m", "--model_peaks",\
	type = "int", dest="n_model_peaks", default='200',\
	help = "number of most highly enriched regions used to estimate peak size: default = 200")

	parser.add_option("-c", "--enrichment_cutoff",\
	type = "float", dest="enrichment_cutoff", default='2',\
	help = "minimum cutoff for fold enrichment: default = 2")
	
	parser.add_option("-s", "--score_cutoff",\
	type = "float", dest="score_cutoff", default='1',\
	help = "minimum cutoff for peak score: default = 1")
	
	parser.add_option("-f", "--fragment_size",\
	type = "int", dest="fragment_size", default='0',\
	help = "manually set fragment size in bp: default = estimate from data")
	
	parser.add_option("-e", "--gaussian",\
	action = "store_false", dest="gaussian", default=True,\
	help = "use empirical model estimate instead of gaussian")
	
	parser.add_option("-p", "--bedpe",\
	action = "store_true", dest="bedpe", default=False,\
	help = "input is paired end and in BEDPE format")
	
	parser.add_option("-l", "--log",\
	type = "str", dest="log", default='log.txt',\
	help = "directory/filename to store log file to: default = log.txt")
	
	parser.add_option("-n", "--negative",\
	action = "store_true", dest="negative", default=False,\
	help = "write negative peaks to negative_peaks.tsv")

	parser.add_option("-q", "--qnorm",\
	action = "store_true", dest="qnorm", default=False,\
	help = "employ quantile normalization instead of size normalization")

	# read arguments and options
	(options, args) = parser.parse_args()
	if len(args) > 2 or len(args) == 0:
		# return help message if argment number is incorrect
		parser.print_help()
		sys.exit(0)
	ip_file = args[0]
	has_control = False
	if len(args) == 2:
		control_file = args[1]
		has_control = True
	
	if options.qnorm and not has_control:
		sys.stderr.write('Quantile normalization requires a control!')
		sys.exit(1)
	
	# load tags
	write_log('Loading tags ...', options.log)
	ip_tags = TagContainer()
	ip_tags(ip_file, options.bedpe)
	control_tags = TagContainer()
	if has_control:
		control_tags(control_file, options.bedpe)
	
	# report tag number
	write_log('Tags in IP: %d' % ip_tags.tag_number, options.log)
	if has_control:
		write_log('Tags in control: %d' % control_tags.tag_number, options.log)

	# model peak size
	if not options.fragment_size:
		if not options.bedpe:
			peak_model = PeakShiftModel(ip_tags, options)
			peak_size = peak_model.peak_size
			write_log('Top %d paired peaks used to estimate peak size' % peak_model.peaks_for_size, options.log)
			write_log('Peak size is %d bp' % peak_size, options.log)
		else:
			write_log('Determine peak size from fragment size ...', options.log)
			peak_size =  2 * ip_tags.fragment_length + 1
			write_log('Peak size is %d bp' % peak_size, options.log)
	else:
		peak_size =  2 * options.fragment_size + 1
		write_log('Fragment size was set manually to %d bp' % options.fragment_size, options.log)
		write_log('Peak size is %d bp' % peak_size, options.log)
		
	# depending on option setting determine model using gaussian or empirically
	if options.gaussian:
		# estimate tag distirbution using gaussian function
		write_log('Estimating tag distribution using gaussian ...', options.log)
		model = generate_ideal_model(peak_size)
		plus_model = model[0]
		minus_model = model[1]
	else:
		# find peaks for modeling
		write_log('Finding peaks for modeling ...', options.log)
		ip_peaks = PeakContainer(ip_tags, control_tags, peak_size, peak_model.plus_model, peak_model.minus_model)
		
		# model tag distribution 
		write_log('Modeling tag distribution ...', options.log)
		plus_model = ip_peaks.model_tag_distribution()[0]
		minus_model = ip_peaks.model_tag_distribution()[1]

	# find peaks using model
	write_log('Finding peaks in ChIP sample ...', options.log)
	ip_peaks = PeakContainer(ip_tags, control_tags, peak_size, plus_model, minus_model)
	
	# find peaks using model in control sample
	if has_control:
		write_log('Finding peaks in control sample ...', options.log)
		control_peaks = PeakContainer(control_tags, ip_tags, peak_size, plus_model, minus_model)
	
	# calculate distribution scores and peak rank
	write_log('Calculating peak scores ...', options.log)
	ip_peaks.determine_distribution_scores(plus_model, minus_model)
	ip_peaks.determine_peak_ranks()
	if has_control:
		control_peaks.determine_distribution_scores(plus_model, minus_model)
		control_peaks.determine_peak_ranks()

	if options.qnorm:
		ip_peaks.qnorm_peaks(control_peaks)
		ip_peaks.determine_peak_ranks()
		control_peaks.determine_peak_ranks()

	# calculate FDR
	if has_control:
		write_log('Calculating FDR ...', options.log)
		ip_peaks.calculate_fdr(control_peaks.peaks)
	else:
		write_log('No FDR calculated as control sample is missing!', options.log)
	
	# write output in a BED like format
	ip_peaks.write_to_stdout(options)
	
	# write peaks in input to file
	if options.negative:
		write_log('Writing input peaks to negative_peaks.tsv', options.log)
		control_peaks.write_artifact_peaks('negative_peaks.tsv')
		
	# run finished successfully
	write_log('Done!', options.log)


def write_log(string, filename):
	# append line of string to log file of filename
	f = open(filename, 'a')
	f.write('%s %s\n' % (strftime("%H:%M:%S", localtime()), string))
	f.close()

def chisquare(f_obs, f_exp):
	# calculates a one-way chi-square for observed versus exprected frequencies
	chisq = 0
	df = len(f_obs)-1
	for i in range(len(f_obs)):
		if f_exp[i] > 0:
			chisq = chisq + (f_obs[i]-f_exp[i]) ** 2 / float(f_exp[i])
	# return chi-sqare value and associated p-value for f_obs == f_exp
	return chisq, chisqprob(chisq, df)

def chisqprob(chisq, df):
	# returns chisquare probability (works only for high degrees of freedom)
	if df < 30:
		raise ValueError('Function does not work for df < 30!')
	if chisq < 15:
		return 1.0
	a = 0.5 * chisq
	y = exp(-a)
	chisq = 0.5 * (df - 1.0)
	if df % 2 == 0:
		e = 1.0
		z = 1.0
	else:
		e = 1.0 / sqrt(pi) / sqrt(a)
		z = 0.5
	c = 0.0
	while (z <= chisq):
		e = e * (a/float(z))
		c = c + e
		z = z + 1.0
	return (c*y)

def median(numlist):
	# calculate median
    s = sorted(numlist)
    l = len(numlist)
    if l == 0:
		return float('nan')
    if l%2 == 0:
        return (s[l/2] + s[l/2-1]) / 2.0
    else:
        return float(s[l/2])

def convolve(signal, filter_width):
	# smooth signal with a flat scanning window of filter_width
	filter_width = float(filter_width)
	overhang = int((filter_width-1) / 2)
	window = deque([])
	result = []
	for i in signal + overhang * [0]:
		window.append(i)
		while len(window) > filter_width:
			window.popleft()
		result.append(sum(window)/ filter_width)
	return result[overhang:]

def generate_ideal_model(peaksize):
	# simulate ideal model
	stdev = peaksize / 5 # appears to fit well with empirical data
	mean_plus_model = (peaksize - 1) / 4
	mean_minus_model = peaksize - (peaksize - 1) / 4
	peak_positions = range(1,peaksize + 1)
	def gauss(x, mu, sigma):
		# gaussian function
		a = 1
		e = math.e
		return(a * e ** (- ((x - mu) ** 2) / (2.0 * sigma ** 2)))
	# generate plus model from gaussian function
	plus_model = []
	for i in peak_positions:
		plus_model.append(gauss(i, mean_plus_model, stdev))
	# generate minus model from gaussian function
	minus_model = []
	for i in peak_positions:
		minus_model.append(gauss(i, mean_minus_model, stdev))
	# normalize model
	norm_factor = (sum(plus_model) + sum(minus_model)) / peaksize
	for i in range(peaksize):
		plus_model[i] = plus_model[i]/norm_factor
		minus_model[i] = minus_model[i]/norm_factor
	return(plus_model, minus_model)


class TagContainer:
	# class for loading, storing and manipulating sequence tags
	def __init__(self,store_mapq=False):
		# intitialize an empty object
		self.tags = {}
		self.mapq = {}
		self.tag_number = 0
		self.fragment_length = 0
		self.store_mapq = store_mapq

	def __call__(self, bed_file, is_pe):
		# when called like a function load bed file and return self
		if not is_pe:
			self.load_bed(bed_file)
		else:
			self.load_bedpe(bed_file)
		self.sort_tags()
		return self
	
	def add_tag(self, chrom, strand, fiveprime, mapq):
		# add tag to dictionary
		if not chrom in self.tags:
			self.tags[chrom] = {}
			# store tags as an array of unsigned integers (4 bytes)
			self.tags[chrom]['+'] = array('i',[])
			self.tags[chrom]['-'] = array('i',[])
			if self.store_mapq:
				self.mapq[chrom] = {}
				self.mapq[chrom]['+'] = array('b',[])
				self.mapq[chrom]['-'] = array('b',[])
		self.tags[chrom][strand].append(fiveprime)
		if self.store_mapq:
			self.mapq[chrom][strand].append(mapq)
		# keep track of total number of tags added
		self.tag_number += 1
		
	def load_bed(self, bed_file):
		# parse a bed file and add contents to self
		for i in csv.reader(open(bed_file), delimiter='\t'):
			try:
				chrom = i[0]
				start = int(i[1])
				end = int(i[2])
				strand = i[5]
				mapq = int(i[4])
				# determine five prime end
				if strand == '+':
					fiveprime = start
				elif strand == '-':
					fiveprime = end
				# add tag to container
				self.add_tag(chrom, strand, fiveprime, mapq)
			except:
				sys.stderr.write("Input file is not in BED format!\n")
				raise
	
	def load_bedpe(self, bed_file):
		# parse a BEDPE file and add unique fragments to self
		uids = set()
		num_fragment_lengths = fragment_length_sum = 0
		try:
			for i in csv.reader(open(bed_file), delimiter='\t'):
				chrom = i[0]
				start = int(i[1])
				end = int(i[5])
				uid = i[0] + i[1] + i[5]
                                mapq = int(i[7])
				if uid not in uids:
					self.add_tag(chrom, '+', start, mapq)
					self.add_tag(chrom, '-', end, mapq)
					uids.add(uid)
					fragment_length_sum += (end-start)
					num_fragment_lengths += 1
			self.fragment_length = fragment_length_sum / num_fragment_lengths
		except:
			sys.stderr.write("Input file is not in BEDPE format!\n")
			raise 
			sys.exit(1)
		
	def sort_tags(self):
		# sort all tags while preserving the array data structure
		if self.store_mapq:
			for chrom in self.mapq.keys():
				self.mapq[chrom]['+'] = array('b', [x for (y,x) in sorted(
					itertools.izip(self.tags[chrom]['+'], self.mapq[chrom]['+']))])
				self.mapq[chrom]['-'] = array('b', [x for (y,x) in sorted(
					itertools.izip(self.tags[chrom]['-'], self.mapq[chrom]['-']))])
		for chrom in self.tags.keys():
			# as sorted returns conversion back to array is required
			self.tags[chrom]['+'] = array('i', sorted(self.tags[chrom]['+']))
			self.tags[chrom]['-'] = array('i', sorted(self.tags[chrom]['-']))
	
	def get_chrom_size(self, chrom):
		# chromosome size to consider for scanning of both strands
		if self.tags[chrom]['+'] and self.tags[chrom]['-']:
			chrom_size = self.tags[chrom]['-'][-1]
			return chrom_size
		else:
			return 0

	def genome_size(self):
		# genome size to consider for scanning of both strands
		genome_size = 0
		for chrom in self.tags.keys():
			genome_size += self.get_chrom_size(chrom)
		return genome_size
			
	def get_tags(self, chrom, strand):
		# return all tags for chromosome
		if chrom in self.tags:
			return self.tags[chrom][strand]
		else:
			return []
	
	def get_mapq(self, chrom, strand):
		if chrom in self.mapq:
			return self.mapq[chrom][strand]
		else:
			return []
			
	def get_tag_number(self, chrom, strand):
		# find out how many tags are mapped to a particular comsomome and strand
		return len(self.tags[chrom][strand])
		
	def get_chrom_names(self):
		# retrieve a sorted list of all chromosome names
		return self.tags.keys()
		

class PeakShiftModel:
	# class for modeling peak size and strand shift
	def __init__(self, tags, options):
		self.tags = tags
		self.window_size = 100
		self.tag_threshold = 10
		self.peak_shifts = []
		self.peak_shift = None
		self.peak_size = None
		self.peak_size_std = None
		self.plus_model = None
		self.minus_model = None
		self.peaks_incorporated = 0
		self.peaks_found = 0
		self.peaks_for_size = 0
		self.peaks = {}
		self.n_model_peaks = options.n_model_peaks
		self.build()

	def build(self):
		# for all chromosomes look for shifted peaks
		for chrom in self.tags.get_chrom_names():
			self.find_simple_peaks(chrom, '+')
			self.find_simple_peaks(chrom, '-')
		for chrom in self.peaks.keys():
			self.determine_shifts(self.peaks[chrom]['+'], self.peaks[chrom]['-'])
		# calculate the median peak_shift of top peaks
		self.peak_shifts = sorted(self.peak_shifts, reverse=True)
		top_shifts = []
		# check whether enough paired peaks have been found, if not reduce peaks for model 
		if self.peaks_incorporated < self.n_model_peaks:
			n_top_peaks = self.peaks_incorporated
		else:
			n_top_peaks = self.n_model_peaks
		self.peaks_for_size = n_top_peaks
		# ceate list of top peaks
		for i in range(n_top_peaks):
			top_shifts.append(self.peak_shifts[i][1])
		self.peak_shift = int(median(top_shifts))
		# peak size is 2 * shift size + 1
		self.peak_size = self.peak_shift * 2 + 1
		self.plus_model = [1] * self.peak_shift  + [0] * (self.peak_shift + 1)
		self.minus_model = [0] * (self.peak_shift + 1) + [1] * self.peak_shift

	def find_simple_peaks(self, chrom, strand):
		# return maxima of tag counts in regions with more tags than threshold
		tags = self.tags.get_tags(chrom, strand)
		window = deque([])
		peak_region = []
		# initiate dicts in case not present
		if not chrom in self.peaks:
			self.peaks[chrom] = {}
			self.peaks[chrom]['+'] = []
			self.peaks[chrom]['-'] = []
		for tag in tags:
			# if we already have same tag in window dont add
			if window and window[-1] == tag:
				pass
			# else add the new tag to window and reposition it
			else:
				window.append(tag)
				window_start = tag - self.window_size
			# get rid of all the tags not fitting in the window
			while window[0] < window_start:
				window.popleft()
			# identify maxima of enriched regions
			tag_count = len(window)
			if tag_count > self.tag_threshold:
				position = tag - self.window_size / 2
				peak_region.append((tag_count, position))
			elif peak_region:
				self.peaks[chrom][strand].append(max(peak_region))
				self.peaks_found += 1
				peak_region = []
	
	def determine_shifts(self, plus_peaks, minus_peaks):
		# looks for minus peaks upstream of plus peaks within fragment size
		minus_peaks = deque(minus_peaks)
		for plus_peak in plus_peaks:
			while minus_peaks:
				minus_peak = minus_peaks[0]
				if minus_peak[1] > plus_peak[1]:
					peak_shift = minus_peak[1] - plus_peak[1]
					if peak_shift < 500:
						self.peak_shifts.append((min(minus_peak[0], plus_peak[0]), peak_shift))
						self.peaks_incorporated += 1
					break
				minus_peaks.popleft()


class Peak:
	# class for peak related infromation and fuctions
	def __init__(self):
		self.size = 0
		self.shift = 0
		self.position = None
		self.tags = ([],[])
		self.mapq = ([],[])
		self.mapq_score = 0
		self.signal = 0
		self.score = 0
		self.background = 0
		self.nrom_signal = 0
		self.norm_background = 0
		self.fold_enrichment = 0
		self.plus_freq_dist = None
		self.minus_freq_dist = None
		self.fdr = 0
		self.dist_score = None
		self.survivals = 0

	def __len__(self):
		# for truth testing and number of tags
		return int(self.score)
	
	def calc_fold_enrichment(self, total_IP, total_control):
		if total_control == 0:
			# avoid zero division if no control sample is available
			total_control = 10**6
		self.fold_enrichment = (self.score / float(total_IP)) / ((self.background + 1) / float(total_control))

	def calc_signal_over_background(self, total_IP, total_control):
		# normalize by median score count in 4 kbp window around peak summit
		if total_control == 0:
			# avoid zero division if no control sample is available
			total_control = 10**6
		self.nrom_signal = self.score * 10**6 / float(total_IP)
		self.norm_background = (self.background + 1) * 10**6 / float(total_control)
		self.signal = (self.nrom_signal - self.norm_background)
	
	def _tag_distribution(self, plus_tags_unadj, minus_tags_unadj, filter_width):
		plus_tags = [tags - self.position + self.shift for tags in plus_tags_unadj]
		minus_tags = [tags - self.position + self.shift for tags in minus_tags_unadj]
		plus_dist = [0] * (self.size)
		minus_dist = [0] * (self.size)
		# project tags to list
		for i in plus_tags:
			plus_dist[i] += 1
		for i in minus_tags:
			minus_dist[i] += 1
		# use a flat moving window to improve S/N ratio
		# smooth by convolution of the singal with the window
		plus_freq_dist = convolve(plus_dist, filter_width)
		minus_freq_dist = convolve(minus_dist, filter_width)
		if (sum(plus_freq_dist) + sum(minus_freq_dist)) > 0:
			# normalize distribution height
			norm_factor = (sum(plus_freq_dist) + sum(minus_freq_dist)) / self.size
			for i in range(self.size):
				plus_freq_dist[i] = plus_freq_dist[i]/norm_factor
				minus_freq_dist[i] = minus_freq_dist[i]/norm_factor

		return (plus_freq_dist,minus_freq_dist)


	def determine_tag_distribution(self, filter_width):
		self.plus_freq_dist, self.minus_freq_dist = self._tag_distribution(
			self.tags[0], self.tags[1], filter_width)

	def determine_tag_distribution_old(self, filter_width):
		# return smoothed frequency distribution position of tags
		# normalize tags for position
		plus_tags = [tags - self.position + self.shift for tags in self.tags[0]]
		minus_tags = [tags - self.position + self.shift for tags in self.tags[1]]
		plus_dist = [0] * (self.size)
		minus_dist = [0] * (self.size)
		# project tags to list
		for i in plus_tags:
			plus_dist[i] += 1
		for i in minus_tags:
			minus_dist[i] += 1
		# use a flat moving window to improve S/N ratio
		# smooth by convolution of the singal with the window
		self.plus_freq_dist = convolve(plus_dist, filter_width)
		self.minus_freq_dist = convolve(minus_dist, filter_width)
		# normalize distribution height
		norm_factor = (sum(self.plus_freq_dist) + sum(self.minus_freq_dist)) / self.size
		for i in range(self.size):
			self.plus_freq_dist[i] = self.plus_freq_dist[i]/norm_factor
			self.minus_freq_dist[i] = self.minus_freq_dist[i]/norm_factor

	def calc_mapq_score(self):
		if len(self.mapq[0]):
			self.mapq_score = (sum(self.mapq[0]) + sum(self.mapq[1])
					)/(len(self.mapq[0]) + len(self.mapq[1]))

	def calc_distribution_score(self, plus_model, minus_model):
		# concatenate plus and minus distributions and models for testing
		model = plus_model[:self.shift] + minus_model[-self.shift:]
		freq_dist = self.plus_freq_dist[:self.shift] + self.minus_freq_dist[-self.shift:]
		# dist score is the p-value returned by the chi-square test
		#print 'calc dist' , freq_dist, model
		self.dist_score = chisquare(freq_dist, model)[1]
	
	def get_score(self):
		# final score is fold enrichment times goodness of fit to model
		return self.signal * self.dist_score


class PeakContainer:
	# a class to identify and classify potential peaks
	def __init__(self, ip_tags, control_tags, peak_size, plus_model, minus_model):
		self.ip_tags = ip_tags
		self.control_tags = control_tags
		self.peak_size = peak_size
		self.peak_shift = (peak_size - 1) / 2
		self.score_threshold = 10
		self.plus_model = plus_model
		self.minus_model = minus_model
		self.peaks = {}
		self.peak_count = 0
		self.plus_window = deque([])
		self.minus_window = deque([])
		self.position = 0
		self.peak2rank = {}
		self.build()

	def build(self):
		# perform main peak finding tasks
		for chrom in self.ip_tags.get_chrom_names():	
			self.find_peaks(chrom)
			self.measure_background(chrom)
			self.determine_fold_enrichment(chrom)
			self.determine_signal_over_background(chrom)

	def qnorm_peaks(self, control_peaks):
		all_peaks = list(itertools.chain(itertools.chain.from_iterable(self.peaks.values()),
			itertools.chain.from_iterable(control_peaks.peaks.values())))
		ip_ranks = sorted(range(len(all_peaks)),
						  lambda i,j: cmp(all_peaks[i].score,
									      all_peaks[j].score))
		control_ranks = sorted(range(len(all_peaks)),
							   lambda i,j: cmp(all_peaks[i].background,
											   all_peaks[j].background))
		#for i in range(len(all_peaks)):
			#print '%.1f\t%.1f\t%.1f\t%.1f' %  (
					#all_peaks[control_ranks[i]].score,
					#all_peaks[control_ranks[i]].background,
					#all_peaks[ip_ranks[i]].score,
					#all_peaks[ip_ranks[i]].background)
		#for chrom in self.ip_tags.get_chrom_names():
			#for peak in self.peaks[chrom]:
				#print peak.score

		#print '--'
		for i in range(len(all_peaks)):
			all_peaks[ip_ranks[i]].score = all_peaks[control_ranks[i]].background
			all_peaks[ip_ranks[i]].signal = all_peaks[ip_ranks[i]].score-all_peaks[ip_ranks[i]].background*10**6/self.control_tags.tag_number
			#all_peaks[ip_ranks[i]].nrom_signal = all_peaks[ip_ranks[i]].nrom_signal
			#all_peaks[ip_ranks[i]].norm_background = all_peaks[ip_ranks[i]].norm_background
			#print '%.1f\t%.1f\t%.1f\t%.1f' %  (
					#all_peaks[control_ranks[i]].score,
					#all_peaks[control_ranks[i]].background,
					#all_peaks[ip_ranks[i]].score,
					#all_peaks[ip_ranks[i]].background)
			
		#for chrom in self.ip_tags.get_chrom_names():
			#for peak in self.peaks[chrom]:
				#print peak.score,peak.get_score()
	def calculate_score(self):
		# calculate score
		score = 0
		tag_shift = self.peak_shift - self.position
		plus = self.plus_model
		minus = self.minus_model
		for tag in self.plus_window:
			score += plus[tag + tag_shift]
		for tag in self.minus_window:
			score += minus[tag + tag_shift]
		return score
	
	def find_peaks(self, chrom):
		# identify peak candidates on chromosome
		self.peaks[chrom] = []
		# convert tag arrays to deque for fast appending and popping
		plus_tags = deque(self.ip_tags.get_tags(chrom, '+'))
		plus_mapq = deque(self.ip_tags.get_mapq(chrom, '+'))
		minus_tags = deque(self.ip_tags.get_tags(chrom, '-'))
		minus_mapq = deque(self.ip_tags.get_mapq(chrom, '-'))
		# initalize windows and stuff
		score_buffer = deque([])
		peak_candidate = Peak()
		# reset scanning windows and position on chromosome
		self.plus_window = deque([])
		self.minus_window = deque([])
		self.plus_mapq = deque([])
		self.minus_mapq = deque([])
		self.position = 0
		while plus_tags and minus_tags:
			if plus_mapq:
				# fill windows
				while plus_tags and plus_tags[0] <= (self.position + self.peak_shift):
						self.plus_window.append(plus_tags.popleft())
						self.plus_mapq.append(plus_mapq.popleft())
				while minus_tags and minus_tags[0] <= (self.position + self.peak_shift):
						self.minus_window.append(minus_tags.popleft())
						self.minus_mapq.append(minus_mapq.popleft())
				# get rid of old tags not fitting in the window any more
				while self.plus_window and self.plus_window[0] < (self.position - self.peak_shift):
						self.plus_window.popleft()
						self.plus_mapq.popleft()
				while self.minus_window and self.minus_window[0] < (self.position - self.peak_shift):
						self.minus_window.popleft()
						self.minus_mapq.popleft()
			else:
				# fill windows
				while plus_tags and plus_tags[0] <= (self.position + self.peak_shift):
						self.plus_window.append(plus_tags.popleft())
				while minus_tags and minus_tags[0] <= (self.position + self.peak_shift):
						self.minus_window.append(minus_tags.popleft())
				# get rid of old tags not fitting in the window any more
				while self.plus_window and self.plus_window[0] < (self.position - self.peak_shift):
						self.plus_window.popleft()
				while self.minus_window and self.minus_window[0] < (self.position - self.peak_shift):
						self.minus_window.popleft()
			# if number of candidates found is high readjust threshold
			if self.peak_count > 35000:
				self.adjust_threshold()
			# add position to region if over threshold
			score = self.calculate_score()
			if score > self.score_threshold:
				# save all scores in buffer
				score_buffer.append(score)
				# get rid of old scores that are outside of the filter
				if len(score_buffer) > self.peak_size:
					score_buffer.popleft()
				# if current score is as big or bigger, consider it instead
				if score >= peak_candidate.score:
					peak_candidate.size = self.peak_size
					peak_candidate.shift = self.peak_shift
					peak_candidate.score = score
					peak_candidate.tags = (list(self.plus_window), list(self.minus_window))
					peak_candidate.mapq = (list(self.plus_mapq), list(self.minus_mapq))
					peak_candidate.survivals = 0
					if self.position >= 0:
						peak_candidate.position = self.position
					else:
						peak_candidate.position = 0
				# candidate survives if current score is smaller
				else:
					peak_candidate.survivals += 1
				# if candidate survives long enough do the expensive lookup
				if peak_candidate.survivals == self.peak_shift:
					# check score buffer to see whether candidate is a maximum
					# candidate is in the middle of the buffer now
					if peak_candidate.score == max(score_buffer):
						self.add_peak(peak_candidate, chrom)
					# consider current score next, reset survivals
					peak_candidate = Peak()
				# while in enriched region move windows in 1 bp steps
				self.position += 1
			else:
				# if we still have a candidate check whether its a max and add
				if peak_candidate:
					if peak_candidate.score == max(score_buffer):
						self.add_peak(peak_candidate, chrom)
					peak_candidate = Peak()
					score_buffer = deque([])
				# determine the next informative position in the genome and move there
				if plus_tags and minus_tags:
					distance_to_next = plus_tags[0] - self.position + 1
					self.position += distance_to_next

	def adjust_threshold(self):
		# allows for dynamic adjustment of peak calling threshold
		# restricts the number of candidate peaks to investigate to 30000
		peak_scores = []
		for chrom in self.peaks.keys():
			for peak in self.peaks[chrom]:
				peak_scores.append(peak.score)
		# set score to 30000th
		self.score_threshold = sorted(peak_scores)[-30000]
		# remove peaks below threshold
		for chrom in self.peaks.keys():
			self.peaks[chrom] = [peak for peak in self.peaks[chrom] if peak.score >= self.score_threshold]
		# recount peaks
		self.peak_count = 0
		for chrom in self.peaks.keys():
			self.peak_count += len(self.peaks[chrom])

	def add_peak(self, peak, chrom):
		# calculate tag distribution frequency and add peak to container
		peak.determine_tag_distribution(11) # 11 is emp optimal window width
		peak.calc_mapq_score()
		self.peaks[chrom].append(peak)
		self.peak_count += 1
	
	def measure_background(self, chrom):
		# for every peak check background level
		plus_tags = deque(self.control_tags.get_tags(chrom, '+'))
		minus_tags = deque(self.control_tags.get_tags(chrom, '-'))
		# convert to deque for super fast and efficient popleft
		self.plus_window = deque([])
		self.minus_window = deque([])
		for peak in self.peaks[chrom]:
			# fill windows
			while plus_tags and plus_tags[0] <= (peak.position + self.peak_shift):
				self.plus_window.append(plus_tags.popleft())
			while minus_tags and minus_tags[0] <= (peak.position + self.peak_shift):
				self.minus_window.append(minus_tags.popleft())
			# get rid of old tags not fitting in the window any more
			while self.plus_window and self.plus_window[0] < (peak.position - self.peak_shift):
				self.plus_window.popleft()
			while self.minus_window and self.minus_window[0] < (peak.position - self.peak_shift):
				self.minus_window.popleft()

			plus_dist, minus_dist = peak._tag_distribution(
				self.plus_window,self.minus_window,11) 

			bg_dist = plus_dist[:peak.shift] + minus_dist[-peak.shift:]
			ip_dist = peak.plus_freq_dist[:peak.shift] + peak.minus_freq_dist[-peak.shift:]

			peak.ip_dist_score = chisquare(bg_dist, ip_dist)[1]

			# calculate normalized background level
			# add position to region if over threshold
			self.position = peak.position
			peak.background = self.calculate_score()

	def determine_fold_enrichment(self, chrom):
		# for evey peak calculate fold enrichment
		for chrom in self.peaks.keys():
			for peak in self.peaks[chrom]:
				peak.calc_fold_enrichment(self.ip_tags.tag_number, self.control_tags.tag_number)
	
	def determine_signal_over_background(self, chrom):
		# for evey peak calculate fold enrichment
		for chrom in self.peaks.keys():
			for peak in self.peaks[chrom]:
				peak.calc_signal_over_background(self.ip_tags.tag_number, self.control_tags.tag_number)
	
	def model_tag_distribution(self):
		# use tags from top 200 peaks to build the distribution model
		ranked_peak_tags = []
		for chrom in self.peaks.keys():
			for peak in self.peaks[chrom]:
				ranked_peak_tags.append((peak.score, (peak.plus_freq_dist, peak.minus_freq_dist)))
		# find the tag count of the 200th largest peak
		tag_threshold = sorted(ranked_peak_tags)[-200][0]
		# add tags from highest peaks to the model
		top_tags = [i[1] for i in ranked_peak_tags if i[0] > tag_threshold]
		plus_model = [0] * self.peak_size
		minus_model = [0] * self.peak_size
		for tags in top_tags:
			plus_tags = tags[0]
			minus_tags = tags[1]
			plus_model = map(add, plus_tags, plus_model)
			minus_model = map(add, minus_tags, minus_model)
		# nromalize model for number of total peaks
		norm_factor = (sum(plus_model) + sum(minus_model)) / self.peak_size
		for i in range(self.peak_size):
			plus_model[i] = plus_model[i]/norm_factor
			minus_model[i] = minus_model[i]/norm_factor
		return (plus_model, minus_model)

	def determine_distribution_scores(self, plus_model, minus_model):
		# calculate distribution similarity of every peak to model
		for chrom in self.peaks.keys():
			for peak in self.peaks[chrom]:
				peak.calc_distribution_score(plus_model, minus_model)
	
	def determine_peak_ranks(self):
		# assign a unique rank to every peak
		all_peaks = []
		for chrom in self.peaks.keys():
			for peak in self.peaks[chrom]:
				all_peaks.append((peak.get_score(), peak.position))
		all_peaks = sorted(all_peaks, reverse=True)
		for i in range(len(all_peaks)):
			self.peak2rank[all_peaks[i]] = str(i + 1)

	def calculate_fdr(self, control_peaks):
		# create a dictionary to correlate scores with FDR values
		score2fdr = {}
		ip_scores = []
		control_scores = []
		for chrom in self.peaks.keys():
			for peak in self.peaks[chrom]:
				ip_scores.append(peak.get_score())
		for chrom in control_peaks.keys():	
			for peak in control_peaks[chrom]:
				control_scores.append(peak.get_score())
		ip_scores = deque(sorted(ip_scores, reverse=True))
		control_scores = deque(sorted(control_scores, reverse=True))
		# calculate FDR at all relevant cutoffs
		ip_count = float(0)
		control_count = float(0)
		while ip_scores:
			ip_score = ip_scores.popleft()
			ip_count += 1
			while control_scores and control_scores[0] >= ip_score:
				control_scores.popleft()
				control_count +=1
			ip_fdr = control_count / (control_count + ip_count) * 100
			score2fdr[str(ip_score)] = ip_fdr
		# add FDR to each peak object
		for chrom in self.peaks.keys():
			for peak in self.peaks[chrom]:
				peak.fdr = score2fdr[str(peak.get_score())]

	def write_to_stdout(self, options):
		# write results to stdout
		sys.stdout.write('#Chromosome\tStart\tEnd\tName\tSummit\tScore\tChIP\tControl\tFoldEnrichment\tDistributionScore\tIPDistributionScore\tFDR\n')
		peak_count = 0
		for chrom in sorted(self.peaks.keys()):
			for peak in self.peaks[chrom]:
				score = peak.get_score()
				if peak.fold_enrichment >= options.enrichment_cutoff and score >= options.score_cutoff:
					peak_count += 1
					summit = peak.position
					start = summit - self.peak_shift
					if start < 0:
						start = 0
					end = summit + self.peak_shift
					name = 'Peak_' + self.peak2rank[(score, summit)]
					signal = peak.nrom_signal
					background = peak.norm_background
					enrichment = peak.fold_enrichment
					dist_score = peak.dist_score
					ip_dist_score = peak.ip_dist_score
					fdr = peak.fdr
					output = (chrom, start, end, name, summit, score, signal, background, enrichment, dist_score, ip_dist_score, fdr)
					sys.stdout.write('%s\t%d\t%d\t%s\t%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n' % output)
		write_log('%d peaks detected' % peak_count, options.log)
		write_log('Enrichment cutoff: %.1f' % options.enrichment_cutoff, options.log)
		write_log('Score cutoff: %.1f' % options.score_cutoff, options.log)

	def write_artifact_peaks(self, control_file_name):
		# write peaks found in input to file
		f = open(control_file_name, 'w')
		f.write('#Chromosome\tStart\tEnd\tName\tSummit\tScore\tChIP\tControl\tFoldEnrichment\tDistributionScore\n')
		f.close()
		peak_count = 0
		for chrom in sorted(self.peaks.keys()):
			for peak in self.peaks[chrom]:
				score = peak.get_score()
				if score > 1:
					peak_count += 1
					summit = peak.position
					start = summit - self.peak_shift
					if start < 0:
						start = 0
					end = summit + self.peak_shift
					name = name = 'Peak_' + self.peak2rank[(score, summit)]
					signal = peak.nrom_signal
					background = peak.norm_background
					enrichment = peak.fold_enrichment
					dist_score = peak.dist_score
					output = (chrom, start, end, name, summit, score, signal, background, enrichment, dist_score)
					f = open(control_file_name, 'a')
					f.write('%s\t%d\t%d\t%s\t%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n' % output)
					f.close()


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("Program canceled by user!\n")
        sys.exit(0)
