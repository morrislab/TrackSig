#!/usr/bin/env python2
# coding=utf-8

# AUTHOR: Jeff Wintersinger, Roujia Li, Yulia Rubanova

#Get copy numbers from CNA files, get VAF from VCF files, correct VAF by
#multiplicity, re-sample the values from Beta distribution for more noise-free predictions

#from create_phylowgs_inputs import *
import argparse
# Requires PyVCF. To install: pip2 install pyvcf
import vcf
import csv
import random
from collections import defaultdict
import numpy as np
from scipy.stats import beta
cutoff = 10000
min_readdepth=0

half_with_highest_depth = False
sample_vaf_from_posterior = True

class ReadCountsUnavailableError(Exception):
  pass

class VariantParser(object):
  def __init__(self):
    # Child classes must give the following variables sensible values in
    # constructor so that list_variants() works subsequently.
    self._cnvs = None
    self._vcf_filename = None

  def list_variants(self):
    variants = self._filter(self._vcf_filename)
    variants_and_reads = []
    for variant in variants:
      try:
        ref_reads, total_reads = self._calc_read_counts(variant)
      except ReadCountsUnavailableError as exc:
        continue
      variants_and_reads.append((variant, ref_reads, total_reads))
    return variants_and_reads

  def _calc_read_counts(self, variant):
    raise Exception('Not implemented -- use child class')

  def _parse_vcf(self, vcf_filename):
    vcfr = vcf.Reader(filename=vcf_filename)
    records = []

    for variant in vcfr:
      variant.CHROM = variant.CHROM.upper()
      # Some VCF dialects prepend "chr", some don't. Remove the prefix to
      # standardize.
      if variant.CHROM.startswith('CHR'):
        variant.CHROM = variant.CHROM[3:]
      records.append(variant)
    return records

  def _does_variant_pass_filters(self, variant):
    if variant.FILTER is None:
      return True
    if len(variant.FILTER) > 0:
      # Variant failed one or more filters.
      return False
    return True

  def _filter(self, vcf_filename):
    variants = []

    all_variants = self._parse_vcf(vcf_filename)

    for variant in all_variants:
      if not is_good_chrom(variant.CHROM):
        continue
      if not self._does_variant_pass_filters(variant):
        continue
      variants.append(variant)
    return variants

  def _get_tumor_index(self, variant, tumor_sample=None):
    """Find the index of the tumor sample.

    Currently hardcodes tumour sample as the last column if name not specified.
    Might not always be true
    """
    if self._tumor_sample:
      tumor_is = [i for i, s in enumerate(variant.samples) if s.sample == tumor_sample]
      assert len(tumor_is) == 1, "Did not find tumor name %s in samples" % tumor_sample
      return tumor_is[0]
    else:
      # Don't make this -1, as some code assumes it will be >= 0.
      return len(variant.samples) - 1


class PcawgConsensusParser(VariantParser):
  def __init__(self, vcf_filename, tumor_sample=None):
    self._vcf_filename = vcf_filename
    self._tumor_sample = tumor_sample

  def _find_ref_and_variant_nt(self, variant):
    assert len(variant.REF) == len(variant.ALT) == 1
    return (str(variant.REF[0]), str(variant.ALT[0]))

  def _calc_read_counts(self, variant):
    if not ('t_alt_count' in variant.INFO and 't_ref_count' in variant.INFO):
      #raise ReadCountsUnavailableError()
      return(None, None)
    #assert len(variant.INFO['t_alt_count']) == len(variant.INFO['t_ref_count']) == 1

    alt_reads = variant.INFO['t_alt_count']
    ref_reads = variant.INFO['t_ref_count']

    if isinstance(alt_reads, list):
    	alt_reads = alt_reads[0]
	ref_reads = ref_reads[0]
    alt_reads = int(alt_reads)
    ref_reads = int(ref_reads)

    total_reads = alt_reads + ref_reads
    # Some variants havezero alt and ref reads.
    #if total_reads == 0:
    #  raise ReadCountsUnavailableError()
    return (ref_reads, total_reads)

  def _get_vaf(self, variant):
     if not ('VAF' in variant.INFO):
       return None

     vaf = variant.INFO['VAF']
     if isinstance(vaf, list):
       vaf = vaf[0]
     vaf = float(vaf)
     return (vaf)


class CnvParser(object):
	def __init__(self, cn_filename):
		self._cn_filename = cn_filename

	def parse(self):
		cn_regions = defaultdict(list)

		with open(self._cn_filename) as cnf:
			reader = csv.DictReader(cnf, delimiter='\t')
			for record in reader:
				chrom = record['chromosome'].upper()
				
				if record["total_cn"] != "NA":
					cn_regions[chrom].append(record)

			# Ensure CN regions are properly sorted, which we later rely on when
			# filtering out regions with multiple abnormal CN states.
			for chrom, regions in cn_regions.items():
				cn_regions[chrom] = sorted(regions, key = lambda r: r['start'])

			return cn_regions

class CnvFormatter(object):
	def __init__(self, cnv_confidence, read_depth, read_length, sampidxs):
		self._cnv_confidence = cnv_confidence
		self._read_depth = read_depth
		self._read_length = read_length
		self._sampidxs = sampidxs

	def _max_reads(self, sampidx):
		return 1e6 * self._read_depth[sampidx]

	def _find_overlapping_variants(self, chrom, cnv, variants):
		overlapping = []
		start = cnv['start']
		end = cnv['end']

		for variant in variants:
			if chrom.upper() == variant['chrom'].upper():
				if start <= variant['pos'] <= end:
					overlapping.append(variant['ssm_id'])
		return overlapping

	def _find_copy_number(self, variant, cnv):
		first = 0
		last = len(cnv)-1
		found = False

		copy_number = 2

		while first<=last and not found:
			midpoint = (first + last)//2
			
			start = int(float(cnv[midpoint]['start']))
			end = int(float(cnv[midpoint]['end']))

			if start <= variant.POS <= end:
				copy_number = int(float(cnv[midpoint]['total_cn']))
				found = True
			else:
				if variant.POS < start:
					last = midpoint-1
				else:
					first = midpoint+1
		return copy_number

def sort_by_list(to_sort, sort_by, reverse=False):
	X = to_sort
	Y = sort_by
	return([x for (y,x) in sorted(zip(Y,X), key=lambda pair: pair[0], reverse=reverse)])

def filter_vcf(vcf_parser, variants):
	variants_filtered = []
	total_counts = []
	for record in variants:
		ref_count, total_count = vcf_parser._calc_read_counts(record)
		if (total_count < min_readdepth):
			continue

		if record.CHROM == "X":
			continue

		variants_filtered.append(record)
		total_counts.append(total_count)

	if (half_with_highest_depth):
		variants_sorted = sort_by_list(variants_filtered, total_counts, reverse=True)
		variants_filtered = variants_sorted[:len(variants_sorted)/2]

	return variants_filtered


def get_correct_vaf(cnv_regions, vcf_file, purity, ouput_file):
	formatter = CnvFormatter(None, None, None, None)
	output = []

	vcf_parser = PcawgConsensusParser(vcf_file)
        variants = vcf_parser._parse_vcf(vcf_file)

	# print("Variants before filtering by depth (" + str(min_readdepth)+ "): " + str(len(variants)))
	variants = filter_vcf(vcf_parser, variants)

	# print("Variants after filtering by depth (" + str(min_readdepth)+ "): " + str(len(variants)))
	if (len(variants) > cutoff):
		variants = random.sample(variants, cutoff)
	# print("Variants after cutoff (" + str(cutoff)+ "): " + str(len(variants)))

	vafs = []
	vafs_new = []
	for record in variants:
		ref_count, total_count = vcf_parser._calc_read_counts(record)
		alt_count = total_count - ref_count

		vaf = vcf_parser._get_vaf(record)
		if vaf is None:
			vaf = round(alt_count / float(total_count),4)
		
		vafs.append(vaf)
		if (sample_vaf_from_posterior):
			vaf_new =  round(beta.rvs(alt_count+1, ref_count+1),4)

			vafs_new.append(vaf_new)
			vaf = vaf_new

		if cnv_regions is not None:
			copy_number = formatter._find_copy_number(record, cnv_regions[record.CHROM])
		else:
			copy_number = 2

		if purity is None:
			purity = 1

		corrected_vaf = (2 + purity * (copy_number-2)) * vaf

		output.append([str(record.CHROM) + "_" + str(record.POS) + "=" + str(corrected_vaf)])

	with open(ouput_file, "w") as ouput_file:
		for record in output:
			ouput_file.write("\t".join([str(x) for x in record]) + "\n")
	
def read_purity(purity_file):
	purities = {}
	with open(purity_file) as f:
		header = f.readline().split()
		for l in f:
			tokens = l.split()
			tumor_id = tokens[header.index('samplename')] 
			purity = float(tokens[header.index('purity')])
			purities[tumor_id] = purity	
	return (purities)

def main():
	parser = argparse.ArgumentParser(description='add copy number to vcf files')
	parser.add_argument('--cnv', dest='cnv', 
	 					help='Path to CNV file')
	parser.add_argument('--vcf', dest='vcf',
	 					help='Path to variants (vcf) file')
	parser.add_argument('--purity', dest='purity_file',
                                                help='Path to purity file')
	parser.add_argument('--output', dest='output',
                                                help='Path to output file')	
	args = parser.parse_args()

	cnv = args.cnv
	vcf = args.vcf
	purity_file = args.purity_file
	output = args.output

	if output is None:
		print("Please provide output file using --output option")
		exit()

	if cnv is not None:
		cnv_parser = CnvParser(cnv)
		cnv_regions = cnv_parser.parse()
	else:
		cnv_regions = None

	start = max(0, (vcf.rfind("/")+1))
	end = vcf[start:].find(".") + start

	tumor_id = vcf[start : end]

	tumor_purity = None
	if purity_file is not None:
		purity = read_purity(purity_file)

		if tumor_id not in purity.keys():
			print("Tumor name not found in purity list")
			exit()
		tumor_purity = purity[tumor_id]

	get_correct_vaf(cnv_regions, vcf, tumor_purity, output)

if __name__ == "__main__":
	main()
