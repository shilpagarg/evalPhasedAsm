#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import vcf
from argparse import ArgumentParser
from collections import defaultdict

def comparison(ph, th, th_offset):

	correct = 0 
	unphasable = 0 
	ambiguous = 0
	switches = 0

	phasings = []  
	ho = 0 
	he = 0 
	gaperrors = 0 
	breaks = 0 
	
	oldx = ('-','-')
	th_pos = th_offset
	ph0_pos = 0
	ph1_pos = 0
	while True:
		if (ph[0][ph0_pos] == '|') and (ph[1][ph1_pos] == '|'):
			breaks += 1
		fb = False
		if ph[0][ph0_pos] == '|':
			fb = True
			ph0_pos += 1
		if ph[1][ph1_pos] == '|':
			fb = True
			ph1_pos += 1
		if fb:
			oldx = ('-','-')
			continue
	
		if (oldx[0] == '-') or (oldx[1] == '-') or (ph[0][ph0_pos] == '-') or (ph[1][ph1_pos] == '-'):
			unphasable += 1
			phasings.append('U')
			if (ph[0][ph0_pos] == '-') or (ph[1][ph1_pos] == '-'): 
				gaperrors += 1
				oldx = ('-','-')
			else:
				oldx = (ph[0][ph0_pos], th[0][th_pos])

		elif (oldx[0] == 'X') or (oldx[1] == 'X') or (ph[0][ph0_pos] == 'X') or (ph[1][ph1_pos] == 'X'):
			ambiguous += 1
			phasings.append('A')
			if (ph[0][ph0_pos] == 'X') or (ph[1][ph1_pos] == 'X'): 
				oldx = ('X','X')
			else:
				oldx = (ph[0][ph0_pos], th[0][th_pos])

		elif (ph[0][ph0_pos] == ph[1][ph1_pos]) and (th[0][th_pos] != th[1][th_pos]): 
			phasings.append('O')
			ho += 1

		elif (ph[0][ph0_pos] != ph[1][ph1_pos]) and (th[0][th_pos] == th[1][th_pos]): 
			phasings.append('E')
			he += 1

		elif ((ph[0][ph0_pos] == oldx[0] and th[0][th_pos] != oldx[1]) or 
		      (ph[0][ph0_pos] != oldx[0] and th[0][th_pos] == oldx[1])):
			switches += 1
			phasings.append('S')

			oldx = (ph[0][ph0_pos], th[0][th_pos])
		else:

			oldx = (ph[0][ph0_pos], th[0][th_pos])
			correct += 1
			phasings.append('C')

		th_pos += 1
		ph0_pos += 1
		ph1_pos += 1

		if ph0_pos == len(ph[0]) : break
		if ph1_pos == len(ph[1]) : break
		if th_pos == len(th[0]) : break

	return switches, ho, he, gaperrors, breaks, unphasable, ambiguous, correct, phasings

def evalphase(phasings):

	flips = 0
	switches = 0
	deserts = 0
	ambiguous = 0
	homo = 0
	hetero = 0
	correct = 0
	oldswitches = 0

	consecutive = 0
	unphasable = True
	for i, phasing in enumerate(phasings):
		if phasing == 'O' or phasing == 'E':
			if phasing == 'O':
				homo += 1
			if phasing == 'E':
				hetero += 1
			continue
		if phasing == 'A':
			ambiguous += 1
			continue
		if phasing == 'C':
			if unphasable:
				flips += (consecutive+1) // 2
				oldswitches += consecutive
			elif not unphasable:
				switches += consecutive % 2
				flips += consecutive // 2
				oldswitches += consecutive
			correct += 1
			consecutive = 0
			unphasable = False

		elif phasing == 'U':
			flips += (consecutive+1) // 2
			oldswitches += consecutive
			deserts += 1
			consecutive = 0
			unphasable = True
		
		elif phasing == 'S':
			consecutive += 1

			
	return flips, switches, homo, hetero, deserts, ambiguous, correct, oldswitches

def read_truth_vcf(filename):

	vr = vcf.Reader(filename)
	samples = vr.samples
	assert len(samples) == 1, 'Expected exactly one sample in truth VCF'
	result = {}
	for record in vr:
		if not record.is_snp:
			continue
		if len(record.ALT) != 1:
			continue
		call = record.samples[0]
		if not call.is_het:
			continue
		genotype = call.data.GT
		assert genotype in ['0|1', '1|0'], 'Expecting all SNPs to be phased (in "|" notation) in truth file'
		result[(record.CHROM, record.POS)] = genotype
	return result


def read_phased_vcf(filename, true_snps = None):

	vr = vcf.Reader(filename)
	samples = vr.samples
	assert len(samples) == 1, 'Expected exactly one sample in phased VCF'
	blocks = defaultdict(list)
	phased = 0
	unphased = 0
	tp_snps = set()
	FP = 0
	FN = 0
	for record in vr:
		if not record.is_snp:
			continue
		if len(record.ALT) != 1:
			continue
		call = record.samples[0]
		if not call.is_het:
			continue
		if true_snps != None:
			if (record.CHROM, record.POS) not in true_snps:
				FP += 1
				continue
		tp_snps.add( (record.CHROM, record.POS) )
		if hasattr(call.data,'PS'):
			phased += 1
			blocks[(record.CHROM, call.data.PS)].append((record.POS,call.data.GT))
		else:
			unphased += 1
	if true_snps != None:
		FN = len(true_snps) - len(tp_snps)
	return phased, unphased, FP, FN, blocks

def main():
	parser = ArgumentParser(prog='evalPhase', description=__doc__)
	parser.add_argument('-v', dest='verbose', action='store_true', default=False,
		help='Be (very) verbose and output statistics on every single phased block.')
	parser.add_argument('truthvcf', metavar='truthvcf',
		help='True phased vcf')
	parser.add_argument('predictedvcf', metavar='predictedvcf',
		help='Predicted phased vcf')

	args = parser.parse_args()

	truth = read_truth_vcf(open(args.truthvcf))
	
	print('Read {} heterozygous SNPs from true phasing'.format(len(truth)))
	
	phased, unphased, FP, FN, phased_blocks = read_phased_vcf(open(args.predictedvcf), truth)
	
	print('Found {} false positive (FP) heterozygous SNPs'.format(FP))
	print('There are {} false negative (FN) heterozygous SNPs'.format(FN))
	print('Retained {} heterozygous SNPs from predicted phasing, out of which:'.format(phased+unphased))
	
	print('  phased SNPs: {}'.format(phased))
	print('  unphased SNPs: {}'.format(unphased))
	print('  number of blocks: {}'.format(len(phased_blocks)))
	print('    --> unphasable due to being first in a block: {}'.format(len(phased_blocks)))
	
	block_list = list(phased_blocks.keys())
	block_list.sort()
	correct_total = 0
	switches_total = 0
	flips_total = 0
	ho_total = 0
	he_total = 0
	blocks_larger_1 = 0
	for (chromosome,block_id) in block_list:
		block = phased_blocks[(chromosome,block_id)]
		if len(block) >= 2:
			blocks_larger_1 += 1
		th = ['','']
		ph = ['','']
		for pos, pred_gt in block:
			ph[0] += pred_gt[0]
			ph[1] += pred_gt[2]
			assert (chromosome, pos) in truth, 'Positions in truth VCF and predictions VCF do not match'
			true_gt = truth[(chromosome, pos)]
			th[0] += true_gt[0]
			th[1] += true_gt[2]
		switches0, ho0, he0, gaperrors0, breaks0, unphasable0, ambiguous0, correct0, phasings = comparison(ph, th, 0)
		flips, switches1, ho1, he1, unphasable1, ambiguous1, correct1, oldswitches = evalphase(phasings)
		assert unphasable1 == 1
		assert ambiguous1 == 0
		if args.verbose:
			print('-'*100)
			print('Block ID:', block_id)
			print('True block:     ', th[0])
			print('Predicted block:', ph[0])
			print("Phasings:", ''.join(phasings))
			print("Correctly phased: ", correct1)
			print("Unphasable (due to gaps or breaks): ", unphasable1)
			print("Switch errors: ", switches1)
			print("Flip errors: ", flips)

		correct_total += correct1
		switches_total += switches1
		flips_total += flips
		ho_total += ho1
		he_total += he1
	print('='*100)
	print('Evaluation of blocks:')
	print('Number of blocks with at least 2 SNPs:', blocks_larger_1)
	print("Correctly phased: ", correct_total)
	print("Switch errors: ", switches_total)
	print("Flip errors: ", flips_total)
	assert ho_total == he_total == 0


if __name__ == '__main__':
	sys.exit(main())
