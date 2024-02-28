import sys
import os
import time
import shutil
import glob
import subprocess
import math
import copy
import pydot
import networkx
import random
import pylab
import datetime
import argparse
import random
import math
import types
from Bio import Phylo
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from io import StringIO

megacc_app = "megacc"
mp_tree_infer_mao = "scripts/pathfinder/infer_NJ_amino_acid.mao"
ancestral_seqs_mao = "ancestral_seqs_ML_protein.mao"
outgroup_file = "outgroup.txt"

print_megacc_cmd = False

parser = argparse.ArgumentParser(description="PathFinder tumor migration path solver.")
parser.add_argument("aln", help="Clone sequence alignment file.", type=str)
parser.add_argument("-o", "--output", help="Output directory to put results in.", type=str, default=".")

args = parser.parse_args()

if args.output != "." and not os.path.exists(args.output):
	os.mkdir(args.output)

def parse_input_aln(aln_file_in, aln_file_out):
	seqs = {}
	with open(aln_file_out, 'w') as out_file:
		with open(aln_file_in, 'r') as file:
			seqname = ""
			out_file.write("#MEGA\n")
			out_file.write("!Title SNVs;\n")
			out_file.write("!Format datatype=Protein;\n")
			for line in file:
				line = line.strip()
				if len(line) == 0: continue
				if line == "#MEGA" or line[0] == "!":
					continue
				if line[0] == "#" or line[0] == ">":
					seqname = line[1:]
				else:
					seqs[seqname] = line
		if 'Normal' not in seqs.keys():
			seqs['Normal'] = 'A' * len(next(iter(seqs.values())))
		seq_len = len(seqs[next(iter(seqs.keys()))])
		target_seq_len = 3
		repeat_count = int(math.ceil(float(target_seq_len) / float(seq_len)))
		primary = ""
		for i in range(seq_len):
			normal_allele = seqs['Normal'][i]
			mut_alleles = list(set([seqs[seqid][i] for seqid in seqs.keys() if seqid != 'Normal']))
			if len(mut_alleles) == 1 and mut_alleles[0] != normal_allele:
				allele = mut_alleles[0]
			else:
				allele = normal_allele
			primary = primary + allele
		for key in seqs.keys(): seqs[key] = seqs[key] * repeat_count
		primary = primary * repeat_count
		for key in seqs.keys():
			out_file.write("#" + key + '\n')
			out_file.write(seqs[key] + '\n')
	return seqs

def read_paths(filename):
	paths = []
	with open(filename, 'r') as file:
		for line in file:
			data = line.strip().split('->')
			paths.append((data[0], data[1]))
	return paths

def infer_mp_tree(mega_aln_filename):
	base_filename = os.path.splitext(os.path.basename(mega_aln_filename))[0]
	tree_filename = os.path.join(scratch_dir, base_filename + ".nwk")
	megacc_cmd = "{} -a {} -d {} -o {}".format(megacc_app, mp_tree_infer_mao, mega_aln_filename, tree_filename)
	if print_megacc_cmd: print(megacc_cmd)
	FNULL = open(os.devnull, 'w')
	return_code = subprocess.call(megacc_cmd, stdout=FNULL, stderr=subprocess.STDOUT, shell=True)
	if return_code != 0:
		raise ValueError('MEGACC returned error code', return_code)
	trees = Phylo.parse(tree_filename, 'newick')
	tree = Phylo.BaseTree.Tree.from_clade(trees.__next__().clade)
	tree.root_with_outgroup({'name': 'Normal'})
	anc_id = 0
	clade_count = 0
	for clade in tree.find_clades():
		clade_count += 1
		if str(clade.name) == "None":
			clade.name = 'anc_node_' + str(anc_id)
			anc_id += 1
	anc_states_files = glob.glob(base_filename + '_ancestral_states_*.txt')
	for file in anc_states_files:
		try:
			os.remove(file)
		except:
			print("Error while deleting ancestral states file : ", file)
	changes_list_files = glob.glob(base_filename + '_changes_list_*.txt')
	for file in changes_list_files:
		try:
			os.remove(file)
		except:
			print("Error while deleting ancestral states file : ", file)
	for clade in tree.find_clades():
		if clade.branch_length is None:
			clade.branch_length = 0.0
	return tree, clade_count

def make_scratch_dir(basename):
	idx = random.randint(1, 1000)
	while os.path.exists("{}_{}".format(basename, idx)):
		idx = random.randint(1, 1000)
	os.mkdir("{}_{}".format(basename, idx))
	return "{}_{}".format(basename, idx)

scratch_dir = make_scratch_dir(os.path.join(args.output, "scratch"))

aln_file_out = os.path.join(scratch_dir, os.path.splitext(os.path.basename(args.aln))[0] + ".meg")
mut_seqs = parse_input_aln(args.aln, aln_file_out)

mega_aln_filename = aln_file_out

initial_tree, node_count = infer_mp_tree(mega_aln_filename)


