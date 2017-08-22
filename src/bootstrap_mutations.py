# AUTHOR: Yulia Rubanova

#!/usr/bin/python

from __future__ import print_function
import sys, random, os
import numpy as np

if __name__ == "__main__":
	if len(sys.argv) < 3:
		print("Please provide more parameters")
		sys.exit()

	file = sys.argv[1]
	n_bootstrap = sys.argv[2]
	if n_bootstrap.isdigit():
		n_bootstrap = int(n_bootstrap)
	
	assert(n_bootstrap > 0)

	tmp_file=str(int(random.uniform(0,10000))) + "1.txt"	
	os.system("wc -l " + file + " > " + tmp_file )
	with open(tmp_file, 'r') as f:
		n_mutations = int(f.readline().split()[0])

	os.system("rm " + tmp_file )

	indices_to_bootstrap = np.random.choice(xrange(n_mutations), n_bootstrap)

	mut_list = []

	sum_phis = 0
	with open(file, 'r') as mut_types:
		for i, line in enumerate(mut_types):
			chr, pos, phi, ref, alt, context = line.split()
			mut_list.append(line.strip())

	for i in indices_to_bootstrap:
		print(mut_list[i])
		





