#Python scripts for the entropy demo

def map_ranks():
  ranks_file = '/Users/ilyanovikov/Documents/rna_et_project/rnasep_example/data/rnasep_sample_ranks.csv'
  infile = open(ranks_file, "r")
  for line in infile:
    residue, raw_score, cov = line.strip().split(',')
    cmd.alter("resi %s and chain_B" % residue, "b=%s" % cov)
  infile.close()
  cmd.spectrum("b", "rainbow_rev", "chain_B", 0,1)
cmd.extend("map_ranks", map_ranks)


def make_pretty():
	cmd.create('chain_B', 'chain B and 3q1q')
	cmd.create('chain_A', 'chain A and 3q1q')
	cmd.create('chain_C', 'chain C and 3q1q')

	cmd.show('surface', 'chain_B')
	cmd.show('cartoon', 'chain_A')
	cmd.show('cartoon', 'chain_C')

	cmd.color('white', 'chain_B')
	cmd.color('yellow', 'chain_A')
	cmd.color('cobalt', 'chain_C')

	cmd.hide('everything', '3q1q')

cmd.extend("make_pretty", make_pretty)

def mark_known_sites():
  fs_file = '/Users/ilyanovikov/Documents/rna_et_project/rnasep_example/input_data/3Q1Q_RF00010_functional_site.txt'
  infile = open(fs_file, "r")
  for line in infile:
    nt, residue = line.strip().split(' ')
    print(residue)
    #cmd.alter("resi %s and chain B" % residue, "b=%s" % rank)
  infile.close()
cmd.extend("mark_known_sites", mark_known_sites)

#spectrum b, rainbow_rev, chain B
#select low_ranked, b > 350