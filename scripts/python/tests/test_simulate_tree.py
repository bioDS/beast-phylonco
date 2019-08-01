import context
from treesim import simulate_tree as tree
from treesim import substitution_model as subs_model
from matplotlib import pyplot as plt
import math
import numpy as np
from scipy.stats import gaussian_kde
from scipy.stats import poisson

def count_diff(a, b):
	min_len = min(len(a), len(b))
	max_len = max(len(a), len(b))
	diff = sum([a[i] != b[i] for i in range(min_len)])
	return diff + max_len - min_len

def plot_kde(data, plot, label='', color='#AAAAFF'):
	kde = gaussian_kde(data)
	f = kde.covariance_factor()
	bw = f * sample.std()
	grid = np.linspace(min(data), max(data), 200)
	plot.fill(grid, kde.evaluate(grid), fc=color, alpha=0.5, label=label)

def example_simulate_tree():
	"""Example usage of tree simulator"""
	print("running example_simulate_tree()")
	lambda_d = 3
	lambda_l = 2
	l = 100
	seed = 777
	newick_tree = '((A:0.5, B:0.5)D:0.5, C:1)E;'
	np.random.seed(seed)		
	model = subs_model.SiFit3(lambda_d, lambda_l)
	seq_map = tree.simulate_tree(newick_tree, model, l)	
	for x in seq_map:
		print("%s: %s\n" % (x, seq_map[x])) 

def test_simulate_tree_1(t=1.0, l=200, N=1000, seed=777):
	"""Test case for tree simulator"""
	print("running test_simulate_tree()")
	lambda_d = 3
	lambda_l = 2	
	np.random.seed(seed)
	model = subs_model.SiFit3(lambda_d, lambda_l)
	freq = model.get_pi()
	ancestor = tree.generate_seq(l, freq)
	diff_child_expm = np.zeros(N)
	diff_child_pois = np.zeros(N)
	mutations = np.zeros(N)
	for i in range(N):
		child_expm = tree.mutate_seq_expm(ancestor, model, t)
		child_pois, mutation_pois = tree.mutate_seq_pois(ancestor, model, t)
		diff_child_expm[i] = count_diff(ancestor, child_expm)
		diff_child_pois[i] = count_diff(ancestor, child_pois)
		mutations[i] = mutation_pois
	
	# expected number of mutations for poisson process
	expected_mean = l * t 
	expected_var = expected_mean
	observed_mean = mutations.mean()
	observed_var = math.pow(mutations.std(), 2)
	print("expected mean of mutations: %f" % expected_mean)
	print("observed mean of mutations: %f" % observed_mean)
	print("expected variance of mutations: %f" % expected_var)
	print("observed variance of mutations: %f" % observed_var)
	
	# plot results
	plt.clf()
	plot_kde(mutations, plt, 'mutate poisson')
	plt.xlabel("Number of mutations (t = 1.0, l = 200)")
	plt.ylabel("Normalized density", rotation=90)
	plt.title("Number of mutations using poisson process (n = %d)" % N)
	x = np.arange(min(mutations), max(mutations))
	plt.plot(x, poisson.pmf(x, expected_mean), c='green', alpha=0.5)
	plt.show()
	# plt.savefig("Fig1_t%f_n%d_%d.pdf" % (t, N, seed))
	
	plt.clf()
	plot_kde(diff_child_expm, plt, 'mutate expm', 'orange')
	plot_kde(diff_child_pois, plt, 'mutate poisson')
	plt.xlabel("Number of different sites (t = 1.0, l = 200)")
	plt.ylabel("Normalized density", rotation=90)
	plt.legend()
	plt.title("Comparison of two mutation implementations (n = %d)" % N)
	plt.show()
	# plt.savefig("Fig2_t%f_n%d_%d.pdf" % (t, N, seed))

def test_simulate_tree_2():
	import numpy as np
	mu_d = -1
	mu_l = 0.5
	l = 200
	seed = 777
	newick_tree = '((A:0.05, B:0.05)D:0.05, C:0.1)E;'
	np.random.seed(seed)		
	model = subs_model.SiFit3(lambda_d, lambda_l)
	file = open("python_lognorm.csv", 'w')
	file.write("tree,lambdaD,lambdaL,node,sequence\n")
	for i in range(N):
		lambda_d = np.random.lognormal(mu_d)
		lambda_l = np.random.lognormal(mu_d)
		seq_map = tree.simulate_tree(newick_tree, model, l)	
		for x in seq_map:
			file.write("%d,%f,%f,%s,%s\n" % (i, lambda_d, lambda_l, x, seq_map[x]))
	file.close()

if __name__ == '__main__':
	example_simulate_tree()	
	test_simulate_tree_1()