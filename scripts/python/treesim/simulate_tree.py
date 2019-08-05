import numpy as np
import newick

def translate_seq(seq, states):
	return ''.join([states[i] for i in seq])

def sample_from(freq):
	cum_freq = np.cumsum(freq)
	r = np.random.rand()
	for i in range(len(cum_freq)):
		if r < cum_freq[i]:
			return i 

def generate_seq(l, freq):
	seq = np.zeros(l, dtype=int)
	for i in range(l):
		seq[i] = sample_from(freq)
	return seq

def mutate_seq_expm(ancestor, subs_model, t):
	"""Mutates sequence using matrix exponentiation"""
	p = subs_model.get_p(t)
	child = np.copy(ancestor)
	for i in range(len(child)):
		state = child[i]
		prob = p[state, :]
		child[i] = sample_from(prob)
	return child

def normalize_mutation(q, state):
	prob = np.copy(q[state, :].A1)
	prob[state] = 0
	scalor = 1 / np.sum(prob)
	return prob * scalor

def mutate_seq_pois(ancestor, subs_model, t):
	"""Mutates sequence using draws from a Poisson distribution"""
	q = subs_model.get_q()
	child = np.copy(ancestor)
	total_mutations = 0
	for i in range(len(child)):
		state = child[i]
		rate = -q[state, state] * t
		mutations = np.random.poisson(rate)
		total_mutations = total_mutations + mutations
		for j in range(mutations):
			prob = normalize_mutation(q, state)
			child[i] = sample_from(prob)
			state = child[i]
	return (child, total_mutations)

def simulate_tree(newick_tree, subs_model, l, mode='expm'):
	"""Simulates sequences on a fixed tree for a substitution model.
	Returns a dictionary of simulated sequences (key-value pairs: node name, sequence).

	Parameters:
	newick_tree -- the tree in newick representation
	subs_model  -- the substitution model
	l           -- the sequence length
	mode        -- the method used to mutate sequences, options are:
	               'expm': mutation by matrix exponentiation
	               'pois': mutation by poisson process
	"""
	seq_map = {}
	tree = newick.loads(newick_tree) 
	root = tree[0]
	freq = subs_model.get_pi()
	seq_map[root.name] = generate_seq(l, freq) # generate ancestor
	nodes = root.descendants.copy()
	while len(nodes) > 0:
		child = nodes.pop()
		t = child.length:
		ancestor = child.ancestor
		ancestor_seq = seq_map[ancestor.name]
		if mode == 'expm':
			child_seq = mutate_seq_expm(ancestor_seq, subs_model, t)
		elif mode == 'pois':
			child_seq, _ = mutate_seq_pois(ancestor_seq, subs_model, t)
		seq_map[child.name] = child_seq
		nodes.extend(child.descendants.copy())
	states = subs_model.get_states()
	for i in seq_map: # translate sequence to states
		seq_map[i] = translate_seq(seq_map[i], states)
	return seq_map
