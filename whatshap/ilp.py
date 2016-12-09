"""
Given a read set, solve MEC by integer linear programming (ILP)
"""

__author__ = "Gunnar Klau, Tobias Marschall"


from .core import ReadSet
from gurobipy import Model, quicksum, GRB, LinExpr
from collections import defaultdict

def phase_by_flipping_ilp(read_set, all_het):
	"""
	document later
	"""
	model = Model('mec')
	# create one variable per read
	x = [ model.addVar(vtype=GRB.BINARY, name='x_'+read.name) for read in read_set ]
	model.update()
	# mapping positions to a list [(read_index, allele), ...]
	columns = defaultdict(list)
	for read_index, read in enumerate(read_set):
		for variant in read:
			columns[variant.position].append((read_index,variant.allele,variant.quality))
	positions = sorted(columns.keys())
	h0 = { p:model.addVar(vtype=GRB.BINARY, name='h0_{}'.format(p)) for p in positions }
	h1 = { p:model.addVar(vtype=GRB.BINARY, name='h1_{}'.format(p)) for p in positions }
	# create flipping variables
	f = {}
	for read_index, read in enumerate(read_set):
		f[read_index] = {}
		for variant in read:
			f[read_index][variant.position] = model.addVar(vtype=GRB.BINARY, obj=variant.quality, name='f_{}_{}'.format(read_index, variant.position))
	model.update()
	for read_index, read in enumerate(read_set):
		for variant in read:
			if variant.allele == 1:
				model.addConstr(f[read_index][variant.position] >= 1 - h0[variant.position] - x[read_index])
				model.addConstr(f[read_index][variant.position] >= 1 - h1[variant.position] - (1 - x[read_index]))
			else:
				model.addConstr(f[read_index][variant.position] >= h0[variant.position] - x[read_index])
				model.addConstr(f[read_index][variant.position] >= h1[variant.position] - (1 - x[read_index]))
	if all_het:
		for p in positions: model.addConstr(h0[p] + h1[p] == 1)
	model.update()
	model.optimize()
	print('model.write')
	model.write('flipping.lp')
	cost = model.ObjVal
	partition = None # [ x[i].X for i,read in enumerate(read_set) ]
	superreads = None # ReadSet()
	return superreads, cost, partition

def phase_by_chen_et_al_ilp(read_set, all_het):
	"""
	The ILP model introduced in [1] and used in [2]

	[1] Chen, Zhi-Zhong, Fei Deng, and Lusheng Wang. 2013. “Exact Algorithms for Haplotype Assembly From Whole-Genome Sequence Data..” Bioinformatics (Oxford, England) 29 (16): 1938–45. doi:10.1093/bioinformatics/btt349.
	[2] Chen, Zhi-Zhong, Fei Deng, Chao Shen, Yiji Wang, and Lusheng Wang. 2016. “Better ILP-Based Approaches to Haplotype Assembly..” Journal of Computational Biology 23 (7): 537–52. doi:10.1089/cmb.2015.0035.
	"""
	assert all_het, 'Model is only for the all-het case'
	model = Model('mec')
	# create haplotype variables
	y = [ model.addVar(vtype=GRB.BINARY, name='y_'+read.name) for read in read_set ]
	model.update()
	columns = defaultdict(list)
	for read_index, read in enumerate(read_set):
		for variant in read:
			columns[variant.position].append((read_index,variant.allele,variant.quality))
	positions = sorted(columns.keys())
	x = { p:model.addVar(vtype=GRB.BINARY, name='x_{}'.format(p)) for p in positions }
	# interpretation: x_j = 1 iff h_j = 1 amd \bar{h}_j = 0
	# create product variables
	z = {}
	for read_index, read in enumerate(read_set):
		z[read_index] = {}
		for variant in read:
			z[read_index][variant.position] = model.addVar(vtype=GRB.BINARY, name='z_{}_{}'.format(read_index, variant.position))
	model.update()
	for read_index, read in enumerate(read_set):
		for variant in read:
			model.addConstr(z[read_index][variant.position] <= x[variant.position])
			model.addConstr(z[read_index][variant.position] <= y[read_index])
			model.addConstr(z[read_index][variant.position] >= x[variant.position] + y[read_index] - 1)
	model.update()
	obj = LinExpr(0.0)
	for read_index, read in enumerate(read_set):
		for variant in read:
			if variant.allele == 0:
				obj.addTerms(2 * variant.quality, z[read_index][variant.position])
				obj.addTerms(-variant.quality, x[variant.position])
				obj.addTerms(-variant.quality, y[read_index])
				obj.addConstant(variant.quality)
			else:
				obj.addTerms(-2 * variant.quality, z[read_index][variant.position])
				obj.addTerms(variant.quality, x[variant.position])
				obj.addTerms(variant.quality, y[read_index])
	model.setObjective(obj, GRB.MINIMIZE)
	model.update()
	model.optimize()
	print('model.write')
	model.write('chen_et_al.lp')
	cost = model.ObjVal
	partition = None # [ x[i].X for i,read in enumerate(read_set) ]
	superreads = None # ReadSet()
	return superreads, cost, partition



def phase_by_bigM_ilp(read_set, all_het):
	"""
	document later
	"""
	assert all_het, 'Not yet implemented'
	model = Model('mec')
	# create one variable per read
	x = [ model.addVar(vtype=GRB.BINARY, name='x_'+read.name) for read in read_set ]
	model.update()
	# mapping positions to a list [(read_index, allele), ...]
	columns = defaultdict(list)
	for read_index, read in enumerate(read_set):
		for variant in read:
			columns[variant.position].append((read_index,variant.allele,variant.quality))
	positions = sorted(columns.keys())
	#for c in columns: print(c, columns[c])
	# Create variables giving the cost of flipping all bits in a given columns
	# to equal the haplotype index (or not).
	# cost of being equal
	e = {p:model.addVar(name='e_{}'.format(p)) for p in positions }
	# cost of being different
	d = {p:model.addVar(name='d_{}'.format(p)) for p in positions }
	# variable giving the minimum in each column
	m = {p:model.addVar(name='m_{}'.format(p)) for p in positions }
	# variable telling whether e or d is the minimum
	y_e = {p:model.addVar(vtype=GRB.BINARY, name='y_e_{}'.format(p)) for p in positions }
	y_d = {p:model.addVar(vtype=GRB.BINARY, name='y_d_{}'.format(p)) for p in positions }
	#add constraints
	for p in positions:
		model.addConstr(e[p] == quicksum(weight*x[read_index] for read_index, allele, weight in columns[p] if allele == 0)
			+ quicksum(weight*(1-x[read_index]) for read_index, allele, weight in columns[p] if allele == 1),
			name='C_e_{}'.format(p))
		model.addConstr(d[p] == quicksum(weight*x[read_index] for read_index, allele, weight in columns[p] if allele == 1)
			+ quicksum(weight*(1-x[read_index]) for read_index, allele, weight in columns[p] if allele == 0),
			name='C_d_{}'.format(p))
		M = sum(weight for _,_,weight in columns[p])
		model.addConstr(m[p] >= e[p] - y_e[p]*M)
		model.addConstr(m[p] >= d[p] - y_d[p]*M)
		model.addConstr(y_e[p] + y_d[p] <= 1)
	model.setObjective(quicksum(m[p] for p in positions), GRB.MINIMIZE)
	model.update()
	model.optimize()
	print('model.write')
	model.write('test.lp')
	cost = model.ObjVal
	partition = None # [ x[i].X for i,read in enumerate(read_set) ]
	superreads = None # ReadSet()
	return superreads, cost, partition
