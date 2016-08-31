"""
Given a read set, solve MEC by integer linear programming (ILP)
"""

from .core import ReadSet
from gurobipy import Model, quicksum, GRB
from collections import defaultdict

def phase_by_ilp(read_set, all_het):
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
	model.update()
	# add constraints
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


#m.model_sense = GRB.MAXIMIZE
#m.addVar(vtype=GRB.BINARY, obj=1, name='x')


