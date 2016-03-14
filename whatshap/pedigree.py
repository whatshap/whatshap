"""
Pedigree-related functions
"""
from collections import namedtuple, Counter


class ParseError(Exception):
	pass


Individual = namedtuple('Individual', ['id', 'mother_id', 'father_id'])


class PedReader:
	"""
	A parser for PED/FAM files as used by PLINK and other tools.

	According to <http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml>:
	The PED file is a white-space (space or tab) delimited file: the first six
	columns are mandatory:
	* Family ID
	* Individual ID
	* Paternal ID
	* Maternal ID
	* Sex (1=male; 2=female; other=unknown)
	* Phenotype

	All fields except the individual, maternal and paternal ID are ignored by
	this class. The entire file is read upon construction.
	"""
	def __init__(self, file):
		if isinstance(file, str):
			with open(file) as f:
				self.individuals = self._parse(f)
		else:
			self.individuals = self._parse(file)

	@staticmethod
	def _parse_record(line):
		"""
		Parse a single non-comment line of a PED or FAM file.
		"""
		fields = line.split()
		if len(fields) < 6:
			raise ParseError("Less than six fields found in PED/FAM file")
		individual_id, paternal_id, maternal_id = fields[1:4]
		if paternal_id == '0':
			paternal_id = None
		if maternal_id == '0':
			maternal_id = None
		return Individual(id=individual_id, mother_id=maternal_id, father_id=paternal_id)

	def _parse(self, file):
		individuals = []
		for line in file:
			if line.startswith('#') or line == '\n':
				continue
			individuals.append(self._parse_record(line))
		self._sanity_check(individuals)
		return individuals

	@staticmethod
	def _sanity_check(individuals):
		"""
		Ensure that each individual occurs only once in the file.
		"""
		individual_ids = [ individual.id for individual in individuals ]
		if not individual_ids:
			return
		id, count = Counter(individual_ids).most_common()[0]
		if count > 1:
			raise ParseError('Individual {!r} occurs more than once in PED file'.format(id))

	def __iter__(self):
		return iter(self.individuals)
