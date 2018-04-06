from whatshap.core import Pedigree, PhredGenotypeLikelihoods, NumericSampleIds

def test_pedigree_no_gls():
	ped = Pedigree(NumericSampleIds())
	genotypes1 = [ 0, 1, 0, 2 ]
	genotypes5 = [ 1, 2, 2, 0 ]
	ped.add_individual('sample1', genotypes1)
	assert len(ped) == 1
	assert ped.variant_count == 4, str(ped.variant_count)
	ped.add_individual('sample5', genotypes5)
	assert len(ped) == 2
	assert ped.variant_count == 4, str(ped.variant_count)
	for i in range(ped.variant_count):
		assert ped.genotype('sample1', i) == genotypes1[i]
		assert ped.genotype_likelihoods('sample1', i) is None
		assert ped.genotype('sample5', i) == genotypes5[i]
		assert ped.genotype_likelihoods('sample5', i) is None


def test_pedigree_with_gls():
	ped = Pedigree(NumericSampleIds())
	genotypes1 = [ 0, 1, 0, 2 ]
	gls1 = [
		PhredGenotypeLikelihoods(0, 1, 2),
		PhredGenotypeLikelihoods(215, 81, 147),
		PhredGenotypeLikelihoods(199, 49, 253),
		PhredGenotypeLikelihoods(167, 200, 163),
	]
	genotypes5 = [ 1, 2, 2, 0 ]
	gls5 = [
		PhredGenotypeLikelihoods(184, 71, 233),
		PhredGenotypeLikelihoods(65, 32, 87),
		PhredGenotypeLikelihoods(28, 215, 131),
		PhredGenotypeLikelihoods(98, 250, 137),
	]
	ped.add_individual('sample1', genotypes1, gls1)
	assert len(ped) == 1
	assert ped.variant_count == 4, str(ped.variant_count)
	ped.add_individual('sample5', genotypes5, gls5)
	assert len(ped) == 2
	assert ped.variant_count == 4, str(ped.variant_count)
	for i in range(ped.variant_count):
		assert ped.genotype('sample1', i) == genotypes1[i]
		assert list(ped.genotype_likelihoods('sample1', i)) == list(gls1[i])
		assert ped.genotype('sample5', i) == genotypes5[i]
		assert list(ped.genotype_likelihoods('sample5', i))  == list(gls5[i])
