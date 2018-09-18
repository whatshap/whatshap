from whatshap.core import Pedigree, GenotypeLikelihoods, NumericSampleIds, Genotype

def test_pedigree_no_gls():
	ped = Pedigree(NumericSampleIds(), 2)
	genotypes1 = [Genotype([0,0]), Genotype([0,1]), Genotype([0,0]), Genotype([1,1])]
	genotypes5 = [Genotype([0,1]), Genotype([1,1]), Genotype([1,1]), Genotype([0,0])]
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
	ped = Pedigree(NumericSampleIds(), 2)
	genotypes1 = [Genotype([0,0]), Genotype([0,1]), Genotype([0,0]), Genotype([1,1])]
	gls1 = [
		GenotypeLikelihoods(2,2,[0, 1, 2]),
		GenotypeLikelihoods(2,2,[215, 81, 147]),
		GenotypeLikelihoods(2,2,[199, 49, 253]),
		GenotypeLikelihoods(2,2,[167, 200, 163]),
	]
	genotypes5 = [Genotype([0,1]), Genotype([1,1]), Genotype([1,1]), Genotype([0,0])]
	gls5 = [
		GenotypeLikelihoods(2,2,[184, 71, 233]),
		GenotypeLikelihoods(2,2,[65, 32, 87]),
		GenotypeLikelihoods(2,2,[28, 215, 131]),
		GenotypeLikelihoods(2,2,[98, 250, 137]),
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
