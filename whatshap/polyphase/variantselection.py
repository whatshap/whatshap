"""
Manages phasable variants for genetic polyploid phasing.

The genetic version of Whatshap polyphase is able to phase certain types of variants only,
depending on the genotype. For a given variant table, computes the set of phasable variants
for the specified parents. Can also filter variants if the coverage and coverage rations
between parents and progenies (derived from the allele depths in the variant table) deviate too
much from the average.
"""

import logging

from typing import List

from whatshap.vcf import VariantTable

logger = logging.getLogger(__name__)


class VariantInfo:
    """
    Stores phasable variants in an easily accessible way. Each variant specifies which allele (as a
    number) is the reference (= majority) and which is the alternative (= minority) allele and how
    often the alt allele occurs on the parent and co-parent sample.
    """

    class ParentVariant:
        __slots__ = ("ref", "alt", "alt_count", "co_alt_count")

        def __init__(self, ref, alt, alt_count, co_alt_count):
            self.ref = ref
            self.alt = alt
            self.alt_count = alt_count
            self.co_alt_count = co_alt_count

    def __init__(self, allowed_types):
        self.allowed_types = allowed_types
        self.phasable = set()
        self.variants = []
        self.node_positions = []
        self.nodes_modified = True

    def __getitem__(self, key):
        if isinstance(key, slice):
            raise NotImplementedError("VariantInfo does not support slices")
        assert isinstance(key, int)
        size = len(self.variants)
        if not (-size <= key < size):
            raise IndexError(f"Index out of bounds: {key}")
        if key < 0:
            key = size + key
        return self.variants[key]

    def __len__(self):
        return len(self.variants)

    def append(self, ref, alt, alt_count, co_alt_count, skip=False):
        self.variants.append(self.ParentVariant(ref, alt, alt_count, co_alt_count))
        if not skip and alt is not None and (alt_count, co_alt_count) in self.allowed_types:
            self.phasable.add(len(self.variants) - 1)
            self.nodes_modified = True

    def correct_type(self, index, alt_count=None, co_alt_count=None):
        old_alt = self.variants[index].alt_count
        old_co_alt = self.variants[index].co_alt_count
        changed = False
        if alt_count is not None and old_alt != alt_count:
            changed = True
            if alt_count < 0:
                raise ValueError(f"Cannot set alt count of variant to {alt_count}")
            self.variants[index].alt_count = alt_count
        if co_alt_count is not None and old_co_alt != co_alt_count:
            changed = True
            if co_alt_count < 0:
                raise ValueError(f"Cannot set alt count of variant to {co_alt_count}")
            self.variants[index].co_alt_count = co_alt_count
        if changed:
            if not self.check_variant_compatibility(old_alt, old_co_alt, alt_count, co_alt_count):
                self.remove_phasable(index)
            self.nodes_modified = True

    def get_phasable(self):
        return sorted(list(self.phasable))

    def remove_phasable(self, pos):
        if pos in self.phasable:
            self.phasable.remove(pos)
            self.nodes_modified = True
        else:
            raise ValueError(f"Marked variant {pos} as unphasable, but it was already before")

    def update_node_positions(self):
        self.node_positions = []
        for p in self.get_phasable():
            for i in range(self.variants[p].alt_count):
                self.node_positions.append(p)
        self.nodes_modified = False

    def node_to_variant(self, node_id):
        if self.nodes_modified:
            self.update_node_positions()
        return self.node_positions[node_id]

    def get_node_positions(self):
        if self.nodes_modified:
            self.update_node_positions()
        return self.node_positions[:]

    @staticmethod
    def check_variant_compatibility(old_alt, old_co_alt, new_alt, new_co_alt):
        if old_alt == 1 and old_co_alt == 0:
            return (new_alt, new_co_alt) in [(1, 0), (1, 1), (2, 0)]
        elif old_alt == 1 and old_co_alt == 1:
            return (new_alt, new_co_alt) in [(1, 1)]
        elif old_alt == 2 and old_co_alt == 0:
            return (new_alt, new_co_alt) in [(1, 0), (1, 1), (2, 0)]
        return False


def compute_phasable_variants(
    variant_table: VariantTable, parent: str, co_parent: str, phasing_param
):
    # determine phasable variants
    if phasing_param.complexity_support == 0:
        allowed_pairs = [(1, 0)]
    elif phasing_param.complexity_support == 1:
        allowed_pairs = [(1, 0), (1, 1)]
    else:
        allowed_pairs = [(1, 0), (2, 0), (1, 1)]
    varinfo = VariantInfo(allowed_pairs)

    gts1 = variant_table.genotypes_of(parent)
    gts2 = variant_table.genotypes_of(co_parent)

    for i, var in enumerate(variant_table.variants):
        gt1 = gts1[i]
        gt2 = gts2[i]
        gt1v = gt1.as_vector()
        gt2v = gt2.as_vector()

        if gt1.is_none() or gt2.is_none():
            varinfo.append(None, None, 0, 0)
            continue

        # check homozygosity
        if gt1.is_homozygous():
            varinfo.append(gt1v[0], None, 0, 0)
            continue

        # count different alleles
        alleles_set = set()
        for gt in [gt1v, gt2v]:
            for a in gt:
                alleles_set.add(a)

        alleles = sorted(list(alleles_set))

        if len(alleles) > 2:
            # genotypes are not bi-allelic
            varinfo.append(None, None, 0, 0)
            continue

        assert len(alleles) == 2

        # determine majority allele and its multiplicity
        gt1v.sort()
        ref = gt1v[int(len(gt1v) / 2 - 1)]
        alt = gt1v[0] if gt1v[0] != ref else gt1v[-1]
        alt_count = sum([1 if a == alt else 0 for a in gt1v])
        co_alt_count = sum([1 if a == alt else 0 for a in gt2v])

        skip = False
        if not phasing_param.allow_deletions:
            if "*" in var.get_alt_allele_list():
                skip = True
        varinfo.append(ref, alt, alt_count, co_alt_count, skip)

    return varinfo


def diff_ratio(ratio):
    if ratio and 0.0 < ratio < 1.0:
        return 1.0 / ratio
    else:
        return ratio


def filter_variants(
    varinfo: VariantInfo,
    parent_cov: List[int],
    co_parent_cov: List[int],
    progeny_cov: List[int],
    cutoff: float,
):
    phasable_indices = varinfo.get_phasable()
    co_parent_ratio = [p / s if s > 0 else 0 for p, s in zip(co_parent_cov, parent_cov)]
    progeny_ratio = [p / s if s > 0 else 0 for p, s in zip(progeny_cov, parent_cov)]

    product_ratio = [progeny_ratio[i] * co_parent_ratio[i] for i in phasable_indices]
    median = sorted(product_ratio)[len(product_ratio) // 2]
    product_ratio = [diff_ratio(x / median) for x in product_ratio]

    for i, n in enumerate(phasable_indices):
        if product_ratio[i] > cutoff:
            varinfo.remove_phasable(n)
