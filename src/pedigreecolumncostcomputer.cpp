#include <cassert>
#include <limits>
#include <utility>
#include <algorithm>
#include <array>
#include <map>
#include <unordered_set>
#include "pedigreecolumncostcomputer.h"

using namespace std;

PedigreeColumnCostComputer::PedigreeColumnCostComputer(const std::vector <const Entry *>&column,const std::vector <unsigned int>& read_marks,
						       unsigned int inheritance_val, std::vector <Pedigree::triple_entry_t > triples, std::vector <unsigned int> id_of_individuals):
column(column),
read_marks(read_marks),
inheritance_val(inheritance_val),
partitioning(0),
num_of_triples(triples.size()),
triples(triples),
id_of_individuals(id_of_individuals)
{
    if (triples.empty()) {
	haps[id_of_individuals[0]] = std::make_pair(0, 1);
    } else {
	roots = compute_roots(triples);
	unsigned int    temp_inheritance = inheritance_val;
      for (auto & triple:triples) {
	    unsigned int    triple_inheritance = temp_inheritance & 3;
	    // mother and father
	    for (int i = 0; i < 2; i++) {
		unsigned int    m = triple[i];
		if (haps.find(m) == haps.end()) {
		    unsigned int    pos = 0;	
		    for (unsigned int j = 0; j < roots.size(); j++) {
			if (roots[j] == m) {
			    pos = j;
			    break;
			}
		    }
		    if (pos >= 0) {
			haps[m] = std::make_pair(2 * pos, 2 * pos + 1);
		    }
		}
	    }
	    // child
	    unsigned int    c = triple[2];
	    unsigned int    haps1;
	    unsigned int    haps2;
	    if ((triple_inheritance / 2 == 0)
		&& (triple_inheritance % 2 == 0)) {
		haps1 = haps[triple[0]].first;
		haps2 = haps[triple[1]].first;
	    } else if ((triple_inheritance / 2 == 0)
		       && (triple_inheritance % 2 == 1)) {
		haps1 = haps[triple[0]].second;
		haps2 = haps[triple[1]].first;
	    } else if ((triple_inheritance / 2 == 1)
		       && (triple_inheritance % 2 == 0)) {
		haps1 = haps[triple[0]].first;
		haps2 = haps[triple[1]].second;
	    } else if ((triple_inheritance / 2 == 1)
		       && (triple_inheritance % 2 == 1)) {
		haps1 = haps[triple[0]].second;
		haps2 = haps[triple[1]].second;
	    }

	    haps[c] = std::make_pair(haps1, haps2);

	    temp_inheritance = (temp_inheritance >> 2);
	}
    }
}

std::vector < unsigned int > PedigreeColumnCostComputer::compute_roots(std::vector <Pedigree::triple_entry_t > triples)
{
    std::set < unsigned int >parents;
    std::set < unsigned int >children;
    std::set < unsigned int >intersect;
    std::vector < unsigned int >roots;

  for (auto & triple:triples) {
	parents.insert(triple[0]);
	parents.insert(triple[1]);
	children.insert(triple[2]);
    }
    std::set_intersection(parents.begin(), parents.end(), children.begin(), children.end(), std::inserter(intersect, intersect.end()));
    std::set_difference(parents.begin(), parents.end(), intersect.begin(), intersect.end(), std::inserter(roots, roots.begin()));
    std::sort(roots.begin(), roots.end());
    return roots;
}

void PedigreeColumnCostComputer::set_partitioning(unsigned int partitioning)
{
    // compute cost from scratch

    unsigned int    border;
    if (triples.empty())
	border = 1;
    else
	border = 2 * roots.size() - 1;
    std::array < unsigned int, 2 > null_array { { 0, 0}};
    for (unsigned int i = 0; i <= border; i++) {
	cost_partition.push_back(null_array);
    }

    partitioning = partitioning;
    for (vector < const Entry * >::const_iterator it = column.begin();
	 it != column.end(); ++it) {
	auto & entry = **it;
	bool  entry_in_partition1 = (partitioning & ((unsigned int) 1)) == 0;
	unsigned int    ind_id = read_marks[entry.get_read_id()];
	switch (entry.get_allele_type()) {

	case Entry::MAJOR_ALLELE:
	    (entry_in_partition1 ? cost_partition[haps[ind_id].first] :cost_partition[haps[ind_id].second])[1] += entry.get_phred_score();
	    break;
	case Entry::MINOR_ALLELE:
	    (entry_in_partition1 ? cost_partition[haps[ind_id].first] :cost_partition[haps[ind_id].second])[0] += entry.get_phred_score();
	    break;
	case Entry::BLANK:
	    break;
	default:
	    assert(false);
	}


	partitioning = partitioning >> 1;
    }
}

void PedigreeColumnCostComputer::update_partitioning(int bit_to_flip)
{
    // update cost based on the changed bit

    const Entry & entry = *column[bit_to_flip];
    partitioning = partitioning ^ (((unsigned int) 1) << bit_to_flip);
    bool entry_in_partition1 = (partitioning & (((unsigned int) 1) << bit_to_flip)) == 0;
    unsigned int ind_id = read_marks[entry.get_read_id()];
    switch (entry.get_allele_type()) {
    case Entry::MAJOR_ALLELE:
	(entry_in_partition1 ? cost_partition[haps[ind_id].second] : cost_partition[haps[ind_id].first])[1] -= entry.get_phred_score();
	(entry_in_partition1 ? cost_partition[haps[ind_id].first] :  cost_partition[haps[ind_id].second])[1] += entry.get_phred_score();
	break;
    case Entry::MINOR_ALLELE:
	(entry_in_partition1 ? cost_partition[haps[ind_id].second] : cost_partition[haps[ind_id].first])[0] -= entry.get_phred_score();
	(entry_in_partition1 ? cost_partition[haps[ind_id].first] :  cost_partition[haps[ind_id].second])[0] += entry.get_phred_score();
	break;
    case Entry::BLANK:
	break;
    default:
	assert(false);
    }
}

unsigned int PedigreeColumnCostComputer::get_cost(std::vector <Pedigree::triple_entry_t > &genotypes)
{
    unsigned int    best_cost = numeric_limits < unsigned int >::max();
    // Enumerate all possible assignments of alleles to haplotypes and 
    // compute costs for those which are compatible with genotypes.
    // TODO: This can be done more efficiently.
    unsigned int    border;
    if (triples.empty())
	border = 4;
    else
	border = std::pow(4, roots.size());


    for (unsigned int i = 0; i < border; ++i) {
	int  flagcounter = 0;

	std::map < unsigned, std::pair < unsigned int, unsigned int  >>enumerate_haps;
	std::vector < unsigned int >duplicate_roots = roots;
	unsigned int    cost = 0;
	unsigned int    counter = 0;
	unsigned int    temp_inheritance = inheritance_val;
	if (triples.empty()) {
	    unsigned int    allele1 = ((i >> counter) & 1);
	    unsigned int    allele2 = ((i >> (counter + 1)) & 1);
	    counter = counter + 2;
	    if (allele1 + allele2 == genotypes[0][0]) {
		flagcounter++;
		enumerate_haps[id_of_individuals[0]] = std::make_pair(allele1, allele2);
	    }

	} else {
	    for (unsigned int j = 0; j < triples.size(); j++) {

		unsigned int    triple_inheritance = temp_inheritance & 3;
		for (unsigned int k = 0; k <= 2; k++) {
		    if (std::find(duplicate_roots.begin(), duplicate_roots.end(), triples[j][k]) != duplicate_roots.end()) {
			duplicate_roots.erase(std::remove(duplicate_roots.begin(), duplicate_roots.end(), triples[j][k]), duplicate_roots.end());
			unsigned int    allele1 = ((i >> counter) & 1);
			unsigned int    allele2 = ((i >> (counter + 1)) & 1);
			counter = counter + 2;
			if (allele1 + allele2 == genotypes[j][k]) {
			    flagcounter++;
			    enumerate_haps[triples[j][k]] = std::make_pair(allele1, allele2);
			}
		    } else {
			if (enumerate_haps.find(triples[j][k]) == enumerate_haps.end()) {
			    unsigned int    allele_c1 =((triple_inheritance & 1) == 0) ? enumerate_haps[triples[j][0]].first : enumerate_haps[triples[j][0]].second;
			    unsigned int    allele_c2 =(((triple_inheritance >> 1) & 1) == 0) ? enumerate_haps[triples[j][1]].first : enumerate_haps[triples[j][1]].second;
			    if (allele_c1 + allele_c2 == genotypes[j][k]) {
				flagcounter++;
				enumerate_haps[triples[j][k]] = std::make_pair(allele_c1, allele_c2);
			    }			   
			}
		    }
		}
		temp_inheritance = (temp_inheritance >> 2);

	    }
	}
	if (flagcounter == id_of_individuals.size()) {	
	    if (triples.empty()) {
		unsigned int    temp1 = enumerate_haps[id_of_individuals[0]].first;
		unsigned int    temp2 = enumerate_haps[id_of_individuals[0]].second;
		cost =  cost + cost_partition[0][temp1] + cost_partition[1][temp2];
	    } else {
		for (unsigned int p = 0; p < roots.size(); p++) {
		    unsigned int    temp1 = enumerate_haps[roots[p]].first;
		    unsigned int    temp2 = enumerate_haps[roots[p]].second;
		    cost += cost_partition[2 * p][temp1];
		    cost += cost_partition[2 * p + 1][temp2];
		}
	    }
	    if (cost < best_cost) {
		best_cost = cost;
	    }
	}
    }

    if (best_cost == numeric_limits < unsigned int >::max()) {
	// throw std::runtime_error("Error: Mendelian conflict");
	best_cost = 10000000;	// hard code here... else overflow
    }

    return best_cost;
}

std::map < unsigned int,std::pair <Entry::allele_t,Entry::allele_t >> PedigreeColumnCostComputer::get_alleles(std::vector <Pedigree::triple_entry_t > &genotypes)
{
    // TODO: avoid code duplication
    unsigned int    best_cost = numeric_limits < unsigned int >::max();
    unsigned int    second_best_cost = numeric_limits < unsigned int >::max();
    std::map < unsigned int, std::pair < Entry::allele_t, Entry::allele_t >> pop_haps;
    unsigned int    border;
    if (triples.empty())
	border = 4;
    else
	border = std::pow(4, roots.size());


    for (unsigned int i = 0; i < border; ++i) {
	int flagcounter = 0;
	std::map < unsigned, std::pair < unsigned int, unsigned int  >>enumerate_haps;
	std::vector < unsigned int >duplicate_roots = roots;
	unsigned int    cost = 0;
	unsigned int    counter = 0;
	unsigned int    temp_inheritance = inheritance_val;
	if (triples.empty()) {
	    unsigned int    allele1 = ((i >> counter) & 1);
	    unsigned int    allele2 = ((i >> (counter + 1)) & 1);

	    counter = counter + 2;
	    if (allele1 + allele2 == genotypes[0][0]) {
		flagcounter++;
		enumerate_haps[id_of_individuals[0]] = std::make_pair(allele1, allele2);
	    }

	} else {
	    for (unsigned int j = 0; j < triples.size(); j++) {
		unsigned int    triple_inheritance = temp_inheritance & 3;
		for (unsigned int k = 0; k <= 2; k++) {
		    if (std::find(duplicate_roots.begin(), duplicate_roots.end(), triples[j][k]) != duplicate_roots.end()) {
			duplicate_roots. erase(std:: remove(duplicate_roots.begin(), duplicate_roots.end(), triples[j][k]), duplicate_roots.end());
			unsigned int    allele1 = ((i >> counter) & 1);
			unsigned int    allele2 =   ((i >> (counter + 1)) & 1);
			counter = counter + 2;
			if (allele1 + allele2 == genotypes[j][k]) {
			    flagcounter++;
			    enumerate_haps[triples[j][k]] = std::make_pair(allele1, allele2);
			}
		    } else {
			if (enumerate_haps.find(triples[j][k]) == enumerate_haps.end()) {
			    unsigned int    allele_c1 =	((triple_inheritance & 1) == 0) ? enumerate_haps[triples[j][0]].first : enumerate_haps[triples[j][0]].second;
			    unsigned int    allele_c2 =(((triple_inheritance >> 1) & 1) == 0) ? enumerate_haps[triples[j][1]].first : enumerate_haps[triples[j][1]].second;
			    if (allele_c1 + allele_c2 == genotypes[j][k]) {
				flagcounter++;
				enumerate_haps[triples[j][k]] = std::make_pair(allele_c1, allele_c2);
			    }
			}
		    }
		}
		temp_inheritance = (temp_inheritance >> 2);

	    }
	}
	if (flagcounter == id_of_individuals.size()) {	
	    if (triples.empty()) {
		unsigned int    temp1 =   enumerate_haps[id_of_individuals[0]].first;
		unsigned int    temp2 =  enumerate_haps[id_of_individuals[0]].second;
		cost = cost + cost_partition[0][temp1] +  cost_partition[1][temp2];
	    } else {
		for (unsigned int p = 0; p < roots.size(); p++) {
		    unsigned int    temp1 = enumerate_haps[roots[p]].first;
		    unsigned int    temp2 = enumerate_haps[roots[p]].second;
		    cost += cost_partition[2 * p][temp1];
		    cost += cost_partition[2 * p + 1][temp2];
		}
	    }
	    if (cost < best_cost) {
		second_best_cost = best_cost;
		best_cost = cost;
	      for (auto & item:enumerate_haps) {
		    unsigned int    key = item.first;
		    unsigned int    temp1 = item.second.first;
		    unsigned int    temp2 = item.second.second;
		    pop_haps[key] = std::make_pair((temp1 ==  0) ? Entry:: MAJOR_ALLELE : Entry:: MINOR_ALLELE, (temp2 == 0) ? Entry:: MAJOR_ALLELE : Entry:: MINOR_ALLELE);
		}

	    }
	}
    }


    if (best_cost == numeric_limits < unsigned int >::max()) {
	throw           std::runtime_error("Error: Mendelian conflict");
    }

    if (second_best_cost == best_cost) {
	// TODO: Having too equal scores does not necissarily imply that
	// the case is
	// undecidable for all three individuals. Be smarter here.

      for (auto & item:pop_haps) {
	    unsigned int    key = item.first;
	    pop_haps[key] =
		std::make_pair(Entry::EQUAL_SCORES, Entry::EQUAL_SCORES);
	}
	return pop_haps;
    }

    return pop_haps;
}

unsigned int PedigreeColumnCostComputer::get_weight(bool second_haplotype)
{
    throw           std::runtime_error("Not yet implemented for trios.");
}