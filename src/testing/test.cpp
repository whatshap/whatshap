#include "../readset.h"
#include "../read.h"
#include "../pedigree.h"
#include "../phredgenotypelikelihoods.h"
#include "../genotypedptable.h"
#include "../pedigree.h"
#include "../genotypecolumncostcomputer.h"
#include "../columniterator.h"
#include "../pedigreepartitions.h"
#include "../backwardcolumniterator.h"
#include "../entry.h"
#include "../transitionprobabilitycomputer.h"
#include "../vector2d.h"

#include <iostream>
#include <string>
#include <vector>
#include <list>
#include <sstream>

#define CATCH_CONFIG_MAIN
#include "catch.hpp"

using namespace std;

size_t popcount(size_t x) {
    unsigned int count = 0;
    for (;x; x >>= 1) {
        count += x & 1;
    }
    return count;
}


ReadSet* string_to_readset(string s, string weights, bool use){
    ReadSet* read_set = new ReadSet;
    stringstream s1(s);
    stringstream s2(weights);
    string line;
    string line_weights;

    if(weights != ""){

    }
    unsigned int index = 0;
    while((std::getline(s1,line,'\n')) && (std::getline(s2,line_weights,'\n'))){
        if(line.length() == 0) continue;
        unsigned int counter = 0;
        Read* read = new Read("read"+std::to_string(index), 50, 0,0);
        for(unsigned int i = 0; i < line.length(); i++){
            counter += 1;
            if(line[i] == ' ') continue;
            unsigned int quality = int(line_weights[i] - '0');
            if(!use){
                read->addVariant((counter)*10,int(line[i] - '0'), quality);
            } else {
                read->addVariant((counter)*10,int(line[i] - '0'), 10);
            }

        }

        read_set->add(read);
        index += 1;
    }
    return read_set;
}

// extract the columns of the matrix as strings
vector<string> get_columns(string matrix, unsigned int col_count){
    stringstream ss(matrix);
    string line;
    vector<string> result(col_count,"");

    unsigned int index = 0;
    while(std::getline(ss,line,'\n')){
        for(unsigned int i = 0; i < col_count; i++){
            if(line[i] != ' ') result[i]+=line[i];
        }
    }

    return result;
}

long double naive_column_cost_computer(string current_column,unsigned int bipartition, unsigned int switch_cost, unsigned int allele1, unsigned int allele2){

    long double result = 1.0L;

    for(unsigned int j = 0; j < current_column.length(); j++){

        // check to which partition the current read belongs
        // check the entry in j-th bit
        bool part2 = ((unsigned int)1 << j) & bipartition;
        if(part2){
            if(int(current_column[j] - '0') == allele2){
                result *= 1-pow(10,-(long double)switch_cost/10.0L);
            } else {
                result *= pow(10,-(long double)switch_cost/10.0L);
            }
        } else {
            if(int(current_column[j] - '0') == allele1){
                result *= 1-pow(10,-(long double)switch_cost/10.0L);
            } else {
                result *= pow(10,-(long double)switch_cost/10.0L);
            }
        }
    }
    return result;
}


// compare vector of entries to string
bool compare_entries(vector<const Entry*> c1, string c2){
    bool result = true;

    //for(unsigned int i = 0; i < c1.size(); i++)
    unsigned int i = 0;
    unsigned int j = 0;
    while((i<c1.size()) && (j<c2.length())){
        switch(c1[i]->get_allele_type()){
        case Entry::REF_ALLELE: if(c2[j] != '0'){result = false;} else {i+=1;j+=1;} break;
        case Entry::ALT_ALLELE: if(c2[j] != '1'){result = false;} else {i+=1;j+=1;} break;
        case Entry::BLANK: i += 1; break;
        default: break;
        }
    }
    return result;
}


TEST_CASE("test forward_backward", "[test forward_backward]")
{
    vector<string> reads = {"11 \n 01","1111111111\n0000011111", "101\n101\n101\n101\n100\n100\n100\n100","01\n11" ,"011011\n110110\n110 10\n110110\n111110\n000 00\n01000 \n000010\n100100"};
    vector<string> weights = {"11 \n 11","1111111111\n1111111111","111\n111\n111\n111\n111\n111\n111\n111","11\n11", "111111\n111111\n111111\n111111\n111111\n111111\n111111\n111111\n111111"};
    vector< vector<unsigned int> > expected_geno = { {1,1,1} , {1,1,1,1,1,2,2,2,2,2}, {2,0,1} ,{1,2}, {1,1,0,1,1,0}};


    SECTION("test single individual1", "[test single individual1]")
    {
        ReadSet* read_set = string_to_readset(reads[0],weights[0], true);
        std::vector<unsigned int>* positions = read_set->get_positions();
        // no genotype likelihoods needed for genotyping mode
        std::vector<PhredGenotypeLikelihoods*> genotype_likelihoods(positions->size(),nullptr);
        // recombination costs set to 1
        std::vector<unsigned int> recombcost(positions->size(), 1);
        // add a single individual to the pedigree
        Pedigree* pedigree = new Pedigree;
        pedigree->addIndividual(0, std::vector<unsigned int >(positions->size(),1), genotype_likelihoods);
        // construct the dp table and perform forward-backward algorithm
        GenotypeDPTable dp_table(read_set, recombcost, pedigree, positions);
        // get the computed likelihoods
        for(unsigned int pos = 0; pos < 3; pos++){
            std::vector<long double> likelihoods = dp_table.get_genotype_likelihoods(0,pos);
        }

        delete read_set;
        delete pedigree;
        delete positions;
    }



    SECTION("test single individual2", "[test single individual2]")
    {
        std::vector<long double> likelihoods;
        for(unsigned int i = 0; i < reads.size(); i++){
            ReadSet* read_set = string_to_readset(reads[i],weights[i], true);
            std::vector<unsigned int>* positions = read_set->get_positions();
            // no genotype likelihoods needed for genotyping mode
            std::vector<PhredGenotypeLikelihoods*> genotype_likelihoods(positions->size(),nullptr);
            // recombination costs set to 1
            std::vector<unsigned int> recombcost(positions->size(), 1);
            // add a single individual to the pedigree
            Pedigree* pedigree = new Pedigree;
            pedigree->addIndividual(0, std::vector<unsigned int >(positions->size(),1), genotype_likelihoods);
            // construct the dp table and perform forward-backward algorithm
            GenotypeDPTable dp_table(read_set, recombcost, pedigree, positions);
            // get the computed likelihoods
            std::vector<long double> likelihoods = dp_table.get_genotype_likelihoods(0,0);

            for(unsigned int pos = 0; pos < positions->size(); pos++){
                std::vector<long double> likelihoods = dp_table.get_genotype_likelihoods(0,pos);
                int pos_max = -1;
                long double cur_max = -1.0L;
                for(unsigned int u = 0; u < likelihoods.size(); u++){
                    if(likelihoods[u] > cur_max){
                        cur_max = likelihoods[u];
                        pos_max = u;
                    }

                }
                REQUIRE(pos_max == expected_geno[i][pos]);
            }

            delete read_set;
            delete pedigree;
            delete positions;

        }

    }

}

TEST_CASE("test transition prob computer", "[test transition prob computer]"){

    SECTION("test simple example", "[test simple example]"){
        TransitionProbabilityComputer trans(10,1,16,4);
        std::vector<long double> expected_cost = {0.9L*0.9L, 0.1L*0.9L, 0.1L*0.1L};
        long double nor = (0.9L*0.9L+2*0.1L*0.9L+0.1L*0.1L)*16;

        for(unsigned int i = 0; i < 4; i++){
            long double row_sum = 0.0L;
            for(unsigned int j = 0; j < 4; j++){
                unsigned int index = popcount(i ^ j);
                REQUIRE((float)(expected_cost[index]/nor) == (float)trans.get(i,j));
                row_sum += trans.get(i,j)*16.0L;
            }
            REQUIRE(float(row_sum) == 1.0);
        }
    }

    SECTION("test for single individual", "[test for single individual]"){
        TransitionProbabilityComputer trans(10,0,4,1);
        REQUIRE(trans.get(0,0) == 0.25);
    }
}


TEST_CASE("test column_cost_computer","[test column_cost_computer]"){

    vector<std::string> reads = {"11\n00", "10\n11", "00\n00", "10\n10", "01\n10"};
    std::string weights = "11\n11";

    for(unsigned int r = 0; r < reads.size(); r++){
        ReadSet* read_set = string_to_readset(reads[r],weights,false);
        std::vector<unsigned int>* positions = read_set->get_positions();
        std::vector<PhredGenotypeLikelihoods*> genotype_likelihoods(positions->size(),nullptr);
        std::vector<unsigned int> recombcost(positions->size(), 1);
        Pedigree* pedigree = new Pedigree;
        pedigree->addIndividual(0, std::vector<unsigned int >(positions->size(),1), genotype_likelihoods);

        // create all pedigree partitions
        std::vector<PedigreePartitions*> pedigree_partitions;
        for(size_t i = 0; i < std::pow(4,pedigree->triple_count()); ++i)
        {
            pedigree_partitions.push_back(new PedigreePartitions(*pedigree,i));
        }

        // translate all individual ids to individual indices
         std::vector<unsigned int> read_sources;
        for(size_t i = 0; i<read_set->size(); ++i)
        {
            read_sources.push_back(pedigree->id_to_index(read_set->get(i)->getSampleID()));
        }
        vector<string> columns = get_columns(reads[r],2);

        ColumnIterator input_column_iterator(*read_set, positions);
        unsigned int col_ind = 0;

        while(input_column_iterator.has_next()){
            unique_ptr<vector<const Entry *> > current_input_column = input_column_iterator.get_next();

            // create column cost computer
            GenotypeColumnCostComputer cost_computer(*current_input_column, col_ind, read_sources, pedigree,*pedigree_partitions[0]);
            cost_computer.set_partitioning(0);

            unsigned int switch_cost = 1;

            // check if costs for initial partition (r1,r2/.) are computed correctly
            REQUIRE(cost_computer.get_cost(0) == naive_column_cost_computer(columns[col_ind],0,switch_cost,0,0));
            REQUIRE(cost_computer.get_cost(1) == naive_column_cost_computer(columns[col_ind],0,switch_cost,0,1));
            REQUIRE(cost_computer.get_cost(2) == naive_column_cost_computer(columns[col_ind],0,switch_cost,1,0));
            REQUIRE(cost_computer.get_cost(3) == naive_column_cost_computer(columns[col_ind],0,switch_cost,1,1));


            // switch first read (r2/r1)
            cost_computer.update_partitioning(0);
            REQUIRE(cost_computer.get_cost(0) == naive_column_cost_computer(columns[col_ind],1,switch_cost,0,0));
            REQUIRE(cost_computer.get_cost(1) == naive_column_cost_computer(columns[col_ind],1,switch_cost,0,1));
            REQUIRE(cost_computer.get_cost(2) == naive_column_cost_computer(columns[col_ind],1,switch_cost,1,0));
            REQUIRE(cost_computer.get_cost(3) == naive_column_cost_computer(columns[col_ind],1,switch_cost,1,1));

            // switch also second read (./r1,r2)
            cost_computer.update_partitioning(1);
            REQUIRE(cost_computer.get_cost(0) == naive_column_cost_computer(columns[col_ind],3,switch_cost,0,0));
            REQUIRE(cost_computer.get_cost(1) == naive_column_cost_computer(columns[col_ind],3,switch_cost,0,1));
            REQUIRE(cost_computer.get_cost(2) == naive_column_cost_computer(columns[col_ind],3,switch_cost,1,0));
            REQUIRE(cost_computer.get_cost(3) == naive_column_cost_computer(columns[col_ind],3,switch_cost,1,1));

            // test partition (r1/r2)
            cost_computer.set_partitioning(2);
            REQUIRE(cost_computer.get_cost(0) == naive_column_cost_computer(columns[col_ind],2,switch_cost,0,0));
            REQUIRE(cost_computer.get_cost(1) == naive_column_cost_computer(columns[col_ind],2,switch_cost,0,1));
            REQUIRE(cost_computer.get_cost(2) == naive_column_cost_computer(columns[col_ind],2,switch_cost,1,0));
            REQUIRE(cost_computer.get_cost(3) == naive_column_cost_computer(columns[col_ind],2,switch_cost,1,1));

            col_ind += 1;
        }

        delete read_set;
        delete positions;
        delete pedigree;
        for(unsigned int i=0; i < pedigree_partitions.size(); i++){
            delete pedigree_partitions[i];
        }
    }
}

TEST_CASE("test BackwardColumnIterator", "[test BackwardColumnIterator]") {

    SECTION("test small examples", "[test small examples]"){
        std::vector<std::string> matrices = {"10 \n010\n000", "01 \n000\n111", "0 1\n1 0\n 11"};
        std::vector<std::string> weights = {"11 \n111\n111", "11 \n111\n111", "1 1\n1 1\n 11"};

        for(unsigned int i = 0; i < matrices.size(); i++){
            ReadSet* read_set = string_to_readset(matrices[i], weights[i],false);
            vector<string> columns = get_columns(matrices[i],3);

            std::vector<unsigned int>* positions = read_set->get_positions();
            BackwardColumnIterator col_it(*read_set, positions);
            REQUIRE(col_it.has_next());

            // iterate backwards from end to start
            for(int j = 2; j >= 0; j--){
                auto col = col_it.get_next();
                REQUIRE(compare_entries(*col,columns[j]));
                if(j>0){
                    REQUIRE(col_it.has_next());
                } else {
                    REQUIRE(!col_it.has_next());
                }
            }

            // use jump column to iterate from end to start
            for(int j = 2; j >= 0; j--){
                col_it.jump_to_column(j);
                auto col = col_it.get_next();
                REQUIRE(compare_entries(*col,columns[j]));
                if(j>0){
                    REQUIRE(col_it.has_next());
                } else {
                    REQUIRE(!col_it.has_next());
                }
            }

            // use jump column to iterate from start to end
            for(int j = 0; j < 3; j++){
                col_it.jump_to_column(j);
                auto col = col_it.get_next();
                REQUIRE(compare_entries(*col,columns[j]));
            }

            delete read_set;
            delete positions;
        }
    }
}

TEST_CASE("test scaling of vector", "[test scaling of vector]"){
    Vector2D<long double> test(2,3,0.8L);
    test.divide_entries_by(0.8L);

    for(unsigned int i = 0; i < test.get_size0(); i++){
        for(unsigned int j = 0; j < test.get_size1(); j++) {
            REQUIRE(test.at(i,j) == 1L);
        }
    }
}
