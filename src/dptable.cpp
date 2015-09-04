#include <stdexcept>
#include <cassert>
#include <limits>
#include <fstream>

#include "columncostcomputer.h"
#include "dptable.h"

//#define DB // db
//#define FF_PARALLEL

#ifdef COLUMN_TIME
#include <sys/time.h>
static inline unsigned long getusec() {
    struct timeval tv;
    gettimeofday(&tv,NULL);
    return (unsigned long)(tv.tv_sec*1e6+tv.tv_usec);
}
#endif

using namespace std;

#ifdef FF_PARALLEL
#include "ff/mapping_utils.hpp"

#ifdef FF_ALLOC
#include "ff/allocator.hpp"
#define ALLOCATOR_INIT()
#define MALLOC(size)   (FFAllocator::instance()->malloc(size))
#define FREE(ptr,unused) (FFAllocator::instance()->free(ptr))
else
#define ALLOCATOR_INIT()
#define MALLOC(size)   (malloc(size))
#define FREE(ptr,unused) (free(ptr))
#endif


DPTable::DPTable(ReadSet* read_set, bool all_heterozygous,int numthreads, long chunksize, unsigned int threshold):  pf(numthreads,true), numthreads(numthreads), chunksize(chunksize), threshold(threshold) {
    this->read_set = read_set;
    this->all_heterozygous = all_heterozygous;
    this->read_count = 0;
    if (this->numthreads<=0) this->numthreads=ff_realNumCores(); // Num of cores //ff_numCores(); - Number of Hyper-Contexts in the platform
    std::cerr << "This is the FastFlow parallel version, max parallelism is set to " << this->numthreads << "\n";
    read_set->reassignReadIds();
    compute_table();
}
#else
DPTable::DPTable(ReadSet* read_set, bool all_heterozygous) {
    this->read_set = read_set;
    this->all_heterozygous = all_heterozygous;
    this->read_count = 0;
    read_set->reassignReadIds();
    compute_table();
}
#endif

DPTable::~DPTable() {
    for(size_t i=0; i<indexers.size(); ++i) {
        delete indexers[i];
    }
    
    for(size_t i=0; i<backtrace_table.size(); ++i) {
        delete backtrace_table[i];
    }
}

shared_ptr<vector<unsigned int> > DPTable::extract_read_ids(const vector<const Entry *>& entries) {
    shared_ptr<vector<unsigned int> > read_ids(new vector<unsigned int>());
    for (size_t i=0; i<entries.size(); ++i) {
        read_ids->push_back(entries[i]->get_read_id());
    }
    return read_ids;
}

// db
#ifdef DB
// helper function to output the bit representation of an unsigned int
string bit_rep(unsigned int a, unsigned int len) {
    
    string s;
    for(int i=0; i< len; ++i) {
        s = ((a&1)?"1":"0") + s;
        a = a >> 1;
    }
    
    return s;
}

void output_vector(const vector<unsigned int> * v) {
    for(int j=v->size()-1; j>= 0; --j) {
        if(v->at(j) == -1) cout << "_ ";
        else cout << v->at(j) << " ";
    }
}

void output_vector_enum(const vector<unsigned int> * v, unsigned int len) {
    for(int j = v->size()-1; j >= 0; --j) {
        cout << j << " [" << bit_rep(j,len) << "] : ";
        if(v->at(j) == -1) cout << "_";
        else cout << v->at(j);
        cout << endl;
    }
}
#endif

void DPTable::compute_table() {
    ColumnIterator column_iterator(*read_set);
    if(!indexers.empty()) { // clear indexers, if present
        for(size_t i=0; i<indexers.size(); ++i) {
            delete indexers[i];
        }
        indexers.resize(0);
    }
    
    if(!backtrace_table.empty()) { // clear backtrace_table, if present
        for(size_t i=0; i<backtrace_table.size(); ++i) {
            delete backtrace_table[i];
        }
        backtrace_table.resize(0);
    }
    
    // empty read-set, nothing to phase, so MEC score is 0
    if (!column_iterator.has_next()) {
        optimal_score = 0;
        optimal_score_index = 0;
        return;
    }
    
    unsigned int n = 0;
    shared_ptr<vector<const Entry *> > current_column;
    shared_ptr<vector<const Entry *> > next_column;
    // get the next column ahead of time
    next_column = column_iterator.get_next();
    shared_ptr<vector<unsigned int> > next_read_ids = extract_read_ids(*next_column);
    ColumnIndexingScheme* next_indexer = new ColumnIndexingScheme(0,*next_read_ids);
    indexers.push_back(next_indexer);
    shared_ptr<vector<unsigned int> > previous_projection_column; // shoul be shared
    shared_ptr<vector<unsigned int> > current_projection_column;
    unsigned int running_optimal_score;
    unsigned int running_optimal_score_index; // optimal score and its index
    double pi = 0.05; // percentage of columns processed
    double pc = pi;
    unsigned int nc = column_iterator.get_column_count();
    
#ifdef FF_PARALLEL
    std::vector<ColumnIndexingIterator*> thiterators(numthreads);
    std::vector<ColumnCostComputer*>     thcostcomputers(numthreads);
    
    for(size_t i=0;i<numthreads;++i) {
        thiterators[i]     = new ColumnIndexingIterator;
        thcostcomputers[i] = new ColumnCostComputer;
    }
#endif

#ifdef DB
    int i = 0;
#endif
    
    long iteration_count = 0;
    #define MAXITER 500
    
    while ((++iteration_count<MAXITER) && (next_indexer != 0)) {
        std::cerr << " Iter " << iteration_count << "\n";
        // move on projection column
        previous_projection_column = std::move(current_projection_column);
        // make former next column the current one
        current_column = std::move(next_column);
        shared_ptr<vector<unsigned int> > current_read_ids = std::move(next_read_ids);
        ColumnIndexingScheme* current_indexer = next_indexer;
        // peek ahead and get the next column
        if (column_iterator.has_next()) {
            next_column = column_iterator.get_next();
            next_read_ids = extract_read_ids(*next_column);
            next_indexer = new ColumnIndexingScheme(current_indexer,*next_read_ids);
            current_indexer->set_next_column(next_indexer);
            indexers.push_back(next_indexer);
        } else {
            assert(next_column.get() == 0);
            assert(next_read_ids.get() == 0);
            read_count = column_iterator.get_read_count();
            next_indexer = 0;
        }
        // reserve memory for the DP column
        vector<unsigned int> dp_column(current_indexer->column_size(),0);
        vector<unsigned int>* backtrace_column = 0;
        // if not last column, reserve memory for forward projections column
        if (next_column.get() != 0) {
#ifdef DB
            cout << i << " : " << endl;
            ++i;
            cout << "allocate current projection / backtrace columns of size : " << current_indexer->forward_projection_size() << endl;
            cout << "forward projection width : " << current_indexer->get_forward_projection_width() << endl << endl;
#endif
            
            
            current_projection_column = shared_ptr<vector<unsigned int> >(new vector<unsigned int>(current_indexer->forward_projection_size(), numeric_limits<unsigned int>::max()));
            // NOTE: forward projection size will always be even
            
            backtrace_column = new vector<unsigned int>(current_indexer->forward_projection_size(), numeric_limits<unsigned int>::max());
            
        }
#ifdef FF_PARALLEL
        {
#if defined(COLUMN_TIME)
            long a = ff::getusec();
#endif
            shared_ptr<ColumnIndexingIterator> iterator = current_indexer->get_iterator();
            {
            
            for(size_t i=0;i<numthreads;++i) {
                thiterators[i]->reset(current_indexer);
                thcostcomputers[i]->reset(&(*current_column), all_heterozygous);
            }
            const bool notlast = (next_column.get() != 0);
            
            ssize_t l = iterator->get_length();
            ssize_t pardegree = 1;
            //if (l>=8) pardegree=numthreads;
            if (l>=262144) pardegree=numthreads; // 18
            else if (l>=65536)  pardegree=std::min(numthreads, 4); // 16
            else if (l>=4096) pardegree=std::min(numthreads, 2); // 12
            //else if (l>=8) pardegree=std::min(numthreads, 2); // 3
            
            long adj_chunksize;
            if (chunksize>(l/pardegree))
                adj_chunksize = l/(pardegree);
            else
                adj_chunksize = l/(2*pardegree);//chunksize;
            std::cerr << "Running with: pardegree " << pardegree << " grain set to " << adj_chunksize << " MAX numthreads " << numthreads << "\n";
        
            
            auto Map = [&thiterators,&thcostcomputers,previous_projection_column,current_indexer](const long start, const long stop, const int thid, ff::ff_buffernode &node) {
                //std::cerr << "PF map body thread" << thid << " (native " << ff_getThreadID() << " ) ";
                //std::cerr << "start " << start << " stop " << stop << "\n";
                //std::cerr.flush();
                
                if (start == stop) return;
                //std::cerr << "PF map body thread" << thid << " (native " << ff_getThreadID() << " )\n";
                
                ColumnIndexingIterator *myiterator                 = thiterators[thid];
                ColumnCostComputer &mycost_computer                = *(thcostcomputers[thid]);
                
                unsigned int cost, forward_index;
                int bit_changed = -1;
                if (start>0) {
                    myiterator->advance_idx(&bit_changed,start-1);
                    mycost_computer.set_partitioning(myiterator->get_partition());
                    //std::cerr << "Partition " << myiterator->get_partition() << "\n";
                }
                
                
                std::vector<ff_task_t> *thcurrent_projection_column = new std::vector<ff_task_t>;
                thcurrent_projection_column->reserve(stop-start);
                
                
                for(long k=start; k<stop; ++k) {
                    myiterator->advance_idx(&bit_changed,k);
                    
                    if (bit_changed >= 0) mycost_computer.update_partitioning(bit_changed);
                    else mycost_computer.set_partitioning(0);
                    
                    cost=0;
                    if (previous_projection_column.get() != nullptr)
                        cost += previous_projection_column->at(myiterator->get_backward_projection());
                    cost += mycost_computer.get_cost();
                    
                    forward_index = myiterator->get_forward_projection();
                    thcurrent_projection_column->push_back(ff_task_t(forward_index,myiterator->get_index(),cost));
                    //printf("thid=%d [%d, %d,%u,%u]\n",thid, bit_changed, cost, myiterator->get_index(), forward_index);
                }
                node.put(thcurrent_projection_column);
                myiterator->reset(current_indexer);
                
            };
            auto Reduce = [&](std::vector<ff_task_t> *v) {
                //std::cerr << "PF reduce body thread native " << ff_getThreadID() << "\n";
                
                const std::vector<ff_task_t> &V = *v;
                for(size_t j=0;j<V.size();++j) {
                    const unsigned int forward_index = V[j].forward_index;
                    const unsigned int cost          = V[j].cost;
                    if (current_projection_column->at(forward_index) > cost) {
                        current_projection_column->at(forward_index) = cost;
                        backtrace_column->at(forward_index) = V[j].index;
                    }
                }
                delete v;
            };
            
            if (notlast) {

                pf.parallel_reduce_idx(0, iterator->get_length(),1,  adj_chunksize, Map, Reduce, pardegree);
                
                
#if defined(COLUMN_TIME)
                long b = ff::getusec();
                printf("(*)lenght=%ld \ttime=%g\n", iterator->get_length(), (double)(b-a)/1000.0);
#endif
            } else {
                
                std::vector<unsigned int> local_idx_min(pardegree); // numthreads to be uniform
                //std::cerr << "LAST\n";
                // Currently forced to sequential, probably the parallel execution has data races
                pf.parallel_for_idx(0, iterator->get_length(),1,0, [&](const long start, const long stop, const int thid, ff::ff_buffernode &) {
                    if (start == stop) return;
                    
                    ColumnIndexingIterator *myiterator   = thiterators[thid];
                    ColumnCostComputer &mycost_computer  = *(thcostcomputers[thid]);
                    
                    unsigned int localidxmin=0;
                    int bit_changed = -1;
                    if (start>0) {
                        myiterator->advance_idx(&bit_changed,start-1);
                        mycost_computer.set_partitioning(myiterator->get_partition());
                        localidxmin = myiterator->get_index();
                    } else localidxmin = 0;
                    for(long k=start; k<stop; ++k) {
                        
                        myiterator->advance_idx(&bit_changed, k);
                        if (bit_changed >= 0)
                            mycost_computer.update_partitioning(bit_changed);
                        else
                            mycost_computer.set_partitioning(0);
                        
                        
                        unsigned int cost = previous_projection_column->at(myiterator->get_backward_projection());
                        cost += mycost_computer.get_cost();
                        dp_column[myiterator->get_index()] = cost;
                        if (dp_column[localidxmin] > cost) localidxmin = myiterator->get_index();
                        //printf("thid=%d [%d, %d,%u,- (%u --> dp_column[localmin]=%u)]\n",thid, bit_changed, cost, myiterator->get_index(), localidxmin, (localidxmin>=0)?dp_column[localidxmin]:-1);
                    }
                    local_idx_min[thid]=localidxmin;
                }, 1 /* pardegree */);
                
                
                /*
                 printf("pardegree %ld local_idx_min[0]=%ld, local_idx_min[1]=%ld\n", pardegree,local_idx_min[0], local_idx_min[1]);
                 for(size_t m=0;m<dp_column.size();++m) {
                 printf(" [%u]  ", dp_column[m]);
                 }
                 printf("\n");
                 */
                
                running_optimal_score_index =local_idx_min[0];
                for(size_t i=1;i<pardegree;++i) {
                    if (dp_column[local_idx_min[i]] < dp_column[running_optimal_score_index])
                        running_optimal_score_index = local_idx_min[i];
                }
#if defined(COLUMN_TIME)
                printf("(*LAST)lenght=%ld, time=%g\n", iterator->get_length(),(double)(ff::getusec()-a)/1000.0);
#endif
            }
            }
        }
#else // FF_PARALLEL
        {
#if defined(COLUMN_TIME)
            long a = getusec();
#endif
            // do the actual compution on current column
            // To be refactored to shared_ptr
            ColumnCostComputer cost_computer(&(*current_column), all_heterozygous);
            shared_ptr<ColumnIndexingIterator> iterator = current_indexer->get_iterator();
            
            // db
#ifdef DB
            cout << "previous projection column (costs) :" << endl << endl;
            if(previous_projection_column.get()!=0) {
                output_vector_enum(previous_projection_column.get(),current_indexer->get_backward_projection_width());
            }
            cout << endl;
            
            cout << "row ids : ";
            output_vector(current_indexer->get_read_ids());
            
            cout << " .. column size : " << current_indexer->column_size() << endl;
            
            cout << "forward projection mask : ";
            if(next_column.get()!=0) {
                output_vector(current_indexer->get_forward_projection_mask());
                cout << " .. width : " << current_indexer->get_forward_projection_width();
            }
            cout << endl;
            
            cout << "------------------" << endl;
#endif
            
            while (iterator->has_next()) {
                int bit_changed = -1;
                iterator->advance(&bit_changed);
                if (bit_changed >= 0) {
                    cost_computer.update_partitioning(bit_changed);
                } else {
                    cost_computer.set_partitioning(iterator->get_partition());
                    if(next_column.get() == 0) { // only if we're at the last column
                        running_optimal_score_index = iterator->get_index(); // default to first
                    }
                }
                
                unsigned int cost = 0;
                if (previous_projection_column.get() != 0) {
                    cost += previous_projection_column->at(iterator->get_backward_projection());
                }
                
                // db
#ifdef DB
                cout << iterator->get_backward_projection() << " [" << bit_rep(iterator->get_backward_projection(), current_indexer->get_backward_projection_width()) << "] -> " << cost;
#endif
                
                cost += cost_computer.get_cost();
                
                // db
#ifdef DB
                cout << " + " << cost_computer.get_cost() << " = " << cost << " -> " << iterator->get_index() << " [" << bit_rep(iterator->get_index(), current_indexer->get_read_ids()->size()) << "]";
                if(next_column.get()!=0) {
                    cout << " -> " << iterator->get_forward_projection() << " [" << bit_rep(iterator->get_forward_projection(), current_indexer->get_forward_projection_width()) << "]";// fpw = " << current_indexer->get_forward_projection_width();
                }
                cout << endl;
#endif
                
                dp_column[iterator->get_index()] = cost;
                // if not last DP column, then update forward projection column and backtrace column
                if (next_column.get() == 0) {
                    // update running optimal score index
                    if (cost < dp_column[running_optimal_score_index]) {
                        running_optimal_score_index = iterator->get_index();
                    }
                } else {
                    unsigned int forward_index = iterator->get_forward_projection();
                    if (current_projection_column->at(forward_index) > cost) {
                        current_projection_column->at(forward_index) = cost;
                        backtrace_column->at(forward_index) = iterator->get_index();
                    }
                }
            }
#if defined(COLUMN_TIME)
            long b = getusec();
            printf("(*)lenght=%ld \ttime=%g\n", iterator->get_length(), (double)(b-a)/1000.0);
#endif

            // db
#ifdef DB
            cout << endl;
#endif
        }
#endif //FF_PARALLEL
        
        if(next_column.get() == 0) { // record optimal score
            running_optimal_score = dp_column[running_optimal_score_index];
        }
        
        // add newly computed backtrace_table column
        backtrace_table.push_back(backtrace_column);
        
        ++n;
        
        // db
#ifdef DB
        for(size_t j=0;j<current_column->size(); ++j) {
            cout << dp_column.at(j) << endl;
        }
        cout << endl;
#endif
        
        // completion percentage output
#ifdef DB
        if((double(n)/double(nc)) > pc) {
            cout << int(pc*100.0) << " %" << endl;
            pc += pi;
        }
#endif
        
    }
    
    // store optimal score for table at end of computation
    optimal_score = running_optimal_score;
    optimal_score_index = running_optimal_score_index;
    
#ifdef FF_PARALLEL
    pf.threadPause();
    for(size_t i=0;i<numthreads;++i) {
        delete thiterators[i];
        delete thcostcomputers[i];
    }
#endif //FF_PARALLEL
}



unsigned int DPTable::get_optimal_score() {
    //if (backtrace_table.empty()) throw runtime_error("Empty backtrace table");
    return optimal_score;
}

shared_ptr<vector<unsigned int> > DPTable::get_index_path() {
    
    shared_ptr<vector<unsigned int> > index_path = shared_ptr<vector<unsigned int> >(new vector<unsigned int>(indexers.size()));
    if (indexers.size() == 0) {
        return index_path;
    }
    unsigned int index = optimal_score_index;
    index_path->at(indexers.size()-1) = index;
    for(size_t i = indexers.size()-1; i > 0; --i) { // backtrack through table
        if(i>0) {
            shared_ptr<ColumnIndexingIterator> iterator = indexers[i]->get_iterator();
            unsigned int backtrace_index = iterator->index_backward_projection(index);
            index = backtrace_table[i-1]->at(backtrace_index);
            index_path->at(i-1) = index;
        }
    }
    
    //db
#ifdef DB
    cout << "index path : " << endl;
    output_vector(index_path.get());
    cout << endl;
#endif
    
    return index_path;
}

void DPTable::get_super_reads(ReadSet* output_read_set) {
    assert(output_read_set != 0);
    
    ColumnIterator column_iterator(*read_set);
    const vector<unsigned int>* positions = column_iterator.get_positions();
    
    Read* r0 = new Read("superread0", -1, 0);
    Read* r1 = new Read("superread1", -1, 0);
    
    if (backtrace_table.empty()) {
        assert(!column_iterator.has_next());
    } else {
        // run through the file again with the column_iterator
        unsigned int i = 0; // column index
        shared_ptr<vector<unsigned int> > index_path = get_index_path();
        while (column_iterator.has_next()) {
            unsigned int index = index_path->at(i);
            shared_ptr<vector<const Entry *> > column = column_iterator.get_next();
            ColumnCostComputer cost_computer(&(*column), all_heterozygous);
            cost_computer.set_partitioning(index);
            
            r0->addVariant(positions->at(i), cost_computer.get_allele(0), cost_computer.get_weight(0));
            r1->addVariant(positions->at(i), cost_computer.get_allele(1), cost_computer.get_weight(1));
            ++i; // next column
        }
    }
    
    output_read_set->add(r0);
    output_read_set->add(r1);
}

vector<bool>* DPTable::get_optimal_partitioning() {
    
    shared_ptr<vector<unsigned int> > index_path = get_index_path();
    vector<bool>* partitioning = new vector<bool>(read_count,false);
    
    for(size_t i=0; i< index_path->size(); ++i) {
        
#ifdef DB
        cout << "index : " << index_path->at(i) << endl;
#endif
        
        unsigned int mask = 1; // mask to pass over the partitioning (i.e., index)
        for(size_t j=0; j< indexers[i]->get_read_ids()->size(); ++j) {
            
#ifdef DB
            cout << indexers[i]->get_read_ids()->at(j) << " : ";
#endif
            
            unsigned int index = index_path->at(i);
            
#ifdef DB
            cout << index << " & " << mask << " = " << (index & mask);
#endif
            
            if((index & mask) == 0) { // id at this index is in p0 (i.e., in the part.)
                partitioning->at(indexers[i]->get_read_ids()->at(j)) = true;
#ifdef DB
                cout << " : true";
#endif
            }
            
#ifdef DB
            cout << endl;
#endif
            mask = mask << 1;
        }
    }
    
    return partitioning;
}
