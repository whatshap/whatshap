
add_executable (whatshap-run
	whatshap-run.cpp
	columnindexingiterator.cpp
	columnindexingiterator.h
	columnindexingscheme.cpp
	columnindexingscheme.h
	columniterator.cpp
	columniterator.h
	entry.cpp
	entry.h
	graycodes.cpp
	graycodes.h
	indexset.cpp
	indexset.h
	pedigreecolumncostcomputer.cpp
	pedigreecolumncostcomputer.h
	pedigree.cpp
	pedigreedptable.cpp
	pedigreedptable.h
	pedigree.h
	pedigreepartitions.cpp
	pedigreepartitions.h
	phredgenotypelikelihoods.cpp
	phredgenotypelikelihoods.h
	read.cpp
	read.h
	readset.cpp
	readset.h
)

# target_link_libraries(whatshap-run ${Boost_LIBRARIES})

# ensure C++11
target_compile_features(whatshap-run PRIVATE cxx_range_for)
