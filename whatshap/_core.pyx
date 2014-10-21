# distutils: language = c++
# distutils: sources = src/activelist.cpp src/activelistdelegator.cpp src/columncostcomputer.cpp src/columnindexingiterator.cpp src/columnindexingscheme.cpp src/columnreader.cpp src/dp.cpp src/dptable.cpp src/entry.cpp src/graycodes.cpp

#print 'hello world'
cdef extern from "../src/dptable.h":
	cdef cppclass DPTable:
		DPTable(bool) except +
