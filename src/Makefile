# ---------------------------------------------------------------------------
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License version 2 as 
#  published by the Free Software Foundation.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
#
#  As a special exception, you may use this file as part of a free software
#  library without restriction.  Specifically, if other files instantiate
#  templates or use macros or inline functions from this file, or you compile
#  this file and link it with other files to produce an executable, this
#  file does not by itself cause the resulting executable to be covered by
#  the GNU General Public License.  This exception does not however
#  invalidate any other reasons why the executable file might be covered by
#  the GNU General Public License.
#
# ---------------------------------------------------------------------------

CC                   = clang
CXX 		     = clang++ -std=c++11
LINK_OPT             =
VERSION              =
OPTIMIZE_FLAGS       = -DNO_PARFOR_SCHEDULER_THREAD -DCOLUMN_TIME -DHAVE_BUILTIN_CLZ -DHAS_CXX11_VARIADIC_TEMPLATES -O3
#OPTIMIZE_FLAGS       = -DNO_PARFOR_SCHEDULER_THREAD -DHAVE_BUILTIN_CLZ -DHAS_CXX11_VARIADIC_TEMPLATES -O3
CXXFLAGS             = -Wall -Wno-deprecated-declarations
CFLAGS               =
LDFLAGS              = -L/opt/local/lib
INCS                 = -I. -I/opt/local/include -I/Users/aldinuc/SVNs/mc-fastflow-code
#LIBS                 = -rdynamic -lboost_program_options-mt -lboost_iostreams-mt -lboost_unit_test_framework-mt
LIBS                 = -rdynamic -lboost_program_options-mt # -lboost_iostreams -lboost_unit_test_framework
OBJ 		     =  columncostcomputer.o columniterator.o graycodes.o  readset.o columnindexingiterator.o dptable.o \
			indexset.o columnindexingscheme.o entry.o read.o
DEPS                 =

TARGET               = dp


all : $(TARGET)


%.o: %.cpp $(DEPS)
	$(CXX) -c -o $@ $< $(CXXFLAGS) $(OPTIMIZE_FLAGS) $(INCS)

dp: $(OBJ)
	$(CXX) $(OPTIMIZE_FLAGS) -o $@ $^ $(LDFLAGS) $(LIBS)


#clean
clean: 
	-rm -fr *.o *~
cleanall: clean
	-rm -fr $(TARGET) *.d
