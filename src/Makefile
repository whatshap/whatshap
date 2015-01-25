CPP = g++
CPPFLAGS = -static -Wall -std=c++11 -O2
LDFLAGS = -static -Wall -std=c++11 -O2

EXEC = dp
CPP_SRC = \
\
	columncostcomputer.cpp \
	columnindexingiterator.cpp \
	columnindexingscheme.cpp \
	columniterator.cpp \
	dptable.cpp \
	entry.cpp \
	graycodes.cpp \
	read.cpp \
	readset.cpp

CPP_OBJ = $(CPP_SRC:.cpp=.o)

all : ${EXEC}

${EXEC} : ${CPP_OBJ}
	${CPP} -o ${@} ${^} dp.cpp ${LDFLAGS}

%.o : %.cpp
	${CPP} -o ${@} -c ${<} ${CPPFLAGS}

clean :
	rm -f *.o ${EXEC}
