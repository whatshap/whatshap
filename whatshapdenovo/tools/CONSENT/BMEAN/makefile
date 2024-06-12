# CC=/usr/bin/g++
CXX= g++
CFLAGS = -O3 -std=c++11 -lpthread -Ispoa/include
LFLAGS = -Ispoa/include/ -Lspoa/build/lib/ -lspoa
EXEC=testLR
OBJS := *.o Complete-Striped-Smith-Waterman-Library/src/*.o

all: $(EXEC)

ifeq ($(prof),1)
 CFLAGS+= -pg
endif
ifeq ($(deb),1)
 CFLAGS+= -O0 -DASSERTS -g
endif

ifeq ($(sani),1)
 CFLAGS= -std=c++11 -lpthread -fsanitize=address -fno-omit-frame-pointer -O1 -fno-optimize-sibling-calls -g
endif



test:
	./testLR

all: $(EXEC)



testLR:  testLR.o bmean.o utils.o
	$(CXX) $(CFLAGS)  $(OBJS) -o $@ $(LFLAGS)

%.o: %.cpp
	$(CXX) -o $@ -c $< $(CFLAGS) $(LFLAGS)


clean:
	rm -rf *.o
	rm testLR

rebuild: clean $(EXEC)
