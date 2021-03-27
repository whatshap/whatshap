SRC_DIR    = src
BIN_DIR    = ../../bin

CC = g++
CPPFLAGS = -Wall -fopenmp -std=c++11 -g -O2
LIBS = -l boost_timer -l boost_system -l boost_program_options

EXEC_NAMES = ViralQuasispecies

EXECS = $(BIN_DIR)/ViralQuasispecies

VPATH = src

.PHONY : all

all: $(EXECS) 

.PHONY : clean

clean:
	rm -rf $(EXECS) $(BIN_DIR)

$(BIN_DIR)/ViralQuasispecies: $(SRC_DIR)/*.cpp $(SRC_DIR)/*.h ${BIN_DIR}
	$(CC) $(CPPFLAGS) $(SRC_DIR)/*.cpp -o $@ ${LIBS}

$(BIN_DIR):
	mkdir -p $(BIN_DIR)

