LIBS := -lamrex
BIN := .
OBJ := .
SRC := .
SRCS := $(wildcard $(SRC)/*.cpp)
OBJS := $(patsubst $(SRC)/%.cpp,$(OBJ)/%.o,$(SRCS))
EXE  := $(BIN)/heat
CFLAGS := -Wall -O2
VFLAGS := --leak-check=full --show-reachable=yes --track-origins=yes

.PHONY: all run clean

all: $(EXE)

$(EXE): $(OBJS) | $(BIN)
	$(CXX) $^ -o $@ $(LIBS) 

$(OBJ)/%.o: $(SRC)/%.cpp | $(OBJ)
	$(CXX) $(CFLAGS) -c $< -o $@

run: $(EXE)
	mpirun -n 2 ./$< input

mem-debug:
	valgrind $(VFLAGS) ./$(EXE) input

pmem-debug:
	mpirun -n 2 valgrind $(VFLAGS) ./$(EXE) input

clean:
	rm -f $(OBJ)/*.o
	rm -f $(EXE)
	rm -rf plt*
	rm -f Back*
