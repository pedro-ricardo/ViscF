# --------------------------------------------------
# Makefile for AMRex Fortran Test examples
# --------------------------------------------------

# Compiled files Location
BIN := .
OBJ := .
EXE  := $(BIN)/solve
# Lib SRC location
MFSim_SRC := ../../src
MFSim_SRCS := $(wildcard $(MFSim_SRC)/*.f08)
MFSim_OBJS := $(patsubst $(MFSim_SRC)/%.f08,$(OBJ)/%.o,$(MFSim_SRCS))
# Example SRC location
SRC := .
SRCS := $(wildcard $(SRC)/*.f08)
OBJS := $(patsubst $(SRC)/%.f08,$(OBJ)/%.o,$(SRCS))
# Join compilation list
OBJS += $(MFSim_OBJS)
# Compilation
INC := -I/opt/apps/jsonfortran/8.2.0-gcc-9.1.1/lib -I$(AMREX_INC)
LIBS := -lamrex -lstdc++ -ljsonfortran
FDEBUG := -fcheck=all 
FCLEAN := -Wall -Wextra -Werror -Wno-unused-dummy-argument
FFLAGS := -g -cpp -O1 -ffree-line-length-none $(FDEBUG) $(FCLEAN)

VFLAGS := --leak-check=full --show-reachable=yes --track-origins=yes

# Unknown
.PHONY: all run clean

# Check for Fortran Modules in file dependencies
%.mod: %.o
	@true

# Build all
all: $(EXE)
# Executable
$(EXE): $(OBJS) | $(BIN)
	$(FC) $^ -o $@ $(LIBS) 
# Files in Lib location
$(OBJ)/%.o: $(SRC)/%.f08 | $(OBJ)
	$(FC) $(FFLAGS) $(INC) -c $< -o $@
# Files in Exemple location
$(OBJ)/%.o: $(MFSim_SRC)/%.f08 | $(OBJ)
	$(FC) $(FFLAGS) $(INC) -c $< -o $@

# Run
run: $(EXE)
	mpirun -n 2 ./$<

# Memory leak tool
mem-debug:
	valgrind $(VFLAGS) ./$(EXE)
pmem-debug:
	mpirun -n 2 valgrind $(VFLAGS) ./$(EXE)

# Clean
clean:
	rm -f $(OBJ)/*.o
	rm -f $(OBJ)/*.mod
	rm -f $(EXE)
	rm -rf plt*
	rm -f time_plt.visit
	

# Example dependency list
main.o: $(SRC)/main.f08 input_reader.mod data_input.mod amr_variables.mod amr_remesh.mod amr_output.mod example_functions.mod
example_functions.o: $(SRC)/example_functions.f08 data_variable.mod data_boxes.mod amr_coords.mod data_input.mod
example_solver.o:  $(SRC)/example_solver.f08 data_variable.mod data_boxes.mod data_mesh.mod

# Lib dependency list
input_reader.o: $(MFSim_SRC)/input_reader.f08 data_input.mod data_mesh.mod
data_input.o: $(MFSim_SRC)/data_input.f08 data_mesh.mod data_remesh.mod data_simulation.mod data_io.mod
data_mesh.o: $(MFSim_SRC)/data_mesh.f08 data_boxes.mod data_variable.mod amr_ghosts.mod
data_remesh.o: $(MFSim_SRC)/data_remesh.f08
data_simulation.o: $(MFSim_SRC)/data_simulation.f08
data_io.o: $(MFSim_SRC)/data_io.f08
data_variable.o: $(MFSim_SRC)/data_variable.f08
data_boxes.o: $(MFSim_SRC)/data_boxes.f08 data_variable.mod
amr_variables.o: $(MFSim_SRC)/amr_variables.f08 data_variable.mod data_mesh.mod amr_startvariables.mod output_list.mod
amr_remesh.o: $(MFSim_SRC)/amr_remesh.f08 amr_variables.mod amr_criteria.mod data_variable.mod data_mesh.mod
amr_startvariables.o: $(MFSim_SRC)/amr_startvariables.f08 data_variable.mod data_boxes.mod amr_coords.mod data_input.mod
amr_criteria.o: $(MFSim_SRC)/amr_criteria.f08 data_input.mod data_variable.mod data_boxes.mod amr_coords.mod
amr_coords.o: $(MFSim_SRC)/amr_coords.f08
output_list.o: $(MFSim_SRC)/output_list.f08 output_names.mod
output_names.o: $(MFSim_SRC)/output_names.f08 data_variable.mod
amr_output.o: $(MFSim_SRC)/amr_output.f08 data_input.mod output_list.mod data_variable.mod data_boxes.mod amr_variables.mod
amr_ghosts.o: $(MFSim_SRC)/amr_ghosts.f08 data_variable.mod
