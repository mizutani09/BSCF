#compiler
FC=$(HOME)/mesasdk/bin/gfortran
FFLAGS=-O3 -fdefault-real-8

F90FILES=bessel.F90 elle.F90 ellf.F90 gammln.F90\
helmadi.F90 main.F90 potential_solver.F90 potsetup.F90 rd.F90 rf.F90 realft.F90 \
setup.F90 sm.F90 tm.F90 tridagr.F90 tridagz.F90 guessrho.F90 poisson_solve.F90\
print2d.F90 print1d.F90 print2default.F90 print1default.F90 testrho.F90 \
findmass.F90 findmom.F90 getinfo.F90 virial.F90 findvol.F90 findj.F90


EXE_DIR := bin/
EXECUTABLE := $(EXE_DIR)bscf
OBJ_DIR := obj/
OBJ_FILES = $(addprefix $(OBJ_DIR),$(notdir $(F90FILES:.F90=.o)))

.PHONY: all dirs clean

all : dirs $(EXECUTABLE)

objs : dirs $(OBJ_FILES)

dirs : $(EXE_DIR) $(OBJ_DIR)

# make folders
$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)
$(EXE_DIR):
	mkdir -p $(EXE_DIR)

# make executable
$(EXECUTABLE): $(OBJ_FILES)
	$(FC) $(OBJ_FILES) -o $(EXECUTABLE)

# make object files
$(OBJ_DIR)%.o : %.F90 runhydro.h
	$(FC) $(FFLAGS) -c $< -o $@

clean:
	rm -rf $(OBJ_DIR)*
	rm -rf $(EXECUTABLE)
