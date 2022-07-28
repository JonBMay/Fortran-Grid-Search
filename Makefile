# Makefile for gridsearch

# ----- Setup -----
# Macros
FC = ifort
OPTFLAGS = 
INCL = -I"lib"
OBJS = sims.o subs.o main.o
PROG = GridSearch.e

# Paths for prerequisites
vpath %.f90 src
vpath %.o lib


# ----- Targets -----
gridsearch: $(PROG)

# clean
clean:
	rm -rf lib/* *.o *.mod $(PROG) *.bak src/*.bak
	clear

cleanall:
	rm -rf lib/  *.o *.mod $(PROG) *.bak src/*.bak AllModels.txt Misfits-m.txt
	clear


# ----- Create Program -----
%.o: %.f90
	@echo "compiling $<"
	@$(FC) -diag-disable 10145 $(INCL) $(OPTFLAGS) -c $^ ||:

$(PROG): $(OBJS)
	@[ -d lib ] || mkdir -p lib
	@echo "creating $(PROG):"
	$(FC) $(INCL) $(OPTFLAGS) -o $@ $^
	@mv *.o *.mod lib 2>/dev/null ||:


# ----- Dependency chains -----
main.o : main.f90 subs.o sims.o
subs.o : subs.f90
sims.o : sims.f90
