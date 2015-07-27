FC := gfortran
AR := ar cr
FPY := f2py 
WFLAGS := -Wall -Wextra -Wconversion -pedantic
DFLAGS := -g -O0 -fbacktrace -fbounds-check -ffpe-trap=zero,overflow,underflow
FFLAGS := -c -cpp $(DFLAGS) $(WFLAGS)
LDLIBS :=
LDFLAGS:= 
SRCDIR := ./src/
OBJDIR := ./include/
BINDIR := ./bin/
SOURCES = $(wildcard $(SRCDIR)*.f90)
OBJECTS = $(addprefix $(OBJDIR), $(notdir $(SOURCES:.f90=.o)))
OUTPUT := $(BINDIR)ns_solver.a
CPPDIR := ../cpp/lib

all: $(OUTPUT)
	cp $(OUTPUT) $(CPPDIR)

$(OUTPUT): $(OBJECTS)
#	$(FPY) -c $^ #-m $(OUTPUT)
	$(AR) $(LDFLAGS) -o $@ $^ $(LDLIBS)
#	$(FC) $(LDFLAGS) -o $@ $^ $(LDLIBS)

$(OBJECTS): $(OBJDIR)%.o: $(SRCDIR)%.f90
	$(FC) -fPIC -J$(OBJDIR) $(FFLAGS) $< -o $@

$(OBJDIR)main.o: $(OBJDIR)precision.o $(OBJDIR)solver.o

clean:
	rm -f $(OBJDIR)*.mod $(OBJECTS)

dist-clean: clean
	rm -f $(OUTPUT)

scratch: clean $(all)

install: $(OUTPUT) clean
