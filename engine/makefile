FC := gfortran
AR := ar cr
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
GUIDIR := ../gui/lib/

all: $(OUTPUT)
	cp $(OUTPUT) $(GUIDIR)

$(OUTPUT): $(OBJECTS)
	$(AR) $(LDFLAGS) -o $@ $^ $(LDLIBS)
#	$(FC) $(LDFLAGS) -o $@ $^ $(LDLIBS)

$(OBJECTS): $(OBJDIR)%.o: $(SRCDIR)%.f90
	$(FC) -J$(OBJDIR) $(FFLAGS) $< -o $@

$(OBJDIR)main.o: $(OBJDIR)precision.o $(OBJDIR)solver.o

clean:
	rm -f $(OBJDIR)*.mod $(OBJDIR)*.o

dist-clean: clean
	rm -f $(OUTPUT)

install: $(OUTPUT) clean
