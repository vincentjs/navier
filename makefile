FC := gfortran
AR := ar cr
CFLAGS := -c -cpp
WFLAGS := -Wall -Wextra -Wconversion -pedantic
DFLAGS := -g -fbacktrace -fbounds-check -ffpe-trap=zero,overflow,underflow
LDLIBS := 
SRCDIR := ./src/
OBJDIR := ./include/
BINDIR := ./bin/
SOURCES = $(wildcard $(SRCDIR)*.f90)
OBJECTS = $(addprefix $(OBJDIR), $(notdir $(SOURCES:.f90=.o)))
OUTPUT := $(BINDIR)ns_solver.a

all: $(OUTPUT)

$(OUTPUT): $(OBJECTS) 
	$(AR) -o $@ $^ $(LDFLAGS)

$(OBJECTS): $(OBJDIR)%.o: $(SRCDIR)%.f90
	$(FC) -J$(OBJDIR) $(CFLAGS) $(WFLAGS) $(DFLAGS) $< -o $@

#$(OBJDIR)kdtree.o: $(OBJDIR)assert.o
#$(OBJDIR)assert.o: $(OBJDIR)error.o
#$(OBJDIR)constants.o: $(OBJDIR)precision.o


clean:
	@rm -f ./include/*.mod ./include/*.o

dist-clean: clean
	@rm -f ./lib/*.a

install: $(OUTPUT) clean
