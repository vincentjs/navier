CXX := gcc
WFLAGS := -Wall -Wextra -Wconversion -pedantic
DFLAGS := -g -O0
CXXFLAGS := -c -cpp $(DFLAGS) $(WFLAGS)
SRCDIR := ./src/
OBJDIR := ./debug/
BINDIR := ./bin/
LDLIBS := ./lib/ns_solver.a
LDFLAGS:= -lgfortran -lm 
INCDIR := -I./src/include -I/usr/include/boost 
SOURCES = $(wildcard $(SRCDIR)*.cpp)
OBJECTS = $(addprefix $(OBJDIR), $(notdir $(SOURCES:.cpp=.o)))
OUTPUT := $(BINDIR)cfd.so
UIDIR := ../ui/lib
GUIDIR := ../gui/lib

all: $(OUTPUT)
	cp $(OUTPUT) $(UIDIR)
	cp $(OUTPUT) $(GUIDIR)

$(OUTPUT): $(OBJECTS)
	$(CXX) -fPIC $(LDFLAGS) $(WFLAGS) -o $@ $^ $(LDLIBS)

$(OBJECTS): $(OBJDIR)%.o: $(SOURCES)
	$(CXX) -J$(OBJDIR) $(CXXFLAGS) $(INCDIR) $< -o $@

clean:
	rm -f $(OBJECTS) $(OUTPUT) 

dist-clean: clean
	rm -f $(OBJECTS) $(OUTPUT) *~ 

