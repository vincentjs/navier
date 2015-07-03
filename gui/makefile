CXX := g++
WFLAGS := -Wall -Wextra -Wconversion -pedantic
DFLAGS := -g -O0
CXXFLAGS := -c -cpp $(DFLAGS) $(WFLAGS)
LDLIBS := ./lib/ns_solver.a
LDFLAGS:= -lgfortran -lm
SRCDIR := ./src/
OBJDIR := ./debug/
BINDIR := ./bin/
INCDIR := ./src/include/
SOURCES = $(wildcard $(SRCDIR)*.cpp)
HEADERS = $(wildcard $(INCDIR)*.h)
OBJECTS = $(addprefix $(OBJDIR), $(notdir $(SOURCES:.cpp=.o)))
OUTPUT := $(BINDIR)cfd

all: $(OUTPUT)

$(OUTPUT): $(OBJECTS)
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)

$(OBJECTS): $(OBJDIR)%.o: $(SOURCES) $(HEADERS)
	$(CXX) -J$(OBJDIR) $(CXXFLAGS) $< -o $@

clean:
	rm -f $(OBJDIR)*.o 

dist-clean: clean
	rm -f $(OUTPUT) *~ 

