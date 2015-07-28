CXX := g++
WFLAGS := -Wall -Wextra -Wconversion -pedantic
DFLAGS := -g -O0
CXXFLAGS := -c -cpp $(DFLAGS) $(WFLAGS)
SRCDIR := ./src/
OBJDIR := ./debug/
BINDIR := ./bin/
LDLIBS := ./lib/cfd.so
LDFLAGS := 
INCDIR := -I./src/include 
SOURCES := $(wildcard $(SRCDIR)*.cpp)
OBJECTS := $(addprefix $(OBJDIR), $(notdir $(SOURCES:.cpp=.o)))
OUTPUT := $(BINDIR)uicfd

all: $(OUTPUT)

$(OUTPUT): $(OBJECTS)
	$(CXX) $(LDFLAGS) $(WFLAGS) -o $@ $^ $(LDLIBS)

$(OBJECTS): $(OBJDIR)%.o: $(SOURCES)
	$(CXX) -J$(OBJDIR) $(CXXFLAGS) $(INCDIR) $< -o $@

clean:
	rm -f $(OBJECTS) 

dist-clean: clean
	rm -f $(OBJECTS) $(OUTPUT) *~
