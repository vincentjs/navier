CXX := g++
WFLAGS := -Wall -Wextra -Wconversion -pedantic
DFLAGS := -g -O0
CXXFLAGS := -c -cpp $(DFLAGS) $(WFLAGS)
SRCDIR := ./src/
HEADIR := ./src/include/
OBJDIR := ./debug/
BINDIR := ./bin/
LDLIBS := ./lib/ns_solver.a
LDFLAGS:= -lgfortran -lm 
INCDIR := -I$(HEADIR)  
SOURCES = $(wildcard $(SRCDIR)*.cpp)
OBJECTS = $(addprefix $(OBJDIR), $(notdir $(SOURCES:.cpp=.o)))
OUTPUT := $(BINDIR)cfd.so
UIDIR := ../ui/lib
UIPYDIR := ../uipy
GUIDIR := ../gui/lib

all: $(OUTPUT)
	cp $(OUTPUT) $(UIDIR)
	cp $(OUTPUT) $(UIPYDIR)
	cp $(OUTPUT) $(GUIDIR)
	cp $(HEADIR)*.h ../ui/src/include

$(OUTPUT): $(OBJECTS)
	$(CXX) $(LDFLAGS) $(WFLAGS) -shared -o $@ $^ $(LDLIBS)

$(OBJDIR)%.o: $(SRCDIR)%.cpp
	$(CXX) -fPIC -J$(OBJDIR) $(CXXFLAGS) $(INCDIR) $(WFLAGS) $< -o $@

clean:
	rm -f $(OBJECTS) $(OUTPUT) 

dist-clean: clean
	rm -f $(OBJECTS) $(OUTPUT) *~ 

