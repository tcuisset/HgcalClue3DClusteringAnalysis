CXX     = g++
OPT_FLAGS = -O3

CXXFLAGS= $(OPT_FLAGS) $(shell root-config --cflags)
LIBS    = $(shell root-config --libs) 

SOURCES = TBNtupleAnalyzer.cc Runclustering.cc 
HEADERS = TBNtupleAnalyzer.h Runclustering.h 
OBJECTS = $(SOURCES:.cc=.o)

EXECUTABLE = runclustering

all: $(SOURCES) $(EXECUTABLE)

%.o: %.cc
	@echo Compiling $<...
	$(CXX) $(CXXFLAGS) -c -o $@ $< 


$(EXECUTABLE): $(OBJECTS)
	@echo "Linking $(PROGRAM) ..."
	@echo "@$(CXX) $(OPT_FLAGS) $(LIBS) $(OBJECTS) -g -o $@"
	@$(CXX) $(OPT_FLAGS) -o $@ $^ $(LIBS)
	@echo "done"


# Specifying the object files as intermediates deletes them automatically after the build process.
.INTERMEDIATE:  $(OBJECTS)

# The default target, which gives instructions, can be called regardless of whether or not files need to be updated.
.PHONY : clean
clean:
	rm -f $(OBJECTS) $(EXECUTABLE)

###
TBNtupleAnalyzer.o: TBNtupleAnalyzer.h
Runclustering.o:TBNtupleAnalyzer.h Runclustering.h
