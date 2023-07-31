CXX=g++
CC=g++

INCLUDE_PATH=/home/tim/.conda/envs/genotyping/include
LIB_PATH=/home/tim/.conda/envs/genotyping/lib
JSON_PATH=./include

BINARY := genotype

CXXFLAGS=-I $(JSON_PATH) -DSEQAN_HAS_ZLIB -lz -lhts -DSEQAN_DISABLE_VERSION_CHECK -std=c++14 -Wall -O2 -fopenmp -lpthread -lboost_filesystem

ifdef INCLUDE_PATH
	CXXFLAGS+=-I $(INCLUDE_PATH)
endif

ifdef LIB_PATH
	CXXFLAGS+=-L $(LIB_PATH)
endif

ODIR:=build
SRCDIR:=src

SOURCES := $(shell find $(SRCDIR) -name '*.cpp')
OBJECTS := $(addprefix $(ODIR)/, $(SOURCES:$(SRCDIR)/%.cpp=%.o))

all: $(BINARY)

$(BINARY): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(OBJECTS) -o $@

$(ODIR)/%.o: $(SRCDIR)/%.cpp
	$(CXX) -c $(CXXFLAGS) $< -o $@

clean:
	if [ -f ./genotype ]; then rm genotype; fi
	rm build/*.o

install:
	if [ -f ./genotype ]; then cp genotype /usr/bin; fi

uninstall:
	if [ -f /usr/bin/genotype ]; then rm /usr/bin/genotype; fi
