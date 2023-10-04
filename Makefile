CXX=g++
CC=g++

#INCLUDE_PATH=/misc/rci-rg/ag_kehr/user/mit16436/anaconda3/envs/genotyping/include
#LIB_PATH=/misc/rci-rg/ag_kehr/user/mit16436/anaconda3/envs/genotyping/lib
JSON_PATH=./include

BINARY := genotype

override CXXFLAGS += -I $(JSON_PATH) -DSEQAN_HAS_ZLIB -lz -lhts -DSEQAN_DISABLE_VERSION_CHECK -std=c++17 -Wall -O2 -DEIGEN_DONT_PARALLELIZE -fopenmp -lpthread -lboost_filesystem

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
	@$(CXX) $(CXXFLAGS) $(OBJECTS) -o $@

$(ODIR)/%.o: $(SRCDIR)/%.cpp
	@mkdir -p $(@D)
	@$(CXX) -c $(CXXFLAGS) $< -o $@

clean:
	if [ -f ./genotype ]; then rm genotype; fi
	rm -r $(ODIR)

install:
	if [ -f ./genotype ]; then cp genotype /usr/bin; fi

uninstall:
	if [ -f /usr/bin/genotype ]; then rm /usr/bin/genotype; fi
