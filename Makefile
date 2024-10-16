CXX=g++
CC=g++

#INCLUDE_PATH=
#LIB_PATH=
JSON_PATH=./include

BINARY := ggtyper

DATE := on $(shell git log --pretty=format:"%cd" --date=iso | cut -f 1,2 -d " " | head -n 1)
VERSION := 0.0.1-$(shell git log --pretty=format:"%h" --date=iso | head -n 1)

override CXXFLAGS += -DDATE=\""$(DATE)"\" -DVERSION=\""$(VERSION)"\"
override CXXFLAGS += -I $(JSON_PATH) -DSEQAN_HAS_ZLIB -lz -lhts -DSEQAN_DISABLE_VERSION_CHECK -DEIGEN_DONT_PARALLELIZE -DEIGEN_DEFAULT_DENSE_INDEX_TYPE=int64_t -std=c++17 -Wall -O2 -fopenmp -lpthread -g 

ifdef INCLUDE_PATH
	override CXXFLAGS+=-I $(INCLUDE_PATH)
endif

ifdef LIB_PATH
	override CXXFLAGS+=-L $(LIB_PATH)
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
	if [ -f $(BINARY) ]; then rm $(BINARY); fi
	if [ -d $(ODIR) ]; then rm -r $(ODIR); fi

install:
	if [ -f $(BINARY) ]; then cp $(BINARY) /usr/bin; fi

uninstall:
	if [ -f /usr/bin/$(BINARY) ]; then rm /usr/bin/$(BINARY); fi
