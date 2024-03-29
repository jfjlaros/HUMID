EXEC := humid
MAIN := humid.cc
LIBS := cluster fastq log ../lib/commandIO/src/error \
  ../lib/commandIO/src/plugins/cli/io ../lib/commandIO/src/plugins/repl/io \
  ../lib/fastp/src/fastqreader ../lib/fastp/src/readpool \
  ../lib/fastp/src/read ../lib/fastp/src/sequence ../lib/fastp/src/writer \
  ../lib/fastp/src/fastareader ../lib/fastp/src/options

PREFIX ?= /usr/local
BINDIR ?= $(PREFIX)/bin

CXX ?= g++
CC_ARGS := -O2 -std=c++20

# Compatibility for MacOS
UNAME_S := $(shell uname -s)

ifeq ($(UNAME_S),Linux)
	CC_ARGS += -fcoroutines
endif

ifeq ($(UNAME_S),Darwin)
	CC_ARGS += -fcoroutines-ts -D_LIBCPP_DISABLE_AVAILABILITY
endif

LDFLAGS += -lisal -ldeflate
LD_STATIC_ARGS := -L../lib/isa-l/.libs/ -static


OBJS := $(addsuffix .o, $(LIBS))

.PHONY: clean distclean doc


all: $(EXEC)

static: $(MAIN) $(OBJS)
	$(CXX) $(CC_ARGS) -o $@ $^ $(LDFLAGS) $(LD_STATIC_ARGS)

$(EXEC): $(MAIN) $(OBJS)
	$(CXX) $(CC_ARGS) -o $@ $^ $(LDFLAGS)

%.o: %.cc
	$(CXX) $(CC_ARGS) -o $@ -c $<

%.o: %.cpp
	$(CXX) $(CC_ARGS) -o $@ -c $<

../docs/cli.rst: ../docs/cli_template.rst $(EXEC)
	cp $< $@ && ./$(EXEC) -h >> $@

clean:
	rm -f $(OBJS)

distclean: clean
	rm -f $(EXEC) static

install:
	install $(EXEC) $(BINDIR)/$(EXEC)

doc: ../docs/cli.rst
