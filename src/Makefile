EXEC := humid
MAIN := humid.cc
LIBS := cluster fastq log ../lib/commandIO/src/error \
  ../lib/commandIO/src/plugins/cli/io ../lib/commandIO/src/plugins/repl/io \
  ../lib/fastp/src/fastqreader ../lib/fastp/src/readpool \
  ../lib/fastp/src/read ../lib/fastp/src/sequence ../lib/fastp/src/writer \
  ../lib/fastp/src/fastareader ../lib/fastp/src/options


CXX ?= g++
CC_ARGS := -O2 -fcoroutines -std=c++20
LD_ARGS := -lisal -ldeflate
LD_STATIC_ARGS := -L../lib/isa-l/.libs/ -static


OBJS := $(addsuffix .o, $(LIBS))

.PHONY: clean distclean doc


all: $(EXEC)

static: $(MAIN) $(OBJS)
	$(CXX) $(CC_ARGS) -o $@ $^ $(LD_ARGS) $(LD_STATIC_ARGS)

$(EXEC): $(MAIN) $(OBJS)
	$(CXX) $(CC_ARGS) -o $@ $^ $(LD_ARGS)

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

doc: ../docs/cli.rst
