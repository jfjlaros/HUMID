EXEC := dedup
MAIN := dedup.cc
LIBS := ngs ../lib/commandIO/src/error \
  ../lib/commandIO/src/plugins/cli/io ../lib/commandIO/src/plugins/repl/io \
  ../lib/fastp/src/fastqreader ../lib/fastp/src/readpool \
  ../lib/fastp/src/read ../lib/fastp/src/sequence ../lib/fastp/src/writer \
  ../lib/fastp/src/fastareader ../lib/fastp/src/options


CC := g++
CC_ARGS := -O2 -fcoroutines -std=c++20
LD_ARGS := -lisal -ldeflate
LD_STATIC_ARGS := -L../lib/isa-l/.libs/ -static


OBJS := $(addsuffix .o, $(LIBS))

.PHONY: clean distclean


all: $(EXEC)

static: $(MAIN) $(OBJS)
	$(CC) $(CC_ARGS) -o $@ $^ $(LD_ARGS) $(LD_STATIC_ARGS)

$(EXEC): $(MAIN) $(OBJS)
	$(CC) $(CC_ARGS) -o $@ $^ $(LD_ARGS)

%.o: %.cc
	$(CC) $(CC_ARGS) -o $@ -c $<

%.o: %.cpp
	$(CC) $(CC_ARGS) -o $@ -c $<

clean:
	rm -f $(OBJS)

distclean: clean
	rm -f $(EXEC) static
