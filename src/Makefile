EXEC := dedup
MAIN := dedup.cc
LIBS := ../lib/ngs/ngs ../lib/commandIO/src/error \
  ../lib/commandIO/src/plugins/cli/io ../lib/commandIO/src/plugins/repl/io \
  ../lib/fastp/src/fastqreader ../lib/fastp/src/readpool \
  ../lib/fastp/src/read ../lib/fastp/src/sequence


CC := g++
CC_ARGS := -O2 -fcoroutines -std=c++20
LD_ARGS := -lisal


OBJS := $(addsuffix .o, $(LIBS))

.PHONY: clean distclean


all: $(EXEC)

$(EXEC): $(MAIN) $(OBJS)
	$(CC) $(CC_ARGS) -o $@ $^ $(LD_ARGS)

%.o: %.c
	$(CC) $(CC_ARGS) -o $@ -c $<

%.o: %.cc
	$(CC) $(CC_ARGS) -o $@ -c $<

clean:
	rm -f $(OBJS)

distclean: clean
	rm -f $(EXEC)
