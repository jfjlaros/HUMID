EXEC := run_tests
MAIN := test_lib
TESTS := test_cluster test_fastq
LIBS := ../src/cluster ../src/fastq \
  ../lib/fastp/src/fastqreader ../lib/fastp/src/readpool \
  ../lib/fastp/src/read ../lib/fastp/src/sequence ../lib/fastp/src/writer \
  ../lib/fastp/src/fastareader ../lib/fastp/src/options
FIXTURES :=


CC := g++
CC_ARGS := -std=c++20 -fcoroutines
LD_ARGS := -lisal -ldeflate


OBJS := $(addsuffix .o, $(TESTS) $(FIXTURES) $(LIBS))

.PHONY: all check clean distclean


all: $(EXEC)

$(EXEC): $(MAIN).cc $(OBJS)
	$(CC) $(CC_ARGS) -o $@ $^ $(LD_ARGS)

%.o: %.cpp
	$(CC) $(CC_ARGS) -o $@ -c $<

%.o: %.cc
	$(CC) $(CC_ARGS) -o $@ -c $<

check: all
	valgrind --error-exitcode=1 --leak-check=full ./$(EXEC) || ((echo "\033[0;31mVALGRIND ERRORS DETECTED\033[0m" && false))

clean:
	rm -f $(OBJS)

distclean: clean
	rm -f $(EXEC)
