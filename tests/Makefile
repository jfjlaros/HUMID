EXEC := run_tests
MAIN := test_lib
TESTS := test_cluster
FIXTURES :=


CC := g++
CC_ARGS := -std=c++20 -fcoroutines -Wno-pmf-conversions


OBJS := $(addsuffix .o, $(TESTS) $(FIXTURES))

.PHONY: all check clean distclean


all: $(EXEC)

$(EXEC): $(MAIN).cc $(OBJS)
	$(CC) $(CC_ARGS) -o $@ $^

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
