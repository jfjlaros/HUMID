EXEC := dedup
MAIN := dedup.cc
LIBS := ../lib/ngs/ngs

CC := g++
CC_ARGS := -O2


OBJS := $(addsuffix .o, $(LIBS))

.PHONY: clean distclean


all: $(EXEC)

$(EXEC): $(MAIN) $(OBJS)
	$(CC) $(CC_ARGS) -o $@ $^

%.o: %.c
	$(CC) $(CC_ARGS) -o $@ -c $<

%.o: %.cc
	$(CC) $(CC_ARGS) -o $@ -c $<

clean:
	rm -f $(OBJS)

distclean: clean
	rm -f $(EXEC)
