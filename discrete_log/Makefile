CC = gcc
CFLAGS = -Wall -Wpointer-arith -Wextra -Og -I. -I../hashmap/include -I../blake2b
LDLIBS = -lgmp

.SUFFIXES: .c .o

PRGS = test

all:\
	$(PRGS)

test:\
	test.o discrete_log.o hashmap.o blake2b.o
	$(CC) $(LDLIBS) test.o discrete_log.o hashmap.o blake2b.o -o $@

hashmap.o:
	$(CC) $(CFLAGS) -c ../hashmap/src/hashmap.c

blake2b.o:
	$(CC) $(CFLAGS) -c ../blake2b/blake2b.c

.c.o: \
	$*.h $*.c
	@echo $*
	$(CC) $(CFLAGS) -c $<

clean:
	find . -name '*.o' -exec rm -f {} ';'
	find . -name '*~' -exec rm -f {} ';'
	find . -name '#*#' -exec rm -f {} ';'
	find . -name '.#*' -exec rm -f {} ';'
	find . -name 'core' -exec rm -f {} ';'
	find . -name '*.core' -exec rm -f {} ';'
	-rm -rf $(PRGS)
