CC = gcc
CFLAGS = -Wall -Wpointer-arith -Wextra -Og -I.
LDLIBS = -lgmp

.SUFFIXES: .c .o

PRGS = test

all:\
	$(PRGS)

test:\
	test.o pollard.o quadratic_sieve.o
	$(CC) $(LDLIBS) test.o pollard.o quadratic_sieve.o -o $@

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
