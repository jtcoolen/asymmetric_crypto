CC = gcc
CFLAGS = -Wall -Wpointer-arith -Wextra -Wpedantic -O3 -I. -I../finite_field

.SUFFIXES: .c .o

PRGS = test

all:\
	$(PRGS)

test:\
	test.o rsa.o helpers.o finite_field.o elgamal.o blake2b.o
	$(CC) test.o rsa.o helpers.o finite_field.o elgamal.o blake2b.o -o $@

finite_field.o:
	$(CC) $(CFLAGS) -c ../finite_field/finite_field.c

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
