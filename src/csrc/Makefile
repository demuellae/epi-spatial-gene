#
#
# Make file for compiling HMM code in this directory.
# Author: Tapas Kanungo
# Date: 23 February 1998
# $Id: Makefile,v 1.3 1998/02/23 08:12:35 kanungo Exp kanungo $
# 
#
CFLAGS= -g
INCS=
# use the following line to "Purify" the code
#CC=purify gcc
CC=gcc
DEPS = hmmTree.h specfunc.h
OBJ = forwardTree.o backwardTree.o BaumWelchTree.o beta_function.o digamma_function.o gamma_function.o ln_gamma_function.o initHmm.o nrutil.o asa121.o readInput.o

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

driver : $(OBJ) driver.o
	gcc -o $@ $^ $(CFLAGS)
	
simulatedTree: $(OBJ) simulatedTree.o
	gcc -o $@ $^ $(CFLAGS)


clean:
	rm *.o driver

