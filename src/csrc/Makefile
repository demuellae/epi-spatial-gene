CC = gcc
CFLAGS = -g
OBJECTS = asa121.o BaumWelchTree.o beta_function.o digamma_function.o forwardTree.o backwardTree.o beta_function.o digamma_function.o gamma_function.o initHmm.o ln_gamma_function.o nrutil.o readInput.o simulatedTree.o
LDFLAGS = -lm
DEPS = hmmTree.h nrutil.h specfunc.h

%.o: %.c $(DEPS) 
	${CC} ${CFLAGS} -c -o $@ $<
	
simulatedTree: $(OBJECTS)
	gcc $(CFLAGS) -o $@ $^ $(LDFLAGS)

asa121.o: specfunc.h

beta_function.o: specfunc.h

digamma_fucntion.o: specfunc.h

forwardTree.o: hmmTree.h

backwardTree.o: hmmTree.h

BaumWelchTree.o: hmmTree.h

nrutil.o: nrutil.h
	
