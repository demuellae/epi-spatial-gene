CC = gcc
CFLAGS = -g
OBJECTS = BaumWelchTree.o forwardTree.o backwardTree.o  initHmm.o nrutil.o readInput.o
LDFLAGS = -lm 
DEPS = hmmTree.h nrutil.h specfunc.h

%.o: %.c
	${CC} ${CFLAGS} -c -o $@ $< $(LDFLAGS)

simulatedTree: $(OBJECTS) specfunc.o simulatedTree.o
	gcc $(CFLAGS) -o $@ $^ $(LDFLAGS)

epiData: $(OBJECTS) epiData.o hmmTree.h specfunc.o
	gcc $(CFLAGS) -o $@ $^ $(LDFLAGS)

specfunc.o: specfunc.h
	gcc -g -I.C:\MinGW\include -c specfunc.c -o specfunc.o -lm

forwardTree.o: hmmTree.h

backwardTree.o: hmmTree.h

BaumWelchTree.o: hmmTree.h

nrutil.o: nrutil.h

clean:
	rm *.o simulatedTree
