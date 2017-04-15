CC=gcc
CXX=g++
override CFLAGS += -Wall -g -std=c99
override CXXFLAGS += -Wall -g -std=c99
override LDFLAGS += -lm

BIN_FILES=lagrange_test.bin lagrange_convergence.bin
OBJECT_FILES=lagrange.o lagrange_test.o lagrange_convergence.o
OUT_FILES = test_full.out\
	    test_fo.out\
	    convergence_slice.out\
	    convergence_norm.out
ORDER_1D=4
ORDER_X=2
ORDER_Y=2
NX=101
NY=101

default: all
all: lagrange_test.bin lagrange_convergence.bin

test: test_fo.out convergence_slice.out 
plot: test_fo.out convergence_slice.out
	python plot_error.py test_fo.out ${NX} ${NY}
	python plot_convergence_1d.py ${ORDER_1D} convergence_slice.out
	python plot_convergence_2d.py convergence_norm.out
convergence_slice.out: lagrange_convergence.bin
	./lagrange_convergence.bin ${ORDER_1D}
test_fo.out: lagrange_test.bin
	./lagrange_test.bin ${ORDER_X} ${ORDER_Y}

lagrange_test.bin: lagrange_test.o lagrange.o
	$(CC) ${CFLAGS} -o $@ $^ ${LDFLAGS}

lagrange_convergence.bin: lagrange_convergence.o lagrange.o
	$(CC) ${CFLAGS} -o $@ $^ ${LDFLAGS}

lagrange_convergence.o: lagrange_convergence.c lagrange.h
lagrange_test.o: lagrange_test.c lagrange.h
lagrange.o: lagrange.c lagrange.h

clean:
	$(RM) ${OBJECT_FILES}
	$(RM) ${BIN_FILES}
	$(RM) ${OUT_FILES}

.PHONY: all test clean default plot
