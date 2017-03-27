CC=gcc
CXX=g++
override CFLAGS += -Wall -g -std=c99
override CXXFLAGS += -Wall -g -std=c99
override LDFLAGS += -lm

BIN_FILES=lagrange_test.bin
OBJECT_FILES=lagrange.o lagrange_test.o
OUT_FILES=test_full.out test_fo.out
ORDER_X=2
ORDER_Y=2
NX=101
NY=101

default: all
all: lagrange_test.bin

test: test_fo.out
plot: test_fo.out
	python plot_error.py test_fo.out ${NX} ${NY}
test_full.out: lagrange_test.bin
	./lagrange_test.bin ${ORDER_X} ${ORDER_Y}

lagrange_test.bin: lagrange_test.o lagrange.o
	$(CC) ${CFLAGS} -o $@ $^ ${LDFLAGS} 

lagrange_test.o: lagrange_test.c lagrange.h
lagrange.o: lagrange.h

clean:
	$(RM) ${OBJECT_FILES}
	$(RM) ${BIN_FILES}
	$(RM) ${OUT_FILES}

.PHONY: all test clean default plot
