CC=g++
CFLAGS=-Wall -g -std=c++11

test_unit: test_unit.o Matrix.o
	$(CC) $(CFLAGS) $^ -o $@

test_unit.o: test_unit.cpp acutest.h
	$(CC) $(CFLAGS) -c $< -o $@

Matrix.o: Matrix.cpp Matrix.hpp
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -rf *.o test_unit