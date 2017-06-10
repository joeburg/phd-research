# Voro++ makefile

# List of executables
EXECUTABLES=main

# Makefile rules
all: $(EXECUTABLES)

main: main.cpp
	g++ -g -std=c++11 -Wall -ansi -pedantic -O3 -I../../src -L../../src -o main main.cpp -lvoro++

clean:
	rm -f $(EXECUTABLES)
