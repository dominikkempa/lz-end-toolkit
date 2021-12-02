SHELL = /bin/sh
CC = g++
CFLAGS = -Wall -Wextra -Wshadow -funroll-loops -DNDEBUG -O3 -march=native -pthread -std=c++0x
CFLAGS_DEBUG = -Wall -Wextra -Wshadow -g2 -pthread -std=c++0x

all: test test_debug

test:
	$(CC) $(CFLAGS) -o test ./src/main.cpp ../../src/utils.cpp

test_debug:
	$(CC) $(CFLAGS_DEBUG) -o test_debug ./src/main.cpp ../../src/utils.cpp

clean:
	/bin/rm -f *.o

nuclear:
	/bin/rm -f test test_debug *.o
