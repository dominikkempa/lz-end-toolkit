SHELL = /bin/sh
CC = g++
CFLAGS = -Wall -Wextra -Wshadow -funroll-loops -DNDEBUG -O3 -march=native -pthread -std=c++0x
#CFLAGS = -Wall -Wextra -Wshadow -g2 -std=c++0x

all: compute_lzend

compute_lzend:
	$(CC) $(CFLAGS) -o compute_lzend ./src/main.cpp ./src/utils.cpp

clean:
	/bin/rm -f *.o

nuclear:
	/bin/rm -f compute_lzend *.o
