SHELL = /bin/sh
CC = g++
CFLAGS = -Wall -Wextra -Wshadow -funroll-loops -DNDEBUG -O3 -march=native -pthread -std=c++0x
#CFLAGS = -Wall -Wextra -Wshadow -g2 -std=c++0x

all: parse verify decode

parse:
	$(CC) $(CFLAGS) -o parse ./src/main_parse.cpp ./src/lz_end_toolkit_src/utils.cpp

verify:
	$(CC) $(CFLAGS) -o verify ./src/main_verify.cpp ./src/lz_end_toolkit_src/utils.cpp

decode:
	$(CC) $(CFLAGS) -o decode ./src/main_decode.cpp ./src/lz_end_toolkit_src/utils.cpp

clean:
	/bin/rm -f *.o

nuclear:
	/bin/rm -f parse verify decode *.o
