CFLAGS = -Wall -g -DINLINE
MFLAGS = -lm
CC = gcc $(CFLAGS)

nrlmsise-test    :     nrlmsise-00.o nrlmsise-00_test.o nrlmsise-00_data.o
	$(CC) -o nrlmsise-test nrlmsise-00.o nrlmsise-00_test.o \
                     nrlmsise-00_data.o $(MFLAGS)

nrlmsise-00.o :		nrlmsise-00.c nrlmsise-00.h
	$(CC) -c  nrlmsise-00.c

nrlmsise-00_test.o :	nrlmsise-00_test.c nrlmsise-00.h
	$(CC) -c  nrlmsise-00_test.c

nrlmsise-00_data.o :	nrlmsise-00_data.c nrlmsise-00.h
	$(CC) -c  nrlmsise-00_data.c

clean   :
	rm -rf  nrlmsise-test nrlmsise-00.o nrlmsise-00_test.o \
                   nrlmsise-00_data.o
