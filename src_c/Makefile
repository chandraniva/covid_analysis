
CC=gcc
CFLAGS=-Wall
OUTFILE=binary

all: fitting

fitting: fitting.o model.o readdata.o
	$(CC) $(CFLAGS) -lm fitting.o model.o readdata.o -o $(OUTFILE)

fitting.o: fitting.c
	$(CC) $(CFLAGS) -c fitting.c

model.o: model.c
	$(CC) $(CFLAGS) -c model.c

readdata.o: readdata.c
	$(CC) $(CFLAGS) -c readdata.c

clean:
	rm *.o $(OUTFILE)

