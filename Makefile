CC = gcc
FLAGS = -fopenmp -ffast-math -Ofast -march=native
INCLUDE = -I /home/obsuser/src/psrdada_ata/ata/src
LIBS =  -L /home/obsuser/src/psrdada_ata/ata/src/mysigproc -l:mysigproc_utils.o -lm

all: clean sumfils

clean:
	touch sumfils; rm sumfils

sumfils:
	$(CC) $(FLAGS) $(INCLUDE) sumfils.c -o sumfils $(LIBS)
