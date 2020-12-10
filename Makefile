CC=gcc
CFLAGS=-I -lm.
DEPS = nrutil.h

Bose: main.c nutil.c odeintM.c rkckM.c rkqsM.c odeintW.c rkckW.c rkqsW.c realft.c four1.c sinft.c
	$(CC) -lm -o Bose main.c nutil.c odeintM.c rkckM.c rkqsM.c odeintW.c rkckW.c rkqsW.c realft.c four1.c sinft.c -I.
