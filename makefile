NAME=pairedtalesf
CC = gcc-9
CFLAGS = -I . -fmax-errors=1 -std=gnu99 -g -O3 -Wall -m64
BCUTILS = -L bcutils -lbcutils
NAMELIB = -L lib -l$(NAME)

all: default

default:
	$(CC) $(CFLAGS) -o lib/lib$(NAME).so $(NAME).c $(BCUTILS) -lm -lz -fopenmp -fPIC -shared -rdynamic

frontend:
	$(CC) $(CFLAGS) -o bin/$(NAME) frontend.c $(BCUTILS) $(NAMELIB) -fopenmp

clean:
	rm -rf *.o *~ *.so *.dSYM bin/* lib/* include/*

install:
	ln -sf ../bcutils/libbcutils.so lib/
	ln -sf $(PWD)/lib/lib$(NAME).so /usr/local/lib/
	ln -sf $(PWD)/lib/libbcutils.so /usr/local/lib/
	ln -sf $(PWD)/bin/$(NAME) /usr/local/bin/
