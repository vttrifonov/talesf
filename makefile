LIB = libpairedtalesf.so
PROG = pairedtalesf

all: default

default:
	gcc -fmax-errors=1 -std=gnu99 -g -O0 -Wall -m64 -o $(LIB) Array.c Hashmap.c pairedtalesf.c -lm -lz -fopenmp -fPIC -shared -rdynamic

frontend:
	gcc -fmax-errors=1 -std=gnu99 -g -O0 -Wall -m64 -I /usr/include/pairedtalesf -o $(PROG) frontend.c -lpairedtalesf -fopenmp

clean:
	rm -f *.o *~ $(LIB) $(PROG)

install:
	install $(LIB) /usr/lib
	mkdir -p /usr/include/pairedtalesf
	cp *.h /usr/include/pairedtalesf
	chmod 644 /usr/include/pairedtalesf/*
	ldconfig