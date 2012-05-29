LIB = libtalesf.so
PROG = talesf

all: default frontend

default:
	gcc -fmax-errors=1 -std=gnu99 -g -O0 -Wall -m64 -o $(LIB) Array.c Hashmap.c talesf.c -lm -lz -fopenmp -fPIC -shared -rdynamic

frontend:
	gcc -fmax-errors=1 -std=gnu99 -g -O0 -Wall -m64 -I /usr/include/talesf -o $(PROG) frontend.c -ltalesf

clean:
	rm -f *.o *~ $(LIB) $(PROG)

install:
	install $(LIB) /usr/lib
	mkdir -p /usr/include/talesf
	cp *.h /usr/include/talesf
	chmod 644 /usr/include/talesf/*
	ldconfig