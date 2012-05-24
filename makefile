LIB = libtalesf.so
PROG = talesf

default:
	gcc -g -O3 -Wall -m64 -o $(LIB) Array.c Hashmap.c talesf.c -lm -lz -fopenmp -fPIC -shared -rdynamic

frontend:
	gcc -g -O3 -Wall -m64 -I /usr/include/talesf -o $(PROG) frontend.c -ltalesf

clean:
	rm -f *.o *~ $(LIB)

install:
	install $(LIB) /usr/lib
	mkdir -p /usr/include/talesf
	cp *.h /usr/include/talesf
	chmod 644 /usr/include/talesf/*
	ldconfig