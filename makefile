LIB = libtalesf.so

default:
	gcc -g -O3 -Wall -m64 -o $(LIB) *.c -lm -lz -fopenmp -fPIC -shared -rdynamic

clean:
	rm -f *.o *~ $(LIB)

install:
	install $(LIB) /usr/lib
	cp talesf.h /usr/include
	chmod 644 /usr/include/talesf.h
	ldconfig
