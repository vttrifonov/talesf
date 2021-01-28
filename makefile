CC = gcc-9
CFLAGS = -I . -fmax-errors=1 -std=gnu99 -g -O3 -Wall -m64
BCUTILS = -L lib -lbcutils
TALESF = -L bcutils -ltalesf

default:
	$(CC) $(CFLAGS) -o lib/libtalesf.so talesf.c $(BCUTILS) -lm -lz -fopenmp -fPIC -shared -rdynamic

frontend:
	$(CC) $(CFLAGS) -o bin/talesf frontend.c $(BCUTILS) $(TALESF) -fopenmp

clean:
	rm -rf *.o *~ *.so *.dSYM bin/* lib/* include/*

install:
	ln -sf ../bcutils/libbcutils.so lib/
	ln -sf $(PWD)/lib/libtalesf.so /usr/local/lib/
	ln -sf $(PWD)/lib/libbcutils.so /usr/local/lib/
	ln -sf $(PWD)/bin/talesf /usr/local/bin/

