LIB = libtalesf.so
PROG = talesf

default:
	gcc -g -O3 -Wall -m64 -o $(LIB) talesf.c -lbcutils -lm -lz -fopenmp -fPIC -shared -rdynamic

frontend:
	gcc -g -O3 -Wall -m64 -I /usr/include/talesf -o $(PROG) frontend.c -lbcutils -ltalesf

clean:
	rm -f *.o *~ $(LIB)

install:
	install $(LIB) /usr/lib
	mkdir -p /usr/include/talesf
	cp *.h /usr/include/talesf
	chmod 644 /usr/include/talesf/*
	ldconfig
