LIB = libpairedtalesf.so
PROG = pairedtalesf

all: default

default:
	gcc -fmax-errors=1 -std=gnu99 -g -O3 -Wall -m64 -o $(LIB) pairedtalesf.c -lbcutils -lm -lz -fopenmp -fPIC -shared -rdynamic

gpu:
	gcc -DPTF_GPU_ENABLED -fmax-errors=1 -std=gnu99 -g -O3 -Wall -m64 -o $(LIB) pairedtalesf.c -lbtfcount -lbcutils -lm -lz -fopenmp -fPIC -shared -rdynamic

frontend:
	gcc -fmax-errors=1 -std=gnu99 -g -O3 -Wall -m64 -I /usr/include/pairedtalesf -o $(PROG) frontend.c -lbcutils -lpairedtalesf -fopenmp

clean:
	rm -f *.o *~ $(LIB) $(PROG)

install:
	install $(LIB) /usr/lib
	mkdir -p /usr/include/pairedtalesf
	cp *.h /usr/include/pairedtalesf
	chmod 644 /usr/include/pairedtalesf/*
	ldconfig
