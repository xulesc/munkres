CC=gcc
CFLAGS=-std=c99 -O3 -ffast-math -lm 
PYCFLAGS=--shared -fPIC -D __PYMOD__ -I/usr/include/python2.7/ -lpython2.7 

munkres.so: munkres.c
	$(CC) -o $@ $< $(CFLAGS) $(PYCFLAGS)
	
exe: munkres.c
	$(CC) $< $(CFLAGS) -D __STANDALONE__ -D __CHATTY__ 

run: munkers.so psalign.py
	python psalign.py

clean:
	rm munkres.so
