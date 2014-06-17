
c = /usr/bin/gcc -ansi -fPIC -Wall -g
libs = yee.so
dependencies = yee.c

$(libs): $(dependencies)
	$(c) $(dependencies) --shared -o $@

test:
	$(c) $(dependencies) -o $@

clean:
	rm -f *.so
	rm -f *.o

rebuild: clean $(libs)
