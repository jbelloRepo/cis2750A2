CC = clang
CFLAGS = -Wall -pedantic -std=c99
# DEPS = mol.h
LIBS=-lm # note: the l means library, m means math

# all: main
all: myprog 

clean:  
	rm -f *.o *.so myprog molecule_wrap.c molecule.py

# Object file for mol.c	
mylib.o:  mol.c mol.h
	$(CC) $(CFLAGS) -c mol.c -fPIC -o mylib.o

# Shared library 
libmol.so: mylib.o
	$(CC) mylib.o -shared -o libmol.so $(LIBS)

# Generates python interface with C-code
molecule_python: 
	swig3.0 -python molecule.i

# Compiles *_wrap.c file 
molecule_wrap.o: molecule_wrap.c
	$(CC) $(FLAGS) -c molecule_wrap.c -I/usr/include/python2.7  -o molecule_wrap.o 

_molecule.so: libmol.so molecule_wrap.o
	$(CC) $(FLAGS) -shared libmol.so molecule_wrap.o -o _molecule.so

# For testing 
main.o:  main.c mol.h
	$(CC) $(FLAGS) -c main.c -o main.o 

myprog:  molecule_python main.o libmol.so _molecule.so
	$(CC) main.o -L. -l mol -o myprog $(LIBS)


# export LD_LIBRARY_PATH=`pwd`
# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$(pwd):/usr/lib/python3.7/config-3.7m-x86_64-linux-gnu/:/usr/include/python3.7m