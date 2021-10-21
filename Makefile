
# Makefile for making libdiversity
# this library contain variuous function for
# get diversity among multiple alignement

# date   : 15 Aug 2002
# author : amikezor
# note   : need libinfoseq (read/transform multiple aligned sequences)

####
# MACHINE DEPEND FLAG
####

# MACHINE=SGI
# MACHINE=LINUX
# MACHINE=MACOSX
SYSTEME = `uname` # Darwin, Linux....
MACHINE = `uname | awk '/Darwin/ {print "MACOSX"} /Linux/ {print "LINUX"}'`


# SGI:    RANLIB=touch
# LINUX:  RANLIB=ranlib
# MACOSX: RANLIB=ranlib
RANLIB=ranlib

####
## MACHINE INDEPENDENT FLAGS
####

AR = ar

SHELL=/bin/sh

EPICSSRC = epics.c
EPOCSSRC = epocs.c

PROGS = epics epocs

LIBTREE = libtree.a
LIBDIR = Library/
INCDIR= Library
CC = gcc

#CFLAGS = -O2 -Wall --pedantic
CFLAGS = -g -Wall --pedantic -I$(INCDIR)
#CFLAGS = -pg -Wall --pedantic
# trouver pourquoi -static ne marche plus sur Darwin !!!

#LIBINCDIR =  ../libinfoseq


all: $(PROGS)


$(LIBDIR)/$(LIBTREE):
	cd Library; make

epics: $(EPICSSRC) $(LIBDIR)/$(LIBTREE)
	$(CC) $(CFLAGS) -D$(MACHINE)  -o $@ $< $(LIBDIR)/$(LIBTREE) -lm;

epocs: $(EPOCSSRC) $(LIBDIR)/$(LIBTREE)
	$(CC) $(CFLAGS) -D$(MACHINE) -o $@ $< $(LIBDIR)/$(LIBTREE) -lm;

%.o: %.c
	$(CC) $(CFLAGS) -D$(MACHINE) -I$(LIBINCDIR) -c -o $@ $< -lm;

clean:
	rm -f *.o $(LIBDIR)/*.o $(LIBDIR)/*.a $(PROGS)

install: all
	cp $(PROGS) ~/bin
