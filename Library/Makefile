
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

# SRC = tree.c   \
#       readwrite.c \
#       partition.c \
#       findnode.c\
#       reroot.c \
#       coevol_EpicsMat.c\
#       coevol_FindML.c\
#       coevol_lnLbranch.c\
#       simul.c\
#       m1qn3.c\
#       events.c\
# 	  coevol_sparseMat.c\
# 	  coevol_FindML_multi.c\
# 	  input_output.c

SRC = tree.c   \
    readwrite.c \
    reroot.c \
	findnode.c\
    coevol_EpicsMat.c\
    coevol_lnLbranch.c\
    events.c\
	coevol_sparseMat.c\
	coevol_FindML_multi.c\
	input_output.c

OBJECTS = $(SRC:.c=.o)

CC = gcc

#CFLAGS = -O2 -Wall --pedantic
CFLAGS = -g -Wall --pedantic
#CFLAGS = -pg -Wall --pedantic
# trouver pourquoi -static ne marche plus sur Darwin !!!

all: libtree.a $(PROGS)

libtree.a: $(OBJECTS)
	if [ -e $@ ] ; then rm $@ ; fi;
	$(AR) -r $@ $+ ;
	$(RANLIB) $@;


%.o: %.c
	$(CC) $(CFLAGS) -D$(MACHINE) -lm  -c -o $@ $<;

clean:
	rm -f *.o *.a *~ $(PROGS)
