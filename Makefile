CC=			gcc
CFLAGS=		-g -Wall -O2
DFLAGS=		-DHAVE_CONFIG_H
OBJS=		error.o sort.o
PROG=		sort_lh3
LIBS=		

.SUFFIXES:.c .o

.c.o:
		$(CC) -c $(CFLAGS) $(DFLAGS) $< -o $@

all:$(PROG)

sort_lh3:$(OBJS)
		$(CC) $(CXXFLAGS) $(DFLAGS) $(OBJS) -o $@ $(LIBS)


clean:
		rm -f *.o a.out *~ *.a $(PROG)
