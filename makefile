
CC       = g++
CFLAGS   = -g -I../lib -ldl -lm -I../lib/qd -lqd -llapack -lblas -std=c++11

all: qpspecial

.cpp.o: 
	$(CC) $(CFLAGS) -c $< -o $@

qpspecial : qpspecial.o
	$(CC) $(CFLAGS) qpspecial.o -o qpspecial -L. -L../lib -lqd -llapack -lblas

clean:
	@rm -f *.o core
	@rm -f *.exe core test
