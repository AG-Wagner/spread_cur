spread_cur_test : spread_cur_test.o spread_cur.o
	cc -o spread_cur_test -lm spread_cur_test.o spread_cur.o

spread_cur.o : spread_cur.c spread_cur.h
	cc -c spread_cur.c
spread_cur_test.o : spread_cur_test.c spread_cur.h
	cc -c spread_cur_test.c

clean :
	rm spread_cur_test spread_cur_test.o spread_cur.o

