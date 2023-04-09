
spkmeans: spkmeans.o spkmeans.h
	gcc -o spkmeans spkmeans.o -lm

spkmeans.o: spkmeans.c
	gcc -ansi -Wall -Wextra -Werror -pedantic-errors -c spkmeans.c -lm

clean:
	rm -f *.o
