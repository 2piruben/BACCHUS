# COMPILATION OF BACPOP

GSL_INCLUDE = /usr/local/include
GSL_LIBRARY= /usr/local/lib
CFLAGS = -DNDEBUG -O -I$(GSL_INCLUDE) -I. -std=c++20 -Wall -msse2
LFLAGS = -lm -lgsl -lgslcblas -L$(GSL_LIBRARY)
GPP = g++
RM = rm -f

# COMPILATION OF LANGIL LIBRARY, MAKES USE OF ITS OWN MAKEFILE


main: main.o
	$(GPP) $(CFLAGS) $(LFLAGS) -o main main.o bacterium.o algebra2d.o population.o molecular.o dish.o diffusible.o

main.o: main.cpp bacterium.o population.o algebra2d.o molecular.o dish.o diffusible.o
	$(GPP) $(CFLAGS) -c main.cpp

dish.o: dish.cpp dish.h population.o
	$(GPP) $(CFLAGS) -c dish.cpp	

population.o: population.cpp population.h bacterium.o algebra2d.o
	$(GPP) $(CFLAGS) -c population.cpp

bacterium.o: bacterium.cpp bacterium.h algebra2d.o molecular.o
	$(GPP) $(CFLAGS) -c bacterium.cpp

molecular.o: molecular.cpp molecular.h
	$(GPP) $(CFLAGS) -c molecular.cpp

diffusible.o: diffusible.cpp diffusible.h
	$(GPP) $(CFLAGS) -c diffusible.cpp

algebra2d.o: algebra2d.cpp algebra2d.h
	$(GPP) $(CFLAGS) -c algebra2d.cpp

all:
	$(MAKE) -C $(LANGILFOLDER)

cleanout:
	$(RM) -r output/**/*.out
	$(RM) -r output/**/*.png
	$(RM) log

clean:
	$(RM) *.o
	$(RM) log

cleanall:
	$(RM) *.o
	$(RM) $(LANGIL_OBJECT) 



