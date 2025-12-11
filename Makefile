# Makefile

MPICXX = /usr/bin/g++
CXXFLAGS = -c -Wall -std=c++11 #-fopenmp

ORG = ./
SRC = ./

mosolid: main.o element.o gibbse.o initcomp.o massb_oxide.o solution.o CGgibbsmin.o
	$(MPICXX) -Wall main.o element.o gibbse.o initcomp.o massb_oxide.o solution.o CGgibbsmin.o -o mosolid
main.o: main.cpp
	$(MPICXX) $(CXXFLAGS) main.cpp
element.o: $(ORG)/element.cpp
	$(MPICXX) $(CXXFLAGS) $(ORG)/element.cpp
gibbse.o: $(SRC)/gibbse.cpp
	$(MPICXX) $(CXXFLAGS) $(SRC)/gibbse.cpp
initcomp.o: $(SRC)/initcomp.cpp
	$(MPICXX) $(CXXFLAGS) $(SRC)/initcomp.cpp
massb_oxide.o: $(SRC)/massb_oxide.cpp
	$(MPICXX) $(CXXFLAGS) $(SRC)/massb_oxide.cpp
solution.o: $(SRC)/solution.cpp
	$(MPICXX) $(CXXFLAGS) $(SRC)/solution.cpp
CGgibbsmin.o: $(SRC)/CGgibbsmin.cpp
	$(MPICXX) $(CXXFLAGS) $(SRC)/CGgibbsmin.cpp
clean:
	rm -f *.o mosolid
