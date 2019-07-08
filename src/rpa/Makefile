CC=g++
CFLAG=-Ofast
LFLAG=-ltrng4

run:./calc_pi.exe
	./calc_pi.exe

OBJECTS=../function/function.o ../function/grid.o
MAIN_OBJECTS=./calculate_pi.o ./w_expansion.o
../function/function.o:../function/function.cpp ../function/function.hpp \
	../function/grid.hpp
	$(CC) -c $(CFLAG) $*.cpp -o $@
../function/grid.o:../function/grid.cpp ../function/grid.hpp
	$(CC) -c $(CFLAG) $*.cpp -o $@
./calculate_pi.o:../function/grid.hpp ../function/function.hpp \
	./calculate_pi.cpp
	$(CC) -c $(CFLAG) $*.cpp -o $@
./w_expansion.o:../function/grid.hpp ../function/function.hpp \
	./calculate_pi.cpp
	$(CC) -c $(CFLAG) $*.cpp -o $@
./calc_pi.exe:$(OBJECTS) ./calculate_pi.o
	$(CC) -o ./calc_pi.exe $^ $(LFLAG)
./pi.txt:./calc_pi.exe
	./calc_pi.exe
./w_exp.exe:$(OBJECTS) ./w_expansion.o
	$(CC) -o ./w_exp.exe $^ $(LFLAG)

compile:$(OBJECTS) $(MAIN_OBJECTS)

link:./calc_pi.exe ./w_exp.exe


clean:
	rm ./*.exe $(OBJECTS) $(MAIN_OBJECTS)

.PHONY: clean compile link run result
