CC=g++
CFLAG=-Ofast
EXE=./calc_pi.exe
OBJ=../function/function.o ../function/grid.o ./calculate_pi.o

run:./calc_pi.exe
	./calc_pi.exe

./calc_pi.exe:../function/function.cpp ../function/function.hpp ../function/grid.hpp ../function/grid.cpp ./calculate_pi.cpp
	$(CC) -o $(EXE) $(CFLAG) ../function/function.cpp ../function/grid.cpp ./calculate_pi.cpp -ltrng4

.PHONY: clean
clean:
	rm ./*.exe
