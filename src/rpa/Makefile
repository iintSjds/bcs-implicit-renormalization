CC=g++
CFLAG=-O3 -w
LFLAG=-ltrng4 -lhdf5 -lhdf5_cpp

result:delta.h5

run:./it_delta.exe w.h5 delta.in
	./it_delta.exe

delta.h5:./it_delta.exe w.h5 delta.in
	./it_delta.exe

w.h5:./w_exp.exe pi.h5 w.in
	./w_exp.exe

pi.h5:./calc_pi.exe pi.in
	./calc_pi.exe

delta.in w.in pi.in:config.py
	python config.py

OBJECTS=../function/function.o ../function/grid.o
MAIN_OBJECTS=./calculate_pi.o ./w_expansion.o ./iterate_delta.o

../function/function.o:../function/function.cpp ../function/function.hpp \
	../function/grid.hpp
	$(CC) -c $(CFLAG) $*.cpp -o $@
../function/grid.o:../function/grid.cpp ../function/grid.hpp
	$(CC) -c $(CFLAG) $*.cpp -o $@
./calculate_pi.o:../function/grid.hpp ../function/function.hpp \
	./calculate_pi.cpp
	$(CC) -c $(CFLAG) $*.cpp -o $@
./w_expansion.o:../function/grid.hpp ../function/function.hpp \
	./w_expansion.cpp
	$(CC) -c $(CFLAG) $*.cpp -o $@
./iterate_delta.o:../function/grid.hpp ../function/function.hpp \
	./iterate_delta.cpp
	$(CC) -c $(CFLAG) $*.cpp -o $@

./calc_pi.exe:$(OBJECTS) ./calculate_pi.o
	$(CC) -o ./calc_pi.exe $^ $(LFLAG)
./w_exp.exe:$(OBJECTS) ./w_expansion.o
	$(CC) -o ./w_exp.exe $^ $(LFLAG)
./it_delta.exe:$(OBJECTS) ./iterate_delta.o
	$(CC) -o ./it_delta.exe $^ $(LFLAG)

compile:$(OBJECTS) $(MAIN_OBJECTS)

link:./calc_pi.exe ./w_exp.exe ./it_delta.exe
clean:
	rm ./*.exe $(OBJECTS) $(MAIN_OBJECTS) *.h5 *.txt *.in

.PHONY: clean compile link run result
