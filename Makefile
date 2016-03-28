CPP=icpc
FLAGS=-O3 -std=c++11 -ltbb -gcc-name=gcc-4.6  -w0
all: main

main:
	$(CPP) $(FLAGS) -DCILKP main.cpp -L/opt/intel/composerxe/lib/intel64  -o ppll

clean:
	rm ppll
