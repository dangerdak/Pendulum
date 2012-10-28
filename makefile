Pendulum: main.o Pendulum.o
	g++ -g -Wall main.o Pendulum.o -o Pendulum

main.o: main.cpp Pendulum.h
	g++ -g -Wall -c main.cpp

Pendulum.o: Pendulum.cpp Pendulum.h
	g++ -g -Wall -c Pendulum.cpp
