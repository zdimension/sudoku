OBJS	= sudoku.o
SOURCE	= sudoku.cpp
OUT	= sudoku
CC	 = g++
FLAGS	 = -g -c -Wall

all: $(OBJS)
	$(CC) -g $(OBJS) -o $(OUT) $(LFLAGS)

sudoku.o: sudoku.cpp
	$(CC) $(FLAGS) sudoku.cpp -std=c++17


clean:
	rm -f $(OBJS) $(OUT)
