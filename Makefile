CC=g++
SOURCE=$(wildcard src/*.cpp src/*/*.cpp)
OBJS=$(SOURCE:src/%.cpp=build/%.o)
FLAGS=-Wall -pedantic -Wextra -O2
OUT_DIRS=$(wildcard src/*/)
DIRS=build $(OUT_DIRS:src/%=build/%)

build/main: $(shell mkdir -p $(DIRS)) build/utils/tinyexpr.o build/main.o $(OBJS)
	$(CC) $^ $(FLAGS) -o $@
	./build/main

build/utils/tinyexpr.o: src/utils/tinyexpr.c src/utils/tinyexpr.h
	gcc -c $< -O2 -o $@

build/main.o: main.cpp
	$(CC) $< -c $(FLAGS) -o $@

build/%.o: src/%.cpp src/%.hpp
	$(CC) $< -c $(FLAGS) -o $@

run:
	./build/main

clean:
	rm build/*.o build/*/*.o build/main
