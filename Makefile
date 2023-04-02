# tem que ter um diretório de build com essa especificação
# \build
#  ---\utils
#  ---\derivatives
#  ---\integration
CC=g++
SOURCE=$(wildcard src/*.cpp src/*/*.cpp)
OBJS=$(SOURCE:src/%.cpp=build/%.o)
FLAGS=-Wall -pedantic -Wextra -O2

build/main: build/utils/tinyexpr.o build/utils/stb.o build/main.o $(OBJS)
	$(CC) $^ $(FLAGS) -o $@
	./build/main

build/utils/tinyexpr.o: src/utils/tinyexpr.c src/utils/tinyexpr.h
	gcc -c $< -O2 -o $@

build/main.o: main.cpp
	$(CC) $< -c $(FLAGS) -o $@

build/%.o: src/%.cpp src/%.hpp
	$(CC) $< -c $(FLAGS) -o $@

build/utils/stb.o: src/utils/stb.cpp
	g++ $< -O2 -c -o $@

run:
	./build/main

clean:
	rm build/*.o build/*/*.o build/main
