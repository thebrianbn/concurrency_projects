CC=g++
CFLAGS= -O3 -mavx -std=c++11
DEPS = 
OBJ = hw4_part2_bbn5024.o

%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

hw4_part2_bbn5024: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)

.PHONY: clean
clean:
	rm -f $(OBJ) hw4_part2_bbn5024
