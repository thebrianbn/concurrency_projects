CC=g++
CFLAGS= -O3 -mavx -pthread -std=c++11
DEPS = 
OBJ = hw3_part1_thread_bbn5024.o

%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

hw3_part1_thread_bbn5024: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)

.PHONY: clean
clean:
	rm -f $(OBJ) hw3_part1_thread_bbn5024
