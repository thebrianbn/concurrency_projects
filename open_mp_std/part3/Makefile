CC=g++
CFLAGS= -O3 -mavx
DEPS = 
OBJ = hw3_part3_omp_bbn5024.o

%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

hw3_part3_omp_bbn5024: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)

.PHONY: clean
clean:
	rm -f $(OBJ) hw3_part3_omp_bbn5024
