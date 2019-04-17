CC=nvcc
CFLAGS=
DEPS = 
OBJ = hw5_part1_bbn5024.o

%.o: %.cu $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

hw5_part1_bbn5024: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)

.PHONY: clean
clean:
	rm -f $(OBJ) hw5_part1_bbn5024
