CC=mpicc
CFLAGS= -O3 -mavx
DEPS = 
OBJ = hw4_part2_bbn5024.o

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

hw4_part2_bbn5024: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)

.PHONY: clean
clean:
	rm -f $(OBJ) hw4_part2_bbn5024
