CC=gcc
CFLAGS= -O3 -mavx -pthread
DEPS = 
OBJ = vector_triad.o

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

vector_triad: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)

.PHONY: clean
clean:
	rm -f $(OBJ) vector_triad