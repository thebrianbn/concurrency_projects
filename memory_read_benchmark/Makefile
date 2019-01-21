CC=g++
CFLAGS= -O3 -mavx
DEPS = 
OBJ = vector_triad.o

%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

vector_triad: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)

.PHONY: clean
clean:
	rm -f $(OBJ) vector_triad