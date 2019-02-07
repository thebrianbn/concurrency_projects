CC=g++
CFLAGS= -O3 -mavx
DEPS = 
OBJ = hw02_slow_code.o

%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

hw02_slow_code: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)

.PHONY: clean
clean:
	rm -f $(OBJ) hw02_slow_code
