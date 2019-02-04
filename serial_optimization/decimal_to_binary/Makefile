CC=g++
CFLAGS= -O3 -mavx
DEPS = 
OBJ = dec2bin.o

%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

dec2bin: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)

.PHONY: clean
clean:
	rm -f $(OBJ) dec2bin 
