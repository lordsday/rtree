SRC = $(wildcard *.c)
OBJ = $(SRC:.c=.o)

CFLAGS = -g -O0
LDFLAGS = -lm

TARGET = test

all: $(TARGET)

test: $(OBJ)
	$(CC) -o $@ $^ $(LDFLAGS)

clean:
	rm -rf *.o $(TARGET)
