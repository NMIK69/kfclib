CC=gcc
CFLAGS=-Wall -pedantic -std=c99 -Wextra -Wmissing-prototypes 
OPTIMIZE=-O3
LDFLAGS=-lm

SRCS=$(wildcard ../kfclib/*.c)
OBJS=$(patsubst ../kfclib/%.c, ./%.o, $(SRCS))

TARGET=example

all: $(TARGET)

debug: CFLAGS+=-ggdb
debug: OPTIMIZE=-O0
debug: $(TARGET)

$(TARGET) : $(OBJS) example.o
	$(CC) $^ $(OPTIMIZE) $(LDFLAGS) -o $@

%.o : %.c
	$(CC) $(CFLAGS) $(OPTIMIZE) -c $< -o $@

%.o : ../kfclib/%.c
	$(CC) $(CFLAGS) $(OPTIMIZE) -c $< -o $@

.PHONY : clean
clean :
	rm -f $(TARGET)
	rm -f *.o

