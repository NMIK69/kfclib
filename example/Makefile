CC=gcc
CFLAGS=-Wall -pedantic -std=c99 -Wextra -Wmissing-prototypes 
OPTIMIZE=-O3
LDFLAGS=-lm

SRCS=$(wildcard ../*.c)
OBJS=$(patsubst ../%.c, ./%.o, $(SRCS))

TARGET=example

all: $(TARGET)

debug: CFLAGS+=-ggdb
debug: OPTIMIZE=-O0
debug: $(TARGET)

.PHONY: debug
.PHONY: clean
.PHONY: example_6dof
.PHONY: example_9dof

example_6dof : $(TARGET)
	./$(TARGET) data_6dof.txt 6
	python3 vis.py out.txt

example_9dof : $(TARGET)
	./$(TARGET) data_9dof.txt 9 
	python3 vis.py out.txt

$(TARGET) : $(OBJS) $(TARGET).o
	$(CC) $^ $(OPTIMIZE) $(LDFLAGS) -o $@

%.o : %.c
	$(CC) $(CFLAGS) $(OPTIMIZE) -c $< -o $@

%.o : ../%.c
	$(CC) $(CFLAGS) $(OPTIMIZE) -c $< -o $@

.PHONY : clean
clean :
	rm -f $(TARGET)
	rm -f *.o

