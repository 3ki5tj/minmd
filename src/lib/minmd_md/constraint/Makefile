SOURCES = $(wildcard *.c)
BINARIES = $(patsubst %.c,%, $(SOURCES))
DEPENDENCIES = $(wildcard *.h) Makefile

CC = gcc
CFLAGS = -O0 -g -Wall
INCLUDES = -I../../../../debug \
	   -I../.. \
	   -I../../minmd_basic \
	   -I../../minmd_basic/def \
	   -I../../minmd_basic/utils \
	   -I../../minmd_math \
	   -I../../minmd_math/rng \
	   -I../../minmd_math/vec \
	   -I../../minmd_math/stat_accum \
	   -I../../minmd_md \
	   -I../../minmd_md/mdutils
LIBS = -lm

all: $(BINARIES)

$(BINARIES): % : %.c $(DEPENDENCIES)
	$(CC) $(CFLAGS) $(INCLUDES) $< -o $@ $(LIBS)

clean:
	rm -f $(BINARIES) *~ .*.un~
