SOURCES = $(wildcard *.c)
BINARIES = $(patsubst %.c,%, $(SOURCES))
DEPENDENCIES = $(wildcard ../*.h) Makefile

CC = gcc
CFLAGS = -O0 -g
INCLUDES = -I../../../../debug \
	   -I../../../lib \
	   -I../../../lib/minmd_basic \
	   -I../../../lib/minmd_basic/def \
	   -I../../../lib/minmd_basic/utils \
	   -I../../../lib/minmd_math \
	   -I../../../lib/minmd_math/rng \
	   -I../../../lib/minmd_math/vec \
	   -I../../../lib/minmd_math/stat_accum \
	   -I../../../lib/minmd_md \
	   -I../../../lib/minmd_md/mdutils \
	   -I../../../lib/minmd_md/thermostat \
	   -I../../../lib/minmd_md/ewald \
	   -I..
LIBS = -lm

all: $(BINARIES)

$(BINARIES): % : %.c $(DEPENDENCIES)
	$(CC) $(CFLAGS) $(INCLUDES) $< -o $@ $(LIBS)

clean:
	rm -f $(BINARIES) *~ .*.un~
