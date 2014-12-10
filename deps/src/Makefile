# Makefile template for shared library

CFLAGS = -fPIC -Wall -Wextra -O2 -g # C flags
LDFLAGS = -shared 
RM = rm -f			# rm command
TARGET_LIB = libpedigree.so	# target lib

.PHONY: all
all: ${TARGET_LIB}

$(TARGET_LIB): inbreeding.o
	$(LD) ${LDFLAGS} -o $@ $^

install: ${TARGET_LIB}
	cp ${TARGET_LIB} ../usr/lib/

.PHONY: clean
clean:
	-${RM} ${TARGET_LIB}
