CC = gcc

CFLAGS = \
	 -g \
	 -Iinclude

CFLAGS_OPT = \
	     -O2 \
	     -march=native

CFLAGS_SANITIZE = \
	 	 -Werror \
		 -Wall \
	 	 -pedantic \
	 	 -g \
	 	 -fsanitize=address \
		 -fno-omit-frame-pointer \
		 -Iinclude

LIB = \
      -lm \
      -lgsl \
      -lgslcblas

OBJ = \
      obj/run.o \
      obj/tools.o
MAIN = \
       obj/main.o

ifeq ($(MAKECMDGOALS),test)
OBJ_TEST = $(patsubst obj/%.o,obj_test/%.o,$(OBJ))
-include unit-test/test.mk
else
CFLAGS += $(CFLAGS_OPT)
endif

E2: obj _E2

_E2: $(MAIN) $(OBJ)
	$(CC) $(CFLAGS_SANITIZE) $^ -o E2 $(LIB)

obj/%.o: src/%.c
	$(CC) -MMD -c $(CFLAGS) $< -o $@ 

obj:
	mkdir -p obj

clean:
	find -iname "*.o" -exec rm {} \;
	find -iname "*.d" -exec rm {} \;
	rm -f program run-test
	rm -rf obj obj_test

.PHONY: clean
