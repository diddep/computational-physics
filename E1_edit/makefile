CC = gcc

CFLAGS = \
	 -pedantic \
	 -Wall \
	 -Werror \
	 -g \
	 -Iinclude

CFLAGS_OPT = \
	     -O2 \
	     -march=native

LIB = \
	  -lm \
	  -lgsl \
      -lgslcblas \

OBJ = \
      obj/tools.o

E1code1 = \
       obj/ht_gen_E1code1.o

E1code3 = \
       obj/pow_gen_E1code3.o


	   

E1code4 = \
       obj/vel_verlet_E1code4.o

ifeq ($(MAKECMDGOALS),test)
-include unit-test/test.mk
else
CFLAGS += $(CFLAGS_OPT)
endif

_ht_gen: $(E1code1) $(OBJ)

	$(CC) $(CFLAGS) $^ -o ht_gen $(LIB)
ht_gen: obj _ht_gen

_pow_gen: $(E1code3) $(OBJ)

	$(CC) $(CFLAGS) $^ -o pow_gen $(LIB)
pow_gen: obj _pow_gen

_vel_verlet: $(E1code4) $(OBJ)

	$(CC) $(CFLAGS) $^ -o vel_verlet $(LIB)
vel_verlet: obj _vel_verlet

obj/%.o: src/%.c
	$(CC) -MMD -c $(CFLAGS) $< -o $@ 

obj:
	mkdir -p obj

clean:
	find -iname "*.o" -exec rm {} \;
	find -iname "*.d" -exec rm {} \;
	rm -f program run-test
	rm -rf obj

.PHONY: clean
