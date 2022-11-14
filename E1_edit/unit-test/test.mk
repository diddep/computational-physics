TEST = \
       obj/test_main.o


CFLAGS += \
	-fsanitize=address \
	-fno-omit-frame-pointer \
	-Iunit-test/include/ \
	-O0
LIB += \
     -lcheck


test: obj run-test

run-test: $(OBJ) $(TEST)
	$(CC) $(CFLAGS) $^ -o $@ $(LIB)

obj/%.o: unit-test/src/%.c
	$(CC) -MMD -c $(CFLAGS) $< -o $@ 
