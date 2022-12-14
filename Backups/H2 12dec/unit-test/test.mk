OBJ_TEST += \
       obj_test/test_main.o


CFLAGS += \
	-fsanitize=address \
	-fno-omit-frame-pointer \
	-Iunit-test/include/ \
	-O0
LIB += \
     -lcheck


test: obj_test run-test

run-test: $(OBJ_TEST)
	$(CC) $(CFLAGS) $^ -o $@ $(LIB)

obj_test/%.o: unit-test/src/%.c
	$(CC) -MMD -c $(CFLAGS) $< -o $@ 

obj_test/%.o: src/%.c
	$(CC) -MMD -c $(CFLAGS) $< -o $@ 

obj_test:
	mkdir -p obj_test
