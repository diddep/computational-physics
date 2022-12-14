#include <check.h>
#define test_setup(arg1, arg2) \
    Suite *s = suite_create((arg1)); \
    TCase *tc_core = tcase_create((arg2));

#define add_test(arg) \
    tcase_add_test((tc_core), (arg))

#define add_loop_test(arg, istart, iend) \
    tcase_add_loop_test((tc_core), (arg), (istart), (iend))

#define test_teardown() \
    suite_add_tcase((s), (tc_core)); \
    SRunner *sr = srunner_create((s)); \
    srunner_run_all((sr), (CK_NORMAL)); \
    srunner_free((sr));

#define check_vectors_equal(v1, v2, len, tol) { \
    for (int __i = 0; __i < (len); __i++) { \
        ck_assert_double_eq_tol((v1)[__i], (v2)[__i], (tol)); \
    } \
}
