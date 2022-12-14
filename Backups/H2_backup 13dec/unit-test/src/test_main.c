#include <stdlib.h>
#include <string.h>

/*****************************************************************************
 * Add your functions that you wanna test here (from, e.g., src/run.c.)
 *
 * Example:
 *
 * extern void get_x(...);
 *
 * ***************************************************************************/

#include "test_main.h"

/* ************************************
 * Here follows the test we want
 * to run
 * ***********************************/
START_TEST(test_what_you_wanna_test)
{
    // empty
}


int
main()
{
    test_setup("testing", "core");
    
    // Tests
    add_test(test_what_you_wanna_test);
    
    test_teardown();
    return 0;
}
