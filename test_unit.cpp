#include "acutest.h"
#include "Matrix.hpp"

void test_Vector_default_constructor(void)
{
    Vector v;
    unsigned int vector_len = 0;
    TEST_CHECK_(v.len() == vector_len, "v.len(): Expected=%d, Actual=%d", vector_len, v.len());
}

void test_Vector_parameterized_constructor(void)
{
    unsigned int vector_len = 10;
    Vector v(vector_len);
    TEST_CHECK_(v.len() == 10, "v.len(): Expected=%d, Actual=%d", vector_len, v.len());
}







// DO NOT MODIFY CODE BELOW
TEST_LIST =
{
    { "void test_Vector_default_constructor(void);", test_Vector_default_constructor },
    { "void test_Vector_parameterized_constructor(void);", test_Vector_parameterized_constructor },
    { NULL }
};