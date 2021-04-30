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

void test_Vector_parameterized_constructor_invalid_arg(void)
{
    unsigned int vector_len = 0;
    TEST_EXCEPTION(Vector{vector_len}, std::invalid_argument);
}

void test_Vector_stdvector_constructor(void)
{
    int num_elements = 10;
    std::vector<double> vec;
    for (int i = 0; i < num_elements; i++)
    {
        vec.push_back(i);
    }

    Vector v(vec);

    TEST_CHECK_(v.len() == num_elements, "v.len(): Expected=%d, Actual=%d", num_elements, v.len());
    TEST_CHECK_(vec.size() == 0, "vec.size(): Expected=%d, Actual=%lu", 0, vec.size());

    for (int i = 0; i < num_elements; i++)
    {
        TEST_CHECK(v[i] == i);
    }
}




// DO NOT MODIFY CODE BELOW
TEST_LIST =
{
    { "void test_Vector_default_constructor(void);", test_Vector_default_constructor },
    { "void test_Vector_parameterized_constructor(void);", test_Vector_parameterized_constructor },
    { "void test_Vector_parameterized_constructor_invalid_arg(void);", test_Vector_parameterized_constructor_invalid_arg },
    { "void test_Vector_stdvector_constructor(void);", test_Vector_stdvector_constructor },
    { NULL }
};