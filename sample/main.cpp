#include <cstdio>

void xorshift_sample();
void mersenne_twister_sample();
void fibonacci_mod_m_sample();
void linear_recurrence_sample();
void linear_recurrence_sample2();
void series_simple_counting_sample();
void series_combinatoric_sample();
void dirichlet_sample();

void test_sample();

int main() {
    //test_sample();
    xorshift_sample();
	mersenne_twister_sample();
	fibonacci_mod_m_sample();
	linear_recurrence_sample();
	linear_recurrence_sample2();
	series_simple_counting_sample();
	series_combinatoric_sample();
	dirichlet_sample();
	return 0;
}
