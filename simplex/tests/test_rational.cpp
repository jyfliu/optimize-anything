#include <iostream>

#include "../src/rational.h"
typedef Simplex::rational<long long> R;

// I probably should've made / used a testing harness
// whatever enjoy this mess

// test constructors and normalize()
void test_constructors() {
#define A(x) assert(x and "test_constructors ")
  R r0;
  A(r0.numerator() == 0);
  A(r0.denominator() == 1);
  R r1(9223372036854775807L);
  A(r1.numerator() == 9223372036854775807L);
  A(r1.denominator() == 1);
  R r2(-9223372036854775807L);
  A(r2.numerator() == -9223372036854775807L);
  A(r2.denominator() == 1);
  R r3(1);
  A(r3.numerator() == 1);
  A(r3.denominator() == 1);

  R r4(1, 4);
  A(r4.numerator() == 1);
  A(r4.denominator() == 4);
  R r5(-3, 1);
  A(r5.numerator() == -3);
  A(r5.denominator() == 1);
  R r6(4, 8);
  A(r6.numerator() == 1);
  A(r6.denominator() == 2);
  R r7(27, 3);
  A(r7.numerator() == 9);
  A(r7.denominator() == 1);
  R r8(-16, 12);
  A(r8.numerator() == -4);
  A(r8.denominator() == 3);
  R r9(1800, -3520);
  A(r9.numerator() == -45);
  A(r9.denominator() == 88);
  R r10(256, -91);
  A(r10.numerator() == -256);
  A(r10.denominator() == 91);
  R r11(-1333, -5000);
  A(r11.numerator() == 1333);
  A(r11.denominator() == 5000);
  R r12(-16, -30);
  A(r12.numerator() == 8);
  A(r12.denominator() == 15);
  R r13(0, 99999);
  A(r13.numerator() == 0);
  A(r13.denominator() == 1);
  try {
    R r14(0, 0);
    A(false and "should have divided by zero");
  } catch (Simplex::division_by_zero e) {}
  try {
    R r15(300, 0);
    A(false and "should have divided by zero");
  } catch (Simplex::division_by_zero e) {}
#undef A
}

// test + - * / unary+ unary- ++ --
void test_arithmetic() {
#define A(x) assert(x and "test_arithmetic ")

#undef A
}

// test < <= == >= > !=
void test_comparisons() {
#define A(x) assert(x and "test_comparisons ")

#undef A
}

int main() {
  test_constructors();
  test_arithmetic();
  test_comparisons();
  std::cout<<"All tests passed." << std::endl;
}

