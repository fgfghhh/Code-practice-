#pragma once
#include<iostream>
using namespace std;

// ¸¡µãÊıÅĞÍ¬
int double_equals(double a, double b)
{
	static const double ZERO = 1e-9;
	return fabs(a - b) < ZERO;
}