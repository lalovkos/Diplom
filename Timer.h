#pragma once
#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <cstdio>

typedef std::chrono::high_resolution_clock hrc;

class Timer {
public:
	Timer() : begin(hrc::now()) {}

	void start() {
		begin = hrc::now();
	}

	double getTime() {
		return std::chrono::duration<double>(hrc::now() - begin).count();
	}

private:
	hrc::time_point begin;
};