#pragma once
#include <string>

class Point {

public:
    double x;
    double y;
    double z;

    Point(const double X, const double Y, const double Z) {
        x = X;
        y = Y;
        z = Z;
    }

    Point(double* points) {
        x = points[0];
        y = points[1];
        z = points[2];
    }

    Point() {
        x = 0;
        y = 0;
        z = 0;
    }

    std::string ToString();
};

