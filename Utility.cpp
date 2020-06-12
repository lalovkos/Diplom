#include "Utility.h"
#include <iostream>

#define precision 0.0001

double side(const Point A, const Point B) {
    return sqrt((B.x - A.x) * (B.x - A.x) + (B.y - A.y) * (B.y - A.y) + (B.z - A.z) * (B.z - A.z));
}

double TriangleSquare(const Point A, const Point B, const Point C) {
    double a = side(B,A);
    double b = side(C,A);
    double c = side(B,C);
    double p = (a + b + c) / 2;
    double s = sqrt(p * (p - a) * (p - b) * (p - c));
    return s;
}

//Проверка на точку внутри треугольника
//TODO: улучшить проверку на isnan, обдумать
bool PointInsideTriangle(const Point A, const Point B, const Point C, const Point point) {
    double TrSquare = TriangleSquare(A, B, C); 
    double tmpSquare = TriangleSquare(A, B, point) + TriangleSquare(A, point, C);
    if (tmpSquare > TrSquare) return false;
    tmpSquare += TriangleSquare(point, B, C);
    //isnan - вырожденный треугольник а значит точка внутри
    if (isnan(tmpSquare) || (tmpSquare - precision < TrSquare && tmpSquare + precision > TrSquare)) return true;
    else return false;
}

Point CreateVector(const Point A, const Point B) {
    return { B.x - A.x, B.y - A.y, B.z - A.z };
}

Point VectorProduct(const Point A, const Point B) {
    return { A.y * B.z - B.y * A.z, A.z * B.x - B.z * A.x, A.x * B.y - B.x * A.y };
}

double DotProduct(const Point A, const Point B) {
    return A.x * B.x + A.y * B.y + A.z * B.z;
}

void Normalize(Point& A) {
    double mlr;
    mlr = sqrt(pow(A.x, 2) + pow(A.y, 2) + pow(A.z, 2));
    A.x = A.x / mlr;
    A.y = A.y / mlr;
    A.z = A.z / mlr;
}

bool PlaneIntersectLine(const Point A, const Point B, const Point C, const Point X, const Point Y, Point* InterPoint) {

    Point rv, N, V, W;
    double e, d;

    bool Intersect = false;
    N = VectorProduct(CreateVector(A, B), CreateVector(A, C));
    Normalize(N);
    V = CreateVector(X, A);

    d = DotProduct(N, V);
    W = CreateVector(X, Y);

    e = DotProduct(N, W);

    if (e != 0) {
        rv.x = X.x + W.x * d / e;
        rv.y = X.y + W.y * d / e;
        rv.z = X.z + W.z * d / e;
        Intersect = true;
    }
    *InterPoint = rv;
    return Intersect;
}

