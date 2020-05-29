#pragma once
#include "math.h"
#include "Point.h"

inline double side(const Point, const Point);

inline double TriangleSquare(const Point, const Point, const Point);

bool PointInsideTriangle(const Point, const Point, const Point, const Point);

Point CreateVector(const Point, const Point);

Point VectorProduct(const Point, const Point);

double DotProduct(const Point, const Point);

void Normalize(Point&);

//Point PlaneIntersectLine(Point A, Point B, Point C, Point X, Point Y)
bool PlaneIntersectLine(const Point, const Point, const Point, const Point, const Point, Point*);