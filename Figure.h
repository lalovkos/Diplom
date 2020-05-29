#include <iostream>
#include <math.h>
#include <vector>
#include <functional>
#include <fstream>
#include "Point.h"
#include "Utility.h"

class Figure {
private:
    std::vector<Point> TopPoints; //Верхние точки фигуры
    std::vector<Point> BottomPoints; //Нижние точки фигуры
    double Phi;                //Пористость
    Point MinPoint, MaxPoint; //Опоясывающий куб
    int Order;
    void UpdateMinMax(void);

public:

    Figure(const std::vector<Point>, const std::vector<Point>, const double, const int); //Конструктор по умолчанию
    std::vector<Point> getBottomPoints();//Возвращает точки фигуры
    std::vector<Point> getTopPoints();
    std::pair<Point, Point> getMinMax();
    double getPhi(); //Возвращает пористость
    int getOrder(); //Возвращает порядок
    void setPoints(const std::vector<Point>, const std::vector<Point>);
    void setPhi(const double);
    void setOrder(int);
    bool PointInside(const Point);

    inline bool operator < (const Figure& figure) const {
        return (Order < figure.Order);
    }

    inline bool operator > (const Figure& figure) const {
        return (Order > figure.Order);
    }

};