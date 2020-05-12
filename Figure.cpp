#include <iostream>
#include <math.h>
#include <vector>
#include <functional>
#include <fstream>


struct Point {
    double x;
    double y;
    double z;
};

class Figure {
    private:
        std::vector<Point> Points; //Точки фигуры
        double Phi;                //Пористость
        Point MinPoint, MaxPoint; //Опоясывающий куб
        void UpdateMinMax();
    
    public:
        //Конструктор по умолчанию
        Figure(const std::vector<Point> points, const double Phi);
        std::vector<Point> getPoints();
        double getPhi();
        void setPoints(const std::vector<Point> points);
        void setPhi(const double Phi);
        bool PointInside(const Point point);
};

//Конструктор по умолчанию
Figure::Figure(const std::vector<Point> points, const double Phi){
    setPoints(points);
    UpdateMinMax();
    setPhi(Phi);
}

//Считаем опоясывающий куб
void Figure::UpdateMinMax() {
    this->MinPoint = this->Points[0];
    this->MaxPoint = this->Points[0];
    for (Point tmp : this->Points) {
        if (tmp.x < this->MinPoint.x) this->MinPoint.x = tmp.x;
        else if (tmp.x > this->MinPoint.x) this->MinPoint.x = tmp.x;

        if (tmp.y < this->MinPoint.y) this->MinPoint.y = tmp.y;
        else if (tmp.y > this->MinPoint.y) this->MinPoint.y = tmp.y;

        if (tmp.z < this->MinPoint.z) this->MinPoint.z = tmp.z;
        else if (tmp.z > this->MinPoint.z) this->MinPoint.z = tmp.z;
    }
}

std::vector<Point> Figure::getPoints() {
    return this->Points;
}

double Figure::getPhi() {
    return this->Phi;
}

void Figure::setPoints(const std::vector<Point> points) {
    this->Points = points;
}

void Figure::setPhi(const double Phi) {
    if (Phi > 0) //Проверям на то, чтобы 
        this->Phi = Phi;
    else throw 2;
}

bool Figure::PointInside(const Point point) {
    bool inside = false;

    return inside;
}
