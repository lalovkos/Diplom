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
        std::vector<Point> TopPoints; //Верхние точки фигуры
        std::vector<Point> BottomPoints; //Нижние точки фигуры
        double Phi;                //Пористость
        Point MinPoint, MaxPoint; //Опоясывающий куб
        void UpdateMinMax();
        int Order; 
    
    public:
       
        Figure(const std::vector<Point> toppoints, const std::vector<Point> bottompoints, const double Phi, const int order); //Конструктор по умолчанию
        std::vector<Point> getTopPoints();
        std::vector<Point> getBottomPoints();//Возвращает точки фигуры
        double getPhi(); //Возвращает пористость
        int getOrder(); //Возвращает порядок
        void setPoints(const std::vector<Point> toppoints, const std::vector<Point> bottompoints);
        void setPhi(const double phi);
        void setOrder(int order);
        bool PointInside(const Point point);

        bool operator < (const Figure& figure) const {
            return (Order < figure.Order);
        }

        bool operator > (const Figure& figure) const {
            return (Order > figure.Order);
        }
};

//Конструктор по умолчанию
Figure::Figure(const std::vector<Point> toppoints, const std::vector<Point> bottompoints, const double Phi,const int order){
    setPoints(toppoints, bottompoints);
    UpdateMinMax();
    setPhi(Phi);
    setOrder(order);
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

std::vector<Point> Figure::getTopPoints() {
    return this->TopPoints;
}

std::vector<Point> Figure::getBottomPoints() {
    return this->BottomPoints;
}

double Figure::getPhi() {
    return this->Phi;
}

int Figure::getOrder() {
    return this->Order;
}

//Выкидывает ошибку при задании неверного количества точек
void Figure::setPoints(const std::vector<Point> toppoints, const std::vector<Point> bottompoints) {
    if (toppoints.size() != bottompoints.size()) throw 22;
    this->BottomPoints = toppoints;
    this->TopPoints = bottompoints;
    //Добавить проверку на скрещивающиеся прямые
}

void Figure::setPhi(const double phi) {
    if (Phi > 0) //Проверям на то, чтобы 
        this->Phi = Phi;
    else throw 2;
}

void Figure::setOrder(const int order) {
    this->Order = order;
}

bool Figure::PointInside(const Point point) {
    bool inside = false;
    int j = this->TopPoints.size() - 1;

    //Проверка на попадание в опоясывающий куб
    if (point.x > MaxPoint.x || point.y > MaxPoint.y || point.z > MaxPoint.z)
        if (point.x < MinPoint.x || point.y < MinPoint.y || point.z < MinPoint.z)
            return false;

    Point InterPoint;
    //Проверка на прохождение луча боковыми гранями
    for (int i = 0; i < TopPoints.size(); i++) {
        int j;
        i + 1 < TopPoints.size() ? j = i + 1 : j = 0;
        if (PlaneIntersectLine(TopPoints[i], TopPoints[j], BottomPoints[i], point, { point.x + 1000, point.y, point.z}, &InterPoint))
            if(InterPoint.x > point.x)
                inside = !inside;
    }

    //Проверка на прохождение луча верхней грани
    if (PlaneIntersectLine (TopPoints[0], TopPoints[1], TopPoints[2], point, { point.x + 1000, point.y, point.z }, &InterPoint))
        if (InterPoint.x > point.x)
            inside = !inside;

    //Проверка на прохождение луча нижней грани
    if (PlaneIntersectLine(BottomPoints[0], BottomPoints[1], BottomPoints[2], point, { point.x + 1000, point.y, point.z }, &InterPoint))
        if (InterPoint.x > point.x)
            inside = !inside;

    return inside;
}

Point CreateVector(Point A, Point B) {
    Point m;
    m.x = B.x - A.x;
    m.y = B.y - A.y;
    m.z = B.z - A.z;
    return m;
}

Point VectorProduct(Point A, Point B) {
    Point VP;
    VP.x = A.y * B.z - B.y * A.z;
    VP.y = A.z * B.x - B.z * A.x;
    VP.z = A.x * B.y - B.x * A.y;
    return VP;
}

double DotProduct(Point A, Point B) {
    double vsp;
    vsp = A.x * B.x + A.y * B.y + A.z * B.z;
    return vsp;
}

void Normalize(Point& A) {
    double R, mlr;
    mlr = sqrt(pow(A.x, 2) + pow(A.y, 2) + pow(A.z, 2));
    A.x = A.x / mlr;
    A.y = A.y / mlr;
    A.z = A.z / mlr;
}

bool PlaneIntersectLine(Point A, Point B, Point C, Point X, Point Y, Point *InterPoint) {
    
    Point rv, N, V, W;
    double e, d;

    bool NotIntersect = true;
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
        NotIntersect = false;
    }
    return NotIntersect;
    *InterPoint = rv;
}
