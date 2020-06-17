#include "Figure.h"

//Считаем опоясывающий куб
void Figure::UpdateMinMax() {
    this->MinPoint.x = 10000000;
    this->MinPoint.y = 10000000;
    this->MinPoint.z = 10000000;
    this->MaxPoint.x = -10000000;
    this->MaxPoint.y = -10000000;
    this->MaxPoint.z = -10000000;

    for (Point tmp : this->TopPoints) {
        if (tmp.x < this->MinPoint.x) this->MinPoint.x = tmp.x;
        else if (tmp.x > this->MaxPoint.x) this->MaxPoint.x = tmp.x;

        if (tmp.y < this->MinPoint.y) this->MinPoint.y = tmp.y;
        else if (tmp.y > this->MaxPoint.y) this->MaxPoint.y = tmp.y;

        if (tmp.z < this->MinPoint.z) this->MinPoint.z = tmp.z;
        else if (tmp.z > this->MaxPoint.z) this->MaxPoint.z = tmp.z;
    }

    for (Point tmp : this->BottomPoints) {
        if (tmp.x < this->MinPoint.x) this->MinPoint.x = tmp.x;
        else if (tmp.x > this->MaxPoint.x) this->MaxPoint.x = tmp.x;

        if (tmp.y < this->MinPoint.y) this->MinPoint.y = tmp.y;
        else if (tmp.y > this->MaxPoint.y) this->MaxPoint.y = tmp.y;

        if (tmp.z < this->MinPoint.z) this->MinPoint.z = tmp.z;
        else if (tmp.z > this->MaxPoint.z) this->MaxPoint.z = tmp.z;
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

std::pair<Point, Point> Figure::getMinMax() {
    return { this->MinPoint, this->MaxPoint };
}

//Выкидывает ошибку при задании неверного количества точек
void Figure::setPoints(std::vector<Point> toppoints, std::vector<Point> bottompoints) {
    if (toppoints.size() != bottompoints.size()) throw 22;
    this->BottomPoints = toppoints;
    this->TopPoints = bottompoints;
    UpdateMinMax();
    //Добавить проверку на скрещивающиеся прямые
}

void Figure::setPhi(const double phi) {
    if (phi > 0 && phi < 1) //Проверям на то, чтобы была правильное значение пористости
        this->Phi = phi;
    else throw 2;
}

void Figure::setOrder(const int order) {
    this->Order = order;
}

bool Figure::PointInside(const Point point) {
    bool inside = false;

    //Проверка на попадание в опоясывающий куб
    if (point.x > MaxPoint.x || point.y > MaxPoint.y || point.z > MaxPoint.z || point.x < MinPoint.x || point.y < MinPoint.y || point.z < MinPoint.z)
            return false;

    Point InterPoint;
    //Алгоритм для прямой призмы
    int j = TopPoints.size() - 1;
    //Проверка на прохождение луча  
    for (int i = 0; i < TopPoints.size(); i++) {
        if ((TopPoints[i].y < point.y && TopPoints[j].y >= point.y || TopPoints[j].y < point.y && TopPoints[i].y >= point.y) &&
            (TopPoints[i].x + (point.y - TopPoints[i].y) / (TopPoints[j].y - TopPoints[i].y) * (TopPoints[j].x - TopPoints[i].x) < point.x))
            inside = !inside;
        j = i;
    }
    return inside;
}



//Конструктор по умолчанию
Figure::Figure(const std::vector<Point> toppoints, const std::vector<Point> bottompoints, const double Phi, const int order) {
    setPoints(toppoints, bottompoints);
    UpdateMinMax();
    setPhi(Phi);
    setOrder(order);
}
