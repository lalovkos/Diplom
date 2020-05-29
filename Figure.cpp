#include "Figure.h"

//—читаем опо€сывающий куб
void Figure::UpdateMinMax() {
    this->MinPoint = this->TopPoints[0];
    this->MaxPoint = this->TopPoints[0];
    for (Point tmp : this->TopPoints) {
        if (tmp.x < this->MinPoint.x) this->MinPoint.x = tmp.x;
        else if (tmp.x > this->MinPoint.x) this->MinPoint.x = tmp.x;

        if (tmp.y < this->MinPoint.y) this->MinPoint.y = tmp.y;
        else if (tmp.y > this->MinPoint.y) this->MinPoint.y = tmp.y;

        if (tmp.z < this->MinPoint.z) this->MinPoint.z = tmp.z;
        else if (tmp.z > this->MinPoint.z) this->MinPoint.z = tmp.z;
    }

    for (Point tmp : this->BottomPoints) {
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

std::pair<Point, Point> Figure::getMinMax() {
    return { this->MinPoint, this->MaxPoint };
}

//¬ыкидывает ошибку при задании неверного количества точек
void Figure::setPoints(std::vector<Point> toppoints, std::vector<Point> bottompoints) {
    if (toppoints.size() != bottompoints.size()) throw 22;
    this->BottomPoints = toppoints;
    this->TopPoints = bottompoints;
    UpdateMinMax();
    //ƒобавить проверку на скрещивающиес€ пр€мые
}

void Figure::setPhi(const double phi) {
    if (phi > 0 && phi < 1) //ѕровер€м на то, чтобы была правильное значение пористости
        this->Phi = phi;
    else throw 2;
}

void Figure::setOrder(const int order) {
    this->Order = order;
}

bool Figure::PointInside(const Point point) {
    //bool inside = false;
    //int j = this->TopPoints.size() - 1;

    ////ѕроверка на попадание в опо€сывающий куб
    //if (point.x > MaxPoint.x || point.y > MaxPoint.y || point.z > MaxPoint.z)
    //    if (point.x < MinPoint.x || point.y < MinPoint.y || point.z < MinPoint.z)
    //        return false;

    //Point InterPoint;
    ////ѕроверка на прохождение луча боковыми гран€ми
    //for (int i = 0; i < TopPoints.size(); i++) {
    //    int j;
    //    i + 1 < TopPoints.size() ? j = i + 1 : j = 0;
    //    //ѕроверка на нахождение провер€емой точки в окрестности грани

    //    if (PlaneIntersectLine(TopPoints[i], TopPoints[j], BottomPoints[i], point, { point.x + 1000, point.y, point.z }, &InterPoint))
    //        if (InterPoint.x > point.x)
    //            inside = !inside;
    //}

    ////ѕроверка на прохождение луча верхней грани
    //if (PlaneIntersectLine(TopPoints[0], TopPoints[1], TopPoints[2], point, { point.x + 1000, point.y, point.z }, &InterPoint))
    //    if (InterPoint.x > point.x)
    //        inside = !inside;

    ////ѕроверка на прохождение луча нижней грани
    //if (PlaneIntersectLine(BottomPoints[0], BottomPoints[1], BottomPoints[2], point, { point.x + 1000, point.y, point.z }, &InterPoint))
    //    if (InterPoint.x > point.x)
    //        inside = !inside;

    //return inside;
    return false;
}

// онструктор по умолчанию
Figure::Figure(const std::vector<Point> toppoints, const std::vector<Point> bottompoints, const double Phi, const int order) {
    setPoints(toppoints, bottompoints);
    UpdateMinMax();
    setPhi(Phi);
    setOrder(order);
}