#include "Figure.h"

//—читаем опо€сывающий куб
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
    bool inside = false;

    //ѕроверка на попадание в опо€сывающий куб
    if (point.x > MaxPoint.x || point.y > MaxPoint.y || point.z > MaxPoint.z)
        if (point.x < MinPoint.x || point.y < MinPoint.y || point.z < MinPoint.z)
            return false;

    Point InterPoint;
    //ѕроверка на прохождение луча боковыми гран€ми
    for (int i = 0; i < TopPoints.size(); i++) {
        int j;
        i + 1 < TopPoints.size() ? j = i + 1 : j = 0;

        //ѕроверка на нахождение провер€емой точки в окрестности грани
        if (!PlaneIntersectLine(TopPoints[i], TopPoints[j], BottomPoints[i], point, { point.x + 1, point.y, point.z }, &InterPoint) && (InterPoint.x >= point.x) && (InterPoint.y >= point.y) && (InterPoint.z >= point.z))
            if ((InterPoint.x + 0.00000001 > TopPoints[j].x) && (InterPoint.x - 0.00000001 < TopPoints[j].x) && (InterPoint.y + 0.00000001 > TopPoints[j].y) && (InterPoint.y - 0.00000001 < TopPoints[j].y)) {
                inside = !inside;
                i++;
                i + 1 < TopPoints.size() ? j = i + 1 : j = 0;
            }
            else if ((PointInsideTriangle(TopPoints[i], TopPoints[j], BottomPoints[i], InterPoint)) || (PointInsideTriangle(TopPoints[j], BottomPoints[i], BottomPoints[j], InterPoint))) {
                inside = !inside;
            }
                   
    }


    return inside;
}

// онструктор по умолчанию
Figure::Figure(const std::vector<Point> toppoints, const std::vector<Point> bottompoints, const double Phi, const int order) {
    setPoints(toppoints, bottompoints);
    UpdateMinMax();
    setPhi(Phi);
    setOrder(order);
}
