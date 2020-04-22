#include <iostream>
#include <fstream>
#include <conio.h>
#include <vector>
#include <functional>
#include <iomanip>
#include <math.h>
#include <clocale>  
#include <time.h>

#define x first
#define y second
#define GRIDNUMBER_X 10
#define GRIDNUMBER_Y 10
using namespace std;

struct FinalElement{
    double Square;
    vector<int> PointsNumbers;
};

//TODO: Убрать из глобальных
bool node[GRIDNUMBER_X + 1][GRIDNUMBER_Y + 1];

//Считаем площадь произвольного многоугольника
double PolygonSquare(const std::vector<std::pair<double, double>> polygon) {
    double r = 0;
    for (int i = 0; i < polygon.size(); i++)
        r += (polygon[(i + 1) % polygon.size()].x - polygon[i].x) * (polygon[i].y + polygon[(i + 1) % polygon.size()].y);
    return abs(r / 2);
}

//Определяет расстояние между точками
double Length(const std::pair<double, double> FirstPoint, const std::pair<double, double> SecondPoint) {
    return sqrt((SecondPoint.x - FirstPoint.x) * (SecondPoint.x - FirstPoint.x) + (SecondPoint.y - FirstPoint.y) * (SecondPoint.y - FirstPoint.y));
}

//Проверка на принадлежность точки нашему многоугольнику
bool PointInsidePolygon(const std::vector<std::pair<double, double>> polygon, std::pair<double, double> point) {
    bool result = false;
    int j = polygon.size() - 1;
    for (int i = 0; i < polygon.size(); i++) {
        if ((polygon[i].y < point.y && polygon[j].y >= point.y || polygon[j].y < point.y && polygon[i].y >= point.y) &&
            (polygon[i].x + (point.y - polygon[i].y) / (polygon[j].y - polygon[i].y) * (polygon[j].x - polygon[i].x) < point.x))
            result = !result;
        j = i;
    }
    return result;
}

//Перегрузка для проверки на принадлежность точки многоугольнику
bool PointInsidePolygon(const std::vector<std::pair<double, double>> polygon, double pointX, double pointY) {
    bool result = false;
    int j = polygon.size() - 1;
    for (int i = 0; i < polygon.size(); i++) {
        if ((polygon[i].y < pointY && polygon[j].y >= pointY || polygon[j].y < pointY && polygon[i].y >= pointY) &&
            (polygon[i].x + (pointY - polygon[i].y) / (polygon[j].y - polygon[i].y) * (polygon[j].x - polygon[i].x) < pointX))
            result = !result;
        j = i;
    }
    return result;
}

//Возможно лучше векторное произведение
//Вычисление площади треугольника
double CalculateTriangleSquare(const std::pair<double, double> FirstPoint, const std::pair<double, double> SecondPoint, const std::pair<double, double> ThirdPoint) {
    double p = (1/2)*Length(FirstPoint, SecondPoint) + Length(FirstPoint,ThirdPoint) + Length(SecondPoint, ThirdPoint);
    return sqrt(p*(p-Length(FirstPoint, SecondPoint)) * (p - Length(FirstPoint, ThirdPoint)) * (p - Length(SecondPoint, ThirdPoint)));
}

//Метод разбивания на ячейки
double CalculateBreakingSquare(const std::vector<std::pair<double, double>> polygon, const int BreakingCoefX, const int BreakingCoefY, const std::pair<double, double> GridStep, const std::pair<double, double> LeftBottomNode){
    double Square = 0; //Вычисляемая площадь    
    std::pair<double, double> CurrentCell; //Текущие координаты ячейки
    for (int i = 0; i < BreakingCoefX; i++) {
        for (int j = 0; j < BreakingCoefY; j++) {
            CurrentCell.x = LeftBottomNode.x + (i + 0.5)*(GridStep.x / BreakingCoefX);
            CurrentCell.y = LeftBottomNode.y + (j + 0.5) * (GridStep.y / BreakingCoefY); 
            if (PointInsidePolygon(polygon,CurrentCell)) {
                Square += (GridStep.x / BreakingCoefX) * (GridStep.y / BreakingCoefY);
            
            }
        }
    }
    return Square;
}

//Определение координат точек пересечения двух отрезков
//точки a и b концы первого отрезка c и d второго
std::pair<double, double> LineSegmentCrossing(std::pair<double, double> const a, std::pair<double, double> const b, std::pair<double, double> const c, std::pair<double, double> const d) { 
    std::pair<double, double> tmp;
    tmp.x = -((a.x * b.y - b.x * a.y) * (d.x - c.x) - (c.x * d.y - d.x * c.y) * (b.x - a.x)) / ((a.y - b.y) * (d.x - c.x) - (c.y - d.y) * (b.x - a.x));
    tmp.y = ((c.y - d.y) * (-tmp.x) - (c.x * d.y - d.x * c.y)) / (d.x - c.x);
    return tmp;
}

//Выбор подходящей угловой точки
std::pair<double, double> InnerCorner(const std::pair<double, double> min, const std::pair<double, double> GridStep, const int i, const int j) {
    std::pair<double, double> Corner;
    
    if (node[i][j]) {
        Corner.x = min.x + i * GridStep.x;
        Corner.y = min.y + j * GridStep.y;
    }
    else  if (node[i + 1][j]) {
        Corner.x = min.x + (i + 1) * GridStep.x;
        Corner.y = min.y + j * GridStep.y;
    }
    else  if (node[i][j + 1]) {
        Corner.x = min.x + (i + 1) * GridStep.x;
        Corner.y = min.y + j * GridStep.y;
    }
    else  if (node[i + 1][j + 1]) {
        Corner.x = min.x + (i + 1) * GridStep.x;
        Corner.y = min.y + (j + 1) * GridStep.y;
    }
    return Corner;
}

//TODO?: Реализовать сортировку точек по местоположению
//Основная программа
int main() {

    std::vector<std::pair<double, double>> points; //вектор для хранения точек

    //TODO: Перенести в отдельную функцию
    ifstream pointsfile("Points.txt");
    if (pointsfile.is_open()) {
        int numberofpoints;
        pointsfile >> numberofpoints;
        if (numberofpoints != 0) {
            for (int i = 0; i < numberofpoints; i++) {
                std::pair<double, double> tmppair;
                pointsfile >> tmppair.x;
                pointsfile >> tmppair.y;
                points.push_back(tmppair);
            }
        }
    }
    else {
        cout << "Problem with file" << endl;
    }
    if (points.size() < 3) {
        std::cout << "Not a figure" << endl;
        return 100;
    }

    //TODO: Перенести в отдельную функцию
    //Вычисление параметров сетки и среднего расстояния между точками
    pair<double, double> max, min, middle;
    pair<double,double> sum;
    min.x = min.y = 0;
    max.x = max.y = 0.1;
    middle.x = middle.y = 0.05; //Подумать

    //TODO: Заменить на краткую запись для упрощения
    //Анализ точек и подсчет основных координат(центр, крайние)
    for (std::pair<double, double> CurPoint : points) {
        if (CurPoint.x < min.x) min.x = CurPoint.x;
        else if (CurPoint.x > max.x) max.x = CurPoint.x;
        if (CurPoint.y < min.y) min.y = CurPoint.y;
        else if (CurPoint.y > max.y) max.y = CurPoint.y;
        sum.x += abs(CurPoint.x);
        sum.y += abs(CurPoint.y);
        middle.x += CurPoint.x;
        middle.y += CurPoint.y;
    }
    middle.x /= points.size();
    middle.y /= points.size();

    //Вычисление шага и инициализация сетки
    std::pair<double, double> GridStep;
    GridStep.x = (max.x - min.x) / GRIDNUMBER_X;
    GridStep.y = (max.y - min.y) / GRIDNUMBER_Y;
    
    FinalElement Grid[GRIDNUMBER_X+1][GRIDNUMBER_Y+1];

    //Обнуляем значения площади на КЭ
    for (int i = 0; i < GRIDNUMBER_X; i++)
        for (int j = 0; j < GRIDNUMBER_Y; j++)
            Grid[i][j].Square= 0;

    //Размечаем узлы на принадлежность многоугольнику
    for (int i = 0; i < GRIDNUMBER_X; i++)
        for (int j = 0; j < GRIDNUMBER_Y; j++)
            node[i][j] = PointInsidePolygon(points, min.x + i*GridStep.x, min.y + j * GridStep.y);
    

    //Разбрасываем точки в их КЭ
    int count = 0;
    for (std::pair<double, double> CurPoint : points) { 
        int i = (int)floor((CurPoint.x - min.x) / GridStep.x);
        int j = (int)floor((CurPoint.y - min.y) / GridStep.y);
        Grid[i][j].PointsNumbers.push_back(count);
        count++;
    }
    
    //TODO: Оптимизировать
    //Проходим все конечные элементы по порядку для расчетов
    //Проверить больше возможностей для некоторых исключений
    std::vector<std::pair<double, double>> workingpolygon;
    for (int i = 0; i < GRIDNUMBER_X+1; i++) {
        for (int j = 0; j < GRIDNUMBER_Y+1; j++){
            if (node[i][j] && node[i + 1][j + 1] && node[i][j + 1] && node[i + 1][j]) //Если все узлы внутри - площадь вся находится в многоугольнике
                Grid[i][j].Square = GridStep.x * GridStep.y;
            else if (!node[i][j] && !node[i + 1][j + 1] && !node[i][j + 1] && !node[i + 1][j]) //Если все узлы снаружи - КЭ пустой
                Grid[i][j].Square = 0;
            else {
                //Сюда попадаем если у нас часть узлов внутри, а часть снаружи
                //Последовательно обходим все границы КЭ проверяя точки пересечения и записывая ее в вектор
                Grid[i][j].Square = CalculateBreakingSquare(points, 100, 100, GridStep, {min.x + i * GridStep.x, min.y + j * GridStep.y});
            }
            workingpolygon.clear();
        }
    }

    std::cout << "_________Points per element_________" << endl;
    for (int i = 0; i < GRIDNUMBER_X+1; i++) {
        for (int j = 0; j < GRIDNUMBER_Y+1; j++) {
            std::cout << Grid[i][j].PointsNumbers.size() << " ";
        }
        std::cout << endl;
    }

    std::cout << "_________Inner Outer Nodes_________" << endl;
    for (int i = 0; i < GRIDNUMBER_X + 1; i++) {
        for (int j = 0; j < GRIDNUMBER_Y + 1; j++) {
            std::cout << node[i][j] << "-";
        }
        std::cout << endl;
    }

    double GeneralSquare = 0;

    std::cout << "_________Square_________" << endl;
    for (int i = 0; i < GRIDNUMBER_X + 1; i++) {
        for (int j = 0; j < GRIDNUMBER_Y + 1; j++) {
            std::cout << std::setw(9) << Grid[i][j].Square;
            GeneralSquare += Grid[i][j].Square;
        }
        std::cout << endl;
    }

    std::cout << "Calculated Square = " << GeneralSquare << "<->";
    GeneralSquare = PolygonSquare(points);
    std::cout << "Real Square = " << GeneralSquare << endl;

    std::cout << "_________Average Phi_________" << endl;
    double porosity;
    for (int i = 0; i < GRIDNUMBER_X + 1; i++) {
        for (int j = 0; j < GRIDNUMBER_Y + 1; j++) {
            porosity = (Grid[i][j].Square * 0.2 + ((GridStep.x * GridStep.y - Grid[i][j].Square) * 0.1)) / GridStep.x * GridStep.y;
            std::cout << std::setw(9) << porosity;
        }
        std::cout << endl;
    }
    std::cout << "Phi in = " << 0.2 << "<->";
    std::cout << "Phi out = " << 0.1;
    
    return 0;
}
