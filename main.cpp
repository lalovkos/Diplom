#include <iostream>
#include <fstream>
#include <conio.h>
#include <vector>
#include <functional>
#include <iomanip>
#include <math.h>
#include <clocale>  
#include "Timer.h"
#include <time.h>

#define x first
#define y second
#define BreakKoeff 1000
#define GlobalPhi 0.4
#define H 10

using namespace std;

struct Grid {
    std::vector<std::pair<double, double>> Setka;
    std::vector<std::vector<int>> Elements;
    double Phi = 0.2;
};

struct Figure {
    std::vector<std::pair<double, double>> points;
    int depth = 0;
    double Phi = 0.1;
    std::pair<double, double> MinPoint = { 0,0 }, MaxPoint = {1,1};
};

struct SquareSt {
    double Square;
    double Phi;
};

//Открытие файла для считывания и его проверки
std::ofstream* OpenFileForWriting(const char* filename) {
    ofstream* stream = new ofstream(filename);
    //Тут будут проверки состояния файла
    return stream;
}

//Считаем площадь произвольного многоугольника
double PolygonSquare(const std::vector<std::pair<double, double>> polygon) {
    double r = 0;
    for (int i = 0; i < polygon.size(); i++)
        r += (polygon[(i + 1) % polygon.size()].x - polygon[i].x) * (polygon[i].y + polygon[(i + 1) % polygon.size()].y);
    return abs(r / 2);
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

//Ввод сетки
Grid InputGrid(const char* filename) {
    Grid tmpGrid;
    std::pair<double, double> tmppoint;
    std::vector<int> tmpvector;
    ifstream gridfile(filename);
    if (gridfile.is_open()) {
        int  nodecount;
        gridfile >> nodecount;
        for (int i = 0; i < nodecount; i++) {
            gridfile >> tmppoint.x;
            gridfile >> tmppoint.y;
            tmpGrid.Setka.push_back(tmppoint);
        }
        int elementscount;
        int pointsincurrent;
        double node;
        gridfile >> elementscount;
        for (int i = 0; i < elementscount; i++) {
            gridfile >> pointsincurrent;
            for (int j = 0; j < pointsincurrent; j++) {
                gridfile >> node;
                tmpvector.push_back(node);
            }
            tmpGrid.Elements.push_back(tmpvector);
            tmpvector.clear();
        }
        gridfile >> tmpGrid.Phi;
    }
    gridfile.close();
    return tmpGrid;
}

//Ввод фигур
std::vector<Figure> InputFigures(const char* filename) {
    ifstream figuresfile(filename);
    char more = '#';
    std::vector<Figure> tmpvector;
    Figure CurrentFigure;
    int numberofpoints;
    if (figuresfile.is_open()) {
        figuresfile >> more;
        while (more == '#') {
            figuresfile >> numberofpoints;
            if (numberofpoints != 0) {
                for (int i = 0; i < numberofpoints; i++) {
                    std::pair<double, double> tmppair;
                    figuresfile >> tmppair.x;
                    figuresfile >> tmppair.y;
                    CurrentFigure.points.push_back(tmppair);

                    if (tmppair.x < CurrentFigure.MinPoint.x) CurrentFigure.MinPoint.x = tmppair.x;
                    else if (tmppair.x > CurrentFigure.MaxPoint.x) CurrentFigure.MaxPoint.x = tmppair.x;

                    if (tmppair.y < CurrentFigure.MinPoint.y) CurrentFigure.MinPoint.y = tmppair.y;
                    else if (tmppair.y > CurrentFigure.MaxPoint.y) CurrentFigure.MaxPoint.y = tmppair.y;
                }
                figuresfile >> CurrentFigure.Phi;
                figuresfile >> CurrentFigure.depth;
            }
            tmpvector.push_back(CurrentFigure);
            CurrentFigure.points.clear();
            figuresfile >> more;
        }
    }
    figuresfile.close();
    return tmpvector;
}


//Основная программа
int main() {
    ofstream* LogFile = OpenFileForWriting("Log.txt");
    try{ 
        Timer TimePassed;

        Grid GeneralGrid = InputGrid("Grid.txt");
        std::vector<Figure> Figures = InputFigures("Points.txt");
        ofstream* StatsFile = OpenFileForWriting("Stats.txt");
        ofstream* PorosityFile = OpenFileForWriting("Porosity.txt");
        ofstream* SquareFile = OpenFileForWriting("Square.txt");
        ofstream* SchemeFile = OpenFileForWriting("GridScheme.txt");

        vector<SquareSt> Square; //Вычисляемая площадь

        //Глобальный Phi
        Square.push_back({ 0,GeneralGrid.Phi });
        //Количество разных площадей потенциально в каждом КЭ с их показателями Phi
        for (int i = 0; i < Figures.size(); i++)
            Square.push_back({ 0,Figures[i].Phi });

        //Главный цикл по КЭ
        for (std::vector<int> CurrentElement : GeneralGrid.Elements) {
            std::vector<std::pair<double, double>> ElementPolygon;
            ElementPolygon.clear();

            //Формируем полигон КЭ, а заодно считаем его опоясывающий прямогоугольник
            std::pair<double, double> MinPoint = { 100000000,100000000 }, MaxPoint = { 0.000000001,0.000000001 }; //Минимальная, максимальная точки элемента
            for (int i : CurrentElement) {
                ElementPolygon.push_back(GeneralGrid.Setka[i]);
                if (GeneralGrid.Setka[i].x < MinPoint.x) MinPoint.x = GeneralGrid.Setka[i].x;
                else if (GeneralGrid.Setka[i].x > MaxPoint.x) MaxPoint.x = GeneralGrid.Setka[i].x;

                if (GeneralGrid.Setka[i].y < MinPoint.y) MinPoint.y = GeneralGrid.Setka[i].y;
                else if (GeneralGrid.Setka[i].y > MaxPoint.y) MaxPoint.y = GeneralGrid.Setka[i].y;
            }


            //Обнуляем суммарные площади фигур на КЭ
            for (int i = 0; i < Figures.size() + 1; i++)
                Square[i].Square = 0;

            std::pair<double, double> CurrentCell; //Текущие координаты ячейки
            std::pair<double, double> Step = { ((MaxPoint.x - MinPoint.x) / BreakKoeff),((MaxPoint.y - MinPoint.y) / BreakKoeff) }; //Размер шага на КЭ

            //Разбиваем на ячейки и опрашиваем фигуры(их необходимо отсортировать)
            for (int i = 0; i < BreakKoeff; i++) {
                for (int j = 0; j < BreakKoeff; j++) {
                    int squarenumber = 0;
                    //Переделать
                    if (PointInsidePolygon(ElementPolygon, { MinPoint.x + (i + 0.5) * Step.x,MinPoint.y + (j + 0.5) * Step.y }))
                    {
                        for (int figurenumber = 0; figurenumber < Figures.size(); figurenumber++) {

                            if (PointInsidePolygon(Figures[figurenumber].points, { MinPoint.x + (i + 0.5) * Step.x,MinPoint.y + (j + 0.5) * Step.y })) {
                                squarenumber = figurenumber + 1;
                                break;
                            }

                        }
                        Square[squarenumber].Square += Step.x * Step.y;
                    }

                }
            }

            //Считаем получившуюся суммарную пористость и площадь КЭ
            double porosity = 0, calcsq = 0;
            for (SquareSt SquarePart : Square) {
                porosity += SquarePart.Square * SquarePart.Phi;
                calcsq += SquarePart.Square;

            }

            porosity /= PolygonSquare(ElementPolygon);

            *PorosityFile << "Average Phi = " << porosity << endl;
            int i = 0;
            for (SquareSt SquarePart : Square) {
                *StatsFile << "Square figure " << i << "=" << SquarePart.Square << endl;
                i++;
            }
            *StatsFile << "CalculatedSquare = " << calcsq << endl;
            *StatsFile << "CalculatedVolume = " << calcsq * H << endl;
            *StatsFile << "RealSquare = " << PolygonSquare(ElementPolygon) << endl;
            *StatsFile << "Square Difference = " << 100 - (calcsq / (PolygonSquare(ElementPolygon) / 100)) << endl;
            *StatsFile << "SplitKoef = " << BreakKoeff << endl;
            *StatsFile << "TimePassed = " << TimePassed.getTime() << endl;
            *StatsFile << "--------------------------------------------------------" << endl;

        }
    }
    catch (int ErrorCode) {
        switch (ErrorCode) {
        case 1: {
            std::cout << "Error in parameters" << std::endl;
            break; 
        }

        case 2: {
            std::cout << "Error in figure parameters" << std::endl;
            break; 
        }

        case 3: {
            std::cout << "File error" << std::endl;
            break;
        }

        default: {std::cout << std::endl << "+++++++++++++++++++++++++++++" << "Unknown exeption = " << ErrorCode << std::endl; break; }
        }
        *LogFile << ErrorCode << "Exception thrown" << std::endl;
        return ErrorCode;
    }
    return 0;
}

//Определяет расстояние между точками
double Length(const std::pair<double, double> FirstPoint, const std::pair<double, double> SecondPoint) {
    return sqrt((SecondPoint.x - FirstPoint.x) * (SecondPoint.x - FirstPoint.x) + (SecondPoint.y - FirstPoint.y) * (SecondPoint.y - FirstPoint.y));
}

//Метод разбивания на ячейки
double CalculateBreakingSquare(const std::vector<std::pair<double, double>> Element, const std::vector<std::pair<double, double>> Figure, const int BreakKoef) {
    double Square = 0; //Вычисляемая площадь    
    std::pair<double, double> CurrentCell; //Текущие координаты ячейки
    std::pair<double, double> MinPoint = { 0,0 }, MaxPoint = { 1,1 }; //Минимальная, максимальная точки элемента
    for (std::pair<double, double> CurrentPoint : Element) {
        if (CurrentPoint.x < MinPoint.x) MinPoint.x = CurrentPoint.x;
        else if (CurrentPoint.x > MaxPoint.x) MaxPoint.x = CurrentPoint.x;

        if (CurrentPoint.y < MinPoint.y) MinPoint.y = CurrentPoint.y;
        else if (CurrentPoint.y > MaxPoint.y) MaxPoint.y = CurrentPoint.y;
    }

    std::pair<double, double> Step = { (MaxPoint.x - MinPoint.x / BreakKoef),(MaxPoint.y - MinPoint.y / BreakKoef) };
    for (int i = 0; i < BreakKoef; i++) {
        for (int j = 0; j < BreakKoef; j++) {
            //0.5 потому что берется центр ячейки для опроса
            if (PointInsidePolygon(Figure, { MinPoint.x + (i + 0.5) * Step.x,MinPoint.y + (j + 0.5) * Step.y }))
                Square += Step.x * Step.y;
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

//Возможно лучше векторное произведение
//Вычисление площади треугольника
double CalculateTriangleSquare(const std::pair<double, double> FirstPoint, const std::pair<double, double> SecondPoint, const std::pair<double, double> ThirdPoint) {
    double p = (1 / 2) * Length(FirstPoint, SecondPoint) + Length(FirstPoint, ThirdPoint) + Length(SecondPoint, ThirdPoint);
    return sqrt(p * (p - Length(FirstPoint, SecondPoint)) * (p - Length(FirstPoint, ThirdPoint)) * (p - Length(SecondPoint, ThirdPoint)));
}


//if (node[i][j] && node[i + 1][j + 1] && node[i][j + 1] && node[i + 1][j]) //Если все узлы внутри - площадь вся находится в многоугольнике
  //    Grid[i][j].Square = GridStep.x * GridStep.y;
  //else if (!node[i][j] && !node[i + 1][j + 1] && !node[i][j + 1] && !node[i + 1][j]) //Если все узлы снаружи - КЭ пустой
  //    Grid[i][j].Square = 0;
  //else {
      //Сюда попадаем если у нас часть узлов внутри, а часть снаружи
      //Последовательно обходим все границы КЭ проверяя точки пересечения и записывая ее в вектор
//            Grid[i][j].Square = CalculateBreakingSquare(points, 100, 100, GridStep, {min.x + i * GridStep.x, min.y + j * GridStep.y});
//        /*}*/
//        workingpolygon.clear();
//    }
//}

//std::cout << "_________Points per element_________" << endl;
//for (int i = 0; i < GRIDNUMBER_X+1; i++) {
//    for (int j = 0; j < GRIDNUMBER_Y+1; j++) {
//        std::cout << Grid[i][j].PointsNumbers.size() << " ";
//    }
//    std::cout << endl;
//}

//std::cout << "_________Inner Outer Nodes_________" << endl;
//for (int i = 0; i < GRIDNUMBER_X + 1; i++) {
//    for (int j = 0; j < GRIDNUMBER_Y + 1; j++) {
//        std::cout << node[i][j] << "-";
//    }
//    std::cout << endl;
//}

//double GeneralSquare = 0;

//std::cout << "_________Square_________" << endl;
//for (int i = 0; i < GRIDNUMBER_X + 1; i++) {
//    for (int j = 0; j < GRIDNUMBER_Y + 1; j++) {
//        std::cout << std::setw(9) << Grid[i][j].Square;
//        GeneralSquare += Grid[i][j].Square;
//    }
//    std::cout << endl;
//}

//std::cout << "Calculated Square = " << GeneralSquare << "<->";
//GeneralSquare = PolygonSquare(points);
//std::cout << "Real Square = " << GeneralSquare << endl;

//std::cout << "_________Average Phi_________" << endl;
//double porosity;
//for (int i = 0; i < GRIDNUMBER_X + 1; i++) {
//    for (int j = 0; j < GRIDNUMBER_Y + 1; j++) {
//        porosity = (Grid[i][j].Square * 0.2 + ((GridStep.x * GridStep.y - Grid[i][j].Square) * 0.1)) / GridStep.x * GridStep.y;
//        std::cout << std::setw(9) << porosity;
//    }
//    std::cout << endl;
//}
//std::cout << "Phi in = " << 0.2 << "<->";
//std::cout << "Phi out = " << 0.1;
//

////Вычисление шага и инициализация сетки
//std::pair<double, double> GridStep;
//GridStep.x = (max.x - min.x) / GRIDNUMBER_X;
//GridStep.y = (max.y - min.y) / GRIDNUMBER_Y;

/* FinalElement Grid[GRIDNUMBER_X+1][GRIDNUMBER_Y+1];*/

 ////Обнуляем значения площади на КЭ
 //for (int i = 0; i < GRIDNUMBER_X; i++)
 //    for (int j = 0; j < GRIDNUMBER_Y; j++)
 //        Grid[i][j].Square= 0;

 ////Размечаем узлы на принадлежность многоугольнику
 //for (int i = 0; i < GRIDNUMBER_X; i++)
 //    for (int j = 0; j < GRIDNUMBER_Y; j++)
 //        node[i][j] = PointInsidePolygon(points, min.x + i*GridStep.x, min.y + j * GridStep.y);
 //

 ////Разбрасываем точки в их КЭ
 //int count = 0;
 //for (std::pair<double, double> CurPoint : points) { 
 //    int i = (int)floor((CurPoint.x - min.x) / GridStep.x);
 //    int j = (int)floor((CurPoint.y - min.y) / GridStep.y);
 //    Grid[i][j].PointsNumbers.push_back(count);
 //    count++;
 //}

