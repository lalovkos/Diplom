#include <iomanip>
#include "Figure.h"
#include <clocale>  
#include "Timer.h"
#include "Point.h"
#include "Utility.h"
#include <time.h>

#define BreakKoeff 10
#define GridInputFileName "Grid.txt"
#define FiguresInputFileName "Points.txt"
#define LogFileName "Log.txt"
#define PorosityFileName "Porosity.txt"
#define ResultsFileName "Stats.txt"

using namespace std;

std::ofstream LogFile; //Глобальный файл логов

struct Grid {
    std::vector<Point> Setka;
    std::vector<std::vector<int>> Elements;
    double Phi;
};

struct VolumeSt {
    std::vector<double> interval;
    double Volume;
    double Phi;
};

//Ввод сетки
Grid InputGrid(const char* filename) {
    Grid tmpGrid;
    Point tmppoint;
    std::vector<int> tmpvector;
    ifstream gridfile(filename);
    if (gridfile.is_open()) {
        int  nodecount;
        if (!(gridfile >> nodecount)) { 
            LogFile << "Mistake in" << filename << " - Nodecount is not a number" << std::endl; 
            throw 1;  
        }
        for (int i = 0; i < nodecount; i++) {
            if (!(gridfile >> tmppoint.x) || !(gridfile >> tmppoint.y) || !(gridfile >> tmppoint.z)) { 
                LogFile << "Mistake in" << filename << " - Node points is not a number" << std::endl;
                throw 1; 
            }
            else tmpGrid.Setka.push_back(tmppoint);
        }
        int elementscount;
        int pointsincurrent;
        double node;
        if (!(gridfile >> elementscount) || (elementscount <= 0)) {
            LogFile << "Mistake in" << filename << " - Elements count is not a number or has unexpected value" << std::endl;
            throw 1;
        }
        for (int i = 0; i < elementscount; i++) {
            if (!(gridfile >> pointsincurrent) || (pointsincurrent <= 0)) {
                LogFile << "Mistake in" << filename << " - Currents elements count is not a number or has unexpected value" << std::endl;
                throw 1;
            }
            for (int j = 0; j < pointsincurrent; j++) {
                if (!(gridfile >> node)) {
                    LogFile << "Mistake in" << filename << " - Elements point is not a number" << std::endl;
                    throw 1;
                }
                tmpvector.push_back(node);
            }
            tmpGrid.Elements.push_back(tmpvector);
            tmpvector.clear();
        }
        if (!(gridfile >> tmpGrid.Phi) || (tmpGrid.Phi < 0) || (tmpGrid.Phi > 1)) {
            LogFile << "Mistake in" << filename << " - Phi is not a number or has unexpected Value" << std::endl;
            throw 1;
        }
    }
    else {
        LogFile << "Cannot open " << filename << std::endl;
        throw 1;
    }
    gridfile.close();
    return tmpGrid;
}

//Ввод фигур
std::vector<Figure> InputFigures(const char* filename) {
    ifstream figuresfile(filename);
   
    std::vector<Figure> tmpvector;
    std::vector<Point> topfigurepoints, bottomfigurepoints;

    if (figuresfile.is_open()) {
        while (!figuresfile.eof()) {
            int numberofpoints;
            double Phi;
            int Order;
            if ((figuresfile >> numberofpoints))
                if (numberofpoints > 0) {
                for (int i = 0; i < numberofpoints; i++) {
                    Point tmppoint;
                    if (!(figuresfile >> tmppoint.x) || !(figuresfile >> tmppoint.y) || !(figuresfile >> tmppoint.z)) {
                        LogFile << "Mistake in " << filename << " - Figures point is not a number" << std::endl;
                        throw 2;
                    }
                    topfigurepoints.push_back(tmppoint);
                }
                for (int i = 0; i < numberofpoints; i++) {
                    Point tmppoint;
                    if (!(figuresfile >> tmppoint.x) || !(figuresfile >> tmppoint.y) || !(figuresfile >> tmppoint.z)) {
                        LogFile << "Mistake in " << filename << " - Figures point is not a number" << std::endl;
                        throw 2;
                    }
                    bottomfigurepoints.push_back(tmppoint);
                }
                if (!(figuresfile >> Phi)) LogFile << "Mistake in " << filename << " - Phi is not a number" << std::endl;
                if (!(figuresfile >> Order)) LogFile << "Mistake in " << filename << " - Order is not a number" << std::endl;

                tmpvector.push_back(Figure(topfigurepoints, bottomfigurepoints, Phi, Order));
                topfigurepoints.clear();
                bottomfigurepoints.clear();
            }
            else {
                LogFile << "Mistake in " << filename << " - Figures number is not a number or has unexpected value" << std::endl;
                throw 2;
            }
            
        }
    }
    else {
        LogFile << "Cannot open" << filename << std::endl;

    }
    figuresfile.close();
    return tmpvector;
}

//Не работает
bool PointInsideVector(std::vector<Point> vect, Point point) {
    for (int i = 0; i < (int)vect.size() - 1; i++)
        if (PointInsideTriangle(vect[i], vect[i + 1], vect[(int)vect.size() - 1], point)) return true;
    return false;
}

//Основная программа
int main() {
    try{ 
        Timer TimePassed;
        Grid GeneralGrid = InputGrid(GridInputFileName);
        LogFile.open(LogFileName);
        std::vector<Figure> Figures = InputFigures(FiguresInputFileName);
        std::sort(Figures.begin(),Figures.end(),greater<Figure>()); //Сортировка по приоритету(по убыванию)

        std::ofstream StatsFile;
        std::ofstream PorosityFile;
        StatsFile.open(ResultsFileName);
        PorosityFile.open(PorosityFileName);

        std::vector<VolumeSt> Volume;

        //Глобальный Phi будет в 0 элементе
        Volume[0].Phi = GeneralGrid.Phi;

        //Количество разных площадей потенциально в каждом КЭ с их показателями Phi
        for (int i = 1; i < Figures.size(); i++)
            Volume[i].Phi = Figures[i].getPhi();

        //Главный цикл по КЭ
        for (std::vector<int> CurrentElement : GeneralGrid.Elements) {
            std::vector<Point> ElementPolygon;

            //Формируем полигон КЭ, а заодно считаем его опоясывающий куб
            Point MinPoint, MaxPoint; //Минимальная, максимальная точки элемента
            for (int i : CurrentElement) {
                i--;
                ElementPolygon.push_back(GeneralGrid.Setka[i]);
                if (GeneralGrid.Setka[i].x < MinPoint.x) MinPoint.x = GeneralGrid.Setka[i].x;
                else if (GeneralGrid.Setka[i].x > MaxPoint.x) MaxPoint.x = GeneralGrid.Setka[i].x;

                if (GeneralGrid.Setka[i].y < MinPoint.y) MinPoint.y = GeneralGrid.Setka[i].y;
                else if (GeneralGrid.Setka[i].y > MaxPoint.y) MaxPoint.y = GeneralGrid.Setka[i].y;

                if (GeneralGrid.Setka[i].z < MinPoint.z) MinPoint.z = GeneralGrid.Setka[i].z;
                else if (GeneralGrid.Setka[i].z > MaxPoint.z) MaxPoint.z = GeneralGrid.Setka[i].z;
            }

            //Обнуляем суммарные площади фигур на КЭ
            for (int i = 0; i < Figures.size() + 1; i++)
                Volume[i].Volume = 0;

            Point Step = { ((MaxPoint.x - MinPoint.x) / BreakKoeff),((MaxPoint.y - MinPoint.y) / BreakKoeff),((MaxPoint.z - MinPoint.z) / BreakKoeff) }; //Размер шага на КЭ
            double StepSquare = Step.y * Step.z;    
            int halfElementsize = (int)(ElementPolygon.size() / 2);

            //Двойная проверка?

            //Разбиваем на "столбы" и работаем с ними
            for (int i = 0; i < BreakKoeff; i++) {
                for (int j = 0; j < BreakKoeff; j++) {
                   
                    //Формируем рабочую прямую параллельную оси x
                    Point StartLine(MinPoint.x, MinPoint.y + Step.y*(i + 0.5), MinPoint.z + Step.z*(j + 0.5));
                    Point EndLine(MinPoint.x + 1, MinPoint.y + Step.y * (i + 0.5), MinPoint.z + Step.z * (j + 0.5));
                    //Точка пересечения
                    Point InterPoint;
                    for (int tfn = 0; tfn < halfElementsize; tfn++) { //TFN - TopFirstNumber //TSN - TopSecondNumber
                        int tsn;
                        tfn + 1 < halfElementsize ? tsn = tfn + 1 : tsn = 0;
                        //Проверяем пересекает ли плоскость грани наш "столб"
                        if (!PlaneIntersectLine(ElementPolygon[tfn], ElementPolygon[tsn], ElementPolygon[tfn + halfElementsize], StartLine, EndLine, &InterPoint)) {
                            //Проверяем, находится ли точка в пределах элемента
                            if (InterPoint.x >= MinPoint.x && InterPoint.x <= MaxPoint.x)
                                //Проверяем, находится ли точка в грани
                                //TODO: Переделать методом триангуляции, перенести в PointInsidePolygon
                                if (PointInsideTriangle(ElementPolygon[tfn], ElementPolygon[tsn], ElementPolygon[tfn+halfElementsize],InterPoint) || (PointInsideTriangle(ElementPolygon[tsn], ElementPolygon[tfn+halfElementsize], ElementPolygon[tsn + halfElementsize],InterPoint))) {
                                    //cout << InterPoint.ToString() << endl;
                                    //cout << ElementPolygon[tfn].ToString() << ElementPolygon[tsn].ToString() << ElementPolygon[tfn + 4].ToString() << ElementPolygon[tsn + 4].ToString() << endl
                                    Volume[0].interval.push_back(InterPoint.x);
                                }
                        }

                    }
                    //Получены главные интервалы внутри элемента, в которых будет считаться обьем
                    std::sort(Volume[0].interval.begin(), Volume[0].interval.end());
                    //Проверяем количество интервалов на четность, если нечетное, добавляем максимальную границу
                    if (Volume[0].interval.size() % 2 != 0) Volume[0].interval.push_back(MaxPoint.x);

                    std::vector<double> interval;
                    int fignum = 0;
                    //Обход фигур по порядку приоритета
                    for (Figure CurFig : Figures) {
                        fignum++;
                        int fnp = (int)CurFig.getTopPoints().size(); //fnp - Figure number of points
                        for (int fn = 0; fn < fnp; fn++) { //FN - FirstNumber //SN - SecondNumber
                            int sn;
                            fn + 1 < fnp ? sn = fn + 1 : sn = 0;
                            //Проверяем пересекает ли плоскость граней фигуры наш "столб"
                            if (!PlaneIntersectLine(CurFig.getTopPoints()[fn], CurFig.getTopPoints()[sn], CurFig.getBottomPoints()[fn],StartLine,EndLine, &InterPoint)) {
                                //Проверяем, находится ли точка в пределах фигуры
                                if (InterPoint.x >= CurFig.getMinMax().first.x && InterPoint.x <= CurFig.getMinMax().second.x)
                                    //Проверяем, находится ли точка в грани
                                    //TODO: Переделать методом триангуляции, перенести в PointInsidePolygon
                                    if (PointInsideTriangle(CurFig.getTopPoints()[fn], CurFig.getTopPoints()[sn], CurFig.getBottomPoints()[fn], InterPoint) || (PointInsideTriangle(CurFig.getTopPoints()[sn], CurFig.getBottomPoints()[fn], CurFig.getBottomPoints()[sn], InterPoint))) {
                                        //Если точка пересечения попадает в главные интервалы добавляем ее в интервал текущей фигуры
                                        for (int num = 0; num = Volume[0].interval.size(); num += 2) {
                                            if (InterPoint.x > Volume[0].interval[num] && InterPoint.x < Volume[0].interval[num + 1])
                                                Volume[fignum].interval.push_back(InterPoint.x);
                                        }
                                    }
                            }

                        }
                        std::sort(Volume[fignum].interval.begin(), Volume[fignum].interval.end());
                        if (Volume[fignum].interval.size() % 2 != 0) Volume[fignum].interval.push_back(MaxPoint.x);
                    }

                    //Добавляем суммарный обьем столба в соответствующем векторе Volume
                    //По каждой фигуре
                    for (int i = 1; i < Volume.size(); i++) {
                        //По интервалам фигуры
                        for (int j = 0; j < Volume[i].interval.size(); j += 2) {
                            //По главному интервалу
                            for (int mn = 0; mn < Volume[0].interval.size(); mn += 2) {
                                //Интервал фигуры внутри
                                if ((Volume[0].interval[mn] <= Volume[i].interval[j]) && (Volume[0].interval[mn + 1] >= Volume[i].interval[j + 1])) {
                                    Volume[i].Volume += StepSquare * (Volume[i].interval[j + 1] - Volume[i].interval[j]);
                                    Volume[0].interval.insert(Volume[0].interval.end() + (mn + 1), { Volume[i].interval[j], Volume[i].interval[j + 1] });
                                    break;
                                }
                                //Интервал имеет границу слева
                                else {
                                    if ((Volume[0].interval[mn] >= Volume[i].interval[j]) && (Volume[0].interval[mn + 1] <= Volume[i].interval[j + 1])) {
                                        Volume[i].Volume += StepSquare * (Volume[i].interval[j + 1] - Volume[0].interval[mn]);
                                        Volume[0].interval[mn] = Volume[i].interval[j + 1];
                                    }
                                    //Интервал имеет границу справа
                                    else {
                                        if ((Volume[0].interval[mn] <= Volume[i].interval[j]) && (Volume[0].interval[mn + 1] <= Volume[i].interval[j + 1])) {
                                            Volume[i].Volume += StepSquare * (Volume[0].interval[mn + 1] - Volume[i].interval[j]);
                                            Volume[0].interval[mn + 1] = Volume[i].interval[j];
                                        }
                                    }
                                }
                            }
                        }
                    }
                    for (int mn = 0; mn < Volume[0].interval.size(); mn += 2)
                        Volume[0].Volume += (Volume[0].interval[mn + 1] - Volume[0].interval[mn]) * StepSquare;
                    for (int i = 0; i < Volume.size(); i++) {
                        std::cout << Volume[i].Volume;
                        Volume[i].interval.clear();
                    }
                    
                }
            }
            //Считаем получившуюся суммарную пористость и площадь КЭ
            double porosity = 0, calcVM = 0;
            for (VolumeSt VolumePart : Volume) {
                porosity += VolumePart.Volume * VolumePart.Phi;
                calcVM += VolumePart.Volume;

            }

            PorosityFile << "Average Phi = " << porosity << endl;
            int i = 0;
            for (VolumeSt VolumePart : Volume) {
                StatsFile << "Square figure " << i << "=" << VolumePart.Volume << endl;
                i++;
            }

            StatsFile << "CalculatedVolume = " << calcVM << endl;
            //StatsFile << "Square Difference = " << 100 - (calcVM / (PolygonSquare(ElementPolygon) / 100)) << endl;
            StatsFile << "SplitKoef = " << BreakKoeff << endl;
            StatsFile << "TimePassed = " << TimePassed.getTime() << endl;
            StatsFile << "--------------------------------------------------------" << endl;

            ElementPolygon.clear();
        }
    }
    catch (int ErrorCode) {
        switch (ErrorCode) {
        case 1: {
            std::cout << "Error in element parameters" << std::endl;
            break;
        }

        case 2: {
            std::cout << "Error in figure parameters" << std::endl;
            break;
        }

        case 404: {
            std::cout << "Can't open file" << std::endl;
            break;
        }


        default: {std::cout << std::endl << "+++++++++++++++++++++++++++++" << "Unknown exeption = " << ErrorCode << std::endl; break; }
        }
        LogFile << ErrorCode << " Exception thrown" << std::endl;
        return ErrorCode;
    }
    return 0;
}


//Проверка на принадлежность точки нашему многоугольнику
//bool PointInsidePolygon(const std::vector<std::pair<double, double>> polygon, std::pair<double, double> point) {
//    bool result = false;
//    int j = polygon.size() - 1;
//    for (int i = 0; i < polygon.size(); i++) {
//        if ((polygon[i].y < point.y && polygon[j].y >= point.y || polygon[j].y < point.y && polygon[i].y >= point.y) &&
//            (polygon[i].x + (point.y - polygon[i].y) / (polygon[j].y - polygon[i].y) * (polygon[j].x - polygon[i].x) < point.x))
//            result = !result;
//        j = i;
//    }
//    return result;
//}

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

