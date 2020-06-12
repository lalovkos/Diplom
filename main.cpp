#include <iomanip>
#include "Figure.h"
#include <clocale>  
#include "Timer.h"
#include "Point.h"
#include <algorithm>
#include "Utility.h"
#include <time.h>

#define BreakKoeff 10
#define CountFileName "C:\\Users\\METAL\\Desktop\\Diplom\\Model\\mesh\\inftry.dat"
#define ElementFileName "C:\\Users\\METAL\\Desktop\\Diplom\\Model\\mesh\\nver.dat"
#define NodesFileName "C:\\Users\\METAL\\Desktop\\Diplom\\Model\\mesh\\xyz.dat"
#define MaterialsFileName "C:\\Users\\METAL\\Desktop\\Diplom\\Model\\properties\\materialsN.txt"
#define Element_MatNumFile "C:\\Users\\METAL\\Desktop\\Diplom\\Model\\mesh\\nvkat.dat"
#define FiguresContoursFile "C:\\Users\\METAL\\Desktop\\Diplom\\Model\\mesh\\ContoursMaterialGeometry"
#define NewMaterialsFileName "C:\\Users\\METAL\\Desktop\\Diplom\\Model\\properties\\materialsN_new.txt"
#define NewElement_MatNumFile "C:\\Users\\METAL\\Desktop\\Diplom\\Model\\mesh\\nvkat_new.dat"
#define FiguresInputFileName "Points.txt"
#define LogFileName "Log.txt"
#define PorosityFileName "Porosity.txt"
#define ResultsFileName "Stats.txt"
#define Depth -5 


#define uint unsigned int


using namespace std;

std::ofstream LogFile; //Глобальный файл логов

struct Material {
    int type;
    double K;
    double Phi;
    double something1;
    double something2;
};

struct Element {
    std::vector<int> points;
    int nmat;
    Element(int* pointmas) {
        for (uint i = 0; i < 8; i++) {
            points.push_back(pointmas[i]);
        }
    }
};

struct Grid {
    std::vector<Point> Nodes;
    std::vector<Element> Elements;
    std::vector<Material> Mats;
};

struct VolumeSt {
    std::vector<double> interval;
    double Volume;
    double Phi;
    double K;
};

//Ввод сетки - countfilename = "\\inftry.dat", 
//nodefilename = "\\xyz.dat", 
//elemfilename = "\\nver.dat", 
//materialsfilename = "\\materialsN"
//elem_matfilename = "\\nvkat.dat" 
Grid* InputGrid(const char* countfilename, const char* nodefilename, const char* elemfilename, const char* materialsfilename, const char *elem_matfilename) {

    Grid *tmpGrid = new Grid;
    Point tmppoint;
    Material mat;
    std::vector<Point> tmppoints;

    ifstream file;
    //Считываем количество узлов и точек
    file.open(countfilename, ios::in);
    if (!file) {
        LogFile << "Cannot open " << countfilename << std::endl;
        throw 1;
    }
    file.ignore(1000, '\n');
    file.ignore(1000, '=');
    uint nodecount;
    if (!(file >> nodecount)) {
        LogFile << "Mistake in " << countfilename << " - Node count is not a number" << std::endl;
        file.close();
        throw 1;
    }
    file.ignore(1000, '=');
    uint elemcount;
    if (!(file >> elemcount)) {
        LogFile << "Mistake in " << countfilename << " - Elements count is not a number" << std::endl;
        file.close();
        throw 1;
    }
    file.close();
    file.clear();

    tmpGrid->Nodes.reserve(nodecount);
    tmpGrid->Elements.reserve(elemcount);
    
    //Считываем точки узлов
    file.open(nodefilename, ios::binary | ios::in);
    if (!file) {
        LogFile << "Cannot open " << nodefilename << std::endl;
        throw 1;
    }
    double* tmpPoints = new double[3 * nodecount];
    if (!file.read((char*)tmpPoints, sizeof(double) * 3 * nodecount)) {
        LogFile << "Mistake in " << nodefilename << " - Node points is not a number" << std::endl;
        file.close();
        throw 1;
    }
    for (int i = 0; i < nodecount; i++)
        tmpGrid->Nodes.push_back(&tmpPoints[i * 3]);
    delete[] tmpPoints;
    file.close();
    file.clear();

    //Считываем номера узлов для элементов
    file.open(elemfilename, ios::binary | ios::in);
    if (!file) {
        LogFile << "Cannot open " << elemfilename << std::endl;
        throw 1;
    }
    //14 это 8 номеров по 4, и 6 байт игнорируем?
    int* tmpElem = new int[14 * elemcount];
    if (!file.read((char*)tmpElem, sizeof(int) * 14 * elemcount)) {
        LogFile << "Mistake in " << elemfilename << " - Element point is not a number" << std::endl;
        file.close();
        throw 1;
    }
    tmpGrid->Elements.push_back(Element(&tmpElem[0]));
    for (int i = 1; i < elemcount; i++)
        tmpGrid->Elements.push_back(Element(&tmpElem[i * 14]));
    delete[]tmpElem;
    file.close();
    file.clear();
    
    for (uint i = 0; i < elemcount; i++)
        for (uint j = 0; j < 8; j++) 
            tmpGrid->Elements[i].points[j]--;

    //Считываем материалы по номерам
    file.open(materialsfilename,ios::in);
    if (!file) {
        LogFile << "Cannot open " << materialsfilename << std::endl;
        throw 1;
    }
    uint matcount;
    if (!(file >> matcount)) {
        LogFile << "Mistake in " << materialsfilename << " - material count not a number" << std::endl;
        file.close();
        throw 1;
    }
    for (uint i = 0; i < matcount; i++) {
        if (!(file >> mat.type >> mat.K  >> mat.Phi >> mat.something1 >> mat.something2)) {
            LogFile << "Mistake in " << materialsfilename << " - Material property is not a number" << std::endl;
            file.close();
            throw 1;
        }
        tmpGrid->Mats.push_back(mat);
    }
    file.close();
    file.clear();

    //Считываем номера материалов для каждого из элементов
    file.open(elem_matfilename, ios::binary | ios::in);
    if (!file) {
        LogFile << "Cannot open " << elem_matfilename << std::endl;
        throw 1;
    }
    int* tmpNumbers = new int[elemcount];
    if (!file.read((char*)tmpNumbers, sizeof(int) * elemcount)) {
        LogFile << "Mistake in " << elem_matfilename << " - Material number is not a number" << std::endl;
        file.close();
        throw 1;
    }
    for (int i = 0; i < elemcount; i++)
        tmpGrid->Elements[i].nmat = --tmpNumbers[i];
    delete[] tmpNumbers;
    file.close();
    file.clear();

    return tmpGrid;
}

//Ввод фигур
std::vector<Figure> InputFigures(const char* figuresfilename) {
    ifstream figuresfile(figuresfilename);

    std::vector<Figure> tmpvector;
    std::vector<Point> topfigurepoints, bottomfigurepoints;

    if (figuresfile.is_open()) {

        int fignum;
        //Временная штука, убрать
        figuresfile >> fignum;
        figuresfile >> fignum;
        figuresfile >> fignum;
        figuresfile >> fignum;
        figuresfile >> fignum;
        //Временная штука, убрать
        
        //Считываем количество фигур
        if (!(figuresfile >> fignum)) {
            LogFile << "Mistake in " << figuresfilename << " - figure number is not a number" << std::endl;
            figuresfile.close();
            throw 2;
        }
        //Для каждой фигуры
        for (uint fn = 0; fn < fignum; fn++) {
            int numberofpoints;
            double Phi;
            int Order;
            if ((figuresfile >> numberofpoints))
                //Число точек
                if (numberofpoints > 0) {
                    for (int i = 0; i < numberofpoints; i++) {
                        Point tmppoint;
                        //Считываем координаты
                        if (!(figuresfile >> tmppoint.x) || !(figuresfile >> tmppoint.y)  /*!(figuresfile >> tmppoint.z)*/) {
                            LogFile << "Mistake in " << figuresfilename << " - figure point is not a number" << std::endl;
                            figuresfile.close();
                            throw 2;
                        }
                        tmppoint.z = 0;
                        topfigurepoints.push_back(tmppoint);
                        //Тут бы по хорошему считывать глубину слоя или типо того
                        tmppoint.z = Depth;
                        bottomfigurepoints.push_back(tmppoint);
                    }

                    //if (!(figuresfile >> Phi)) LogFile << "Mistake in " << figuresfilename << " - Phi is not a number" << std::endl;
                    Phi = 0.6;
                    Order = fignum - fn;
                    //if (!(figuresfile >> Order)) LogFile << "Mistake in " << figuresfilename << " - Order is not a number" << std::endl;

                    tmpvector.push_back(Figure(topfigurepoints, bottomfigurepoints, Phi, Order));
                    topfigurepoints.clear();
                    bottomfigurepoints.clear();
                }
                else {
                    LogFile << "Mistake in " << figuresfilename << " - Figures number is not a number or has unexpected value" << std::endl;
                    figuresfile.close();
                    throw 2;
                }
        }
    }
    else {
        LogFile << "Cannot open" << figuresfilename << std::endl;

    }
    figuresfile.close();
    return tmpvector;
}

//TODO: Добавить и привыкнуть к vector.reserve
//TODO: Поменять вывод материалов
//TODO: Добавить НОРМАЛЬНОЕ MinMax начальные значения
//TODO: Реализовать HeshMap для вывода материалов
//TODO: Добавить LogFile везде
//TODO: Исправить поведение по с слоями

//Основная программа
int main() {
    try{ 
        bool firstmeth = false;
        double tick;
        LogFile.open(LogFileName);
        std::ofstream PorosityFile;
        PorosityFile.open(PorosityFileName);
        Timer TimePassed;
        Grid *GeneralGrid = InputGrid(CountFileName, NodesFileName,ElementFileName, MaterialsFileName, Element_MatNumFile);
        std::vector<Figure> Figures = InputFigures(FiguresContoursFile);
        std::sort(Figures.begin(),Figures.end(),greater<Figure>()); //Сортировка по приоритету(по убыванию)
        
        //Открываем для записи нужные файлы
        ofstream newelem_numfile;
        newelem_numfile.open(NewElement_MatNumFile, ios::binary | ios::out);
        ofstream newmaterialsfile;
        newmaterialsfile.open(NewMaterialsFileName, ios::out);
        newmaterialsfile << GeneralGrid->Elements.size() << " " << std::endl;

        LogFile << "InputTime " << TimePassed.getTime() << std::endl;

        //Вектор обьемов по фигурам внутри КЭ
        std::vector<VolumeSt> Volume;
        
        //Очищаем фоновую(которая была в КЭ )
        Volume.reserve(Figures.size() + 1);
        Volume.push_back({ {}, 0, 0});

        //Количество разных площадей потенциально в каждом КЭ
        for (int i = 1; i < Figures.size()+1; i++)
            Volume.push_back({ {}, 0 , 0});
        tick = TimePassed.getTime();
        int pos = 0;

        //Главный цикл по КЭ
        for (Element CurrentElement : GeneralGrid->Elements) {
            std::vector<Point> ElementPolygon;
            Material CurMat = GeneralGrid->Mats[CurrentElement.nmat];
            pos++;
            
            //Формируем полигон КЭ, а заодно считаем его опоясывающий куб
            Point MinPoint(100000000,100000000, 100000000), MaxPoint(-100000000, -100000000, -100000000); //Минимальная, максимальная точки элемента
            for (int i : CurrentElement.points) {
                ElementPolygon.push_back(GeneralGrid->Nodes[i]);
                if (GeneralGrid->Nodes[i].x < MinPoint.x) MinPoint.x = GeneralGrid->Nodes[i].x;
                else if (GeneralGrid->Nodes[i].x > MaxPoint.x) MaxPoint.x = GeneralGrid->Nodes[i].x;

                if (GeneralGrid->Nodes[i].y < MinPoint.y) MinPoint.y = GeneralGrid->Nodes[i].y;
                else if (GeneralGrid->Nodes[i].y > MaxPoint.y) MaxPoint.y = GeneralGrid->Nodes[i].y;

                if (GeneralGrid->Nodes[i].z < MinPoint.z) MinPoint.z = GeneralGrid->Nodes[i].z;
                else if (GeneralGrid->Nodes[i].z > MaxPoint.z) MaxPoint.z = GeneralGrid->Nodes[i].z;
            }
            
            //Меняем 3 с 4 и 7 с 8 точки, для соблюдения обхода граней //Переделать в swap
            std::swap(ElementPolygon[2], ElementPolygon[3]);
            std::swap(ElementPolygon[6], ElementPolygon[7]);
       
            //Делаем для фоновой Phi значение равное значению до вычисления
            Volume[0].Phi = GeneralGrid->Mats[CurrentElement.nmat].Phi;
            
            //Обнуляем суммарные площади фигур на КЭ
            for (int i = 0; i < Figures.size() + 1; i++)
                Volume[i].Volume = 0;

            Point Step = { ((MaxPoint.x - MinPoint.x) / BreakKoeff),((MaxPoint.y - MinPoint.y) / BreakKoeff),((MaxPoint.z - MinPoint.z) / BreakKoeff) }; //Размер шага на КЭ
            double StepSquare = Step.y * Step.z;
            int halfElementsize = (int)(ElementPolygon.size() / 2);
            
            //Выбор между первым и следующим методом 
            if (firstmeth) {
                for (int i = 0; i < BreakKoeff; i++) {
                    for (int j = 0; j < BreakKoeff; j++) {
                        for (int k = 0; k < BreakKoeff; k++) {
                            int p = 0;
                            for (Figure CurFig : Figures) {
                                p++;
                                if (CurFig.PointInside({ MinPoint.x + Step.x * (i + 0.5) ,MinPoint.y + Step.y * (j + 0.5), MinPoint.z + Step.z * (k + 0.5) })) {
                                    Volume[p].Volume += StepSquare * Step.x;
                                    p = -1;
                                    break;
                                }
                            }
                            if (p != -1) Volume[0].Volume += StepSquare * Step.x;
                        }
                    }
                }
            }
            else {
                //Разбиваем на "столбы" и работаем с ними
                for (int i = 0; i < BreakKoeff; i++) {
                    for (int j = 0; j < BreakKoeff; j++) {

                        //Формируем рабочую прямую параллельную оси x
                        Point StartLine(MinPoint.x, MinPoint.y + Step.y * (i + 0.5), MinPoint.z + Step.z * (j + 0.5));
                        Point EndLine(MinPoint.x + 1, MinPoint.y + Step.y * (i + 0.5), MinPoint.z + Step.z * (j + 0.5));
                        //Точка пересечения
                        Point InterPoint;
                        for (int tfn = 0; tfn < halfElementsize; tfn++) { //TFN - TopFirstNumber //TSN - TopSecondNumber
                            int tsn;
                            tfn + 1 < halfElementsize ? tsn = tfn + 1 : tsn = 0;
                            //Проверяем пересекает ли плоскость грани наш "столб"
                            if (PlaneIntersectLine(ElementPolygon[tfn], ElementPolygon[tsn], ElementPolygon[tfn + halfElementsize], StartLine, EndLine, &InterPoint)) {
                                //Проверяем, находится ли точка в пределах элемента
                                if (InterPoint.x >= MinPoint.x && InterPoint.x <= MaxPoint.x)
                                    //Проверяем, находится ли точка в грани
                                    //TODO: Переделать методом триангуляции, перенести в PointInsidePolygon
                                    if (PointInsideTriangle(ElementPolygon[tfn], ElementPolygon[tsn], ElementPolygon[tfn + halfElementsize], InterPoint) || (PointInsideTriangle(ElementPolygon[tsn], ElementPolygon[tfn + halfElementsize], ElementPolygon[tsn + halfElementsize], InterPoint))) {
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
                            //Если точка внутри, добавляем первый начальный интервал
                            if (CurFig.PointInside({ StartLine.x + (Step.x * 0.001), StartLine.y, StartLine.z }))
                                Volume[fignum].interval.push_back(MinPoint.x);
                            int fnp = (int)CurFig.getTopPoints().size(); //fnp - Figure number of points
                            for (int fn = 0; fn < fnp; fn++) { //FN - FirstNumber //SN - SecondNumber
                                int sn;
                                fn + 1 < fnp ? sn = fn + 1 : sn = 0;
                                //Проверяем пересекает ли плоскость граней фигуры наш "столб"
                                if (PlaneIntersectLine(CurFig.getTopPoints()[fn], CurFig.getTopPoints()[sn], CurFig.getBottomPoints()[fn], StartLine, EndLine, &InterPoint)) {
                                    //Проверяем, находится ли точка в пределах фигуры
                                    if (InterPoint.x >= CurFig.getMinMax().first.x && InterPoint.x <= CurFig.getMinMax().second.x)
                                        //Проверяем, находится ли точка в грани
                                        //TODO: Переделать методом триангуляции, перенести в PointInsidePolygon
                                        if (PointInsideTriangle(CurFig.getTopPoints()[fn], CurFig.getTopPoints()[sn], CurFig.getBottomPoints()[fn], InterPoint) || (PointInsideTriangle(CurFig.getTopPoints()[sn], CurFig.getBottomPoints()[fn], CurFig.getBottomPoints()[sn], InterPoint))) {
                                            //Если точка пересечения попадает в главные интервалы добавляем ее в интервал текущей фигуры
                                            for (int num = 0; num < Volume[0].interval.size(); num += 2) {
                                                if (InterPoint.x >= Volume[0].interval[num] && InterPoint.x <= Volume[0].interval[num + 1])
                                                    Volume[fignum].interval.push_back(InterPoint.x);
                                            }
                                        }
                                }
                            }
                            std::sort(Volume[fignum].interval.begin(), Volume[fignum].interval.end());
                            if (Volume[fignum].interval.size() % 2 != 0) Volume[fignum].interval.push_back(MaxPoint.x);
                        }

                        //Перенести выше
                        //Добавляем суммарный обьем столба в соответствующем векторе Volume
                        //По каждой фигуре
                        for (int i = 1; i < Volume.size(); i++) {
                            //По интервалам фигуры
                            for (int j = 0; j < Volume[i].interval.size(); j += 2) {
                                //По главному интервалу
                                for (int mn = 0; mn < Volume[0].interval.size(); mn += 2) {
                                    bool leftinter = false;
                                    bool rightinter = false;
                                    //Левый интервал внутри
                                    if ((Volume[0].interval[mn] < Volume[i].interval[j]) && (Volume[0].interval[mn + 1] > Volume[i].interval[j]))
                                        leftinter = true;
                                    //Правый интервал внутри
                                    if ((Volume[0].interval[mn] < Volume[i].interval[j + 1]) && (Volume[0].interval[mn + 1] > Volume[i].interval[j + 1]))
                                        rightinter = true;

                                    //Интервал фигуры внутри
                                    if (leftinter && rightinter) {
                                        if ((Volume[0].interval[mn] != Volume[i].interval[j]) && (Volume[0].interval[mn + 1] != Volume[i].interval[j + 1]))
                                            Volume[0].interval.insert(Volume[0].interval.begin() + (mn + 1), { Volume[i].interval[j], Volume[i].interval[j + 1] });
                                        //Полностью равно интервалу
                                        else {
                                            //Два раза, потому что элемент смещается
                                            Volume[0].interval.erase(Volume[0].interval.begin() + (mn));
                                            Volume[0].interval.erase(Volume[0].interval.begin() + (mn));
                                        }
                                        Volume[i].Volume += StepSquare * (Volume[i].interval[j + 1] - Volume[i].interval[j]);
                                        break;
                                    }

                                    //Интервал наложен на правую границу
                                    else if (leftinter && !rightinter) {
                                        Volume[i].Volume += StepSquare * (Volume[0].interval[mn + 1] - Volume[i].interval[j]);
                                        //Меняем правую границу
                                        Volume[0].interval[mn + 1] = Volume[i].interval[j];
                                        j -= 2;
                                        break;
                                    }

                                    //Интервал наложен на левую границу
                                    else if (!leftinter && rightinter) {
                                        Volume[i].Volume += StepSquare * (Volume[i].interval[j + 1] - Volume[0].interval[mn]);
                                        //Меняем левую границу
                                        Volume[0].interval[mn] = Volume[i].interval[j + 1];
                                        j -= 2;
                                        break;
                                    }

                                    //Интервал наложен целиком, либо находится вне
                                    else {
                                        //Наложен целиком
                                        if ((Volume[0].interval[mn] >= Volume[i].interval[j]) && (Volume[0].interval[mn + 1] <= Volume[i].interval[j + 1])) {
                                            Volume[i].Volume += StepSquare * (Volume[0].interval[mn + 1] - Volume[0].interval[mn]);
                                            //Два раза, потому что элемент смещается
                                            Volume[0].interval.erase(Volume[0].interval.begin() + (mn));
                                            Volume[0].interval.erase(Volume[0].interval.begin() + (mn));
                                            j -= 2;
                                            break;
                                        }
                                    }
                                }
                            }
                        }
                        //Считаем оставшийся фоновый обьем
                        for (int mn = 0; mn < Volume[0].interval.size(); mn += 2)
                            Volume[0].Volume += (Volume[0].interval[mn + 1] - Volume[0].interval[mn]) * StepSquare;
                        for (int i = 0; i < Volume.size(); i++)
                            Volume[i].interval.clear();
                    }
                }
            }

            //Считаем получившуюся суммарную пористость и площадь КЭ
            double porosity = 0, calcVM = 0;
            for (VolumeSt VolumePart : Volume) {
                porosity += VolumePart.Volume * CurMat.Phi;
                calcVM += VolumePart.Volume;
            }
            porosity /= calcVM;

            //Записываем в файлы результаты для отладки
            PorosityFile << "El N = " << pos << " Phi = " << porosity << " V = " << calcVM << " XYZlb " << ElementPolygon[0].ToString() << std::endl;

            //Записываем новый материал
            Material tmpMat =  GeneralGrid->Mats[CurrentElement.nmat];
            newmaterialsfile << tmpMat.type << " " << tmpMat.K << " " << tmpMat.Phi << " " << tmpMat.something1 << " " << tmpMat.something2 << " " << std::endl;;
            
            //Привязываем этот материал к КЭ
            newelem_numfile.write((char*)&pos,sizeof(int));

            //Очищение полигона и переход к новому элементу 
            ElementPolygon.clear();

        }
        LogFile << "SumTime = " << TimePassed.getTime() << std::endl;
        LogFile << "Success";
        newelem_numfile.close();
        newmaterialsfile.close();
        PorosityFile.close();
        LogFile.close();
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
        LogFile.close();
        return ErrorCode;
    }
    
    LogFile.close();
    return 0;
}
