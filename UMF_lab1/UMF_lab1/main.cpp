#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>
using namespace std;

enum Area
{
    Sham = 0,
    Bound = 1,
    Inner = 2
};

class Matrix
{
private:
    vector<vector<double>> Mat;
    int n{ 0 };
    int M{ 0 };
    vector<int> Shift;
public:
    Matrix()
    {

    }
    Matrix(int n, int m, vector<int>& shift)
    {
        this->n = n;
        M = m;
        Mat.resize(M, vector<double>(n));
        for (int i{ 0 }; i < shift.size(); i++)
            Shift.push_back(shift[i]);
    }
    double GetElement(int i, int j)
    {
        int numDiag = j - i;
        int k{ 0 };
        if (numDiag == 0)
        {
            while (Shift[k++] < 0);
            return Mat[k - 1][i];
        }
        int index{ -1 };
        for (int i{ 0 }; i < Shift.size(); i++)
            if (numDiag == Shift[i])
            {
                index = (numDiag < 0) ? i : i + 1;
                break;
            }
        if (index == -1) return 0;
        return Mat[index][i];
    }
    int size()
    {
        return n;
    }
    vector<int> GetShift()
    {
        return Shift;
    }
    void SetElement(int i, int j, double value)
    {
        int numDiag = j - i;
        int k{ 0 };
        if (numDiag == 0)
        {
            while (Shift[k++] < 0);
            Mat[k - 1][i] = value;
            return;
        }
        int index{ -1 };
        for (int i{ 0 }; i < Shift.size(); i++)
            if (numDiag == Shift[i])
            {
                index = (numDiag < 0) ? i : i + 1;
                break;
            }
        if (index == -1) return;//throw new string("Diag not exist");
        Mat[index][i] = value;
    }
};

struct Point
{
    double X;
    double Y;

    Point(double pointX, double pointY)
    {
        X = pointX;
        Y = pointY;
    }
    Point()
    {
        X = 0;
        Y = 0;
    }
};

class Name
{
private:
    vector<vector<Point>> Grid;
    vector<Point> Area = vector<Point>();
    vector<double> F = vector<double>();
    Matrix A;
    vector<double> U = vector<double>();
    vector<double> Utr = vector<double>();
public:
    Name()
    {

    }

    double GetUtr(double x, double y)
    {
        return pow(x, 4) + pow(y, 4);
    }

    double GetF(double x, double y)
    {
        return -12 * x * x + -12 * y * y + x * x * x * x + y * y * y * y;
        //return x * x * x * x + y * y * y * y;
        //return -4 + x * x + y * y;
        //return 0;
    }

    void sort(vector<double>& v)
    {
        for (int i{ 0 }; i < v.size(); i++)
            for (int j{0};j < v.size();j++)
                if (v[i] < v[j])
                {
                    auto tmp = v[i];
                    v[i] = v[j];
                    v[j] = tmp;
                }
    }

    double GetGamma(double x, double y)
    {
        return 1;
    }

    double GetLambda(double x, double y)
    {
        return 1;
    }

    double GetBound(double x, double y, int edge)
    {
       // return 5;
        /*switch (edge)
        {
        case 1:
            return y * y * y * y;
            break;
        case 2: 
            return x * x * x * x + 1;
            break;
        case 3:
            return 1 + y * y * y * y;
            break;
        case 4:
            return x * x * x * x;
            break;
        }*/
        switch (edge)
        {
        case 1:
            return 0 + y * y * y * y;
            break;
        case 2:
            return x * x * x * x + 1;
            break;
        case 3:
            return 1 + y * y * y * y;
            break;
        case 4:
            return x * x * x * x + 0.5 * 0.5 * 0.5 * 0.5;
            break;
        case 5:
            return 2 * 2 * 2 * 2 + y * y * y * y;
            break;
        case 6:
            return x * x * x * x + 1;
            break;
        case 7:
            return 3 * 3 * 3 * 3 + y * y * y * y;
            break;
        case 8:
            return x * x * x * x + 0.5 * 0.5 * 0.5 * 0.5;
            break;
        case 9:
            return 4 * 4 * 4 * 4 + y * y * y * y;
            break;
        case 10:
            return x * x * x * x + 1;
            break;
        case 11:
            return 5 * 5 * 5 * 5 + y * y * y * y;
            break;
        case 12:
            return x * x * x * x;
            break;
        }
    }

    void GetArea(string fileName)
    {
        ifstream inputFile(fileName);
        if (!inputFile.is_open()) throw string("File not found");
        double x, y;
        while (inputFile >> x >> y)
        {
            Area.push_back(Point(x, y));
        }
    }

    void DevGrid(double hx, double hy, bool gridType = true)
    {
        double xmax{ 0 }, ymax{ 0 }, xmin{ Area[0].X }, ymin{ Area[0].Y };
        vector<double> xi, yi;
        for (int i{ 0 }; i < Area.size(); i++)
        {
            if (xi.empty()) xi.push_back(Area[i].X);
            else
            {
                bool isContain = false;
                for (int k{ 0 }; k < xi.size(); k++)
                {
                    if (xi[k] == Area[i].X) { isContain = true; break; }
                }
                if (!isContain) xi.push_back(Area[i].X);
            }
            if (yi.empty()) yi.push_back(Area[i].Y);
            else
            {
                bool isContain = false;
                for (int k{ 0 }; k < yi.size(); k++)
                {
                    if (yi[k] == Area[i].Y) { isContain = true; break; }
                }
                if (!isContain) yi.push_back(Area[i].Y);
            }
            if (Area[i].X > xmax) xmax = Area[i].X;
            if (Area[i].X < xmin) xmin = Area[i].X;
            if (Area[i].Y > ymax) ymax = Area[i].Y;
            if (Area[i].Y < ymin) ymin = Area[i].Y;
        }
        if (gridType)
        {
            double l{ xmax - xmin };
            double m{ ymax - ymin };
            int j{ 0 };
            int i{ 0 };
            while (ymin + (double)i / hy < ymax)
            {
                vector<Point> tmp;
                j = 0;
                while (xmin + (double)j / hx < xmax)
                {
                    tmp.push_back(Point(xmin + (double)j / hx, ymin + (double)i / hy));
                    j++;
                }
                tmp.push_back(Point(xmax, ymin + (double)i  / hy));
                Grid.insert(Grid.begin(), tmp);
                i++;
            }
            vector<Point> tmp;
            j = 0;
            while (xmin + (double)j / hx < xmax)
            {
                tmp.push_back(Point(xmin + (double)j / hx, ymax));
                j++;
            }
            tmp.push_back(Point(xmax, ymax));
            Grid.insert(Grid.begin(), tmp);
        }
        else
        {
            int nx = hx, ny = hy;
            sort(xi);
            sort(yi);
            double l{ xmax - xmin };
            double m{ ymax - ymin };
            const double pi{ 3.1415926535897932 };
            int shifty = yi.size() - 2;
            for (int i{ ny }; i >= 0; i--)
            {
                vector<Point> tmp;
                double y = (i == ny) ? ymax : m * (1 - cos(pi * i / (2 * ny)));
                double yt = (i == ny) ? ymax : m * (1 - cos(pi * (i - 1) / (2 * ny)));
                int shiftx = 1;
                if (y > yi[shifty] && yt < yi[shifty])
                {
                    for (int j{ 0 }; j <= nx; j++)
                    {
                        double x = (j == nx) ? xmax : l * (1 - cos(pi * j / (2 * nx)));
                        double xt{ (j == nx) ? xmax : l * (1 - cos(pi * (j + 1) / (2 * nx))) };
                        if (x < xi[shiftx] && xt > xi[shiftx])
                        {
                            tmp.push_back(Point(x, y));
                            tmp.push_back(Point(xi[shiftx++], y));
                        }
                        else
                            tmp.push_back(Point(x, y));
                    }
                    Grid.push_back(tmp);
                    tmp.resize(0);
                    y = yi[shifty--];
                    shiftx = 1;
                    for (int j{ 0 }; j <= nx; j++)
                    {
                        double x = (j == nx) ? xmax : l * (1 - cos(pi * j / (2 * nx)));
                        double xt{ (j == nx) ? xmax : l * (1 - cos(pi * (j + 1) / (2 * nx))) };
                        if (x < xi[shiftx] && xt > xi[shiftx])
                        {
                            tmp.push_back(Point(x, y));
                            tmp.push_back(Point(xi[shiftx++], y));
                        }
                        else
                            tmp.push_back(Point(x, y));
                    }
                    Grid.push_back(tmp);

                }
                else
                {
                    for (int j{ 0 }; j <= nx; j++)
                    {
                        double x = (j == nx) ? xmax : l * (1 - cos(pi * j / (2 * nx)));
                        double xt{ (j == nx) ? xmax : l * (1 - cos(pi * (j + 1) / (2 * nx))) };
                        if (x < xi[shiftx] && xt > xi[shiftx])
                        {
                            tmp.push_back(Point(x, y));
                            tmp.push_back(Point(xi[shiftx++], y));
                        }
                        else
                            tmp.push_back(Point(x, y));
                    }
                    Grid.push_back(tmp);
                }

            }
        }
    }

    /*void GetGrid(int nx, int ny, bool gridType = true)
    {
        double xmax{ 0 }, ymax{ 0 }, xmin{ Area[0].X }, ymin{ Area[0].Y };
        vector<double> xi, yi;
        for (int i{ 0 }; i < Area.size(); i++)
        {
            if (xi.empty()) xi.push_back(Area[i].X);
            else
            {
                bool isContain = false;
                for (int k{ 0 }; k < xi.size(); k++)
                {
                    if (xi[k] == Area[i].X) { isContain = true; break; }
                }
                if (!isContain) xi.push_back(Area[i].X);
            }
            if (yi.empty()) yi.push_back(Area[i].Y);
            else
            {
                bool isContain = false;
                for (int k{ 0 }; k < yi.size(); k++)
                {
                    if (yi[k] == Area[i].Y) { isContain = true; break; }
                }
                if (!isContain) yi.push_back(Area[i].Y);
            }
            if (Area[i].X > xmax) xmax = Area[i].X;
            if (Area[i].X < xmin) xmin = Area[i].X;
            if (Area[i].Y > ymax) ymax = Area[i].Y;
            if (Area[i].Y < ymin) ymin = Area[i].Y;
        }
        if (gridType)
        {
            double l{ xmax - xmin };
            double m{ ymax - ymin };
            double dx{ l / nx };
            double dy{ m / ny };
            for (int i{ ny }; i >= 0; i--)
            {
                vector<Point> tmp;
                for (int j{ 0 }; j <= nx; j++)
                {
                    tmp.push_back(Point(xmin + j * dx, ymin + i * dy));
                }
                Grid.push_back(tmp);
            }

        }
        else
        {
            sort(xi);
            sort(yi);
            double l{ xmax - xmin };
            double m{ ymax - ymin };
            const double pi{ 3.1415926535897932 };
            int shifty = yi.size() - 2;
            for (int i{ ny }; i >= 0; i--)
            {
                vector<Point> tmp;
                double y = (i == ny) ? ymax : m * (1 - cos(pi * i / (2 * ny)));
                double yt = (i == ny) ? ymax : m * (1 - cos(pi * (i - 1) / (2 * ny)));
                int shiftx = 1;
                if (y > yi[shifty] && yt < yi[shifty])
                {
                    for (int j{ 0 }; j <= nx; j++)
                    {
                        double x = (j == nx) ? xmax : l * (1 - cos(pi * j / (2 * nx)));
                        double xt{ (j == nx) ? xmax : l * (1 - cos(pi * (j + 1) / (2 * nx))) };
                        if (x < xi[shiftx] && xt > xi[shiftx])
                        {
                            tmp.push_back(Point(x, y));
                            tmp.push_back(Point(xi[shiftx++], y));
                        }
                        else
                            tmp.push_back(Point(x, y));
                    }
                    Grid.push_back(tmp);
                    tmp.resize(0);
                    y = yi[shifty--];
                    shiftx = 1;
                    for (int j{ 0 }; j <= nx; j++)
                    {
                        double x = (j == nx) ? xmax : l * (1 - cos(pi * j / (2 * nx)));
                        double xt{ (j == nx) ? xmax : l * (1 - cos(pi * (j + 1) / (2 * nx))) };
                        if (x < xi[shiftx] && xt > xi[shiftx])
                        {
                            tmp.push_back(Point(x, y));
                            tmp.push_back(Point(xi[shiftx++], y));
                        }
                        else
                            tmp.push_back(Point(x, y));
                    }
                    Grid.push_back(tmp);

                }
                else
                {
                    for (int j{ 0 }; j <= nx; j++)
                    {
                        double x = (j == nx) ? xmax : l * (1 - cos(pi * j / (2 * nx)));
                        double xt{ (j == nx) ? xmax : l * (1 - cos(pi * (j + 1) / (2 * nx))) };
                        if (x < xi[shiftx] && xt > xi[shiftx])
                        {
                            tmp.push_back(Point(x, y));
                            tmp.push_back(Point(xi[shiftx++], y));
                        }
                        else
                            tmp.push_back(Point(x, y));
                    }
                    Grid.push_back(tmp);
                }
                
            }
        }
        for (int i{ 0 }; i < Grid.size(); i++)
        {
            cout << Grid[i][0].Y << ": ";
            for (int j{ 0 }; j < Grid[i].size(); j++)
                cout << Grid[i][j].X << "  ";
            cout << endl;
        }
    }*/

    int GetNum(int i, int j, int count)
    {
        return (i - 1) * (count - 2) + j - 1;
    }

    int GetType(const Point& point, int& west)
    {
        int total{ 0 };
        for (int i{ 0 }; i < Area.size(); i++)
        {
            Point* a = &Area[i];
            Point* b = (i == Area.size() - 1) ? &Area[0] : &Area[i + 1];
            if (a->X > b->X || a->Y > b->Y)
            {
                Point* tmp = a;
                a = b;
                b = tmp;
            }

            if (a->X != b->X && a->Y != b->Y)
            {
                if ((point.X - a->X) / (b->X - a->X) == (point.Y - a->Y) / (b->Y - a->Y))
                {
                    west = i;
                    return Bound;
                }

            }
            else
                if (a->X == b->X)
                {
                    if (point.X == a->X && point.Y >= a->Y && point.Y <= b->Y)
                    {
                        west = i;
                        return Bound;
                    }
                }
                else
                    if (a->Y == b->Y)
                    {
                        if (point.Y == a->Y && point.X >= a->X && point.X <= b->X)
                        {
                            west = i;
                            return Bound;
                        }
                    }

            if (point.X <= a->X && point.Y > a->Y && point.Y < b->Y)
                total++;
        }
        west = -1;
        return (total % 2 == 0) ? Sham : Inner;
    }

    void GetMatrix()
    {
        int n = Grid.size();
        int level = 0;
        int count = (n - 2) * (Grid[0].size() - 2);
        vector<int> shift = vector<int>{ -(int)Grid[0].size() + 2, -1, 1, (int)Grid[0].size() - 2 };
        A = Matrix(count, 5, shift);
        F.resize(count);
        U.resize(count);
        int west{ 0 };
        for (int i{ 1 }; i < n - 1; i++)
        {
            int m = Grid[i].size();
            for (int j{ 1 }; j < m - 1; j++)
            {
                // ???????????? ????????????????????:
                double hx_j{ Grid[i][j + 1].X - Grid[i][j].X }, hx_jp{ Grid[i][j].X - Grid[i][j - 1].X };
                double hy_i{ Grid[i - 1][j].Y - Grid[i][j].Y }, hy_ip{ Grid[i][j].Y - Grid[i + 1][j].Y };
                //i,j
                int index = GetNum(i, j, m);
                cout << i << " " << j << endl;
                A.SetElement(level, index, (2 / (hx_j * hx_jp) + 2 / (hy_i * hy_ip)) + GetGamma(Grid[i][j].X, Grid[i][j].Y));
                
                //j - 1
                index = GetNum(i, j - 1, m);
                if (index > GetNum(i, 0, m))
                    A.SetElement(level, index, -2 / (hx_jp * (hx_j + hx_jp)));
                else
                {
                    if (GetType(Grid[i][j - 1], west) == Bound)
                        F[level] += 2 * GetBound(Grid[i][j - 1].X, Grid[i][j - 1].Y, west + 1) / (hx_jp * (hx_j + hx_jp));
                }

                // j + 1
                index = GetNum(i, j + 1, m);
                if (index < GetNum(i, m - 1, m))
                    A.SetElement(level, index, -2 / (hx_j * (hx_j + hx_jp)));
                else
                {
                    if (GetType(Grid[i][j + 1], west) == Bound)
                        F[level] += 2 * GetBound(Grid[i][j + 1].X, Grid[i][j + 1].Y, west + 1) / (hx_j * (hx_j + hx_jp));
                }
                // i + 1
                index = GetNum(i + 1, j, m);
                if (index < count)
                    A.SetElement(level, index, -2 / (hy_ip * (hy_i + hy_ip)));
                else
                {
                    if (GetType(Grid[i + 1][j], west) == Bound)
                        F[level] += 2 * GetBound(Grid[i + 1][j].X, Grid[i + 1][j].Y, west + 1) / (hy_ip * (hy_i + hy_ip));
                }
                //i-1
                index = GetNum(i - 1, j, m);
                if (index >= 0)
                    A.SetElement(level, index, -2 / (hy_i * (hy_i + hy_ip)));
                else
                {
                    if (GetType(Grid[i - 1][j], west) == Bound)
                        F[level] += 2 * GetBound(Grid[i - 1][j].X, Grid[i - 1][j].Y, west + 1) / (hy_i * (hy_i + hy_ip));
                }

                F[level]+= GetF(Grid[i][j].X, Grid[i][j].Y);
                level++;
            }
        }
    }
   
    void SetBoundOne()
    {
        int west{ 0 };
        int n = Grid.size();
        int level = 0;
        int count = A.size();
        for (int i{ 1 }; i < n - 1; i++)
        {
            int m = Grid[i].size();
            for (int j{ 1 }; j < m - 1; j++)
            {
                int type = GetType(Grid[i][j], west);

                if (type == Bound)
                {
                    for (int i{ 0 }; i < A.size(); i++)
                        A.SetElement(level, i, 0);
                    A.SetElement(level, level, 1);
                    F[level] = GetBound(Grid[i][j].X, Grid[i][j].Y, west + 1);
                    // ???????????????? ?????? ???????????? level, f[level] = value
                }
                if (type == Sham)
                {
                    for (int i{ 0 }; i < A.size(); i++)
                        A.SetElement(level, i, 0);
                    A.SetElement(level, level, 1);
                    F[level] = 0;
                }
                level++;
            }
        }


        ofstream out("out.txt");
        for (int i{ 0 }; i < A.size(); i++)
        {
            for (int j{ 0 }; j < A.size(); j++)
                out << A.GetElement(i,j) << "  ";
            out << endl;
        }
        out.close();
        ofstream out1("f.txt");
        for (int i{ 0 }; i < F.size(); i++)
            out1 << scientific << setprecision(16) << F[i] << endl;
    }

    template <typename T>
    double get_norm(const vector<T>& v)
    {
        double result{ 0 };
        int size{ (int)v.size() };
        for (int i{ 0 }; i < size; i++)
            result += v[i] * v[i];
        return sqrt(result);

    }

    double iteration(Matrix& a, vector<double>& y0, vector<double>& y, const vector<double>& b, double w)
    {
        int n = a.size();
        double sum{ 0 }, sumdiscrepancy{ 0 }, discrepancy{ 0 };
        vector<int> shift = A.GetShift();
        for (int i{ 0 }; i < n; i++)
        {
            double sum{ 0 };
            for (int j0{ 0 }; j0 < 5; j0++)
            {
                int j = (j0 == 2) ? 0 : shift[(j0 < 2) ? j0 : j0 - 1];
                sum += a.GetElement(i, i + j) * y0[((i + j) >= 0 && i + j < n) ? i + j : 0];
            }
            discrepancy = b[i] - sum;
            y[i] = y0[i] + (w / a.GetElement(i,i)) * discrepancy;
            sumdiscrepancy += discrepancy * discrepancy;
        }
        sumdiscrepancy = sqrt(sumdiscrepancy);
        return sumdiscrepancy;
    }

    int jacobi(Matrix& a, vector<double>& x_start, vector<double>& x, const vector<double>& b, double w, int max_it, double eps)
    {
        int n = a.size();
        double eps_b{ 0 };
        vector<double> x_t(n);
        for (int i{ 0 }; i < n; i++)
            x_t[i] = x_start[i];
        int i{ 1 };
        double norm{ get_norm(b) };
        double discrepancy{ 0 };
        do
        {
            discrepancy = iteration(a, x_t, x, b, w);
            for (int k{ 0 }; k < n; k++)
                x_t[k] = x[k];

            eps_b = discrepancy / norm;
            //cout << "iteration: " << i << " discrepancy: " << eps_b << endl;
            i++;
        } while (i <= max_it && eps_b > eps);
        if (i > max_it) cout << "max iteration reached" << endl;
        else cout << "Ja:discrepancy " << discrepancy << endl;
        //for (int i{ 0 }; i < 10; i++)
            //cout << endl;
        return i;
    }

    void SolveSLAE(double w)
    {
        vector<double> xStart = vector<double>(A.size(), 1);
        jacobi(A, xStart, U, F, w, 10000, 1e-14);
        
       
    }

    void PrintResult(string fileName)
    {
        ofstream outputFile(fileName);
        int n = Grid.size();
        Utr.resize(A.size());
        int level{ 0 };
        for (int i{ 1 }; i < n - 1; i++)
        {
            int m = Grid[i].size();
            for (int j{ 1 }; j < m - 1; j++)
            {
                Utr[level] = GetUtr(Grid[i][j].X, Grid[i][j].Y);
                outputFile << scientific << setprecision(16) << Grid[i][j].X << " " << Grid[i][j].Y << " " << U[level] << endl;
                level++;
            }
        }
        ofstream out("f_true.txt");
        for (int i{ 0 }; i < A.size(); i++)
        {
            double sum = 0;
            for (int j{ 0 }; j < A.size(); j++)
            {
                sum += A.GetElement(i, j) * Utr[j];
            }
            out << scientific << setprecision(16) << sum << endl;
        }
    }
};

int main()
{
    Name a;
    //
    a.GetArea("area.txt");
    a.DevGrid(24, 24, false);
    //a.GetGrid(10,4, false);
    a.GetMatrix();
    a.SetBoundOne();
    a.SolveSLAE(1);
    a.PrintResult("result.txt");
    return 0;
}