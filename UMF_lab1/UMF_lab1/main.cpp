#include <iostream>
#include <vector>
#include <fstream>

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
    Matrix(int n, int m, vector<int>& shift)
    {
        n = n;
        M = m;
        Mat.resize(n, vector<double>(M));
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
            return Mat[i][k - 1];
        }
        int index{ -1 };
        for (int i{ 0 }; i < Shift.size(); i++)
            if (numDiag == Shift[i])
            {
                index = i;
                break;
            }
        if (index == -1) throw new string("Diag not exist");
        return Mat[i][index];
    }

    void SetElement(int i, int j, double value)
    {
        int numDiag = j - i;
        int k{ 0 };
        if (numDiag == 0)
        {
            while (Shift[k++] < 0);
            Mat[i][k - 1] = value;
        }
        int index{ -1 };
        for (int i{ 0 }; i < Shift.size(); i++)
            if (numDiag == Shift[i])
            {
                index = i;
                break;
            }
        if (index == -1) throw new string("Diag not exist");
        Mat[i][index] = value;
    }
};

struct Point
{
    double X;
    double Y;
    short int Type;

    Point(double pointX, double pointY, short int type = Sham)
    {
        X = pointX;
        Y = pointY;
        Type = type;
    }
    Point()
    {
        X = 0;
        Y = 0;
        Type = Sham;
    }
};

class Name
{
private:
    vector<vector<Point>> Grid; /*{ {Point(0,4),Point(1,4),Point(2,4),Point(3,4),Point(4,4)},
                                   {Point(0,3),Point(1,3),Point(2,3),Point(3,3),Point(4,3)},
                                   {Point(0,2),Point(1,2),Point(2,2),Point(3,2),Point(4,2)},
                                   {Point(0,1),Point(1,1),Point(2,1),Point(3,1),Point(4,1)},
                                   {Point(0,0),Point(1,0),Point(2,0),Point(3,0),Point(4,0)}};*/
    vector<Point> Area = vector<Point>();
    vector<double> F = vector<double>();
    //Matrix A;
    vector<vector<double>> A = vector<vector<double>>();
    vector<double> U = vector<double>();
public:
    double GetF(double x, double y)
    {
        return x;
    }

    double GetGamma(double x, double y)
    {
        return 1;
    }

    double GetLambda(double x, double y)
    {
        return 0;
    }
    double GetBound(double x, double y, int edge)
    {
        switch (edge)
        {
        case 1:
            return 0;
            break;
        case 2:
            return x;
            break;
        case 3:
            return 1;
            break;
        case 4:
            return x;
            break;
        case 5:
            return 2;
            break;
        case 6:
            return x;
            break;
        case 7:
            return 3;
            break;
        case 8:
            return x;
            break;
        case 9:
            return 4;
            break;
        case 10:
            return x;
            break;
        case 11:
            return 5;
            break;
        case 12:
            return x;
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

    void GetGrid(int nx, int ny, bool gridType = true)
    {
        double xmax{ 0 }, ymax{ 0 }, xmin{ Area[0].X }, ymin{ Area[0].Y };
        for (int i{ 0 }; i < Area.size(); i++)
        {
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

        }
        for (int i{ 0 }; i < Grid.size(); i++)
        {
            cout << Grid[i][0].Y << ": ";
            for (int j{ 0 }; j < Grid[i].size(); j++)
                cout << Grid[i][j].X << "  ";
            cout << endl;
        }
    }

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

    // Если перепрыгивает строку, думает, что все ок
    void GetMatrix()
    {
        int n = Grid.size();
        int level = 0;
        int count = (n - 2) * (Grid[0].size() - 2);
        A.resize(count);
        for (int i{ 0 }; i < A.size(); i++)
            A[i].resize(count);
        F.resize(count);
        U.resize(count);
        for (int i{ 1 }; i < n - 1; i++)
        {
            int m = Grid[i].size();
            for (int j{ 1 }; j < m - 1; j++)
            {
                // Введем переменные:
                double hx_j{ Grid[i][j + 1].X - Grid[i][j].X }, hx_jp{ Grid[i][j].X - Grid[i][j - 1].X };
                double hy_i{ Grid[i - 1][j].Y - Grid[i][j].Y }, hy_ip{ Grid[i][j].Y - Grid[i + 1][j].Y };

                int index = GetNum(i, j, m);
                cout << i << " " << j << endl;
                A[level][index] = -(2 / (hx_j * hx_jp) + 2 / (hy_i * hy_ip)) + GetGamma(Grid[i][j].X, Grid[i][j].Y);

                index = GetNum(i, j - 1, m);
                if (index > GetNum(i, 0, m))
                    A[level][index] = 2 / (hx_jp * (hx_j + hx_jp));

                index = GetNum(i, j + 1, m);
                if (index < GetNum(i, m - 1, m))
                    A[level][index] = 2 / (hx_j * (hx_j + hx_jp));

                index = GetNum(i + 1, j, m);
                if (index < count)
                    A[level][index] = 2 / (hy_ip * (hy_i + hy_ip));

                index = GetNum(i - 1, j, m);
                if (index >= 0)
                    A[level][index] = 2 / (hy_i * (hy_i + hy_ip));


                F[level] = GetF(Grid[i][j].X, Grid[i][j].Y);
                level++;
            }
        }
        ofstream out("out.txt");
        for (int i{ 0 }; i < A.size(); i++)
        {
            for (int j{ 0 }; j < A[i].size(); j++)
                out << A[i][j] << "  ";
            out << endl;
        }

    }

    void SetBoundOne(string fileName)
    {
        ifstream inputFile(fileName);
        if (!inputFile.is_open()) throw string("File not found");
        int number{ 0 };
        //int start{ 0 }, end{ 0 };
        double value{ 0 };
        int west{ 0 };

        int n = Grid.size();
        int level = 0;
        int count = A.size();
        for (int i{ 1 }; i < n - 1; i++)
        {
            int m = Grid[i].size();
            for (int j{ 1 }; j < m - 1; j++)
            {
                double hx_j{ Grid[i][j + 1].X - Grid[i][j].X }, hx_jp{ Grid[i][j].X - Grid[i][j - 1].X };
                double hy_i{ Grid[i - 1][j].Y - Grid[i][j].Y }, hy_ip{ Grid[i][j].Y - Grid[i + 1][j].Y };


                // Check i-1, i+1, j-1, j+1

                int index = GetNum(i, j - 1, m);
                int type = GetType(Grid[i][j - 1], west);
                if (index <= GetNum(i, 0, m) && type == Bound)
                    F[level] += -2 * GetBound(Grid[i][j - 1].X, Grid[i][j - 1].Y, west + 1) / (hx_jp * (hx_j + hx_jp));

                index = GetNum(i, j + 1, m);
                type = GetType(Grid[i][j + 1], west);
                if (index >= GetNum(i, m - 1, m) && type == Bound)
                    F[level] += -2 * GetBound(Grid[i][j + 1].X, Grid[i][j + 1].Y, west + 1) / (hx_j * (hx_j + hx_jp));

                index = GetNum(i - 1, j, m);
                type = GetType(Grid[i - 1][j], west);
                if (index <= 0 && type == Bound)
                    F[level] += -2 * GetBound(Grid[i - 1][j].X, Grid[i - 1][j].Y, west + 1) / (hy_i * (hy_i + hy_ip));

                index = GetNum(i + 1, j, m);
                type = GetType(Grid[i + 1][j], west);
                if (index >= count && type == Bound)
                    F[level] += -2 * GetBound(Grid[i + 1][j].X, Grid[i + 1][j].Y, west + 1) / (hy_ip * (hy_i + hy_ip));

                index = GetNum(i, j, m);
                type = GetType(Grid[i][j], west);

                if (type == Bound)
                {
                    for (int i{ 0 }; i < A[level].size(); i++)
                        A[level][i] = 0;
                    A[level][level] = 1;
                    F[level] = GetBound(Grid[i][j].X, Grid[i][j].Y, west + 1);
                    // Занулить всю строку level, f[level] = value
                }
                if (type == Sham)
                {
                    for (int i{ 0 }; i < A[level].size(); i++)
                        A[level][i] = 0;
                    A[level][level] = 1;
                    F[level] = 0;
                }
                level++;
            }
        }


        ofstream out("out.txt");
        for (int i{ 0 }; i < A.size(); i++)
        {
            for (int j{ 0 }; j < A[i].size(); j++)
                out << A[i][j] << "  ";
            out << endl;
        }
        out.close();
        ofstream out1("f.txt");
        for (int i{ 0 }; i < F.size(); i++)
            out1 << F[i] << endl;

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

    template <typename T>
    double iteration(vector<vector<T>>& a, vector<T>& y0, vector<T>& y, const vector<T>& b, double w)
    {
        int n = a.size();
        double sum{ 0 }, sumdiscrepancy{ 0 }, discrepancy{ 0 };
        for (int i{ 0 }; i < n; i++)
        {
            double sum{ 0 };
            for (int j0{ 0 }; j0 < n; j0++)
            {
                sum += a[i][j0] * y0[j0];
            }
            discrepancy = b[i] - sum;
            y[i] = y0[i] + (w / a[i][i]) * discrepancy;
            sumdiscrepancy += discrepancy * discrepancy;
        }
        sumdiscrepancy = sqrt(sumdiscrepancy);
        return sumdiscrepancy;
    }

    template <typename T>
    int jacobi(vector<vector<T>>& a, vector<double>& x_start, vector<double>& x, const vector<double>& b, double w, int max_it, double eps)
    {
        int n = a.size();
        double eps_b{ 0 };
        vector<T> x_t(n);
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
        for (int i{ 0 }; i < 10; i++)
            cout << endl;
        return i;
    }

    void SolveSLAE(double w)
    {
        vector<double> xStart = vector<double>(A.size(), 0);
        jacobi(A, xStart, U, F, w, 10000, 1e-6);
    }

};

int main()
{
    Name a;
    //
    a.GetArea("area.txt");
    a.GetGrid(5, 10);
    a.GetMatrix();
    a.SetBoundOne("bound1.txt");
    a.SolveSLAE(1);
    return 0;
}