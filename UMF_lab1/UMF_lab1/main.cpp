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
    vector<vector<short int>> Area = vector<vector<short int>>();
    vector<double> F = vector<double>();
    //Matrix A;
    vector<vector<double>> A = vector<vector<double>>();
public:
    double GetF(double x, double y)
    {
        return 0;
    }
    void GetArea(string fileName)
    {
        int n{ 0 }, m{ 0 };
        ifstream inputFile(fileName);
        if (!inputFile.is_open()) throw string("File not found");
        inputFile >> n >> m;
        Area.resize(n);
        for (int i{ 0 }; i < n; i++)
        {
            Area[i].resize(m);
            for (int j{ 0 }; j < m; j++)
                inputFile >> Area[i][j];
        }
    }

    void GetGrid(int n, bool gridType = true)
    {
        if (gridType)
        {
            for (int i{ 0 }; i < Area.size(); i++)
            {
                vector<Point> tmp;
                ///double dx = Area[i].size() < n ? (double)Area[i].size() / n : 1;
                double dx{ 1.0 / n };
                double step{ 0 };
                //int count = (Area[i].size() % 2 == 0) ?  : Area[i].size() + 1;
                for (int j{ 0 }; step < Area[i].size() - 1; j++)
                {
                    step = j * dx;
                    tmp.push_back(Point(step, Area.size() - i - 1, Area[i][(int)step]));
                }
                Grid.push_back(tmp);
            }


        }
        else
        {

        }

    }

    int GetNum(int i, int j, int count)
    {
        return (i - 1) * (count - 1) + j - 1;
    }

    void GetMatrix()
    {
        int n = Grid.size();
        int level = 0;
        int count = (n - 2) * (Grid[0].size() - 2);
        A.resize(count);
        for (int i{ 0 }; i < A.size(); i++)
            A[i].resize(count);
        F.resize(count);
        for (int i{ 1 }; i < n - 1; i++)
        {
            int m = Grid[i].size();
            for (int j{ 1 }; j < m - 1; j++)
            {
                // ¬ведем переменные:
                double hx_j{ Grid[i][j + 1].X - Grid[i][j].X }, hx_jp{ Grid[i][j].X - Grid[i][j - 1].X };
                double hy_i{ Grid[i - 1][j].Y - Grid[i][j].Y }, hy_ip{ Grid[i][j].Y - Grid[i + 1][j].Y };
                
                int index = GetNum(i, j, n);
                cout << i << " " << j << endl;
                A[level][index] = -(2 / (hx_j * hx_jp) + 2 / (hy_i * hy_ip));
                index = GetNum(i, j - 1, n);
                if(index >= 0)
                    A[level][index] = 2 / (hx_jp * (hx_j + hx_jp));
                index = GetNum(i, j + 1, n);
                if (index  < count)
                    A[level][index] = 2 / (hx_j * (hx_j + hx_jp));
                index = GetNum(i + 1, j, n);
                if(index < count)
                    A[level][index] = 2 / (hy_ip * (hy_i + hy_ip));
                index = GetNum(i - 1, j, n);
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
        char fixed{ '\0' };
        int lock{ 0 }, start{ 0 }, end{ 0 };
        double value{ 0 };
        while (inputFile >> fixed >> lock >> start >> end >> value)
        {
            if (fixed == 'y')
            {
                int n = Grid.size();
                int level = 0;
                int count = A.size();
                for (int i{ 1 }; i < n - 1; i++)
                {
                    int m = Grid[i].size();
                    for (int j{ 1 }; j < m - 1; j++)
                    {
                        if (!((int)Grid[i][j].X < start || (int)Grid[i][j].X > end || (int)Grid[i][0].Y != lock))
                        {

                            //Centre
                            if (Grid[i][j].Type == Sham)
                            {
                                for (int k{ 0 }; k < count; k++)
                                    A[level][k] = 0;
                                A[level][GetNum(i, j, n)] = 1;
                                F[level] = 0;
                            }
                            if (Grid[i][j].Type == Bound)
                            {
                                for (int k{ 0 }; k < count; k++)
                                    A[level][k] = 0;
                                A[level][GetNum(i, j, n)] = 1;
                                F[level] = value;
                            }
                            if (Grid[i][j].Type == Inner)
                            {
                                int index = GetNum(i, j - 1, n);
                                if (index >= 0)
                                    A[level][index];
                                index = GetNum(i, j + 1, n);
                                if (index < count)
                                    A[level][index];
                                index = GetNum(i + 1, j, n);
                                if (index < count)
                                    A[level][index];
                                index = GetNum(i - 1, j, n);
                                if (index >= 0)
                                    A[level][index];
                            }
                        }
                        level++;
                    }
                }
            }
            else
            {

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

    void SolveSLAE()
    {

    }

};

int main()
{
    Name a;
    //
    a.GetArea("area.txt");
    a.GetGrid(1);
    a.GetMatrix();
    a.SetBoundOne("bound1.txt");
    return 0;
}