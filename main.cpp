#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matplotlibcpp.h"
#include <iostream>
#include <cmath>
#include <vector>

namespace plt = matplotlibcpp;
// namespace pyplt = matplotlibcpp::

const int K = 7;
const int K1 = 8;
const int N = 11;

// void gaussMethod(double c1[K1][K1], double d1[K1], double b[K1])
void gaussMethod(std::vector<std::vector<double>> &c, std::vector<double> &d, std::vector<double> &b)
{
    double q[K1][K1];
    int indD[K1];
    int cold[K1];

    for (int i = 0; i < K1; ++i)
    {
        indD[i] = i;
    }

    // reducing coefficient matrix to a triangular one
    for (int i = 0; i < K; ++i)
    {
        int i0 = i;
        int iflag = -1;
        int j1 = i + 1;
        double cmax = c[i][i];

        // массив не меняется
        // ищется максимум (cmax)
        // i0 становится равно j, а j = i + 1
        // iflag тогглится

        // НИХУЯ НЕ ДЕЛАЕТ
        for (int j = j1; j < K1; ++j)
        {
            if (fabs(cmax) < fabs(c[j][i]))
            {
                i0 = j;
                cmax = c[j][i];
                iflag = 1;
            }
        }

        // подмена какая-то
        // которая никогда не вызывается
        if (iflag == 1)
        {
            printf("iflag test\n");
            for (int j = 0; j < K1; ++j)
            {
                cold[j] = c[i][j];
                c[i][j] = c[i0][j];
                c[i0][j] = cold[j];
            }
            double dold = d[i];
            d[i] = d[i0];
            d[i0] = dold;

            double indDold = indD[i];
            indD[i] = indD[i0];
            indD[i0] = indDold;
        }

        // хер пойми что
        // массив меняется здесь
        for (int j = j1; j < K1; ++j)
        {
            q[j][i] = c[j][i] / c[i][i];
            d[j] -= q[j][i] * d[i];
            for (int j2 = j1; j2 < K1; ++j2)
            {
                c[j][j2] -= q[j][i] * c[i][j2];
            }
        }
    }

    b[K1 - 1] = d[K1 - 1] / c[K1 - 1][K1 - 1];

    for (int i = 0; i < K; ++i)
    {
        int IM = K - i - 1;
        double aux = d[IM];
        for (int j = IM; j < K; ++j)
        {
            aux -= c[IM][j + 1] * b[j + 1];
        }

        b[IM] = aux / c[IM][IM];
    }
}

int main()
{
    int K0 = 0;

    std::vector<double> coefVec(K1);
    std::vector<double> D(K1);
    std::vector<int> NI(K);
    std::vector<double> sumdy(K);

    // double C[K1][K1]
    // double YC[K][N];
    // double dyc[K][N];
    // double X[K][N];
    // double Y[K][N];
    std::vector<std::vector<double>> C(K1, std::vector<double>(K1));
    std::vector<std::vector<double>> YC(K, std::vector<double>(N));
    std::vector<std::vector<double>> dyc(K, std::vector<double>(N));
    std::vector<std::vector<double>> X(K, std::vector<double>(N));
    std::vector<std::vector<double>> Y(K, std::vector<double>(N));

    FILE *inp_file;
    if ((inp_file = fopen("OTa_q_gg_inp.dat", "r")) == NULL)
    {
        printf("Error opening a file");
        exit(1);
    }

    fscanf(inp_file, "%d", &K0);

    int NII = 0;
    

    if (K0 <= K)
    {
        for (int i = 0; i < K0; ++i)
        {
            fscanf(inp_file, "%d", &NI[i]);
        }

        std::cout<< "123" << std::endl;

        for (int i = 0; i < K0; ++i)
        {
            NII = NI[i];

            char str[100];
            for (int i1 = 0; i1 < 3; ++i1)
            {
                fscanf(inp_file, "%s", str);
            }

            // here X[i][j] is the Q_gg value of isotope j, element i
            // and Y[i][j] is the value of cross section's logarithm (called "Sigma") of isotope j, element i
            for (int j = 0; j < NII; ++j)
            {
                fscanf(inp_file, "%lf %lf", &X[i][j], &Y[i][j]);
            }
        }
    }
    else
    {
        fclose(inp_file);
        exit(1);
    }

    

    fclose(inp_file);

    // With given elements' values Q_gg and Sigma (cross-section log)
    // we calculate matrix elements values that are found by using
    // Least Squares method to our data
    for (int i = 1; i < K1; ++i)
    {
        NII = NI[i - 1];
        C[i][i] = (double)NII;
        for (int j = 0; j < NII; ++j)
        {
            C[0][0] += X[i - 1][j] * X[i - 1][j];
            C[0][i] += X[i - 1][j];
            C[i][0] += X[i - 1][j];
            D[0] += X[i - 1][j] * Y[i - 1][j];
            D[i] += Y[i - 1][j];
        }
    }

    printf("Matrix before gauss method: \n");

    for (int i = 0; i < K1; ++i)
    {
        for (int j = 0; j < K1; ++j)
        {
            printf("%f\t", C[i][j]);
        }
        printf("\n");
    }

    gaussMethod(C, D, coefVec);

    FILE *out_file;
    if ((out_file = fopen("OTa_q_gg_out.dat", "w")) == NULL)
    {
        printf("Error opening a file");
        exit(1);
    }

    for (int i = 0; i < K1; ++i)
    {
        fprintf(out_file, " \t %d\t %.17g\n", i + 1, coefVec[i]);
    }

    fprintf(out_file, "\n");

    double a = coefVec[0];

    // solving the y = ax + b equation

    // решаем уравнение y = ax + b
    for (int i = 0; i < K0; ++i)
    {
        double b = coefVec[i + 1];
        sumdy[i] = 0.0;
        for (int j = 0; j < NI[i]; ++j)
        {
            YC[i][j] = a * X[i][j] + b;
            dyc[i][j] = YC[i][j] - Y[i][j];
            if (fabs(YC[i][j]) > 1.0)
            {
                sumdy[i] += (dyc[i][j] * dyc[i][j] / fabs(YC[i][j]));
            }
        }
    }


    for (int i = 0; i < K0; ++i)
    {
        fprintf(out_file, "sigma\t %d\t %d\t %.17g\n", i + 1, NI[i], sumdy[i]);
    }

    fclose(out_file);

    std::vector<std::string> colorVec {"b", "r", "g", "c", "m", "y", "k"};

    for (int i = 0; i < X.size(); ++i) {
        for (int j = 0; j < X.size(); ++j) {
            plt::scatter(X[i], Y[i], 40);
        }
    }

    plt::xlim(-55, 5);
    plt::ylim(-20, 15);

    plt::plot(std::vector<int> {2000, -2000}, std::vector<double> {622.226, -605.471}, "b"); // 8.37
    plt::plot(std::vector<int> {2000, -2000}, std::vector<double> {628.855, -598.842}, "r"); // 15.00
    plt::plot(std::vector<int> {2000, -2000}, std::vector<double> {619.517, -608.179}, "g"); // 5.67
    plt::plot(std::vector<int> {2000, -2000}, std::vector<double> {618.458, -609.238}, "c"); // 4.61
    plt::plot(std::vector<int> {2000, -2000}, std::vector<double> {616.405, -611.291}, "m"); // 2.56
    plt::plot(std::vector<int> {2000, -2000}, std::vector<double> {615.13, -612.566}, "y");  // 1.28
    plt::plot(std::vector<int> {2000, -2000}, std::vector<double> {613.301, -614.395}, "k"); // -0.55
    
    plt::show();

    printf("we are back\n");

    return 0;
}
