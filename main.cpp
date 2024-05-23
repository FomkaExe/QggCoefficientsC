#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matplotlibcpp.h"
#include <vector>

namespace plt = matplotlibcpp;

const int K = 7;
const int K1 = 8;
const int N = 11;

void gaussMethod(double c1[K1][K1], double d1[K1], double b[K1])
{
    double q[K1][K1];
    int indD[K1];
    int cold[K1];

    // зануление
    for (int i = 0; i < K1; ++i)
    {
        indD[i] = i;
    }

    // printf("\nD1!!!!!!!\n");
    // for (int kok = 0; kok < K1; ++kok) {
    //     printf("%f\t", d1[kok]);
    // }
    // printf("\n\n");
    // reducing coefficient matrix to a triangular one

    // Приводим матрицу коэффициентов к треугольной
    for (int i = 0; i < K; ++i)
    {
        int i0 = i;
        int iflag = -1;
        int j1 = i + 1;
        double cmax = c1[i][i];

        // массив не меняется
        // ищется максимум (cmax)
        // i0 становится равно j, а j = i + 1
        // iflag тогглится

        // НИХУЯ НЕ ДЕЛАЕТ
        for (int j = j1; j < K1; ++j)
        {
            if (fabs(cmax) < fabs(c1[j][i]))
            {
                i0 = j;
                cmax = c1[j][i];
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
                cold[j] = c1[i][j];
                c1[i][j] = c1[i0][j];
                c1[i0][j] = cold[j];
            }
            double dold = d1[i];
            d1[i] = d1[i0];
            d1[i0] = dold;

            double indDold = indD[i];
            indD[i] = indD[i0];
            indD[i0] = indDold;
        }

        // хер пойми что
        // массив меняется здесь
        for (int j = j1; j < K1; ++j)
        {
            q[j][i] = c1[j][i] / c1[i][i];
            d1[j] -= q[j][i] * d1[i];
            for (int j2 = j1; j2 < K1; ++j2)
            {
                c1[j][j2] -= q[j][i] * c1[i][j2];
            }
            // printf("\nMATRIX IN MNK\ni = %d, j = %d\n", i, j);
            // for (int hui = 0; hui < K1; ++hui)
            // {
            //     for (int jopa = 0; jopa < K1; ++jopa)
            //     {
            //         printf("%f\t", c1[hui][jopa]);
            //     }
            //     printf("\n");
            // }
        }
    }

    // printf("\nD1!!!!!!!\n");
    // for (int kok = 0; kok < K1; ++kok) {
    //     printf("%f\t", d1[kok]);
    // }
    // printf("\n");


    b[K1 - 1] = d1[K1 - 1] / c1[K1 - 1][K1 - 1];

    for (int i = 0; i < K; ++i)
    {
        int IM = K - i - 1;
        double aux = d1[IM];
        for (int j = IM; j < K; ++j)
        {
            aux -= c1[IM][j + 1] * b[j + 1];
        }

        b[IM] = aux / c1[IM][IM];
    }
}

int main()
{
    int K0 = 0;

    double Z[K1];
    double D[K1];
    double C[K1][K1];

    int NI[K];
    double sumdy[K];

    double YC[K][N];
    double dyc[K][N];
    double X[K][N];
    double Y[K][N];

    for (int i = 0; i < K; ++i)
    {
        NI[i] = 0;
        sumdy[i] = 0;
        for (int j = 0; j < N; ++j)
        {
            X[i][j] = 0;
            Y[i][j] = 0;
            YC[i][j] = 0;
            dyc[i][j] = 0;
        }
    }
    for (int i = 0; i < K1; ++i)
    {
        D[i] = 0;
        Z[i] = 0;
        for (int j = 0; j < K1; ++j)
        {
            C[i][j] = 0;
        }
    }

    FILE *inp_file;
    if ((inp_file = fopen("OTa_q_gg_inp.dat", "r")) == NULL)
    {
        printf("Error opening a file");
        exit(1);
    }

    fscanf(inp_file, "%d", &K0);

    std::vector<std::vector<double>> iks(K, std::vector<double>(N));
    std::vector<std::vector<double>> igrek(K, std::vector<double>(N));

    // for (int i = 0; i < K; ++i) {
    //     for (int j = 0; j < N; ++j) {
    //         iks[i][j] = X[i][j];
    //         igrek[i][j] = Y[i][j];
    //     }
    // }

    int NII = 0;

    if (K0 <= K)
    {
        for (int i = 0; i < K0; ++i)
        {
            fscanf(inp_file, "%d", &NI[i]);
        }

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

            // Здесь X[i][j] - значение величины Q_gg j-ого изотопа i-ого элемента
            // и Y[i][j] - значение логарифма сечения "Сигма" j-ого изотопа i-ого элемента
            for (int j = 0; j < NII; ++j)
            {
                fscanf(inp_file, "%lf %lf", &X[i][j], &Y[i][j]);
                iks[i][j] = X[i][j];
                igrek[i][j] = Y[i][j];
            }
        }
    }
    else
    {
        fclose(inp_file);
        exit(1);
    }

    fclose(inp_file);

    plt::plot(iks.at(0).at(0));
    plt::show();

    // With given elements' values Q_gg and Sigma (cross section log)
    // we calculate matrix elements values that are found by using
    // Least Squares method to our data

    // По известным значениям величин Q_gg и Сигма
    // находим элементы матрицы которые получаются применением
    // метода наименьших квадратов к исходным данным
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

    gaussMethod(C, D, Z);

    FILE *out_file;
    if ((out_file = fopen("OTa_q_gg_out.dat", "w")) == NULL)
    {
        printf("Error opening a file");
        exit(1);
    }

    for (int i = 0; i < K1; ++i)
    {
        fprintf(out_file, " \t %d\t %.17g\n", i + 1, Z[i]);
    }

    fprintf(out_file, "\n");

    double a = Z[0];

    // solving the y = ax + b equation

    // решаем уравнение y = ax + b
    for (int i = 0; i < K0; ++i)
    {
        double b = Z[i + 1];
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

    printf("\nwe are back\n");

    return 0;
}
