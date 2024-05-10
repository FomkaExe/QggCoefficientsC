// Calculation of Qgg coefficients (Q - thermal energy, g - ground state, 
// Qgg - difference of ground state energies of the colliding ions, also called as elements lines)
// for each element produced in the nuclear reactions cross-sections of different isotopes of each element 'S' form a line in Q-S plane
// our task is to find the parallel lines y=a*x+b that descibes the cross-sections of different elements with the best accuracy
// for this we use least-square root method
// after applying it to our problem we obtain the system of linear equations
// to solve it we use Gauss method
// the last step is to find the accuracy of our solution

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

const int K = 7;
const int K1 = 8;
const int N = 11;

void gaussMethod(double c1[K1][K1], double d1[K1], double b[K1])
{
    double cnew[K1][K1];
    double q[K1][K1];
    int indD[K1];
    int cold[K1];

    for (int i = 0; i < K1; ++i)
    {
        indD[i] = i;
        for (int j = 0; j < K1; ++j)
        {       
            cnew[i][j] = 0.0;
        }
    }

    // reducing coefficient matrix to a triangular one

    // Приводим матрицу коэффициентов к треугольной
    for (int i = 0; i < K; ++i) {
        int i0 = i;
        int iflag = -1;
        int j1 = i + 1;
        double cmax = c1[i][i];

        for (int j = j1; j < K1; ++j) {
            if (fabs(cmax) < fabs(c1[j][i])) {
                i0 = j;
                cmax = c1[j][i];
                iflag = 1;
            }
        }

        if (iflag == 1) {
            for (int j = 0; j < K1; ++j) {
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

        for (int j = j1; j < K1; ++j) {
            q[j][i] = c1[j][i]/c1[i][i];
            cnew[j][i] = q[j][i];
            d1[j] -= q[j][i]*d1[i];
            for (int j2 = j1; j2 < K1; ++j2) {
                c1[j][j2] -= q[j][i]*c1[i][j2];
            }
        }
    }

    // finding the results

    // Находим результаты
    b[K1 - 1] = d1[K1 - 1]/c1[K1 - 1][K1 - 1];

    for (int i = 0; i < K; ++i) {
        int IM = K - i - 1;
        double aux = d1[IM];
        for (int j = IM; j < K; ++j) {
            aux -= c1[IM][j + 1]*b[j + 1];
        }

        b[IM] = aux/c1[IM][IM];
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


    // zero initialising arrays

    // инициализируем массивы нулями
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

    // opening input file "inp_file"
    // if input file cannot open we terminate the program

    // открываем файл входных данных "inp_file"
    // если файл не может открыться мы заканчиваем работу программы
    FILE *inp_file;
    if ((inp_file = fopen("OTa_q_gg_inp.dat", "r")) == NULL)
    {
        printf("Error opening a file");
        exit(1);
    }
    // reading the first line of the file and assigning the value to K0
    // K0 is the amount of different elements in inp_file

    // считываем первую строку файла и записываем ее значение в K0
    // K0 - кол-во разных элементов в inp_file
    fscanf(inp_file, "%d", &K0);

    int NII = 0;

    // If amount of elements in file (K0) is bigger than the K constant,
    // we terminate the program
    // else we read data from inp_file

    // Если кол-во элементов в файле (K0) больше чем константа K,
    // то мы заканчиваем работу программу
    // в ином случае мы читаем данные из inp_file
    if (K0 <= K)
    {
        // reading elements' isotope count into NI array

        // считываем кол-во изотопов элементов и помещаем в массив NI
        for (int i = 0; i < K0; ++i)
        {
            fscanf(inp_file, "%d", &NI[i]);
        }


        for (int i = 0; i < K0; ++i)
        {
            NII = NI[i];


            // skipping the line in inp_file

            // пропускаем строчку в inp_file
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
            }
        }
    }
    else
    {
        fclose(inp_file);
        exit(1);
    }

    fclose(inp_file);

    // With given elements' values Q_gg and Sigma (cross section log)
    // we calculate matrix elements values that are found by using 
    // Least Squares method to our data

    // По известным значениям величин Q_gg и Сигма
    // находим элементы матрицы которые получаются применением
    // метода наименьших квадратов к исходным данным
    for (int i = 1; i < K1; ++i)
    {
        // On the main matrix diagonal (sans the C[0][0])
        // amount of isotopes of each element will be written

        // на главной диагонали матрицы (Пропуская элемент C[0][0])
        // будет записано количество изотопов каждого элемента
        NII = NI[i - 1];
        C[i][i] = (double)NII;
        for (int j = 0; j < NII; ++j)
        {
            C[0][0] += X[i - 1][j] * X[i - 1][j];
            C[0][i] += X[i - 1][j];
            C[i][0] += X[i - 1][j];
            D[0]    += X[i - 1][j] * Y[i - 1][j];
            D[i]    += Y[i - 1][j];
        }
    }

    gaussMethod(C, D, Z);

    // Opening output data file "out_file"
    // if the file does not open we terminate the program

    // Открываем файл выходных данных "out_file"
    // Если файл не открывается завершаем работу программы
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


    // writing to file

    // записываем в файл
    for (int i = 0; i < K0; ++i)
    {
        fprintf(out_file, "sigma\t %d\t %d\t %.17g\n", i + 1, NI[i], sumdy[i]);
    }

    fclose(out_file);

    return 0;
}
