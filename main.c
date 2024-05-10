// calculation of Qgg coefficients (Q-thermal energy, g-ground state, Qgg -defference of ground state energies of the colliding ions)
// also called as elements lines
// for each element produced in the nuclear reactions cross-sections of different isiotopes of each element 'S' form a line in Q-S plane
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

    // +
    for (int i = 0; i < K1; ++i)
    {
        indD[i] = i;
        for (int j = 0; j < K1; ++j)
        {       
            cnew[i][j] = 0.0;
        }
    }

    // 
    for (int j = 0; j < K; ++j) {
        int i0 = j;
        int iflag = -1;
        int j1 = j + 1;
        double cmax = c1[j][j];

        for (int is = j1; is < K1; ++is) {
            if (fabs(cmax) < fabs(c1[is][j])) {
                i0 = is;
                cmax = c1[is][j];
                iflag = 1;
            }
        } //+
        if (iflag == 1) {
            for (int jc = 0; jc < K1; ++jc) {
                cold[jc] = c1[j][jc];
                c1[j][jc] = c1[i0][jc];
                c1[i0][jc] = cold[jc];
            }
            double dold = d1[j];
            d1[j] = d1[i0];
            d1[i0] = dold;

            double indDold = indD[j];
            indD[j] = indD[i0];
            indD[i0] = indDold;
        }
        for (int i = j1; i < K1; ++i) {
            q[i][j] = c1[i][j]/c1[j][j];
            cnew[i][j] = q[i][j];
            d1[i] -= q[i][j]*d1[j];
            for (int j2 = j1; j2 < K1; ++j2) {
                c1[i][j2] -= q[i][j]*c1[j][j2];
            }
        }
    }

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
    // const int K = 7;
    // const int K1 = 8;
    // const int N = 11;
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

    // zero initialising arrays

    // инициализируем массивы нулями
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

    // we should always close the file after using it

    // после использования файла его всегда нужно закрывать
    fclose(inp_file);

    // with given elements' values Q_gg and Sigma (cross section log)
    // we calculate matrix elements values that are found by using 
    // Least Squares method to our data

    // По известным значениям величин Q_gg и Сигма
    // находим элементы матрицы которые получаются применением
    // метода наименьших квадратов к исходным данным
    for (int i = 1; i < K1; ++i)
    {
        // на главной диагонали матрицы (С элемента C[1][1])
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
    // Calling gauss method

    // Вызываем метод гаусса
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
        fprintf(out_file, "%d %lf\n", i + 1, Z[i]);
    }


    double a = Z[0];

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
        fprintf(out_file, "sigma %d %d %lf\n", i, NI[i], sumdy[i]);
    }

    return 0;
}
