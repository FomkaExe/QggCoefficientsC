#include <iostream>
#include <fstream>
#include <stdexcept>
#include <vector>
#include <cmath>

// helpers --------------------------------------------------------------------
// generic aliases
template <typename T = double>
using matrix1D = std::vector<T>;
template <typename T = double>
using matrix2D = std::vector<std::vector<T>>;

using size_type = std::size_t;

// factory function for one-dimensional matrices
template <typename T = double>
matrix1D<T> make_matrix_1D(size_type _sz, T _default = T{}) {
    return matrix1D<T>(_sz, _default);
}

// factory function for two-dimensional matrices
template <typename T = double>
matrix2D<T> make_matrix_2D(size_type _szm, size_type _szn, T _default = T{}) {
    return matrix2D<T>(_szm, matrix1D<T>(_szn, _default));
}
// ----------------------------------------------------------------------------

// i/o ------------------------------------------------------------------------
constexpr std::string_view INPUT{"../input/OTa_qgg.dat"};
constexpr std::string_view OUTPUT{"../output/OTa_qgg_out.dat"};
// ----------------------------------------------------------------------------

// _C and _D are acquired by value
template <typename T = double>
matrix1D<T> gauss_method(matrix2D<T> _C, matrix1D<T> _D) {
    const int n = _C.size();

    // Forward elimination with partial pivoting
    for (int col = 0; col < n; ++col) {
        // Find pivot row
        int pivot = col;
        for (int row = col + 1; row < n; ++row) {
            if (std::abs(_C[row][col]) > std::abs(_C[pivot][col])) {
                pivot = row;
            }
        }

        // Check singularity
        if (std::abs(_C[pivot][col]) < 1e-12) {
            std::cout << "ERROR";
            return matrix1D<T>{};
        }

        // Swap rows
        if (pivot != col) {
            std::swap(_C[pivot], _C[col]);
            std::swap(_D[pivot], _D[col]);
        }

        // Eliminate below
        for (int row = col + 1; row < n; ++row) {
            T factor = _C[row][col] / _C[col][col];
            for (int j = col; j < n; ++j) {
                _C[row][j] -= factor * _C[col][j];
            }
            _D[row] -= factor * _D[col];
        }
    }

    matrix1D<T> R(n);
    // Back substitution
    for (int i = n - 1; i >= 0; --i) {
        T sum = _D[i];
        for (int j = i + 1; j < n; ++j) {
            sum -= _C[i][j] * R[j];
        }
        R[i] = sum / _C[i][i];
    }

    return R;
}

/* ERROR CODES:
 * 0 - success
 * 1 - array error
 * 2 - file open error
 * 3 - file read/write error
 * 4 - algebraic operation error
 */
auto main() -> int {
    std::ifstream input(INPUT.data());
    if (!input.is_open()) {
        std::cerr << "Ошибка при открытии файла";
        return 2;
    }

    // K0 is amount of elements in dataset (7 elements in O+Ta)
    size_type K0{};
    input >> K0;

    // array of number of rows of each elements' points
    auto NI{make_matrix_1D(K0, 0)};
    // array containing data in the first column of input file
    auto X{make_matrix_2D(K0, 0)};
    // array containing data in the second column of the input file
    auto Y{make_matrix_2D(K0, 0)};
    // Matrix of errors
    auto p{make_matrix_2D(K0, 0)};

    for (size_type i = 0; i < K0; ++i) {
        input >> NI[i];
        // Reshape the X, Y and p matrices
        X.at(i) = make_matrix_1D(NI[i]);
        Y.at(i) = make_matrix_1D(NI[i]);
        p.at(i) = make_matrix_1D(NI[i]);
    }

    // fill the X, Y and p matrices with dataset contents
    for (size_type i = 0; i < K0; ++i) {
        [[maybe_unused]] std::string tag, name;
        input >> tag >> name;
        for (size_type j = 0; j < NI[i]; ++j) {
            [[maybe_unused]] size_type ik, in;
            double xij, yij;

            input >> ik >> in >> xij >> yij;

            X.at(i).at(j) = xij;
            Y.at(i).at(j) = yij;

            p.at(i).at(j) = (i == 0)
                ? std::exp(yij + 6.0L)
                : std::exp(yij + 2.3L * (i - 1));

            // left for legacy reasons
            // p.at(i).at(j) = 1.0;
            // p.at(i).at(j) = std::exp(yij + 6.0);
        }
    }

    input.close();

    // K1 is the size of matrices used in for solving equations system
    size_type K1 = K0 + 1;
    // rightmost matrix in the system of linear equations
    auto D{make_matrix_1D(K1)};
    // leftmost matrix in the system of linear equations
    auto C{make_matrix_2D(K1, K1)};
    // result matrix of a, b1, b2, ..., bK
    auto Z{make_matrix_1D(K1)};

    // matrix population ---------------------------------------------
    for (size_type i = 1; i < K1; ++i) {
        size_type currentElem = NI[i - 1];
        for (size_type j = 0; j < currentElem; ++j) {
            C[0][0] += p[i - 1][j] * (X[i - 1][j] * X[i - 1][j]);
            C[0][i] += p[i - 1][j] * X[i - 1][j];
            C[i][0] += p[i - 1][j] * X[i - 1][j];
            C[i][i] += p[i - 1][j];

            D[0]    += p[i - 1][j] * X[i - 1][j] * Y[i - 1][j];
            D[i]    += p[i - 1][j] * Y[i - 1][j];
        }
    }
    // ---------------------------------------------------------------

    auto R = gauss_method(C, D);

    if (R.empty()) {
        std::cout << "Singular or nearly singular matrix" << std::endl;
        return 4;
    }

    // y = a * x + b
    long double a{R.front()};
    auto sigma{make_matrix_1D(K0)};

    auto sumP{0.0L};
    auto sumErr{0.0L};
    auto nev{make_matrix_1D(K0)};

    for (size_t i = 0; i < K0; ++i) {
        long double b = R.at(i + 1);

        auto err{0.0L};
        auto sumPikin{0.0L};
        for (size_t j = 0; j < NI[i]; ++j) {
            long double yc{a * X.at(i).at(j) + b};
            long double dyc{yc - Y.at(i).at(j)};

            sumPikin += p.at(i).at(j);

            if (std::abs(dyc) > 1.0) {
                sigma.at(i) += dyc*dyc / std::abs(yc);
            }

            if (std::abs(yc) > 1.0e-6L) {
                err += (p.at(i).at(j) * std::abs(dyc / Y[i][j]));
            }
        }
        sumP += sumPikin;
        sumErr += err;
        nev.at(i) = (err / K0) / sumPikin;
    }

    auto sumNev{(sumErr / K0) / sumP};

    std::ofstream output(OUTPUT.data());
    if (!output.is_open()) {
        std::cout << "Ошибка записи результатов" << std::endl;
        return 2;
    }

    output << "a\t\t\t = " << R.front() << '\n';
    for (size_type i = 1; i < R.size(); ++i) {
        output << "b(" << i << ") \t\t = "
               << R.at(i) << '\n';

        output << "sigma(" << i << ") \t = "
               << sigma.at(i - 1) << '\n';

        output << "невязка(" << i << ") \t = "
               << nev.at(i - 1) << '\n';
    }

    output << "невязка общая \t= " << sumNev;

    return 0;
}
