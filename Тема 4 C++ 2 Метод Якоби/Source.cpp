#include <iostream>
#include <vector>
#include <cmath>
#include <stdexcept>

using namespace std;

vector<double> jacobi(vector<vector<double>> a, vector<double> b, double tol = 1e-2, int max_iterations = 100) {
    int n = b.size();
    vector<double> x(n, 0.0);
    vector<double> x_new(n, 0.0);

    for (int k = 0; k < max_iterations; ++k) {
        for (int i = 0; i < n; ++i) {
            double s1 = 0.0;
            for (int j = 0; j < n; ++j) {
                if (j != i) {
                    s1 += a[i][j] * x[j];
                }
            }
            x_new[i] = (b[i] - s1) / a[i][i];
        }

        bool convergence = true;
        for (int i = 0; i < n; ++i) {
            if (fabs(x_new[i] - x[i]) > tol) {
                convergence = false;
                break;
            }
        }

        if (convergence) {
            break;
        }

        x = x_new;
    }

    return x;
}

int main() {
    vector<vector<double>> a = { {5, 3, -2},
                                {2, 1, -1},
                                {3, -2, -3} };

    vector<double> b = { -1, 0, 2 };

    vector<double> solution = jacobi(a, b);

    cout << "–ешение методом якоби: ";
    for (double x : solution) {
        cout << x << " ";
    }
    cout << endl;

    return 0;
}
