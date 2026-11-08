// SEIQR COVID-19 model solved with the fourth-order Runge-Kutta method (RK4)
#include <iostream>
#include <fstream>
using namespace std;

// parameters (Case 1: R0 > 1). For Case 2 change beta to 1e-5.
const double A = 100, mu = 0.02, beta = 4e-5;
const double alpha = 0.1, gamma_ = 0.05, sigma = 0.08, theta = 0.07;
const double k1 = alpha + mu + gamma_, k2 = sigma + mu, k3 = theta + mu;

// right-hand side of the five equations: dy/dt = f(y)
void f(const double y[5], double dy[5]) {
    double S = y[0], E = y[1], I = y[2], Q = y[3], R = y[4];
    double inf = beta * S * (E + I);
    dy[0] = A - mu * S - inf;
    dy[1] = inf - k1 * E;
    dy[2] = alpha * E - k2 * I;
    dy[3] = gamma_ * E + sigma * I - k3 * Q;
    dy[4] = theta * Q - mu * R;
}

int main() {
    double N = A / mu;
    double R0 = beta * N * (1.0 / k1 + alpha / (k1 * k2));
    cout << "R0 = " << R0 << endl;

    double y[5] = {4990, 10, 0, 0, 0};
    double h = 0.1, T = 365;
    int steps = (int)(T / h + 0.5);

    ofstream out("seiqr_rk4_output.csv");
    out << "t,S,E,I,Q,R\n";
    double peakI = 0, dayPeak = 0;

    for (int n = 0; n <= steps; n++) {
        double t = n * h;
        out << t << "," << y[0] << "," << y[1] << "," << y[2] << "," << y[3] << "," << y[4] << "\n";
        if (y[2] > peakI) { peakI = y[2]; dayPeak = t; }

        // one RK4 step: four slopes a, b, c, d
        double a[5], b[5], c[5], d[5], tmp[5];
        f(y, a);
        for (int i = 0; i < 5; i++) tmp[i] = y[i] + 0.5 * h * a[i];
        f(tmp, b);
        for (int i = 0; i < 5; i++) tmp[i] = y[i] + 0.5 * h * b[i];
        f(tmp, c);
        for (int i = 0; i < 5; i++) tmp[i] = y[i] + h * c[i];
        f(tmp, d);
        for (int i = 0; i < 5; i++) y[i] += h / 6.0 * (a[i] + 2 * b[i] + 2 * c[i] + d[i]);
    }
    cout << "Peak infected = " << peakI << " on day " << dayPeak << endl;
    out.close();
    return 0;
}
