// SEIQR COVID-19 model solved with Euler's method (basic C++)
#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

int main() {
    // parameters (Case 1: R0 > 1). For Case 2 change beta to 1e-5.
    double A = 100, mu = 0.02, beta = 4e-5;
    double alpha = 0.1, gamma = 0.05, sigma = 0.08, theta = 0.07;
    double k1 = alpha + mu + gamma, k2 = sigma + mu, k3 = theta + mu;

    // basic reproduction number
    double N = A / mu;
    double R0 = beta * N * (1.0 / k1 + alpha / (k1 * k2));
    cout << "R0 = " << R0 << endl;

    // initial values and step size
    double S = 4990, E = 10, I = 0, Q = 0, R = 0;
    double h = 0.1, T = 365;
    int steps = (int)(T / h);

    ofstream out("seiqr_output.csv");
    out << "t,S,E,I,Q,R\n";
    double peakI = 0, dayPeak = 0;

    for (int n = 0; n <= steps; n++) {
        double t = n * h;
        out << t << "," << S << "," << E << "," << I << "," << Q << "," << R << "\n";
        if (I > peakI) { peakI = I; dayPeak = t; }

        // rates of change
        double inf = beta * S * (E + I);
        double dS = A - mu * S - inf;
        double dE = inf - k1 * E;
        double dI = alpha * E - k2 * I;
        double dQ = gamma * E + sigma * I - k3 * Q;
        double dR = theta * Q - mu * R;

        // Euler update
        S += h * dS; E += h * dE; I += h * dI; Q += h * dQ; R += h * dR;
    }
    cout << "Peak infected = " << peakI << " on day " << dayPeak << endl;
    out.close();
    return 0;
}
