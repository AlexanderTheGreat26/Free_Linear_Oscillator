#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <string>
#include <array>
#include <sstream>


const double l = 0.1;
const double w_0 = 1.0;
const double step = 0.1;
const double t_0 = 0;
const double t_N = 4*acos(-1) / std::sqrt(std::pow(w_0, 2) - std::pow(l, 2)); // two periods
const double x_0 = 0;
const double v_0 = 1;


typedef std::array<double, 4> coeff;
typedef std::vector<double> grid;


void Runge_Kutta_method (grid & velocities, grid & coordinates, const double & h);

void Euler_Cauchy_method (grid & velocities, grid & coordinates, const double & h);

void method_implementation (grid & times, grid coordinates, grid velocities, const double & mesh_step,
                            const std::string & method_name, void method(grid & velocities, grid & coordinates, const double & h));

grid mesh (double left_border, const double & right_border, const double & mesh_step);


int main() {
    grid coordinates {x_0}, velocities {v_0};
    grid times = std::move(mesh(t_0, t_N, step));
    method_implementation(times, coordinates, velocities, step, "RKM", Runge_Kutta_method);
    method_implementation(times, coordinates, velocities, step, "PCM", Euler_Cauchy_method);
    return 0;
}


double velocity (double x, double v) {
    return -2.0*l*v - std::pow(w_0, 2) * x;
}


double coordinate (double x, double v) {
    return v;
}

void coefficients (const double & h, double & x_n, double & v_n, coeff & k, coeff & q,
                   double f (double x, double y), double g (double x, double y)) {
    k[0] = g(x_n, v_n);
    q[0] = f(x_n, v_n);
    k[1] = g(x_n + h * q[0] / 2.0, v_n + h * k[0] / 2.0);
    q[1] = f(x_n + h * q[0] / 2.0 ,v_n + h * k[0] / 2.0);
    k[2] = g(x_n + h * q[1] / 2.0, v_n + h * k[1] / 2.0);
    q[2] = f(x_n + h * q[1] / 2.0, v_n + h * k[1] / 2.0);
    k[3] = g(x_n + h * q[2], v_n + h * k[2]);
    q[3] = f(x_n + h * q[2], v_n + h * k[2]);
}


void Runge_Kutta_method (grid & velocities, grid & coordinates, const double & h) {
        coeff k, q;
        long n = velocities.size()-1;
        double x_n = coordinates[n];
        double v_n = velocities[n];
        coefficients(h, x_n, v_n, k, q, coordinate, velocity);
        coordinates.emplace_back(x_n + h/6.0*(q[0] + 2.0*q[1] + 2.0*q[2] + q[3]));
        velocities.emplace_back(v_n + h/6.0*(k[0] + 2.0*k[1] + 2.0*k[2] + k[3]));
}


double predictor (double & x_n, double & y_n, const double & h, double f(double x, double y)) {
    return y_n + h * f(x_n, y_n);
}

double corrector (double & x_n, double & y_n, double & y_pred, const double & h, double f(double x, double y)) {
    return y_n + h/2.0*(f(x_n, y_n) + f(x_n, y_pred));
}


void Euler_Cauchy_method (grid & velocities, grid & coordinates, const double & h) {
    long n = velocities.size()-1;
    double x_n = coordinates[n];
    double v_n = velocities[n];
    double v_pred = predictor(x_n, v_n, h, velocity);
    double v_new = corrector (x_n, v_n, v_pred, h, velocity);
    double x_new = x_n + h * v_new;
    velocities.emplace_back(v_new);
    coordinates.emplace_back(x_new);
}


//There only technical functions below.
grid mesh (double left_border, const double & right_border, const double & mesh_step) {
    std::vector <double> xx ((right_border-left_border) / mesh_step + 1);
    xx[0] = left_border;
    std::generate(xx.begin()+1, xx.end(), [&] {left_border += step; return left_border;});
    return xx;
}


template <typename T>
std::string toString (T val) {
    std::ostringstream oss;
    oss << val;
    return oss.str();
}


// Creates text file (.txt) from with coordinates (x, y, z) for plotting.
void data_file_creation (const std::string & name, grid & x, grid & y, grid & z) {
    std::ofstream fout;
    fout.open(name + ".txt", std::ios::out | std::ios::trunc);
    for (int i = 0; i < x.size(); ++i)
        fout << toString(x[i]) << '\t' << toString(y[i]) << '\t' << toString(z[i]) << '\n';
    fout.close();
}


void plot (const std::string & name, const std::string & title, const std::string & xlabel, const std::string & ylabel,
           const std::string & x_column, const std::string & y_column) {
    FILE *gp = popen("gnuplot -persist", "w");
    if (!gp) throw std::runtime_error("Error opening pipe to GNUplot.");
    std::vector<std::string> stuff = {"set term jpeg size 1920, 1080 font \"Helvetica,30\"",
                                      "set output \'" + title + ".jpg\'",
                                      "set title \'" + title + "\'",
                                      "set grid xtics ytics",
                                      "set autoscale xfix",
                                      "set xlabel \'" + xlabel + "\'",
                                      "set ylabel \'" + ylabel + "\'",
                                      "set key off",
                                      "set ticslevel 0",
                                      "set border 4095",
                                      "plot \'" + name + ".txt\' using " + x_column + ':' + y_column + " w lines lw 4",
                                      "set tics font \'Helvetica,20\'",
                                      "set terminal pop",
                                      "set output",
                                      "replot", "q"};
    for (const auto & it : stuff)
        fprintf(gp, "%s\n", it.c_str());
    pclose(gp);
}


void method_implementation (grid & times, grid coordinates, grid velocities, const double & mesh_step,
                            const std::string & method_name, void method(grid & velocities, grid & coordinates, const double & h)) {
    for (int i = 0; i < times.size(); ++i)
        method(velocities, coordinates, mesh_step);
    data_file_creation(method_name, times, coordinates, velocities);
    plot(method_name, "Solution " + method_name, "t", "x", "1", "2");
    plot(method_name, "Phase portrait " + method_name, "x", "v", "2", "3");
}