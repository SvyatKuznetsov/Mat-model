#include <iostream>
#include <cmath>
#include <string>
#include <iomanip>

const int N = 16;
const int INIT_M = 8;
const int hi = 1;
const double T = 10;
const double L = 0;
const double R = 3;

double u(double x, double t, int type)
{
    switch (type)
    {
        case 0:
            return 1;
        case 1:
            return x + 1;
        case 2:
            return x * x;
        case 3:
            return x + 1;
        case 4:
            return x * x;
        case 5:
            return x * std::exp(-3 * t);
        default:
            return -1;
    }
}

double k(double x, double t, int type)
{
    switch (type)
    {
        case 0:
            return 1;
        case 1:
            return x + 2;
        case 2:
            return x + 1;
        case 3:
            return x + 2;
        case 4:
            return x + 1;
        case 5:
            return (x + 1) * std::exp(-3 * t);
        default:
            return -1;
    }
}

double q(double x, double t, int type)
{
    switch (type)
    {
        case 0:
            return 1;
        case 1:
            return 2 * x;
        case 2:
            return 2 * x;
        case 3:
            return 2 * x;
        case 4:
            return 2 * x;
        case 5:
            return 2 * x + t;
        default:
            return -1;
    }
}

double f(double x, double t, int type)
{
    switch (type)
    {
        case 0:
            return 1;
        case 1:
            return 2 * x * x + 2 * x - 1;
        case 2:
            return 2 * x * x * x - 4 * x - 2;
        case 3:
            return 2 * x * x + 2 * x - 1;
        case 4:
            return 2 * x * x * x - 4 * x - 2;
        case 5:
            return std::exp(-3 * t) * (-3 * x - std::exp(-3 * t) + 2 * x * x + t * x);
        default:
            return -1;
    }
}

double phi(double x, double t, int type)
{
    switch (type)
    {
        case 0:
            return 1;
        case 1:
            return x + 1;
        case 2:
            return x * x;
        case 3:
            return 1;
        case 4:
            return 1;
        case 5:
            return x * std::exp(-3 * t);
        default:
            return -1;
    }
}

double nu1(double t, int type)
{
    switch (type)
    {
        case 0:
            return 1;
        case 1:
            return -1;
        case 2:
            return 0;
        case 3:
            return -1;
        case 4:
            return 0;
        case 5:
            return -std::exp(-6 * t);
        default:
            return 0;
    }
}

double nu2(double t, int type)
{
    switch (type)
    {
        case 0:
            return 1;
        case 1:
            return 4;
        case 2:
            return 9;
        case 3:
            return 4;
        case 4:
            return 9;
        case 5:
            return 3 * std::exp(-3 * t);
        default:
            return -1;
    }
}

void splitX(double *x, int N, double xL, double h)
{
    for (int i = 0; i < N + 1; ++i)
    {
        x[i] = xL + i * h;
    }
}

void splitXHelp(double *xHelp, const double *x, int N)
{
    for (int i = 0; i < N; ++i)
    {
        xHelp[i] = (x[i + 1] + x[i]) / 2;
    }
}

void Explicit(const double *prev_v, double *v, double *x, double *xHelp, double t, double tay, double h, int type)
{
    v[0] = prev_v[0] + tay * (k(xHelp[0], t, type) * (prev_v[1] - prev_v[0]) / h / (h / 2)
            - (hi * prev_v[0] - nu1(t, type)) / (h / 2) - q(x[0], t, type) * prev_v[0] + f(x[0], t, type));
    for (int i = 1; i < N; ++i)
    {
        v[i] = prev_v[i] + tay * (k(xHelp[i], t, type) * (prev_v[i + 1] - prev_v[i]) / (h * h)
                - k(xHelp[i - 1], t, type) * (prev_v[i] - prev_v[i - 1]) / (h * h)
                - q(x[i], t, type) * prev_v[i] + f(x[i], t, type));
    }
    v[N] = nu2(t, type);
}

void run(const double *a, const double *c, const double *b, const double *g, int n, double *v)
{
    double alpha[n + 1];
    double beta[n + 1];
    alpha[1] = -b[0] / c[0];
    beta[1] = g[0] / c[0];
    for (int i = 1; i < n; ++i)
    {
        alpha[i + 1] = -b[i] / (a[i] * alpha[i] + c[i]);
        beta[i + 1] = (g[i] - a[i] * beta[i]) / (a[i] * alpha[i] + c[i]);
    }
    v[n] = (g[n] - a[n] * beta[n]) / (a[n] * alpha[n] + c[n]);
    for (int i = n - 1; i >= 0; --i)
    {
        v[i] = alpha[i + 1] * v[i + 1] + beta[i + 1];
    }
}

void Implicit(const double *prev_v, double *v, const double *x, const double *xHelp, double t, double tay, double h, int type)
{
    double *a = new double[N + 1];
    double *c = new double[N + 1];
    double *b = new double[N + 1];
    double *F = new double[N + 1];
    double halfH = h / 2;

    t = t + tay;
    F[0] = halfH * f(x[0], t, type) + nu1(t, type) + halfH / tay * prev_v[0];
    for (int i = 1; i < N; ++i)
    {
        F[i] = h * f(x[i], t, type) + h / tay * prev_v[i];
    }
    F[N] = nu2(t, type);
    a[0] = 0;
    c[0] = halfH / tay + k(xHelp[0], t, type) / h + hi + halfH * q(x[0], t, type);
    b[0] = -k(xHelp[0], t, type) / h;
    for (int i = 1; i < N; ++i)
    {
        a[i] = -k(xHelp[i - 1], t, type) / h;
        c[i] = h / tay + k(xHelp[i], t, type) / h
                + k(xHelp[i - 1], t, type) / h + h * q(x[i], t, type);
        b[i] = -k(xHelp[i], t, type) / h;
    }
    a[N] = 0;
    c[N] = 1;
    b[N] = 0;
    run(a, c, b, F, N, v);
    delete[] F;
    delete[] b;
    delete[] c;
    delete[] a;
}

double maxErr(double *v, int n, int m, double *x, double tay, int type)
{
    double mErr = std::abs(u(x[0], m * tay, type) - v[0]);
    for (int i = 1; i < n; ++i)
    {
        if (mErr < std::abs(u(x[i], m * tay, type) - v[i]))
        {
            mErr = std::abs(u(x[i], m * tay, type) - v[i]);
        }
    }
    return mErr;
}

int main()
{
    int method = 1;
    for (int j = 0; j < 2; ++j)
    {
        if (method == 1)
        {
            std::cout << "\nExplicit\n";
        }
        else
        {
            std::cout << "\nImplicit\n";
        }
        const double h = (R - L) / N;
        const int type = 1;
        double *x = new double[N + 1];
        double *xHelp = new double[N];
        splitX(x, N, L, h);
        splitXHelp(xHelp, x, N);
        int bound = std::pow(2, 15);
        for (int M = INIT_M; M <= bound; M *= 2)
        {
            const double tay = T / M;
            double v[N + 1];
            double vNext[N + 1];
            for (int i = 0; i < N + 1; ++i)
            {
                v[i] = phi(i * h, 0, type);
            }
            for (int m = 0; m < M; ++m)
            {
                if (method == 1)
                {
                    Explicit(v, vNext, x, xHelp, tay * m, tay, h, type);
                }
                else
                {
                    Implicit(v, vNext, x, xHelp, tay * m, tay, h, type);
                }
                std::copy(vNext, vNext + N + 1, v);
            }
            double mErr = maxErr(v, N, M, x, tay, type);
            std::cout << std::setw(10) << "M = " << M << "\t";
            if (_isnan(mErr))
            {
                std::cout << std::setw(20) << "inf" << "\n";
            }
            else
            {
                std::cout << std::setw(20) << mErr << "\n";
            }
        }
        ++method;
    }
    return 0;
}