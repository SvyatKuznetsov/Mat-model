#include <iostream>

double u(double x, int type)
{
    switch (type)
    {
        case 0:
            return 1;
        case 1:
            return x + 1;
        case 2:
            return x * x;
        default:
            return 0;
    }
}

double k(double x, int type)
{
    switch (type)
    {
        case 0:
            return 1;
        case 1:
            return x + 2;
        case 2:
            return x + 1;
        default:
            return 0;
    }
}

double hi1()
{
    return 0;
}

double q(double x, int type)
{
    switch (type)
    {
        case 0:
            return 1;
        case 1:
            return 1;
        case 2:
            return 1;
        default:
            return 0;
    }
}

double fx(double x, int type)
{
    switch (type)
    {
        case 0:
            return 1;
        case 1:
            return x;
        case 2:
            return x * x - 4 * x - 2;
        default:
            return 0;
    }
}

double du_dx(double x, int type)
{
    switch (type)
    {
        case 0:
            return 0;
        case 1:
            return 1;
        case 2:
            return 2 * x;
        default:
            return 0;
    }
}

double nu1(double x, int type)
{
    return -k(x, type) * du_dx(x, type) + hi1() * u(x, type);
}

double nu2(double x, int type)
{
    return u(x, type);
}

int main()
{
    const double left = 0;
    const double right = 2;
    int n = 4;
    const int type = 0;
    double arr4[5], arr8[9], arr16[17], arr32[33], arr64[65], arr128[129], arr256[257], arr512[513];

    for (n = 4; n <= 512; n *= 2)
    {
        double h = (right - left) / n;
        double h_ = h / 2;
        double x[n + 1], x_[n], a[n + 1], b[n + 1], c[n + 1], u_[n + 1], v[n + 1], alpha[n], betta[n], g[n + 1], f[n + 1], fArr;
        for (int i = 0; i < n + 1; ++i)
        {
            x[i] = left + i * h;
        }

        for (int i = 0; i < n; ++i)
        {
            x_[i] = (x[i] + x[i + 1]) / 2;
        }

        for (int i = 0; i < n + 1; ++i)
        {
            u_[i] = u(x[i], type);
        }

        a[0] = 0;
        b[0] = -k(x_[0], type) / h;
        c[0] = -b[0] + hi1() + h_ * q(x[0], type);

        for (int i = 1; i < n; ++i)
        {
            a[i] = -k(x_[i - 1], type) / h;
            b[i] = -k(x_[i], type) / h;
            c[i] = -a[i] - b[i] + h * q(x[i], type);
        }

        for (int i = 0; i < n + 1; ++i)
        {
            f[i] = fx(x[i], type);
        }

        a[n] = 0;
        b[n] = 0;
        c[n] = 1;

        g[0] = h_ * f[0] + nu1(x[0], type);
        for (int i = 1; i < n; ++i)
        {
            g[i] = h * f[i];
        }
        g[n] = nu2(x[n], type);

        alpha[1] = -b[0] / c[0];
        betta[1] = g[0] / c[0];

        for (int i = 1; i < n; ++i)
        {
            alpha[i + 1] = -b[i] / (a[i] * alpha[i] + c[i]);
            betta[i + 1] = (g[i] - a[i] * betta[i]) / (a[i] * alpha[i] + c[i]);
        }

        v[n] = (g[n] - a[n] * betta[n]) / (a[n] * alpha[n] + c[n]);
        for (int i = n - 1; i >= 0; --i)
        {
            v[i] = alpha[i + 1] * v[i + 1] + betta[i + 1];
        }

        fArr = std::abs(u_[0] - v[0]);
        std::cout << "n: " << n;
        for (int i = 1; i < n + 1; ++i)
        {
            if (std::abs(u_[i] - v[i]) > fArr)
            {
                fArr = std::abs(u_[i] - v[i]);
            }
        }
        std::cout << "\nmaxErr: " << fArr << "\n";

        std::cout << "\n";
    }

    return 0;
}