#include <iostream>
#include <iomanip>
#include <cmath>




double I(double theta, double a, double b, double c);
double midpoint_rule(double lb, double ub, double stepsize, double theta, double a, double b, double c);
double theta_ast(double a, double b, double c);
int main()
{
  double theta, a, b, c, lb, ub, stepsize, se;
  theta = 5;
  std::cout << "This program uses the midpoint numerical analysis method to \n"
            << "integrate from x0 to x1 of the Information function.\n";

  std::cout << "\na: ";
  std::cin >> a;
  std::cout << "\nb: ";
  std::cin >> b;
  std::cout << "\nc: ";
  std::cin >> c;
  std::cout << "What step size to use for integration? The smaller the more accurate\n"
            << "Step size: ";
  std::cin >> stepsize;
  std::cout << "\nstandard error: ";
  std::cin >> se;
  lb = theta - 2*se;
  ub = theta + 2*se;
  std::cout << "Theta asterisk: " << theta_ast(a, b, c) << "\n";
   std::cout << "Weight of theta(" ") " <<  std::setprecision(15) << midpoint_rule(lb, ub, stepsize, theta, a, b, c) << "\n";
  return 0;
}


double I(double theta, double a, double b, double c)
{
  double numer = (2.89 * std::pow(a, 2)) * (1 - c);
  double denom1 = c + std::exp( (1.7*a)*(theta - b));
  double denom2 = std::pow(1 + std::exp((-1.7*a)*(theta-b)), 2);
  return (numer)/(denom1 * denom2);
}

// lb = lower bound
// ub = upper bound
// step size is the size of delta x
double midpoint_rule(double lb, double ub, double count, double theta, double a, double b, double c)
{
  long double integral = 0;
  long double step = (ub - lb) / count;
  for (unsigned i = 1; i <= count; ++i) {
    integral += step * I(lb + (i-1) * step, a, b, c);
  }
  return  integral;
}

double theta_ast(double a, double b, double c)
{
    return (b + (1/1.7) * std::log( ( 1 + std::sqrt(1 + 8*c))/2));
}
