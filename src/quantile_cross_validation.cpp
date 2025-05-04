#include <Rcpp.h>
#include <cmath>
#include <cstring>
#include <math.h>
#include <string.h>
#include <stdio.h>

using namespace Rcpp;

double simpsons_rule(const std::vector<double>& x, const std::vector<double>& y) {
  std::vector<double> x_clean, y_clean;

  for (size_t i = 0; i < x.size(); ++i) {
    if (std::isfinite(x[i]) && std::isfinite(y[i])) {
      x_clean.push_back(x[i]);
      y_clean.push_back(y[i]);
    }
  }

  int n = x_clean.size();
  if (n < 2) return 0.0;
  if (n % 2 == 0) n -= 1; // usamos solo hasta el punto impar más cercano

  double h = (x_clean[n - 1] - x_clean[0]) / (n - 1);
  double result = y_clean[0] + y_clean[n - 1];

  for (int i = 1; i < n - 1; i += 2) result += 4 * y_clean[i];
  for (int i = 2; i < n - 2; i += 2) result += 2 * y_clean[i];

  return result * h / 3.0;
}

double kernel_function_density(double u, const char *kernel) {
    if (strcmp(kernel, "gaussian") == 0) {
        return exp(-0.5 * u * u) / sqrt(2 * M_PI);
    } else if (strcmp(kernel, "epanechnikov") == 0) {
        return (fabs(u) < 1) ? 0.75 * (1 - u * u) : 0.0;
    } else if (strcmp(kernel, "rectangular") == 0) {
        return (fabs(u) < 1) ? 0.5 : 0.0;
    } else if (strcmp(kernel, "triangular") == 0) {
        return (fabs(u) < 1) ? 1.0 - fabs(u) : 0.0;
    } else if (strcmp(kernel, "biweight") == 0) {
        return (fabs(u) < 1) ? (15.0 / 16.0) * pow(1 - u * u, 2) : 0.0;
    } else if (strcmp(kernel, "cosine") == 0) {
        return (fabs(u) < 1) ? (1 + cos(M_PI * u)) / 2.0 : 0.0;
    } else if (strcmp(kernel, "optcosine") == 0) {
        return (fabs(u) < 1) ? (M_PI / 4.0) * cos(M_PI * u / 2.0) : 0.0;
    } else {
        Rcpp::Rcerr << "unknown kernel" << std::endl;
        return NAN;
    }
}

//' Quantile cross-validation
//'
//' Performs quantile-based cross-validation for bandwidth selection.
//'
//' @param Tis Numeric matrix.
//' @param vals_minus_i Numeric matrix.
//' @param ti Numeric vector.
//' @param vals Numeric vector.
//' @param hs Numeric vector of bandwidths.
//' @param kernel Character, kernel name.
//' @return Numeric vector with cross-validation scores.
//' @export
// [[Rcpp::export]]
NumericVector quantile_cross_validation(NumericMatrix Tis,
                                        NumericMatrix vals_minus_i,
                                        NumericVector ti,
                                        NumericVector vals,
                                        NumericVector hs,
                                        std::string kernel) {
    int nh = hs.size();
    int n = Tis.nrow();
    int segs = 25;
    NumericVector cvh(nh);

    for (int k = 0; k < nh; k++) {
        cvh[k] = 0.0;

        for (int i = 0; i < n; i++) {

            // We create the grid of points for the integration
            NumericVector ti_row = Tis(i, _); // ti_row[j] = T(j) when (i) is not in the sample
            int t_len = segs * (n-1);
            NumericVector t(t_len);

            // First segment [0, T(1)]
            for (int j = 0; j < segs; j++) {
              t[j] = 0 + (ti_row[0] - 0) * j / (segs - 1);
              //Rcpp::Rcout << "t: " << t[j] << std::endl;


            }
            int idx_t = segs;
            // The other segments [T(1), 1]
            for (int m = 0; m < n - 2; m++) {
              for (int j = 0; j < segs; j++) {
                t[idx_t] = ti_row[m] + (ti_row[m + 1] - ti_row[m]) * j / (segs - 1);
                //Rcpp::Rcout << "t: " << t[idx_t] << std::endl;

                idx_t++;

              }

            }

            NumericVector aux(t_len);
            for (int j = 0; j < t_len; j++) {
                double u = (ti[i] - t[j]) / hs[k];
                //Rcpp::Rcout << "u: " << u << std::endl;
                aux[j] = kernel_function_density(u, kernel.c_str()) / hs[k];
                //Rcpp::Rcout << "aux: " << aux[j] << std::endl;
            }

            NumericVector integral_value(n-1);
            double est = 0.0;

            int n_blocks = t_len / segs;
            for (int b = 0; b < n_blocks; ++b) {
                int m = b * segs;
                std::vector<double> t_sub(t.begin() + m, t.begin() + m + segs);
                std::vector<double> aux_sub(aux.begin() + m, aux.begin() + m + segs);

                double val = simpsons_rule(t_sub, aux_sub);
                val = std::isnan(val) ? 0.0 : val;

                integral_value[b] = val;
                //Rcpp::Rcout << "integral_value: " << val << std::endl;

                est += vals_minus_i(i, b) * val;
            }


            cvh[k] += pow(vals[i] - est, 2)/n;
        }
    }
    return cvh;
}
