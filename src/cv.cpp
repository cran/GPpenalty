#include <Rcpp.h>
#include <random> // for std::default_random_engine
#include <algorithm> // for std::shuffle
// #include <armadillo>
#include <cmath>

using namespace std;
using namespace Rcpp;
// using namespace arma;


// in R code, "fold <- which(random.folds==j)" part
IntegerVector which_equal(IntegerVector vec, int j) {
  IntegerVector indices;
  for(int i=0; i<vec.size(); i++) {
    if(vec[i] == j) {
      indices.push_back(i);
    }
  }
  return indices;
}


// row means of matrix
NumericVector row_means(const NumericMatrix& mat) {
  int n_rows = mat.nrow(), n_cols = mat.ncol();
  NumericVector means(n_rows);
  for(int i=0; i<n_rows; i++) {
    double total = 0;
    for(int j=0; j<n_cols; j++) {
      total = total + mat(i, j);
    }
    means[i] = total / n_cols;
  }
  return means;
}




// [[Rcpp::export]]
List cv(NumericVector y, NumericMatrix x, NumericVector lambda, NumericMatrix initialvals, int n=-1, int k=-1, int p=-1, int dim=-1,
        bool profile=true, bool mu=false, bool g=false, double fixed_g=-1, bool scad=false, bool penalty=true, int ncores=1) {

  // R functions
  Rcpp::Function multi("multi", Rcpp::Environment::namespace_env("GPpenalty"));
  Rcpp::Function kriging("kriging", Rcpp::Environment::namespace_env("GPpenalty"));
  Function score("score");


  if (k==-1) { // if k is not specified leave one out cv
    k = y.size();
  }
  if (n==-1) {
    n = y.size();
  }
  if (p==-1) {
    p = x.ncol();
  }
  // folds
  IntegerVector folds(n);
  for(int i=0; i<n; i++) {
    folds[i] = (i % k) + 1;
  }
  // random folds
  random_device rd;// seed
  mt19937 gen(rd());
  shuffle(folds.begin(), folds.end(), gen);

  // for score and mse results
  NumericMatrix score_matrix(lambda.size(), k);
  NumericMatrix mse_matrix(lambda.size(), k);

  for(int j=1; j<=k; j++) {
    IntegerVector valid_indice = which_equal(folds, j); // validation data index
    IntegerVector train_indice;

    for(int i=0; i<n; i++) { // training data index
      if(folds[i] !=j) {
        train_indice.push_back(i);
      }
    }
    // training and validation data vectors and matrices
    NumericVector y_valid(valid_indice.size());
    NumericVector y_train(train_indice.size());
    NumericMatrix x_valid(valid_indice.size(), p);
    NumericMatrix x_train(train_indice.size(), p);

    // train
    for(int i=0; i<train_indice.size(); i++) {
      y_train[i] = y[train_indice[i]];
      for(int col=0; col<p; col++) {
        x_train(i, col) = x(train_indice[i], col);
      }
    }
    // valid
    for(int i=0; i<valid_indice.size(); i++) {
      y_valid[i] = y[valid_indice[i]];
      for(int col=0; col<p; col++) {
        x_valid(i, col) = x(valid_indice[i], col);
      }
    }

    List multi_list = List::create(
      Named("y")=y_train,
      Named("x")=x_train,
      Named("dim")=dim,
      Named("mu")=mu,
      Named("g")=g,
      Named("fixed_g")=fixed_g,
      Named("profile")=profile,
      Named("scad")=scad,
      Named("initialvals")= initialvals,
      Named("lambda")=lambda,
      Named("ncores")=ncores
    );
    // estimation
    List out = multi(multi_list, Named("penalty")=true);
    // estimated values
    NumericMatrix theta = out["theta"];
    NumericVector s2 = out["s2"];
    NumericVector mu_vec;
    NumericVector g_vec;

    if (mu) mu_vec = out["mu"];
    if (g) g_vec = out["g"];


    for(int i=0; i<lambda.size(); i++) {
      if (i >= theta.nrow()) stop("Index i out of bounds for theta");
      if (i >= s2.size()) stop("Index i out of bounds for s2");
      List pred;
      if (mu) {
        if (g) {
          pred = kriging(Named("y")=y_train,
                         Named("x")=x_train,
                         Named("xx")=x,
                         Named("sigma2")=s2[i],
                         Named("theta")=theta(i,_),
                         Named("mu")=mu_vec[i],
                         Named("g")=g_vec[i]);
        } else {
          pred = kriging(Named("y")=y_train,
                         Named("x")=x_train,
                         Named("xx")=x,
                         Named("sigma2")=s2[i],
                         Named("theta")=theta(i,_),
                         Named("mu")=mu_vec[i]);
        }
      } else {
        if (g) {
          pred = kriging(Named("y")=y_train,
                         Named("x")=x_train,
                         Named("xx")=x,
                         Named("sigma2")=s2[i],
                         Named("theta")=theta(i,_),
                         Named("g")=g_vec[i]);
        } else {
          pred = kriging(Named("y")=y_train,
                         Named("x")=x_train,
                         Named("xx")=x,
                         Named("sigma2")=s2[i],
                         Named("theta")=theta(i,_));
        }
      }

      NumericVector pred_error = y-as<NumericVector>(pred["mup"]);
      double mse_value = sum(pow(pred_error, 2.0));

      // mse and score
      mse_matrix(i, j-1) = mse_value;

      if(k != n) {
        score_matrix(i, j-1) = as<double>(
          score(Named("y")=y,
                Named("mu")=pred["mup"],
                Named("sigma")=pred["Sigmap"])
        );
      }

    }
  }

  // max score, 1se score, min mse, min 1se
  NumericVector score_means, mse_means;
  int score_max, score_1se, mse_min, mse_1se;

  if (k != n) { // if k fold, return both score and mse
    // score
    score_means = row_means(score_matrix);
    score_max = max_element(score_means.begin(), score_means.end()) - score_means.begin(); // need +1

    // standard error
    double score_sd = 0.0;
    for(int j=0; j<k; j++) {
      double score_val = score_matrix(score_max, j);
      score_sd += pow(score_val - score_means[score_max], 2);
    }
    score_sd = sqrt(score_sd/k); // standard deviation
    double one_se = score_means[score_max] - score_sd / sqrt(k);

    // find score_1se index
    NumericVector se_indices;
    for(int i=0; i<score_means.size(); i++) {
      if(score_means[i] >= one_se) {
        se_indices.push_back(i);
      }
    }
    // grab the largest index value
    if(se_indices.size() > 0) {
      score_1se = se_indices[se_indices.size()-1]; // need +1
    } else {
      score_1se = score_max;
    }

    // mse
    mse_means = row_means(mse_matrix);
    mse_min = min_element(mse_means.begin(), mse_means.end()) - mse_means.begin(); // need +1

    // standard deviation
    double mse_sd = 0.0;
    for(int j=0; j<k; j++) {
      double mse_val = mse_matrix(mse_min, j);
      mse_sd += pow(mse_val - mse_means[mse_min], 2);
    }
    mse_sd = sqrt(mse_sd/k); // standard deviation
    double mse_one_se = mse_means[mse_min] + mse_sd / sqrt(k);

    // find mse_1se index
    NumericVector mse_se_indices;
    for(int i=0; i<mse_means.size(); i++) {
      if(mse_means[i] <= mse_one_se) {
        mse_se_indices.push_back(i);
      }
    }
    // grab the largest index value
    if(mse_se_indices.size() > 0) {
      mse_1se = mse_se_indices[mse_se_indices.size()-1]; // need +1
    } else {
      mse_1se = mse_min;
    }
    return List::create(
      _["folds"] = folds,
      _["score"] = score_means,
      _["mse"] = mse_means,
      _["score_max"] = score_max + 1,
      _["mse_min"] = mse_min + 1,
      _["score_1se"] = score_1se + 1,
      _["mse_1se"] = mse_1se + 1
    );
  } else { // if leave one out, return only mse
    // mse
    mse_means = row_means(mse_matrix);
    mse_min = min_element(mse_means.begin(), mse_means.end()) - mse_means.begin(); // need +1

    // standard deviation
    double mse_sd = 0.0;
    for(int j=0; j<k; j++) {
      double mse_val = mse_matrix(mse_min, j);
      mse_sd += pow(mse_val - mse_means[mse_min], 2);
    }
    mse_sd = sqrt(mse_sd/k); // standard deviation
    double mse_one_se = mse_means[mse_min] + mse_sd / sqrt(k);

    // find mse_1se index
    NumericVector mse_se_indices;
    for(int i=0; i<mse_means.size(); i++) {
      if(mse_means[i] <= mse_one_se) {
        mse_se_indices.push_back(i);
      }
    }
    // grab the largest index value
    if(mse_se_indices.size() > 0) {
      mse_1se = mse_se_indices[mse_se_indices.size()-1]; // need +1
    } else {
      mse_1se = mse_min;
    }
    return List::create(
      _["folds"] = folds,
      //_["score"] = score_means,
      _["mse"] = mse_means,
      //_["score_max"] = score_max + 1,
      _["mse_min"] = mse_min + 1,
      //_["score_1se"] = score_1se + 1,
      _["mse_1se"] = mse_1se + 1
    );
  }
}
