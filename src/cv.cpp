//#include <Rcpp.h>
#include <random> // for std::default_random_engine
#include <algorithm> // for std::shuffle

#include <cmath>
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace arma;
using namespace std;
using namespace Rcpp;



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

// how many are inside the range of upper and lower bounds for md
NumericVector md_sums_range(const NumericMatrix& mat, double md_upper, double md_lower) {
  int n_rows = mat.nrow(), n_cols = mat.ncol();
  NumericVector md_num(n_rows);
  for(int i=0; i<n_rows; i++) {
    int count = 0;
    for(int j=0; j<n_cols; j++) {
      double md_value = mat(i, j);
      if(md_value >= md_lower && md_value <= md_upper) {
        count++;
      }
    }
    md_num[i] = count;
  }
  return(md_num);
}

// [[Rcpp::export]]
List cv(NumericVector y, NumericMatrix x, NumericVector lambda, NumericMatrix initialvals, int n=-1, int k=-1, int p=-1, int d=-1,
        bool profile=true, bool mu=false, bool g=false, double fixed_g=-1, bool scad=false, double theta_upper=1000, double theta_lower=0.001,
        bool penalty=true, int ncores=1) {

  // R functions
  Rcpp::Function multi("multi", Rcpp::Environment::namespace_env("GPpenalty"));
  Rcpp::Function kriging("kriging", Rcpp::Environment::namespace_env("GPpenalty"));
  Function score("score");
  Function dpe("dpe");


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

  // for score, mse, dpe, and md results
  NumericMatrix score_matrix(lambda.size(), k);
  NumericMatrix mse_matrix(lambda.size(), k);
  NumericMatrix dpe_matrix(lambda.size(), k);
  NumericMatrix md_matrix(lambda.size(), k);

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
      Named("d")=d,
      Named("mu")=mu,
      Named("g")=g,
      Named("fixed_g")=fixed_g,
      Named("profile")=profile,
      Named("scad")=scad,
      Named("initialvals")= Rcpp::clone(initialvals),
      Named("lambda")=lambda,
      Named("theta_upper")=theta_upper,
      Named("theta_lower")=theta_lower,
      Named("ncores")=ncores
    );
    // estimation
    List out = multi(multi_list, Named("penalty")=true);
    // estimated values
    arma::mat theta = as<arma::mat>(out["theta"]);
    arma::vec s2 = as<arma::vec>(out["s2"]);
    arma::vec mu_vec = as<arma::vec>(out["mu"]);
    arma::vec g_vec = as<arma::vec>(out["g"]);

    // Rcout << "theta.nrow(): " << theta.n_rows << ", s2.size(): " << s2.n_elem << ", lambda.size(): " << lambda.size() << std::endl;


    for(int i=0; i<lambda.size(); i++) {
      if (i >= static_cast<int>(theta.n_rows)) Rcpp::stop("Index i out of bounds for theta");
      if (i >= static_cast<int>(s2.n_elem)) Rcpp::stop("Index i out of bounds for s2");
      List pred;
      pred = kriging(Named("y")=y_train,
                     Named("x")=x_train,
                     Named("xx")=x_valid,
                     Named("sigma2")=s2[i],
                     Named("theta")=theta.row(i),
                     Named("mu")=mu_vec[i],
                     Named("g")=g_vec[i]);

      NumericVector pred_error = y_valid-as<NumericVector>(pred["mup"]);
      double mse_value = sum(pow(pred_error, 2.0));

      // mse, score, md, and dpe
      mse_matrix(i, j-1) = mse_value;

      if(k != n) {
        List score_and_md = score(Named("y")=y_valid,
                             Named("mu")=pred["mup"],
                             Named("sigma")=pred["Sigmap"],
                             Named("md")=true);

        score_matrix(i, j-1) = score_and_md["score"];
        md_matrix(i, j-1) = score_and_md["md"];
        dpe_matrix(i, j-1) = as<double>(
          dpe(Named("y")=y_valid,
              Named("mu")=pred["mup"],
              Named("R")=pred["R"])
        );
      }

    }
  }

  // max score, 1se score, min mse, min 1se
  NumericVector score_means, mse_means, dpe_means, md_means, md_count;
  int score_max, score_1se, mse_min, mse_1se, dpe_min, dpe_1se, md_min, md_1se;

  if (k != n) { // if k fold, return both score and mse
    // score
    score_means = row_means(score_matrix);
    score_max = max_element(score_means.begin(), score_means.end()) - score_means.begin(); // need +1

    // dpe
    dpe_means = row_means(dpe_matrix);
    dpe_min = min_element(dpe_means.begin(), dpe_means.end()) - dpe_means.begin(); // need +1

    // mse
    mse_means = row_means(mse_matrix);
    mse_min = min_element(mse_means.begin(), mse_means.end()) - mse_means.begin(); // need +1

    // md
    md_means = row_means(md_matrix);
    md_min = min_element(md_means.begin(), md_means.end()) - md_means.begin(); // need +1


    // standard error
    double score_var = 0.0;
    double mse_var = 0.0;
    double dpe_var = 0.0;
    double md_var = 0.0;
    for(int j=0; j<k; j++) {
      score_var += pow(score_matrix(score_max, j) - score_means[score_max], 2);
      mse_var += pow(mse_matrix(mse_min, j) - mse_means[mse_min], 2);
      dpe_var += pow(dpe_matrix(dpe_min, j) - dpe_means[dpe_min], 2);
      md_var += pow(md_matrix(md_min, j) - md_means[md_min], 2);
    }
    // score
    double score_sd = sqrt(score_var/(k-1)); // standard deviation
    double score_one_se = score_means[score_max] - score_sd / sqrt(k);
    // mse
    double mse_sd = sqrt(mse_var/(k-1)); // standard deviation
    double mse_one_se = mse_means[mse_min] + mse_sd / sqrt(k);
    // dpe
    double dpe_sd = sqrt(dpe_var/(k-1));
    double dpe_one_se = dpe_means[dpe_min] + dpe_sd / sqrt(k);
    // md
    double md_sd = sqrt(md_var/(k-1));
    double md_one_se = md_means[md_min] + md_sd / sqrt(k);

    // find 1se index
    NumericVector score_se_indices;
    NumericVector mse_se_indices;
    NumericVector dpe_se_indices;
    NumericVector md_se_indices;
    for(int i=0; i<score_means.size(); i++) {
      if(score_means[i] >= score_one_se) {
        score_se_indices.push_back(i);
      }
      if(mse_means[i] <= mse_one_se) {
        mse_se_indices.push_back(i);
      }
      if(dpe_means[i] <= dpe_one_se) {
        dpe_se_indices.push_back(i);
      }
      if(md_means[i] <= md_one_se) {
        md_se_indices.push_back(i);
      }
    }
    // grab the largest index value (score)
    if(score_se_indices.size() > 0) {
      score_1se = score_se_indices[score_se_indices.size()-1]; // need +1
    } else {
      score_1se = score_max;
    }

    // grab the largest index value (mse)
    if(mse_se_indices.size() > 0) {
      mse_1se = mse_se_indices[mse_se_indices.size()-1]; // need +1
    } else {
      mse_1se = mse_min;
    }

    // grab the largest index value (dpe)
    if(dpe_se_indices.size() > 0) {
      dpe_1se = dpe_se_indices[dpe_se_indices.size()-1]; // need +1
    } else {
      dpe_1se = dpe_min;
    }

    // grab the largest index value (md)
    if(md_se_indices.size() > 0) {
      md_1se = md_se_indices[md_se_indices.size()-1]; // need +1
    } else {
      md_1se = md_min;
    }

    // md
    // md_count = md_sums_range(md_matrix, md_upper, md_lower);
    // find max index where md_count==k
    // int md_opt=1;
    // for(int i=md_count.size()-1; i>=0; i--) {
    //  if(md_count[i]==k) {
    //    md_opt = i; // need +1
    //    break;
    //  }
    //}

    return List::create(
      _["folds"] = folds,
      _["score"] = score_means,
      _["mse"] = mse_means,
      _["dpe"] = dpe_means,
      _["md_matrix"] = md_matrix,
      _["score_max"] = score_max + 1,
      _["score_1se"] = score_1se + 1,
      _["mse_min"] = mse_min + 1,
      _["mse_1se"] = mse_1se + 1,
      _["dpe_min"] = dpe_min + 1,
      _["dpe_1se"] = dpe_1se + 1,
      _["md_min"] = md_min + 1,
      _["md_1se"] = md_1se + 1
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
      _["mse"] = mse_means,
      _["mse_min"] = mse_min + 1,
      _["mse_1se"] = mse_1se + 1
    );
  }
}
