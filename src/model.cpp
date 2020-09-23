// [[Rcpp::depends(RcppParallel, RcppArmadillo)]]

#include <stdio.h>
#include <RcppParallel.h>
// Note: RcppArmadillo is thread-safe
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace RcppParallel;

//' Fast extraction of unique rows of a numeric matrix
//' @param expected: expected number of unique rows to pre-reserve memory and avoid costly
//' resizing.
// [[Rcpp::export]]
arma::mat fastNumericUnique(
  arma::mat x,
  size_t expected
) {
  // Works best when there are a lot of duplicates and only a few of unique rows
  std::vector<unsigned int> idx;
  idx.reserve(expected);

  for (int i = 0; i < x.n_rows;) {
    for (int j = 0; j < idx.size(); ++j) {
      if (all(x.row(i) == x.row(idx[j]))) {
        goto next;
      }
    }
    idx.push_back(i);
next:
    ++i;
  }

  return  x.rows(arma::uvec(idx));
}

struct Plan {
  arma::uvec rowSubset;
  int m;
  int file;
  int iiMax;
};

//' The `Worker` paradigm comes from `RcppParallel`.
struct DecayModelExecutor : public Worker {
  arma::vec &ii;
  arma::mat &allAlphas;
  arma::mat &allBetas;
  double timeinterval;
  arma::vec &prompts;
  arma::vec &promptStarts;
  arma::vec &promptEnds;
  std::size_t planLength;
  Plan *plans;
  arma::vec &output;
  arma::mat *outputGrad;

  DecayModelExecutor(
    arma::vec &ii,
    arma::mat &allAlphas,
    arma::mat &allBetas,
    double timeinterval,
    arma::vec &prompts,
    arma::vec &promptStarts,
    arma::vec &promptEnds,
    List plans,
    arma::vec &output,
    arma::mat *outputGrad
  )
    : ii(ii), allAlphas(allAlphas), allBetas(allBetas),
      timeinterval(timeinterval),
      prompts(prompts), promptStarts(promptStarts), promptEnds(promptEnds),
      planLength(plans.length()),
      output(output),
      outputGrad(outputGrad) {
    this->plans = new Plan[plans.length()];
    for (int i = 0; i < plans.length(); ++i) {
      List curPlan = plans[i];
      this->plans[i].rowSubset = as<arma::uvec>(curPlan["rowSubset"]);
      this->plans[i].m = as<int>(curPlan["m"]);
      this->plans[i].file = as<int>(curPlan["file"]);
      this->plans[i].iiMax = as<int>(curPlan["iiMax"]);
    }
  }

  ~DecayModelExecutor() {
    delete[] plans;
  }

  void modelOne(int i) {
    Plan &curPlan = plans[i];

    arma::vec thetas = allAlphas.col(curPlan.m - 1);
    arma::vec alphas = exp(thetas);
    arma::vec Bs = allBetas.col(curPlan.m - 1);
    int BsLength = Bs.size();
    int iiMax = curPlan.iiMax;
    int iiLength = curPlan.rowSubset.size();
    int file = plans[i].file;
    arma::vec prompt = prompts.subvec(
      promptStarts[file - 1] - 1,
      promptEnds[file - 1] - 1
    );

    arma::vec expB = exp(-alphas * timeinterval);

    arma::vec ans(iiLength, arma::fill::zeros);
    int idxMax = iiMax;
    for (int k = 0; k < BsLength; ++k) {
      arma::vec powExpB(idxMax);
      for (int idx = 0; idx < idxMax; ++idx) {
        powExpB(idx) = pow(expB(k), idx);
      }

      arma::vec ansFFT = real(ifft(fft(prompt) % fft(powExpB, prompt.size())));
      arma::vec curAns = Bs[k] * ansFFT;
      ans += curAns;

      if (this->outputGrad != NULL) {
        // Compute gradient (Jacobian) with respect to thetas
        arma::vec cst(idxMax);
        for (int idx = 0; idx < idxMax; ++idx) {
          cst(idx) = powExpB(idx) * timeinterval * idx;
        }
        arma::vec thetasGrad =
          - (Bs[k] * real(ifft(fft(prompt) % fft(cst, prompt.size())))) * alphas[k];
        arma::uvec uk1(1);
        uk1 << k * 2;
        outputGrad->submat(curPlan.rowSubset, uk1) = thetasGrad;

        // Compute gradient (Jacobian) with respect to Bs
        arma::uvec uk2(1);
        uk2 << k * 2 + 1;
        outputGrad->submat(curPlan.rowSubset, uk2) = ansFFT;
      }
    }

    output.elem(curPlan.rowSubset) = ans;
  }

  void operator()(std::size_t begin, std::size_t end) {
    for (int i = begin; i < end; ++i) {
      modelOne(i);
    }
  }
};

// [[Rcpp::export]]
NumericVector decayModelExecute(
  int outputLength,
  arma::vec ii,
  arma::mat allAlphas,
  arma::mat allBetas,
  double timeinterval,
  arma::vec prompts,
  arma::vec promptStarts,
  arma::vec promptEnds,
  List plans,
  bool grad,
  bool multithread
) {
  if (grad) {
    arma::vec ans(outputLength);
    arma::mat ansGrad(outputLength, allAlphas.n_rows * 2);
    DecayModelExecutor executor(
      ii,
      allAlphas,
      allBetas,
      timeinterval,
      prompts,
      promptStarts,
      promptEnds,
      plans,
      ans,
      &ansGrad
    );
    if (multithread) {
      parallelFor(0, plans.length(), executor);
    } else {
      executor(0, plans.length());
    }
    NumericVector nvAns(ans.begin(), ans.end());
    NumericMatrix nmAnsGrad = wrap(ansGrad);
    CharacterVector gradColnames(nmAnsGrad.ncol());
    for (int i = 0; i < nmAnsGrad.ncol() / 2; ++i) {
      std::ostringstream stheta;
      stheta << "theta" << (i + 1);
      gradColnames[i * 2] = stheta.str();

      std::ostringstream sB;
      sB << "B" << (i + 1);
      gradColnames[i * 2 + 1] = sB.str();
    }
    colnames(nmAnsGrad) = gradColnames;
    nvAns.attr("gradient") = nmAnsGrad;
    return nvAns;
  } else {
    arma::vec ans(outputLength);
    DecayModelExecutor executor(
        ii,
        allAlphas,
        allBetas,
        timeinterval,
        prompts,
        promptStarts,
        promptEnds,
        plans,
        ans,
        NULL
    );
    if (multithread) {
      parallelFor(0, plans.length(), executor);
    } else {
      executor(0, plans.length());
    }
    return NumericVector(ans.begin(), ans.end());
  }
}


