#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerMatrix gibbsSampler(IntegerMatrix z, 
                           double alpha,
                           double beta,
                           int ncolors,
                           NumericMatrix seeds,         // Cada linha: (x, y) de uma seed
                           IntegerVector seedColors,    // Cor associada a cada seed
                           NumericVector seedDeltas,    // Parâmetro escalar de cada seed
                           int steps) {
  int nrows = z.nrow();
  int ncols = z.ncol();
  int n = nrows * ncols;
  int nseeds = seeds.nrow();
  
  NumericMatrix distances(n, nseeds);
  
  for(int pos = 0; pos < n; pos++) {
    int pos_x = pos % nrows;
    int pos_y = pos / nrows;
    
    for(int s = 0; s < nseeds; s++) {
      double seed_x = seeds(s, 0) - 1;
      double seed_y = seeds(s, 1) - 1;
      
      double dx = pos_x - seed_x;
      double dy = pos_y - seed_y;
      distances(pos, s) = std::sqrt(dx * dx + dy * dy);
      
    }
  }
  
  IntegerVector positions = seq(0, n - 1);
  
  for(int iter = 0; iter < steps; iter++) {
    IntegerVector updateOrder = sample(positions, n, false);
    
    for(int idx = 0; idx < n; idx++) {
      int pos = updateOrder[idx];
      int pos_x = pos % nrows;
      int pos_y = pos / nrows;
      
      std::vector<int> neighbors;
      neighbors.reserve(4);
      
      if(pos_x > 0) neighbors.push_back(z(pos_x - 1, pos_y));
      if(pos_x < nrows - 1) neighbors.push_back(z(pos_x + 1, pos_y));
      if(pos_y > 0) neighbors.push_back(z(pos_x, pos_y - 1));
      if(pos_y < ncols - 1) neighbors.push_back(z(pos_x, pos_y + 1));
      
      NumericVector energies(ncolors);
      double max_energy = R_NegInf;
      
      for(int v = 0; v < ncolors; v++) {
        double contrib = 0.0;
        for(size_t j = 0; j < neighbors.size(); j++) {
          contrib += (v == neighbors[j]) ? alpha : 0.0;
        }
        
        double energy = contrib + (v == 0 ? beta : 0.0);
        
        for(int s = 0; s < nseeds; s++) {
          if(v == seedColors[s]) {
            energy += seedDeltas[s] / (1 + distances(pos, s));
          }
        }
        
        energies[v] = energy;
        if(energies[v] > max_energy) max_energy = energies[v];
      }
      
      NumericVector probs(ncolors);
      double sum_probs = 0.0;
      
      for(int v = 0; v < ncolors; v++) {
        probs[v] = exp(energies[v] - max_energy);
        sum_probs += probs[v];
      }
      
      for(int v = 0; v < ncolors; v++) {
        probs[v] /= sum_probs;
      }
      
      double u = R::runif(0, 1);
      double cumsum = 0.0;
      int selected = 0;
      
      for(int v = 0; v < ncolors; v++) {
        cumsum += probs[v];
        if(u <= cumsum) {
          selected = v;
          break;
        }
      }
      
      z(pos_x, pos_y) = selected;
    }
  }
  
  return z;
}

// [[Rcpp::export]]
IntegerMatrix conditionalGibbsSampler(NumericMatrix y,
                                      IntegerMatrix z, 
                                      double alpha,
                                      double beta,
                                      int ncolors,
                                      NumericMatrix seeds,
                                      IntegerVector seedColors,
                                      NumericVector seedDeltas,
                                      NumericVector mus,
                                      NumericVector sigmas,
                                      int steps) {
  int nrows = z.nrow();
  int ncols = z.ncol();
  int n = nrows * ncols;
  int nseeds = seeds.nrow();
  
  NumericMatrix distances(n, nseeds);
  
  for(int pos = 0; pos < n; pos++) {
    int pos_x = pos % nrows;
    int pos_y = pos / nrows;
    
    for(int s = 0; s < nseeds; s++) {
      double seed_x = seeds(s, 0) - 1;
      double seed_y = seeds(s, 1) - 1;
      
      double dx = pos_x - seed_x;
      double dy = pos_y - seed_y;
      distances(pos, s) = std::sqrt(dx * dx + dy * dy);
      
    }
  }
  
  IntegerVector positions = seq(0, n - 1);
  
  for(int iter = 0; iter < steps; iter++) {
    IntegerVector updateOrder = sample(positions, n, false);
    
    for(int idx = 0; idx < n; idx++) {
      int pos = updateOrder[idx];
      int pos_x = pos % nrows;
      int pos_y = pos / nrows;
      
      std::vector<int> neighbors;
      neighbors.reserve(4);
      
      if(pos_x > 0) neighbors.push_back(z(pos_x - 1, pos_y));
      if(pos_x < nrows - 1) neighbors.push_back(z(pos_x + 1, pos_y));
      if(pos_y > 0) neighbors.push_back(z(pos_x, pos_y - 1));
      if(pos_y < ncols - 1) neighbors.push_back(z(pos_x, pos_y + 1));
      
      NumericVector energies(ncolors);
      double max_energy = R_NegInf;
      
      for(int v = 0; v < ncolors; v++) {
        double contrib = 0.0;
        for(size_t j = 0; j < neighbors.size(); j++) {
          contrib += (v == neighbors[j]) ? alpha : 0.0;
        }
        
        double energy = contrib + (v == 0 ? beta : 0.0);
        
        for(int s = 0; s < nseeds; s++) {
          if(v == seedColors[s]) {
            energy += seedDeltas[s] / (1 + distances(pos, s));
          }
        }
        
        energy += R::dnorm(y(pos_x, pos_y), mus[v], sigmas[v], true);
        
        energies[v] = energy;
        if(energies[v] > max_energy) max_energy = energies[v];
      }
      
      NumericVector probs(ncolors);
      double sum_probs = 0.0;
      
      for(int v = 0; v < ncolors; v++) {
        probs[v] = exp(energies[v] - max_energy);
        sum_probs += probs[v];
      }
      
      for(int v = 0; v < ncolors; v++) {
        probs[v] /= sum_probs;
      }
      
      double u = R::runif(0, 1);
      double cumsum = 0.0;
      int selected = 0;
      
      for(int v = 0; v < ncolors; v++) {
        cumsum += probs[v];
        if(u <= cumsum) {
          selected = v;
          break;
        }
      }
      
      z(pos_x, pos_y) = selected;
    }
  }
  
  return z;
}