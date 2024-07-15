//
//  spreads.hpp
//  
//
//  Created by Céline Nguyen on 14/04/2020.
//

#ifndef spreads_hpp
#define spreads_hpp

#include <stdio.h>
#include <random>
#include <algorithm>
#include <cmath>
#include <utility>
#include <string>
#include <fstream>
#include <vector>
#include <string>
#include <stdexcept>
#include <unordered_set>
#include <iostream>

class Stats{ // Keeps in storage all simulations of a Monte Carlo method in order to plot convergence 
protected:
    std::vector<double> Values;
    std::vector<double> Mean;
    std::vector<double> Cumsum;
    std::vector<double> Var;
public:
    Stats(double x=0.): Values(), Mean(), Cumsum(), Var() {};
    std::vector<double> get_Values() const;
    std::vector<double> get_Mean() const;
    std::vector<double> get_Var() const;
    std::vector<double> get_uconfidence() const;
    std::vector<double> get_lconfidence() const;
    double get_MonteCarlo() const;
    friend Stats & operator+=(Stats &, double);
    friend std::ostream & operator << (std::ostream &, const Stats &);
    friend void Export(std::string, const Stats &);
};

// Monte Carlo function
template <class Statistique, class Measurement, class RNG>
void MonteCarlo(Statistique & res, const Measurement & f, RNG & G, long unsigned int n){
    for (long unsigned i=0;i<n; i++){
        res += f(G);
    }
};



class Spread {
protected:
    double x1;
    double x2;
    double sigma1;
    double sigma2;
    double rho;
    double r;
    double K;
    int T;
public:
    Spread(double x1_, double x2_, double sigma1_, double sigma2_, double rho_, double r_, double K_, int T_): x1(x1_), x2(x2_), sigma1(sigma1_), sigma2(sigma2_), rho(rho_), r(r_), K(K_), T(T_) {};
    Stats MonteCarlo_price(std::mt19937 & G, long unsigned int n, double t=0.) const;
    double payoff(double U, double V) const;
    double Bachelier_price() const;
    double Margrabe_price() const;
    double Kirk_price() const;
};



#endif /* spreads_hpp */
