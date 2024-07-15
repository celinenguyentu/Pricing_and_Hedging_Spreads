//
//  spreads.cpp
//  
//
//  Created by Céline Nguyen on 14/04/2020.
//
#include <stdio.h>
#include "spreads.hpp"
#include <random>
#include <iostream>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <string>
#include <math.h>

// STATS
std::vector<double> Stats::get_Values() const {return Values;};

std::vector<double> Stats::get_Mean() const {return Mean;};

std::vector<double> Stats::get_Var() const {return Var;};

std::vector<double> Stats::get_uconfidence() const {
    std::vector<double> Uconfidence(Values.size());
    for (int i=0; i<Values.size(); i++){
        Uconfidence[i] = Mean[i]+(1.96*sqrt(Var[i]/(i+1)));
    }
    return Uconfidence;
};

std::vector<double> Stats::get_lconfidence() const {
    std::vector<double> Lconfidence(Values.size());
    for (int i=0; i<Values.size(); i++){
        Lconfidence[i] = Mean[i]-(1.96*sqrt(Var[i]/(i+1)));
    }
    return Lconfidence;
};

double Stats::get_MonteCarlo() const {return Mean.back();};

Stats & operator += (Stats & stat, double x){
    (stat.Values).push_back(x);
    int n = (stat.Values).size();
    if (n == 1){
        (stat.Mean).push_back(x);
        (stat.Cumsum).push_back(x*x);
        (stat.Var).push_back(0);
    }
    else{
        double newMean = ((stat.Mean).back()*(stat.Mean).size()+x)/n;
        (stat.Mean).push_back(newMean);
        double newCumsum = (stat.Cumsum).back()+x*x;
        (stat.Cumsum).push_back(newCumsum);
        double newVar = (stat.Cumsum).back()/n - (stat.Mean).back()*(stat.Mean).back();
        (stat.Var).push_back(newVar);
    }
    return stat;
}

std::ostream & operator << (std::ostream & flux, const Stats & stat){
    std::vector<double> uconf = stat.get_uconfidence();
    std::vector<double> lconf = stat.get_lconfidence();
    for (unsigned k=0; k<(stat.Mean).size(); k++){
        flux << k+1 << " " << (stat.Values)[k] << " " << (stat.Mean)[k] << " " << uconf[k] << " " << lconf[k] <<  std::endl;
    }
    return flux;
}

void Export(std::string s, const Stats & stat){
    std::ofstream fichier(s);
    std::vector<double> values = stat.get_Values();
    std::vector<double> mean = stat.get_Mean();
    std::vector<double> uconf = stat.get_uconfidence();
    std::vector<double> lconf = stat.get_lconfidence();
    for (int i=0; i<values.size(); i++){
        fichier << i+1 << " " << values[i] << " " << mean[i] << " " << uconf[i] << " " << lconf[i] << std::endl;
    }
    fichier.close();
};


// SPREAD

// Density of distribution N(0,1)
double normalDensity(double x)
{
    return exp(-x*x/2)/sqrt(2*M_PI);
}

// Cumulative distribution function of N(0,1)
double normalCDF(double x) // Phi(-∞, x) aka N(x)
{
    return std::erfc(-x/std::sqrt(2))/2;
}

// Computes a Monte Carlo estimator
Stats Spread::MonteCarlo_price(std::mt19937 & G, long unsigned int n, double t) const {
    auto Function_to_evaluate = [=](std::mt19937 & G){
        std::normal_distribution<double> N(0,1);
        double U = N(G);
        double V = N(G);
        return exp(-r*T)*payoff(U,V);};
    Stats stat;
    MonteCarlo(stat, Function_to_evaluate, G, n);
    return stat;
};

// Returns payoff at time T of spread option using samples of two independant variables distributed normally
double Spread::payoff(double U, double V) const {
    double S1T = x1*exp((r-sigma1*sigma1/2)*T+sigma1*rho*sqrt(T)*U+sigma2*sqrt(1-rho*rho)*sqrt(T)*V);
    double S2T = x2*exp((r-sigma2*sigma2/2)*T +sigma2*sqrt(T)*U);
    if ((S2T-S1T-K)>0) return S2T-S1T-K;
    else return 0;
}

// Computes Bachelier approximation
double Spread::Bachelier_price() const {
    double m = (x2-x1);
    double s2 = x1*x1*(exp(sigma1*sigma1*T)-1)-2*x1*x2*(exp(rho*sigma1*sigma2*T)-1)+x2*x2*(exp(sigma2*sigma2*T)-1);
    double temp = m-K*exp(-r*T);
    return temp*normalCDF(temp/sqrt(s2))+sqrt(s2)*normalDensity(temp/sqrt(s2));
}

// Computes Margrabe formula when K!=0
double Spread::Margrabe_price() const {
    try {
        if (K!=0) throw std::invalid_argument("Invalid K");
        else {
            double sigma = sqrt(sigma1*sigma1+sigma2*sigma2-2*rho*sigma1*sigma2);
            double d1 = log(x2/x1)/(sigma*sqrt(T))+sigma*sqrt(T)/2;
            double d2 = log(x2/x1)/(sigma*sqrt(T))-sigma*sqrt(T)/2;
            return x2*normalCDF(d1)-x1*normalCDF(d2);
        }
    }
    catch (const std::invalid_argument& msg) {
        std::cout << std::endl;
        std::cerr << msg.what() << std::endl;
        std::exit( EXIT_FAILURE );
    }
}

// Computes Kirk approximation
double Spread::Kirk_price() const {
    double Y0 = x1 - K*exp(-r*T);
    double sigma = sqrt(sigma1*sigma1*(x1/Y0)*(x1/Y0)+sigma2*sigma2 - 2*rho*sigma1*sigma2*(x1/Y0));
    double d1 = log(x2/Y0)/(sigma*sqrt(T))+sigma*T/2;
    double d2 = log(x2/Y0)/(sigma*sqrt(T))-sigma*T/2;
    return x2*normalCDF(d1)-Y0*normalCDF(d2);
}
