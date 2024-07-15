//
//  spreads_examples.cpp
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

int main(){
    std::mt19937 G(time(NULL));
    {
        std::cout << "SPREAD OPTION" << std::endl;
        // Parameters :
        int T = 1;
        double S10 = 50;
        double S20 = 80;
        double sigma1 = 0.3;
        double sigma2 = 0.7;
        double rho = 0.2;
        double r = 0.05;
        double K = 20;
        int n = 10000;
        Spread S(S10, S20, sigma1, sigma2, rho, r, K, T);
        
        // Monte Carlo approximation :
        Stats MonteCarlo = S.MonteCarlo_price(G,n);
        std::cout << "Monte Carlo : " << MonteCarlo.get_MonteCarlo() << std::endl;
        // Export data:
        Export("MonteCarloPrice.dat", MonteCarlo);
        
        // Bachelier approximation :
        double Bachelier = S.Bachelier_price();
        std::cout << "Bachelier : " << Bachelier << std::endl;
        
        // Margrabe price :
        if (K==0){
            double Margrabe = S.Margrabe_price();
            std::cout << "Margrabe : " << Margrabe << std::endl;
        }
        // Kirk price :
        double Kirk = S.Kirk_price();
        std::cout << "Kirk : " << Kirk << std::endl;
        
    }
    return 0;

}
