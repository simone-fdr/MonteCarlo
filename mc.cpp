#include <random>
#include <iostream>
#include <chrono>
#include <array>
#include <functional>
#include <math.h>
#include "V.hpp"
#include "RandomVector.hpp"

// Everything is done with doubles for simplicity, but can be done with templates

//! @param f is the function to integrate
//! @param domain of integration
//! @param sample_size iterations
double MonteCarlo_integral(std::function<double (V<double>)> const &f, Geometry const & domain ,int sample_size){
    double volume = domain.getVolume();
    double value;
    for(int i = 0; i < sample_size; i++){
        value += f(domain.generateRandomVector());
    }
    value = value * volume / sample_size;
    return value;
}

// ==================== MAIN ==================
// Define the function
// Define the domain of integration
// Call the MonterCarlo_integral function

int main(int argc, char ** argv){

    if(argc != 2){
        std::cout << "./integration #sample" << std::endl;
        return 0;
    }
    int sample_size = std::atoi(argv[1]);
    std::random_device rd;
    auto seed = rd();
    std::default_random_engine eng(seed);

    // Function: x1*x2*...*xn
    auto f = [](V<double> x){
        double r = 1;
        for(int i = 0; i < x.size(); i++){
            r*=x[i];
        }
        return r;
    };

    /*V<double> center{3, 4, 5};
    double radius = 2;
    
    auto start = std::chrono::high_resolution_clock::now();
    double integral = MonteCarlo_integral_hb(f, center, radius, sample_size, eng);
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);  
    std::cout << "The result is " << integral << std::endl;
    std::cout << "In " << duration.count()  << " milliseconds" << std::endl;*/


    return 0;
}
