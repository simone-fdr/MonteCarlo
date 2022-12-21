#include "mc.hpp"
#include <random>
#include <iostream>
#include <chrono>
#include <array>
#include <functional>
#include <math.h>


//! @param f is the function to integrate
//! @param domain of integration
//! @param sample_size iterations
double MonteCarlo_integral(std::function<double (V<double>)> const &f, Geometry& domain ,int sample_size){
    double volume = domain.getVolume();
    double value;
    for(int i = 0; i < sample_size; i++){
        value += f(domain.generateRandomVector());
    }
    value = value * volume / sample_size;
    return value;
}

void test_rectangle(int sample_size, std::default_random_engine eng){
    // Function: x1+x2+...+xn
    auto f = [](V<double> x){
        double r = 0;
        for(int i = 0; i < x.size(); i++){
            r+=x[i];
        }
        return r;
    };

    std::vector<double> boundaries{0, 1, 0, 1, 0, 1, 0, 1};

    HyperRectangle rect(boundaries, eng);

    auto start = std::chrono::high_resolution_clock::now();
    double integral = MonteCarlo_integral(f, rect, sample_size);
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);  
    std::cout << "The result is " << integral << std::endl;
    std::cout << "Expected was " << boundaries.size()/4.0 << std::endl;
    std::cout << "In " << duration.count()  << " milliseconds" << std::endl;
}

void test_ball(int sample_size, std::default_random_engine eng){
    auto normal = [](V<double> x){
        double exponent = 0;
        for(int i = 0; i < x.size(); i++){
            exponent -= x[i]*x[i];
        }
        return std::exp(exponent);
    };
    V<double> center{0, 0};
    double radius = 3;
    
    HyperBall sphere(center,radius, eng);

    auto start = std::chrono::high_resolution_clock::now();
    double integral = MonteCarlo_integral(normal, sphere, sample_size);
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);  
    std::cout << "The result is " << integral << std::endl;
    std::cout << "Expected was " << std::sqrt(std::pow(M_PI, center.size()))<< std::endl;
    std::cout << "In " << duration.count()  << " milliseconds" << std::endl;
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

    test_ball(sample_size, eng);
    test_rectangle(sample_size, eng);
    return 0;
}