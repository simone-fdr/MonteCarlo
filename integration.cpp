#include <random>
#include <iostream>
#include <chrono>
#include <array>
#include<functional>

// Create an hyper rectangle distribution
//! @param boundaries in pairs [a1,b1,a2,b2,...,an,bn]
//! @return vector of distributions of size n
std::vector<std::uniform_real_distribution<double>> generate_hr_dist(std::vector<double> &boundaries){
    std::vector<std::uniform_real_distribution<double>> hr_dis; // Vettore di distribuzioni grande n
    for(int i = 0; i < boundaries.size(); i+=2){
        hr_dis.emplace_back(std::uniform_real_distribution<double>(boundaries[i], boundaries[i+1]));
    }
    return hr_dis;
}

// Returns a single random vector within an hyper rectangle
//! @param hr_dist vector of distributions
//! @return random vector of size n inside
std::vector<double> random_vector(std::vector<std::uniform_real_distribution<double>> &hr_dis, std::default_random_engine &eng){
    std::vector<double> rv; //random vector
    for(int i = 0; i < hr_dis.size(); i++){
        rv.emplace_back(hr_dis[i](eng));
    }
    return rv;
}

double hr_hv(std::vector<double> &boundaries){
    double hv = 1;
    for(int i = 0; i < boundaries.size(); i+=2)
        hv *= (boundaries[i+1] - boundaries[i]);
    return hv;
}

//! @param f is the function to integrate
//! @param boundaries of the hyper rectangle
//! @param sample_size iterations
// TODO togli e metti le reference per vedere l'efficenza, non dovrebbe cambiare molto perché sono vettori n-dimensionali
// TODO Non è bello che eng passi fra così tante funzioni
double MonteCarlo_integral_hr(std::function<double (std::vector<double>)> const &f, std::vector<double> &boundaries,
                             int sample_size, std::default_random_engine &eng){
    std::vector<double> rand;
    std::vector<std::uniform_real_distribution<double>> hr_dis = generate_hr_dist(boundaries);
    double area = hr_hv(boundaries);
    double value;
    for(int i = 0; i < sample_size; i++){
        value += f(random_vector(hr_dis, eng));
    }
    value = value * area / sample_size;
    return value;
}

void print_vector(std::vector<double> &x){
    for(auto el : x){
        std::cout << " " << el;
    }
    std::cout << std::endl;
}

int main(int argc, char ** argv){

    if(argc != 2){
        std::cout << "./integration #sample" << std::endl;
        return 0;
    }
    int sample_size = std::atoi(argv[1]);
    std::random_device rd;
    auto seed = rd();
    std::default_random_engine eng(seed);

    // Funzione: x1*x2*...*xn
    auto produttoria = [](std::vector<double> x){
        double r = 1;
        for(int i = 0; i < x.size(); i++){
            r*=x[i];
        }
        return r;
    };

    std::vector<double> boundaries{ 1,2,
                                    2,4,
                                    2,3 };

    
    auto start = std::chrono::high_resolution_clock::now();
    double integral = MonteCarlo_integral_hr(produttoria, boundaries, sample_size, eng);
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);  
    std::cout << "L'integrale vale " << integral << std::endl;
    std::cout << "In " << duration.count()  << " millisecondi" << std::endl;

    return 0;
}