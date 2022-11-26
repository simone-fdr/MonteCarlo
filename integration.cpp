#include <random>
#include <iostream>
#include <chrono>
#include <array>
#include <functional>
#include <math.h>

void print_vector(std::vector<double> x){
    for(auto el : x){
        std::cout << " " << el;
    }
    std::cout << std::endl;
}

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

double norm_vector(std::vector<double> v){
    double norm = 0;
    for(int vi:v){
        norm += (vi*vi);
    }
    return norm;    
}

// TODO Study between spherical coordinates vs trial-error
//! @param c center
//! @param r radius
//! @return random vector within hypersphere
std::vector<double> random_vector_hs(std::vector<double> c, double r, std::default_random_engine &eng){
    std::vector<double> rv;
    std::uniform_real_distribution<double> dis(0,1);
    int TMP_COUNT = 0;
    do{
        ////////////////////////rv.clear();
        for(int i = 0; i < c.size(); i++){
            rv.emplace_back(dis(eng));
        }
        if(TMP_COUNT) print_vector(rv);
        TMP_COUNT++;
        // TODO Check this do-while()
    }while(norm_vector(rv) > 1.);
    for(int i = 0; i < c.size(); i++){
        rv[i] = rv[i] * r + c[i];
    }
    return rv;
}

// Returns a single random vector within an hyper rectangle
//! @param hr_dist vector of distributions
//! @return random vector of size n inside
std::vector<double> random_vector_hr(std::vector<std::uniform_real_distribution<double>> &hr_dis, std::default_random_engine &eng){
    std::vector<double> rv; //random vector
    for(int i = 0; i < hr_dis.size(); i++){
        rv.emplace_back(hr_dis[i](eng));
    }
    return rv;
}

//! @return (n-1)!
double gamma(double n){
    if(n == 1.) return 1;
    if(n == 0.5) return std::sqrt(M_PI);
    return gamma(n-1)*(n-1);
}

double hs_hv(double r, double n){
    // TODO choose a better threshold
    if(n > 20){
        return (1./std::sqrt(n*M_PI))*std::pow(2*M_PI*M_E/n,n/2.0)*std::pow(r,n);
    }
    return std::pow(r,n)*std::pow(M_PI, n/2.)/gamma(n/2. + 1);
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
        value += f(random_vector_hr(hr_dis, eng));
    }
    value = value * area / sample_size;
    return value;
}

int main(int argc, char ** argv){

    /*if(argc != 2){
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
    std::cout << "In " << duration.count()  << " millisecondi" << std::endl;*/


    std::random_device rd;
    auto seed = rd();
    std::default_random_engine eng(seed);

    std::vector<double> center{0,0,0};
    for(double i = 0; i < 100; i+=1){
        random_vector_hs(center, 1, eng);
    }

    return 0;
}