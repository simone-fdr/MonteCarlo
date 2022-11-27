#include <random>
#include <iostream>
#include <chrono>
#include <array>
#include <functional>
#include <math.h>


// ================ VECTOR UTILS =======================
void print_vector(std::vector<double> x){
    for(auto el : x){
        std::cout << " " << el;
    }
    std::cout << std::endl;
}


double norm_vector(std::vector<double> v){
    double norm = 0;
    for(double vi:v){
        norm += (vi*vi);
    }
    return norm;    
}

// ============== HYPER RECTANGLE ===================

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
std::vector<double> random_vector_hr(std::vector<std::uniform_real_distribution<double>> &hr_dis, std::default_random_engine &eng){
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

// ================== HYPER SPHERE =======================

//! @param c center
//! @param r radius
//! @return random vector on the hypersphere
// Monte Carlo integration does not work on hyperspheres since they are volume equals to 0
std::vector<double> random_vector_hs(std::vector<double> c, double r, std::default_random_engine &eng){
    std::vector<double> rv;
    std::uniform_real_distribution<double> dis(-10,10); //-10 to 10 in order to avoid floating point errors while normalizing dividing by a small value
    for(int i = 0; i < c.size(); i++){
        rv.emplace_back(dis(eng));
    }
    double norm = norm_vector(rv);
    for(int i = 0; i < c.size(); i++){
        rv[i] = rv[i] * r / norm + c[i];
    }
    return rv;
}

// http://compneuro.uwaterloo.ca/files/publications/voelker.2017.pdf
//! @param c center
//! @param r radius
//! @return random vector within hyperball
std::vector<double> random_vector_hb(std::vector<double> c, double r, std::default_random_engine &eng){
    std::vector<double> rv;
    std::vector<double> origin(c.size()+2, 0.);
    rv = random_vector_hs(origin,1,eng);
    rv.pop_back();
    rv.pop_back();
    for(int i = 0; i < c.size(); i++){
        rv[i] = rv[i] * r + c[i];
    }
    return rv;
}

//! @param c center
//! @param r radius
//! @return random vector within hyperball
std::vector<double> random_vector_hb_inefficent(std::vector<double> c, double r, std::default_random_engine &eng){
    std::vector<double> rv;
    std::uniform_real_distribution<double> dis(-1,1);
    do{
        rv.clear();
        for(int i = 0; i < c.size(); i++){
            rv.emplace_back(dis(eng));
        }
    }while(norm_vector(rv) > 1.);
    for(int i = 0; i < c.size(); i++){
        rv[i] = rv[i] * r + c[i];
    }
    return rv;
}


//! @return (n-1)!
double gamma(double n){
    if(n == 1.) return 1;
    if(n == 0.5) return std::sqrt(M_PI);
    return gamma(n-1)*(n-1);
}

double hb_hv(double r, double n){
    // TODO choose a better threshold
    if(n > 20){
        return (1./std::sqrt(n*M_PI))*std::pow(2*M_PI*M_E/n,n/2.0)*std::pow(r,n);
    }
    return std::pow(r,n)*std::pow(M_PI, n/2.)/gamma(n/2. + 1);
}

// ======================== MONTECARLO INTEGRALS ===================

//! @param f is the function to integrate
//! @param boundaries of the hyper rectangle
//! @param sample_size iterations
// TODO togli e metti le reference per vedere l'efficenza, non dovrebbe cambiare molto perche' sono vettori n-dimensionali
// TODO Non e' bello che eng passi fra cosi' tante funzioni
double MonteCarlo_integral_hr(std::function<double (std::vector<double>)> const &f, std::vector<double> &boundaries,
                             int sample_size, std::default_random_engine &eng){
    std::vector<std::uniform_real_distribution<double>> hr_dis = generate_hr_dist(boundaries);
    double area = hr_hv(boundaries);
    double value;
    for(int i = 0; i < sample_size; i++){
        value += f(random_vector_hr(hr_dis, eng));
    }
    value = value * area / sample_size;
    return value;
}

double MonteCarlo_integral_hb(std::function<double (std::vector<double>)> const &f, std::vector<double> &c, double r,
                             int sample_size, std::default_random_engine &eng){
    double area = hb_hv(r, c.size());
    double value;
    for(int i = 0; i < sample_size; i++){
        value += f(random_vector_hb_inefficent(c,r,eng));
    }
    value = value * area / sample_size;
    return value;
}

// ==================== MAIN ==================

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

    std::vector<double> center{3, 4, 5};
    double radius = 2;
    
    auto start = std::chrono::high_resolution_clock::now();
    double integral = MonteCarlo_integral_hb(produttoria, center, radius, sample_size, eng);
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);  
    std::cout << "L'integrale vale " << integral << std::endl;
    std::cout << "In " << duration.count()  << " millisecondi" << std::endl;

    return 0;
}