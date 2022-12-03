#ifndef RANDOM_VECTOR_HPP
#define RANDOM_VECTOR_HPP

#include <random>
#include <vector>
#include "V.hpp"

class HyperRectangle{
    public:
    // Create an hyper rectangle distribution and fills distributions
    //! @param seed to initialize engine
    //! @param boundaries in pairs [a1,b1,a2,b2,...,an,bn]
    HyperRectangle(std::vector<double> boundaries_, std::default_random_engine engine_);
    // Returns a single random vector within an hyper rectangle
    //! @param hr_dist vector of distributions
    //! @return random vector of size n inside
    V<double> generateRandomVector();
    //! @param boundaries of the hyperrectangle
    //! @return hypervolume of hyperrectangle
    double getVolume();
    private:
    // Collection of values representing the boundaries
    std::vector<double> boundaries;
    // Collection of distributions
    std::vector<std::uniform_real_distribution<double>> distribuitions;
    // Engine
    std::default_random_engine engine;
};

class HyperBall{
    public:
    // http://compneuro.uwaterloo.ca/files/publications/voelker.2017.pdf
    //! @param c center
    //! @param r radius
    //! @return random vector within hyperball
    V<double> generateRandomVector(V<double> c, double r, std::default_random_engine &eng);
    // It's very inefficient for n>>0
    //! @param c center
    //! @param r radius
    //! @return random vector within hyperball
    V<double> generateRandomVectorInefficient(V<double> c, double r, std::default_random_engine &eng);
    //For low dimension it uses recursive definition, otherwise it use Stirling approximation
    //! @param r radius
    //! @param n dimension
    //! @return hypervolume of hyperball
    double getVolume(double r, double n);
    private:
    //! @param c center
    //! @param r radius
    //! @return random vector on the hypersphere
    // Monte Carlo integration does not work on hyperspheres since they have volume equals to 0
    V<double> generateRandomVectorSphere(V<double> c, double r, std::default_random_engine &eng);
    //! @param n
    //! @return (n-1)!
    double gamma(double n);
};

class HyperSimplex{
    //! @param vs collection of vertices
    //! @param eng engine for random numbers
    //! @return random vector inside the simplex
    V<double> generateRandomVector(std::vector<V<double>> vs, std::default_random_engine &eng);
    // https://en.wikipedia.org/wiki/Simplex#Volume
    //! @param vs collection of vertices
    //! @return hypervolume of simplex
    double getVolume(std::vector<V<double>> vs);
    //! @param m matrix
    //! @return determinant of matrix
    double det(V<V<double>> m);
};

class HyperPolytope{
    //! @param vs collection of vertices
    //! @param eng engin for random numbers
    //! @return rnadom vector inside Polytopes using representation theorem
    V<double> generateRandomVector(std::vector<V<double>> vs, std::default_random_engine &eng);
};

#endif