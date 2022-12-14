#ifndef RANDOM_VECTOR_HPP
#define RANDOM_VECTOR_HPP

#include <random>
#include <vector>
#include "V.hpp"

class Geometry{
    public Geometry(){}
    V<double> generateRandomVector() = 0;
    double getVolume() = 0;
}

class HyperRectangle : Geometry{
    public:
    // Create an hyper rectangle distribution and fills distributions
    //! @param engine_ to initialize engine
    //! @param boundaries_ in pairs [a1,b1,a2,b2,...,an,bn]
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

class HyperBall : Geometry{
    public:
    // Create an HyperBall via center and radius
    HyperBall(V<double> center_, double radius_, std::default_random_engine engine_);
    // http://compneuro.uwaterloo.ca/files/publications/voelker.2017.pdf
    //! @param c center
    //! @param r radius
    //! @return random vector within hyperball
    V<double> generateRandomVector();
    // It's very inefficient for n>>0
    //! @param c center
    //! @param r radius
    //! @return random vector within hyperball
    V<double> generateRandomVectorInefficient();
    //For low dimension it uses recursive definition, otherwise it use Stirling approximation
    //! @param r radius
    //! @param n dimension
    //! @return hypervolume of hyperball
    double getVolume();
    private:
    //! @param c center
    //! @param r radius
    //! @return random vector on the hypersphere
    // Monte Carlo integration does not work on hyperspheres since they have volume equals to 0
    V<double> generateRandomVectorSphere();
    //! @param n
    //! @return (n-1)!
    double gamma(double n);
    V<double> center;
    double radius;
    std::default_random_engine engine;
};

class Simplex : Geometry{
    Simplex(std::vector<V<double>> vertices_, std::default_random_engine engine_);
    //! @param vs collection of vertices
    //! @param eng engine for random numbers
    //! @return random vector inside the simplex
    V<double> generateRandomVector();
    // https://en.wikipedia.org/wiki/Simplex#Volume
    //! @param vs collection of vertices
    //! @return hypervolume of simplex
    double getVolume();
    private:
    //! @param m matrix
    //! @return determinant of matrix
    double det(V<V<double>> m);
    std::vector<V<double>> vertices;
    std::default_random_engine engine;
};

class Polytope : Geometry{
    // It must be convex
    Polytope(std::vector<V<double>> vertices_, std::default_random_engine engine_);
    //! @param vs collection of vertices
    //! @param eng engin for random numbers
    //! @return rnadom vector inside Polytopes using representation theorem
    V<double> generateRandomVector();
    double getVolume();
    private:
    std::vector<V<double>> vertices;
    std::default_random_engine engine;
};

#endif