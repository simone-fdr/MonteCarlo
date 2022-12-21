#include "RandomVector.hpp"

// RECTANGLE ===========================================

// TODO Check whether is right or not to pass the engine 
HyperRectangle::HyperRectangle(std::vector<double> const & boundaries_, std::default_random_engine engine_): boundaries{boundaries_}, engine{engine_}{
    for(int i = 0; i < boundaries.size(); i+=2){
        distribuitions.emplace_back(std::uniform_real_distribution<double>(boundaries[i], boundaries[i+1]));
    }
}

double HyperRectangle::getVolume() const{
    double hv = 1;
    for(int i = 0; i < boundaries.size(); i+=2)
        hv *= (boundaries[i+1] - boundaries[i]);
    return hv;
}

V<double> HyperRectangle::generateRandomVector(){
    V<double> rv; //random vector
    for(int i = 0; i < distribuitions.size(); i++){
        rv.emplace_back(distribuitions[i](engine));
    }
    return rv;
}

// BALL ==========================================================

HyperBall::HyperBall(V<double> center_, double radius_, std::default_random_engine engine_) :
        center(center_), radius(radius_), engine(engine_) {
            n = center.size();
        }

V<double> HyperBall::generateRandomVectorInefficient() {
    V<double> rv;
    rv = generateRandomVectorSphere();
    rv.pop_back();
    rv.pop_back();
    for(int i = 0; i < n; i++){
        rv[i] = rv[i] * radius + center[i];
    }
    return rv;
}

V<double> HyperBall::generateRandomVector() {
    V<double> rv;
    std::uniform_real_distribution<double> dis(-1,1);
    do{
        rv.clear();
        for(int i = 0; i < n; i++){
            rv.emplace_back(dis(engine));
        }
    }while(rv.norm() > 1.);
    for(int i = 0; i < n; i++){
        rv[i] = rv[i] * radius + center[i];
    }
    return rv;
}

double HyperBall::getVolume() const{
    if(n > 20){
        return (1./std::sqrt(n*M_PI))*std::pow(2*M_PI*M_E/n,n/2.0)*std::pow(radius,n);
    }
    return std::pow(radius,n)*std::pow(M_PI, n/2.)/gamma(n/2. + 1);
}

V<double> HyperBall::generateRandomVectorSphere() {
    V<double> rv;
    std::uniform_real_distribution<double> dis(-10,10); //-10 to 10 in order to avoid floating point errors while normalizing dividing by a small value
    for(int i = 0; i < n + 2; i++){
        rv.emplace_back(dis(engine));
    }
    double norm = rv.norm();
    for(int i = 0; i < n; i++){
        rv[i] = rv[i] / norm;
    }
    return rv;
}

double Geometry::gamma(double n) const{
    if(n == 1.) return 1;
    if(n == 0.5) return std::sqrt(M_PI);
    return gamma(n-1)*(n-1);
}

// SIMPLEX =============================================================

Simplex::Simplex(std::vector<V<double>> vertices_, std::default_random_engine engine_):
        vertices(vertices_), engine(engine_) {}

V<double> Simplex::generateRandomVector() {
    int i,j;
    //Inefficent check
    if(vertices.size() != (vertices[0].size() + 1)){
        std::cout << "Wrong dimension for HyperTriangle" << std::endl;
        return vertices[0].v;
    }
    V<double> rv = vertices[0];
    std::uniform_real_distribution<double> dis(0,1);
    std::vector<double> as;
    rv = vertices[0];
    // a has size n-1
    for(i = 1; i < vertices.size(); i++){
        as.emplace_back(dis(engine));
    }
    // Check in inside triangle or not
    double sum = 0;
    for(auto a:as){
        sum += a;
    }
    // If the conic combination exceed the convex combination
    if(sum > 1){
        // Reflection by the plane a1+a2+...+an=1
        double t = 1;
        for(auto a:as){
            t-=a;
        }
        t /= as.size();
        for(auto &a:as){
            a = a + 2*t;
        }
    }
    for(i = 1; i < vertices.size(); i++){
        rv += (vertices[i] - vertices[0]) * as[i-1];
    }
    return rv.v;
}

double Simplex::getVolume() const{
    V<V<double>> m{vertices};
    V<double> v0 = m[0];
    m.v.erase(vertices.begin());
    m -= v0;
    return 1./gamma(v0.size() + 1) * det(m); //https://en.wikipedia.org/wiki/Simplex#Volume
}

double Simplex::det(V<V<double>> m) const{
    double det = 1;
    for (int c = 0; c < m.size(); c++)
    {
         det = det*m[c][c];
        for (int r = c + 1; r < m.size(); r++)
        {
            double ratio = m[r][c] / m[c][c];
            for (int k = c; k < m.size(); k++)
            {
                m[r][k] = m[r][k] - ratio * m[c][k];
            }
        }
    }
    return det;
}

// CONVEX POLYTOPE ======================================

Polytope::Polytope(std::vector<V<double>> vertices_, std::default_random_engine engine_) :
        vertices(vertices_), engine(engine_) {}

V<double> Polytope::generateRandomVector() {
    // Convex combination
    V<double> rv(vertices[0].size(), 0);
    std::uniform_real_distribution<double> dis{0,1};
    V<double> lambda;
    for(int i = 0; i < vertices.size(); i++){
        lambda.emplace_back(dis(engine));
    }
    double sum = lambda.sum();
    for(int i = 0; i < vertices.size(); i++){
        lambda[i] /= sum;
        rv += (vertices[i] * lambda[i]);
    }
    return rv;
}

// It is very hard to implement the volume of n-dimensional polytopes
double Polytope::getVolume() const{
    return 1;
}