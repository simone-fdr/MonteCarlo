#include "RandomVector.hpp"

// TODO Check whether is right or not to pass the engine 
HyperRectangle::HyperRectangle(std::vector<double> boundaries_, std::default_random_engine engine_): boundaries{boundaries_}, engine{engine_}{
    for(int i = 0; i < boundaries.size(); i+=2){
        distribuitions.emplace_back(std::uniform_real_distribution<double>(boundaries[i], boundaries[i+1]));
    }
}

double HyperRectangle::getVolume(){
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