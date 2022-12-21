#ifndef V_HPP
#define V_HPP

#include <vector>
#include <iostream>

//Template classes cannot be initialized outside of header

template<typename T>
class V{
    public: 
    std::vector<T> v;
    
    V<T>()=default;

    V<T>(std::vector<T> _v): v{_v} {}

    V<T>(std::initializer_list<T> const & _v): v{_v} {}

    V<T>(int size, T value){
        v = std::vector<T>(size,value);
    }

    V<T>& operator += (const V<T>& rhs){
        for(int i = 0; i< v.size(); i++){
            this->v[i] += rhs.v[i];
        }
        return *this;
    } 

    V<T> operator + (V<T> const &r){
        V<T> tmp = *this;
        tmp+=r;
        return tmp;
    }

    V<T>& operator += (T a){
        for(int i = 0; i< v.size(); i++){
            this->v[i] += a;
        }
        return *this;
    } 

    V<T> operator+ (T a){
        V<T> tmp = *this;
        tmp+=a;
        return tmp;
    }

    V<T>& operator -= (const V<T>& rhs){
        for(int i = 0; i< v.size(); i++){
            this->v[i] -= rhs.v[i];
        }
        return *this;
    } 

    V<T> operator - (V<T> const &r){
        V<T> tmp = *this;
        tmp-=r;
        return tmp;
    }

    V<T>& operator-= (T a){
        for(int i = 0; i< v.size(); i++){
            this->v[i] -= a;
        }
        return *this;
    }

    V<T> operator - (T a){
        V<T> tmp = *this;
        tmp-=a;
        return tmp;
    }

    V<T> operator* (T a){
        V<T> tmp = *this;
        for(int i = 0; i< v.size(); i++){
            tmp.v[i] *= a;
        }
        return tmp;
    }

    V<T> operator/ (T a){
        V<T> tmp = *this;
        for(int i = 0; i< v.size(); i++){
            tmp.v[i] /= a;
        }
        return tmp;
    }

    T & operator[] (int i){
        return v[i];
    }

    T operator[] (int i) const{
        return v[i];
    }

    int size(){
        return v.size();
    }

    auto begin(){
        return v.begin();
    }

    auto end(){
        return v.end();
    }

    template <class... Args>
    void emplace_back(Args&&... args){
        v.emplace_back(args ...);
    }

    void clear(){
        v.clear();
    }

    void pop_back(){
        v.pop_back();
    }

    void print(){
        for(T el : v){
            std::cout << " " << el;
        }
        std::cout << std::endl;
    }

    double sum(){
        double sum = 0;
        for(int i = 0; i < size(); i++){
            sum += v[i];
        }
        return sum;
    }

    double norm(){
        double norm = 0;
        for(int i = 0; i < size(); i++){
            norm += v[i] * v[i];
        }
        return norm;
    }
    // I might generalize the former 2 methods with "norm_p(int p)", but it is not useful
};

#endif