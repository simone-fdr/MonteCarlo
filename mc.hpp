#include "V.hpp"
#include "RandomVector.hpp"
#include <functional>


double MonteCarlo_integral(std::function<double (V<double>)> const &f, Geometry const & domain ,int sample_size);