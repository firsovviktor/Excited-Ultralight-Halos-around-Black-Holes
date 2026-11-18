#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <complex>
#include <cmath>
#include <vector>
#include <string>
#include <numeric>
#include <algorithm>
#include <random>
#include <optional>
#include <thread>
#include <chrono>

/*#include <matplot/matplot.h>
#include <indicators/progress_bar.hpp>
indicators::ProgressBar bar{
    indicators::option::BarWidth{50},
    indicators::option::Start{"["}, indicators::option::End{"]"},
    indicators::option::Fill{"="}, indicators::option::Lead{">"},
    indicators::option::ShowPercentage{true},
    indicators::option::ShowElapsedTime{true},
    indicators::option::ShowRemainingTime{true},
};*/

#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <boost/math/special_functions/laguerre.hpp>
#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/special_functions/bessel.hpp>

#include "coulomb.hpp"

double CONST_ALPHA = 0.0008;
double CONST_XI = 2;
const double CONST_K_MAX = 10;
thread_local std::mt19937_64 RNG;
int TYPE_OF_UNBOUND_FUNCTION = 0;

struct Vec3d{
    double mag, theta, phi;
    double x, y, z;

    Vec3d(double mag=0, double theta=0, double phi=0){
        this->mag = mag;
        this->theta = theta;
        this->phi = phi;
        this->x = mag * std::sin(theta) * std::cos(phi);
        this->y = mag * std::sin(theta) * std::sin(phi);
        this->z = mag * std::cos(theta);
    }
    Vec3d(double mag, double theta, double phi, double x, double y, double z){
        this->mag = mag;
        this->theta = theta;
        this->phi = phi;
        this->x = x;
        this->y = y;
        this->z = z;
    }
    Vec3d normalize() const{
        //a zero-length vector has no direction, so pick z-hat: the choice is arbitrary
        //but keeps the postcondition (a genuine unit vector, no NaN in any field -
        //note Vec3d_from_Cartesian(0,0,0) leaves theta = acos(0/0) = NaN)
        if (!(mag > 0.0))
            return Vec3d(1.0, 0.0, 0.0, 0.0, 0.0, 1.0);
        return Vec3d(1.0, theta, phi, x/mag, y/mag, z/mag);
    }
    void print() const{
        //std::cout << "r: " << mag << ", theta: " << theta << ", phi: " << phi << std::endl;
        //std::cout << "x: " << x << ", y: " << y << ", z: " << z << std::endl;
    }

    double operator*(const Vec3d& other) const{
        return x*other.x + y*other.y + z*other.z;
    }
    Vec3d operator*(double scalar) const{
        return Vec3d(mag * scalar, theta, phi, x*scalar, y*scalar, z*scalar);
    }
};

Vec3d Vec3d_from_Cartesian(double x, double y, double z){
    double r = std::sqrt(x*x + y*y + z*z);
    double theta = std::acos(z/r);
    double phi = std::atan2(y, x);
    return Vec3d(r, theta, phi, x, y, z);
}

double cos_between(const Vec3d& a, const Vec3d& b){
    //a zero-length vector has no direction; returning 1 makes the k*r*(1-khat.rhat)
    //argument of the scattering wavefunction vanish, which is the correct k->0 and
    //r->0 limit. The !(x>0) form also catches a NaN magnitude.
    if (!(a.mag > 0.0) or !(b.mag > 0.0))
        return 1.0;
    return std::clamp((a*b)/(a.mag*b.mag), -1.0, 1.0);
}

Vec3d generate_random_on_shell_momentum(){
    std::normal_distribution<double> dist_xy(0.0, 1.0/std::sqrt(2.0));
    std::uniform_real_distribution<double> dist_phi_xy(0.0, 2*M_PI);
    std::normal_distribution<double> dist_z(1.0, 1.0/std::sqrt(2.0));

    double xy = std::abs(dist_xy(RNG));
    double phi_xy = dist_phi_xy(RNG);
    double z = dist_z(RNG);

    return Vec3d_from_Cartesian(xy*std::cos(phi_xy), xy*std::sin(phi_xy), z);
}

Vec3d generate_random_off_shell_momentum(){
    std::uniform_real_distribution<double> dist_r(0, CONST_K_MAX);
    std::normal_distribution<double> dist_x(0.0, 1.0);
    std::normal_distribution<double> dist_y(0.0, 1.0);
    std::normal_distribution<double> dist_z(0.0, 1.0);

    double mag = dist_r(RNG);
    double x = dist_x(RNG);
    double y = dist_y(RNG);
    double z = dist_z(RNG);

    return Vec3d_from_Cartesian(x, y, z).normalize()*mag;
}

Vec3d generate_random_fixed_energy_momentum(double energy){
    std::normal_distribution<double> dist_x(0.0, 1.0);
    std::normal_distribution<double> dist_y(0.0, 1.0);
    std::normal_distribution<double> dist_z(0.0, 1.0);

    double mag = std::sqrt(2.0*energy)*CONST_XI;
    double x = dist_x(RNG);
    double y = dist_y(RNG);
    double z = dist_z(RNG);

    return Vec3d_from_Cartesian(x, y, z).normalize()*mag;   
}

std::complex<double> hydrogen_bound_wavefunction(int n, int l, int m, Vec3d r_vec){
    // r measured in Bohr radii, theta and phi in radians
    // factor of 1/r_Bohr^{3/2} is omitted
    if (n==0){
        return std::complex<double>(1.0, 0.0);
    }
    double radial_part = std::sqrt(std::pow(2.0/n, 3)) * std::pow(2.0*r_vec.mag/n, l) * std::exp(-r_vec.mag/n) *
                         std::sqrt(boost::math::factorial<double>(n - l - 1) / (2.0 * n * boost::math::factorial<double>(n + l))) *
                         boost::math::laguerre(n-l-1, 2*l+1, 2.0*r_vec.mag/n);

    std::complex<double> angular_part;
    if (l == 0){
        angular_part = std::complex<double>(0.5/std::sqrt(M_PI), 0);
    }else if (l==1){
        if (m == -1)
            angular_part = 0.5*std::sqrt(1.5/M_PI)*std::sin(r_vec.theta)*std::exp(std::complex<double>(0, -r_vec.phi));  
        if (m == 0)
            angular_part = 0.5*std::sqrt(3/M_PI)*std::cos(r_vec.theta);
        if (m == 1)
            angular_part =-0.5*std::sqrt(1.5/M_PI)*std::sin(r_vec.theta)*std::exp(std::complex<double>(0, r_vec.phi));         
    }else if (l==2){
    	if (m == -2)
	    angular_part = 0.25*std::sqrt(7.5/M_PI)*std::sin(r_vec.theta)*std::sin(r_vec.theta)*std::exp(std::complex<double>(0, -2*r_vec.phi));
	if (m == -1)
	    angular_part = 0.5 *std::sqrt(7.5/M_PI)*std::sin(r_vec.theta)*std::cos(r_vec.theta)*std::exp(std::complex<double>(0, -r_vec.phi));
	if (m == 0)
	    angular_part = 0.25*std::sqrt(5/M_PI)*(3*std::cos(r_vec.theta)*std::cos(r_vec.theta)-1);
	if (m == 1)
            angular_part =-0.5 *std::sqrt(7.5/M_PI)*std::sin(r_vec.theta)*std::cos(r_vec.theta)*std::exp(std::complex<double>(0, r_vec.phi));
	if (m == 2)
            angular_part = 0.25*std::sqrt(7.5/M_PI)*std::sin(r_vec.theta)*std::sin(r_vec.theta)*std::exp(std::complex<double>(0, 2*r_vec.phi));
    }else{
        angular_part = boost::math::spherical_harmonic(l, m, r_vec.theta, r_vec.phi);
    }

    return radial_part * angular_part;
}

std::complex<double> approximate_hydrogen_unbound_wavefunction(Vec3d k_vec, Vec3d r_vec){
    // r_vec.mag measured in Bohr radii, theta and phi in radians
    // k_vec.mag measured in m*v_dm, theta and phi in radians
    // factor of 1/sqrt{|k|*r_Bohr} omitted
    // DOESN'T WORK FOR k_vec.mag>=xi AS APPROXIMATION BREAKS DOWN
    double k = k_vec.mag/CONST_XI; //now measured in 1/r_Bohr
    double phi = (std::log(std::max(1e-12, k))+1)/k - M_PI/4 + k*r_vec.mag*(1+cos_between(k_vec, r_vec))/2;
    std::complex<double> wavefunction = std::exp(std::complex<double>(0, phi)) * std::sqrt(2*M_PI) *
                                        boost::math::cyl_bessel_j(0, 2.0*std::sqrt(r_vec.mag*(1 - cos_between(k_vec, r_vec))));

    return wavefunction;
}

std::complex<double> hydrogen_unbound_wavefunction(Vec3d k_vec, Vec3d r_vec){
    // exact scattering wavefunction, eq. (154), NOT the focusing approximation:
    // psi = e^{i k.r} Gamma(1 - i/(k r_b)) e^{pi/(2 k r_b)} M(i/(k r_b), 1, i k r (1-khat.rhat))
    // r_vec.mag measured in Bohr radii, theta and phi in radians
    // k_vec.mag measured in m*v_dm, theta and phi in radians
    // factor of 1/sqrt{|k|*r_Bohr} omitted, same convention as the approximation
    // above, so this reduces to it for k_vec.mag<<xi and to a plane wave for
    // k_vec.mag>>xi - no restriction on k any more
    double k = k_vec.mag/CONST_XI; //now measured in 1/r_Bohr
    return coulomb::ScatteringState(std::max(1e-12, k)).psi(r_vec.mag, cos_between(k_vec, r_vec));
}

struct State{
    bool is_bound;
    bool is_on_shell;
    std::optional<int> n, l, m;
    std::optional<Vec3d> momentum;
    std::optional<coulomb::ScatteringState> scattering; //per-momentum cache, see evaluate()

    State(bool is_bound=true, bool is_on_shell=true){
        this->is_bound = is_bound;
        this->is_on_shell = is_on_shell;
    }
    State(bool is_bound, bool is_on_shell, int n, int l, int m){
        this->is_bound = is_bound;
        this->is_on_shell = is_on_shell;
        this->n = n;
        this->l = l;
        this->m = m;
    }
    State(bool is_bound, bool is_on_shell, Vec3d momentum){
        this->is_bound = is_bound;
        this->is_on_shell = is_on_shell;
        this->momentum = momentum;
    }
    std::complex<double> evaluate(Vec3d r_vec){
        if (is_bound) {
            return hydrogen_bound_wavefunction(n.value(), l.value(), m.value(), r_vec);
        } else {
            //the Gamma function and prefactor of eq. (154) depend on |k| only, so they
            //are built once per momentum rather than once per sample point. Keyed on k
            //so a reassigned momentum (or a changed CONST_XI) can't leave it stale.
            double k = std::max(1e-12, momentum.value().mag/CONST_XI); //in 1/r_Bohr
            if (TYPE_OF_UNBOUND_FUNCTION == 0){
                if ((!scattering) or (scattering->k != k))
                        scattering.emplace(k);
                    return scattering->psi(r_vec.mag, cos_between(momentum.value(), r_vec));
            }else if (TYPE_OF_UNBOUND_FUNCTION == 1){
                if (k>1.5){
                    if ((!scattering) or (scattering->k != k))
                        scattering.emplace(k);
                    return scattering->psi(r_vec.mag, cos_between(momentum.value(), r_vec));
                }else{
                    return approximate_hydrogen_unbound_wavefunction(momentum.value(), r_vec);
                }
            }else if (TYPE_OF_UNBOUND_FUNCTION == 2){
                return approximate_hydrogen_unbound_wavefunction(momentum.value(), r_vec);
            }else{
                std::cerr << "Unknwon unbound function type" << std::endl;
                return std::complex<double>(0.0, 0.0);
            }
        }
    }
    double energy(){
        //without particle mass, measured in units of m\alpha^2
        if (is_bound){
            //INCORRECT, fine splitting /n^3 not /n^4
            return -0.5/std::pow(n.value(), 2) + 15.0/8.0*std::pow(CONST_ALPHA, 2)/std::pow(n.value(), 4) - 6.0/(2.0*l.value()+1.0)*std::pow(CONST_ALPHA, 2)/std::pow(n.value(), 3);
        }else{
        //without particle mass, measured in units of m\alpha^2
        //assuming momentum in m*v_dm
            return std::pow(momentum.value().mag, 2)/2.0/std::pow(CONST_XI, 2);
        }
    }
    double weight_factor(){
        //weight factor for j<>1 states
        if (is_on_shell){
            return std::pow(M_PI, 2)*std::sqrt(std::pow(momentum.value().x, 2)+std::pow(momentum.value().y, 2))/momentum.value().mag;
        }else{
            return 4.0*M_PI*CONST_K_MAX*momentum.value().mag;
        }
    }
    double fixed_energy_weight_factor(){
        //weight factor for j=1 states
        if (is_on_shell){
            double k = momentum.value().mag;
            return std::exp(-(k*k - 2.0*momentum.value().z + 1.0));
        }else{
            return 1.0;
        }
    }
    void print(){
        if (is_bound){
            //std::cout << "bound   ";
        }else{
            //std::cout << "unbound ";
        }
        if (is_on_shell){
            //std::cout << "on shell ";
        }else{
            //std::cout << "dashed   ";
        }
        if (is_bound){
            //std::cout << n.value() << " " << l.value() << " " << m.value() << std::endl;
        }else{
            momentum.value().print();
        }
    }
};

struct Matrix_element{
    State in_state_1, in_state_2, out_state_1, out_state_2;
    double weight;

    Matrix_element(){
        this->in_state_1 = State();
        this->in_state_2 = State();
        this->out_state_1 = State();
        this->out_state_2 = State();
        weight = 1.0;
    }

    Matrix_element(const State& in_state_1, const State& in_state_2, 
                   const State& out_state_1, const State& out_state_2,
                   double weight)
                 : in_state_1(in_state_1), 
                   in_state_2(in_state_2), 
                   out_state_1(out_state_1), 
                   out_state_2(out_state_2),
                   weight(weight) {}

    int get_lowest_n(){
        int scale = -1;
        if ((in_state_1.is_bound)and(in_state_1.n > 0)){
            if ((scale == -1) or (in_state_1.n < scale)){
                scale = in_state_1.n.value();
            }
        }
        if ((in_state_2.is_bound)and(in_state_2.n > 0)){
            if ((scale == -1) or (in_state_2.n < scale)){
                scale = in_state_2.n.value();
            }
        }
        if ((out_state_1.is_bound)and(out_state_1.n > 0)){
            if ((scale == -1) or (out_state_1.n < scale)){
                scale = out_state_1.n.value();
            }
        }
        if ((out_state_2.is_bound)and(out_state_2.n > 0)){
            if ((scale == -1) or (out_state_2.n < scale)){
                scale = out_state_2.n.value();
            }
        }
        return scale;
    }

    std::pair<std::complex<double>, std::complex<double>> evaluate(int N_points_per_matrix_element, double range){
        double total_real_part = 0.0;
        double total_imag_part = 0.0;
        double total_square_real_part = 0.0;
        double total_square_imag_part = 0.0;
        int lowest_n = get_lowest_n();

        std::uniform_real_distribution<double> dist_uniform(-range, range);
        std::exponential_distribution<double> dist_magnitude(1.0/lowest_n);
        std::uniform_real_distribution<double> dist_cos_theta(-1.0, 1.0);
        std::uniform_real_distribution<double> dist_phi(0.0, 2.0 * M_PI);

        Vec3d r_sample(0, 0, 0);
        
        for (int i=0; i<N_points_per_matrix_element; i++){
            if (lowest_n < 0.5){
                double x, y, z;
                x = dist_uniform(RNG);
                y = dist_uniform(RNG);
                z = dist_uniform(RNG);
                r_sample = Vec3d_from_Cartesian(x, y, z);
            }else{
                double mag, theta, phi;
                mag = dist_magnitude(RNG);
                theta = std::acos(dist_cos_theta(RNG));
                phi = dist_phi(RNG);
                r_sample = Vec3d(mag, theta, phi);
            }
            
            std::complex<double> sample_value = in_state_1.evaluate(r_sample)*
                                                in_state_2.evaluate(r_sample)*
                                                std::conj(out_state_1.evaluate(r_sample))*
                                                std::conj(out_state_2.evaluate(r_sample));

            if (lowest_n > 0){
                sample_value *= 4*M_PI*lowest_n* std::pow(r_sample.mag, 2) * std::exp(r_sample.mag/lowest_n);
            }
            total_real_part += sample_value.real();
            total_imag_part += sample_value.imag();
            total_square_real_part += sample_value.real() * sample_value.real();
            total_square_imag_part += sample_value.imag() * sample_value.imag();
        }
        std::complex<double> mean_value(total_real_part/N_points_per_matrix_element, total_imag_part/N_points_per_matrix_element);
        double real_variance = total_square_real_part/N_points_per_matrix_element - mean_value.real()*mean_value.real();
        double imag_variance = total_square_imag_part/N_points_per_matrix_element - mean_value.imag()*mean_value.imag();
        std::complex<double> error_value(std::sqrt(std::max(0.0, real_variance))/std::sqrt(N_points_per_matrix_element), 
                                         std::sqrt(std::max(0.0, imag_variance))/std::sqrt(N_points_per_matrix_element));
        if (lowest_n < 1){
            double volume = std::pow(2*range, 3);
            return std::make_pair(mean_value*volume, error_value*volume);
        }else{
            return std::make_pair(mean_value, error_value);
        }
    }

    void fill_one_momentum(){
        double energy_dif;
        energy_dif = in_state_1.energy()+in_state_2.energy()-out_state_1.energy()-out_state_2.energy();
        energy_dif = std::abs(energy_dif);

        if (!in_state_1.is_bound){
            in_state_1.momentum = generate_random_fixed_energy_momentum(energy_dif);
            weight *= in_state_1.fixed_energy_weight_factor();
        }

        if (!in_state_2.is_bound){
            in_state_2.momentum = generate_random_fixed_energy_momentum(energy_dif);
            weight *= in_state_2.fixed_energy_weight_factor();
        }

        if (!out_state_1.is_bound){
            out_state_1.momentum = generate_random_fixed_energy_momentum(energy_dif);
            weight *= out_state_1.fixed_energy_weight_factor();
        }

        if (!out_state_2.is_bound){
            out_state_2.momentum = generate_random_fixed_energy_momentum(energy_dif);
            weight *= out_state_2.fixed_energy_weight_factor();
        }
    }

    void fill_two_momenta(){
        double energy_dif;
        energy_dif = in_state_1.energy()+in_state_2.energy()-out_state_1.energy()-out_state_2.energy();

        if (energy_dif >= 0){
            if(!in_state_1.is_bound){
                in_state_1.momentum = generate_random_on_shell_momentum();
                energy_dif += in_state_1.energy();
                weight *= in_state_1.weight_factor();
            }
            if(!in_state_2.is_bound){
                in_state_2.momentum = generate_random_on_shell_momentum();
                energy_dif += in_state_2.energy();
                weight *= in_state_2.weight_factor();
            }
            if(!out_state_1.is_bound){
                out_state_1.momentum = generate_random_fixed_energy_momentum(energy_dif);
                weight *= out_state_1.fixed_energy_weight_factor();
            }
            if(!out_state_2.is_bound){
                out_state_2.momentum = generate_random_fixed_energy_momentum(energy_dif);
                weight *= out_state_2.fixed_energy_weight_factor();
            }
        }else{
            if(!out_state_1.is_bound){
                out_state_1.momentum = generate_random_on_shell_momentum();
                energy_dif -= out_state_1.energy();
                weight *= out_state_1.weight_factor();
            }
            if(!out_state_2.is_bound){
                out_state_2.momentum = generate_random_off_shell_momentum();
                energy_dif -= out_state_2.energy();
                weight *= out_state_2.weight_factor();
            }
            if(!in_state_1.is_bound){
                in_state_1.momentum = generate_random_fixed_energy_momentum(-1.0*energy_dif);
                weight *= in_state_1.fixed_energy_weight_factor();
            }
            if(!in_state_2.is_bound){
                in_state_2.momentum = generate_random_fixed_energy_momentum(-1.0*energy_dif);
                weight *= in_state_2.fixed_energy_weight_factor();
            }
        }
    }

    void fill_three_momenta(){
        if ((in_state_1.is_bound)or(in_state_2.is_bound)){
            out_state_1.momentum = generate_random_on_shell_momentum();
            weight *= out_state_1.weight_factor();
            out_state_2.momentum = generate_random_off_shell_momentum();
            weight *= out_state_2.weight_factor();
            double energy_dif = in_state_1.energy() + in_state_2.energy() - out_state_1.energy() - out_state_2.energy();
            if (!in_state_1.is_bound){
                in_state_1.momentum = generate_random_fixed_energy_momentum(-1.0*energy_dif);
                weight *= in_state_1.fixed_energy_weight_factor();
            }else{
                in_state_2.momentum = generate_random_fixed_energy_momentum(-1.0*energy_dif);
                weight *= in_state_2.fixed_energy_weight_factor();
            }
        }
        if ((out_state_1.is_bound)or(out_state_2.is_bound)){
            in_state_1.momentum = generate_random_on_shell_momentum();
            weight *= in_state_1.weight_factor();
            in_state_2.momentum = generate_random_on_shell_momentum();
            weight *= in_state_2.weight_factor();
            double energy_dif = in_state_1.energy() + in_state_2.energy() - out_state_1.energy() - out_state_2.energy();
            if (!out_state_1.is_bound){
                out_state_1.momentum = generate_random_fixed_energy_momentum(energy_dif);
                weight *= out_state_1.fixed_energy_weight_factor();
            }else{
                out_state_2.momentum = generate_random_fixed_energy_momentum(energy_dif);
                weight *= out_state_2.fixed_energy_weight_factor();
            }
        }
    }

    void fill(){
        int number_of_unbound_pairs = 0;
        if (!in_state_1.is_bound) 
            number_of_unbound_pairs++;
        if (!in_state_2.is_bound) 
            number_of_unbound_pairs++;
        if (!out_state_1.is_bound) 
            number_of_unbound_pairs++;
        if (!out_state_2.is_bound) 
            number_of_unbound_pairs++;

        if (number_of_unbound_pairs == 0){
            return;
        }else{
            weight = 4.0*M_PI;
            if (number_of_unbound_pairs == 1){
                fill_one_momentum();
            }else if (number_of_unbound_pairs == 2){
                fill_two_momenta();
            }else if (number_of_unbound_pairs == 3){
                fill_three_momenta();
            }else{
                std::cerr << "undefined behaviour for Matrix_element.fill(int)" << std::endl;
                std::cerr << "number of unbound states should be 0..3" << std::endl;
            }
        }
    }

    void print(){
        in_state_1.print();
        in_state_2.print();
        out_state_1.print();
        out_state_2.print();
        //std::cout << in_state_1.energy()+in_state_2.energy()-out_state_1.energy()-out_state_2.energy() << " = ";
        //std::cout << in_state_1.energy() << " + " << in_state_2.energy() << " - ";
        //std::cout << out_state_1.energy() << " - " << out_state_2.energy() << std::endl;
        //std::cout << std::endl;
    }
};

std::pair<Matrix_element, Matrix_element> generate_diagram(std::string input_string){
    // Assumes n,l,m<10
    std::string diagram_template = input_string;
    int delimiter = 0;

    delimiter = diagram_template.find("+");
    std::string first_in_state_1_string = diagram_template.substr(0, delimiter);
    diagram_template.erase(0, delimiter+1);

    delimiter = diagram_template.find("-");
    std::string first_in_state_2_string = diagram_template.substr(0, delimiter);
    diagram_template.erase(0, delimiter+1);

    delimiter = diagram_template.find("=");
    std::string first_out_state_1_string = diagram_template.substr(0, delimiter);
    diagram_template.erase(0, delimiter+1);

    delimiter = diagram_template.find("=");
    std::string first_out_state_2_string = diagram_template.substr(0, delimiter);
    std::string second_in_state_2_string = first_out_state_2_string;
    diagram_template.erase(0, delimiter+1);

    delimiter = diagram_template.find("+");
    std::string second_out_state_1_string = diagram_template.substr(0, delimiter);
    diagram_template.erase(0, delimiter+1);
    
    delimiter = diagram_template.find("-");
    std::string second_out_state_2_string = diagram_template.substr(0, delimiter);
    diagram_template.erase(0, delimiter+1);

    std::string second_in_state_1_string = diagram_template;


    Matrix_element first_matrix_element;
    Matrix_element second_matrix_element;

    if (first_in_state_1_string[0] == 'k'){
        first_matrix_element.in_state_1 = State(false, true, Vec3d());
    }else{
        first_matrix_element.in_state_1 = State(true, true, int(first_in_state_1_string[0])-int('0'),
                                                            int(first_in_state_1_string[1])-int('0'),
                                                            int(first_in_state_1_string[2])-int('0'));
    } 
    if (first_in_state_2_string[0] == 'k'){
        first_matrix_element.in_state_2 = State(false, true, Vec3d());
    }else{
        first_matrix_element.in_state_2 = State(true, true, int(first_in_state_2_string[0])-int('0'),
                                                            int(first_in_state_2_string[1])-int('0'),
                                                            int(first_in_state_2_string[2])-int('0'));
    }
    if (first_out_state_1_string[0] == 'k'){
        first_matrix_element.out_state_1 = State(false, true, Vec3d());
    }else{
        first_matrix_element.out_state_1 = State(true, true, int(first_out_state_1_string[0])-int('0'),
                                                            int(first_out_state_1_string[1])-int('0'),
                                                            int(first_out_state_1_string[2])-int('0'));
    } 
    if (first_out_state_2_string[0] == 'k'){
        first_matrix_element.out_state_2 = State(false, false, Vec3d());
    }else{
        first_matrix_element.out_state_2 = State(true, false, int(first_out_state_2_string[0])-int('0'),
                                                            int(first_out_state_2_string[1])-int('0'),
                                                            int(first_out_state_2_string[2])-int('0'));
    } 

    if (second_in_state_1_string[0] == 'k'){
        second_matrix_element.in_state_1 = State(false, true, Vec3d());
    }else{
        second_matrix_element.in_state_1 = State(true, true, int(second_in_state_1_string[0])-int('0'),
                                                            int(second_in_state_1_string[1])-int('0'),
                                                            int(second_in_state_1_string[2])-int('0'));
    } 
    if (second_in_state_2_string[0] == 'k'){
        second_matrix_element.in_state_2 = State(false, false, Vec3d());
    }else{
        second_matrix_element.in_state_2 = State(true, false, int(second_in_state_2_string[0])-int('0'),
                                                            int(second_in_state_2_string[1])-int('0'),
                                                            int(second_in_state_2_string[2])-int('0'));
    }
    if (second_out_state_1_string[0] == 'k'){
        second_matrix_element.out_state_1 = State(false, true, Vec3d());
    }else{
        second_matrix_element.out_state_1 = State(true, true, int(second_out_state_1_string[0])-int('0'),
                                                            int(second_out_state_1_string[1])-int('0'),
                                                            int(second_out_state_1_string[2])-int('0'));
    } 
    if (second_out_state_2_string[0] == 'k'){
        second_matrix_element.out_state_2 = State(false, true, Vec3d());
    }else{
        second_matrix_element.out_state_2 = State(true, true, int(second_out_state_2_string[0])-int('0'),
                                                            int(second_out_state_2_string[1])-int('0'),
                                                            int(second_out_state_2_string[2])-int('0'));
    }


    first_matrix_element.fill();

    Vec3d momenta[4];
    //save momenta of the first vertex according to their indices in the template
    if (!first_matrix_element.in_state_1.is_bound) 
        momenta[int(first_in_state_1_string[1])-int('1')] = first_matrix_element.in_state_1.momentum.value();
    if (!first_matrix_element.in_state_2.is_bound) 
        momenta[int(first_in_state_2_string[1])-int('1')] = first_matrix_element.in_state_2.momentum.value();
    if (!first_matrix_element.out_state_1.is_bound) 
        momenta[int(first_out_state_1_string[1])-int('1')] = first_matrix_element.out_state_1.momentum.value();
    if (!first_matrix_element.out_state_2.is_bound) 
        momenta[int(first_out_state_2_string[1])-int('1')] = first_matrix_element.out_state_2.momentum.value();
    //assign the second vertex momenta in accordance to the template
    if (!second_matrix_element.in_state_1.is_bound) 
        second_matrix_element.in_state_1.momentum = momenta[int(second_in_state_1_string[1])-int('1')];
    if (!second_matrix_element.in_state_2.is_bound) 
        second_matrix_element.in_state_2.momentum = momenta[int(second_in_state_2_string[1])-int('1')];
    if (!second_matrix_element.out_state_1.is_bound) 
        second_matrix_element.out_state_1.momentum = momenta[int(second_out_state_1_string[1])-int('1')];
    if (!second_matrix_element.out_state_2.is_bound) 
        second_matrix_element.out_state_2.momentum = momenta[int(second_out_state_2_string[1])-int('1')];    


    return std::make_pair(first_matrix_element, second_matrix_element);
}

std::string generate_timestamped_filename(const std::string& prefix, const std::string& ext) {
    auto now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    std::stringstream ss;
    
    // Format timestamp
    ss << prefix << "_" 
       << std::put_time(std::localtime(&now), "%Y-%m-%d_%H-%M-%S") 
       << ext;
       
    return ss.str();
}

int main(){
    // Test for the hydrogen bound wavefunction
    /*
    double max_R = 10.0;
    double step = 0.1;

    std::vector<double> x_values;
    std::vector<double> y_values;

    for (double x = -max_R; x <= max_R; x += step){
        x_values.push_back(x);
        std::complex<double> psi = hydrogen_bound_wavefunction(2, 1, 0, Vec3d(std::abs(x), 0, 0));
        y_values.push_back(std::abs(psi));
    }
    auto f = matplot::gcf();
    f->backend()->run_command("unset warnings");
    matplot::plot(x_values, y_values);
    matplot::show();
    */
    
    // Test 1 for the hydrogen unbound wavefunction
    /*
    double max_R = 2.0;
    double step = 0.04;
    Vec3d momentum = Vec3d_from_Cartesian(-1, 0, 0);

    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> c;
    for (double cur_x = -max_R; cur_x <= max_R; cur_x += step) {
        for (double cur_y = -max_R; cur_y <= max_R; cur_y += step) {
            if (std::sqrt(cur_x*cur_x + cur_y*cur_y) < step*step) {
                continue;
            }
            Vec3d coordinate = Vec3d_from_Cartesian(cur_x, cur_y, 0);

            coordinate.print();
            //std::cout << hydrogen_unbound_wavefunction(momentum, coordinate) << std::endl;

            x.push_back(coordinate.x);
            y.push_back(coordinate.y);
            c.push_back(std::abs(hydrogen_unbound_wavefunction(momentum, coordinate)));
        }
    }
    auto f = matplot::gcf();
    f->backend()->run_command("unset warnings");
    auto l = matplot::scatter(x, y, 5, c);
    l->marker_face(true);
    matplot::colorbar();
    matplot::colormap(matplot::palette::viridis());
    matplot::show();
    */
    
    // Test 2 for the hydrogen unbound wavefunction
    /*
    double max_R = 10;
    double step = 0.04;
    Vec3d momentum = Vec3d_from_Cartesian(-1, 0, 0);

    std::vector<double> x;
    std::vector<std::complex<double>> exact;
    std::vector<std::complex<double>> approx;

    for (double cur_x = -max_R; cur_x<max_R; cur_x += step){
        Vec3d coordinate = Vec3d_from_Cartesian(cur_x, 0, 0);

        x.push_back(cur_x);
        exact.push_back(hydrogen_unbound_wavefunction(momentum, coordinate));
        approx.push_back(approximate_hydrogen_unbound_wavefunction(momentum, coordinate));
    }
    
    std::vector<double> y_values_1;
    std::vector<double> y_values_2;
    std::vector<double> z_values_1;
    std::vector<double> z_values_2;
    for (int i=0; i<x.size(); i++){
        y_values_1.push_back(exact[i].real());
        y_values_2.push_back(approx[i].real());
        z_values_1.push_back(exact[i].imag());
        z_values_2.push_back(approx[i].imag());
    }


    auto f = matplot::gcf(); 
    // Set the width and height (in pixels)
    f->width(1200);  // Makes the window wider for side-by-side plots
    f->height(600);

    // Create a 1x2 grid, and set the active plot to the 1st position (left)
    matplot::subplot(1, 2, 1);  
    matplot::plot(x, y_values_1, x, y_values_2);
    auto l = matplot::legend({"Exact re", "Approx re"});
    l->font_size(18);

    // Set the active plot to the 2nd position (right)
    matplot::subplot(1, 2, 2);
    matplot::plot(x, z_values_1, x, z_values_2);
    auto m = matplot::legend({"Exact im", "Approx im"});
    m->font_size(18);

    // Render both plots in the same window
    matplot::show();
    */
   
    // Test for the matrix element computation
    /*
    int N_points_per_matrix_element = 10000000;
    double range = 20.0; //20 is enough up to (and including) n=3

    std::string diagram_template = "k1+k2-100=k3=k1+k2-100";
    auto [old_first_matrix_element, old_second_matrix_element] = generate_diagram(diagram_template);
    auto [first_matrix_element, second_matrix_element] = generate_diagram(diagram_template); 

    for (int i=0; i<10; i++) {
        auto [result, error] = second_matrix_element.evaluate(N_points_per_matrix_element, range);
        //std::cout << "///" << std::endl;
        second_matrix_element.print();
        //std::cout << i << " Matrix element: " << result << ", error: " << error << std::endl;
        //std::cout << "///" << std::endl;
        //std::cout << std::endl;
        //std::cout << std::endl;
    }
    */

    //Programme parameters
    auto start = std::chrono::steady_clock::now();
    const char* LSF_idx = std::getenv("LSB_JOBINDEX");
    int LSF_task = LSF_idx ? std::atoi(LSF_idx) : 0;

    std::string log_dir = "/home/labs/gperez/viktorf/logs";
    std::string filename = generate_timestamped_filename("simulation", ".log");
    std::ofstream log_file(log_dir+"/"+filename+"_"+std::to_string(LSF_task), std::ios::app);
    if (!log_file.is_open()){
        std::cerr << "Cannot open " << log_dir << "/" << filename << " - aborting" << std::endl;
        return 1;
    }

    bool easy_mode = true;
    if (const char* env_easy = std::getenv("EASY_MODE")){
        easy_mode = (std::atoi(env_easy) != 0);
    }

    std::vector<double> vector_critical_xis = {0.2, 1.5, 4.5, 1e9};
    std::vector<int> vector_N_matrix_elements = {(int) 3e7, (int) 1e7, (int) 5e6, (int)1e6};
    //vector_N_matrix_elements = {(int) 6e6, (int) 2e6, (int) 1e6, (int)1e6};
    std::vector<int> vector_N_points_per_matrix_element = {(int) 1e3, (int) 1e3, (int) 5e2, (int)5e2};
    int N_matrix_elements;
    int N_points_per_matrix_element;

    double range = 20.0; //20 is enough up to (and including) n=3
    double matrix_element_tolerance = 0.05;
    //to denote negative m use special symbols
    // -1='/', -2='.' lower doesn't work
    //a+b-c=d=a+b-c
    
    std::vector<std::string> vector_diagram_template_2p = {"211+211-k1=100=211+211-k1", "210+210-k1=100=210+210-k1", "21/+21/-k1=100=21/+21/-k1",
                                                        "211+210-k1=100=211+210-k1", "211+21/-k1=100=211+21/-k1", "210+21/-k1=100=210+21/-k1",
                                                        "211+21/-k1=100=210+210-k1",

                                                        "211+211-100=k1=211+211-100", "210+210-100=k1=210+210-100", "21/+21/-100=k1=21/+21/-100",
                                                        "211+210-100=k1=211+210-100", "211+21/-100=k1=211+21/-100", "210+21/-100=k1=210+21/-100",
                                                        "211+21/-100=k1=210+210-100",

                                                        "100+k1-211=211=100+k1-211", "100+k1-210=211=100+k1-210", "100+k1-21/=211=100+k1-21/",
                                                        "100+k1-211=210=100+k1-211", "100+k1-210=210=100+k1-210", "100+k1-21/=210=100+k1-21/",
                                                        "100+k1-211=21/=100+k1-211", "100+k1-210=21/=100+k1-210", "100+k1-21/=21/=100+k1-21/",
                                                    
                                                        "100+k1-211=200=100+k1-211", "100+k1-210=200=100+k1-210", "100+k1-21/=200=100+k1-21/",

                                                        "200+k1-211=k2=200+k1-211", "200+k1-210=k2=200+k1-210", "200+k1-21/=k2=200+k1-21/",

                                                        "k1+k2-211=k3=k1+k2-211", "k1+k2-210=k3=k1+k2-210", "k1+k2-21/=k3=k1+k2-21/",

                                                        "211+200-100=k1=211+200-100", "210+200-100=k1=210+200-100", "21/+200-100=k1=21/+200-100",

                                                        "211+200-k1=100=211+200-k1", "210+200-k1=100=210+200-k1", "21/+200-k1=100=21/+200-k1",

                                                        "211+k1-100=k2=211+k1-100", "210+k1-100=k2=210+k1-100", "21/+k1-100=k2=21/+k1-100",
                                                        "211+k1-200=k2=211+k1-200", "210+k1-200=k2=210+k1-200", "21/+k1-200=k2=21/+k1-200",
                                                        "211+k1-k3=k2=211+k1-k3", "210+k1-k3=k2=210+k1-k3", "21/+k1-k3=k2=21/+k1-k3",

                                                        "211+k1-k2=100=211+k1-k2", "210+k1-k2=100=210+k1-k2", "21/+k1-k2=100=21/+k1-k2", 
                                                        "211+k1-k2=200=211+k1-k2", "210+k1-k2=200=210+k1-k2", "21/+k1-k2=200=21/+k1-k2",

                                                        "100+k1-200=211=100+k1-200", "100+k1-200=210=100+k1-200", "100+k1-200=21/=100+k1-200",

                                                        "100+k1-k2=211=100+k1-k2", "100+k1-k2=210=100+k1-k2", "100+k1-k2=21/=100+k1-k2",
                                                        "200+k1-k2=211=200+k1-k2", "200+k1-k2=210=200+k1-k2", "200+k1-k2=21/=200+k1-k2",
                                                        "k3+k1-k2=211=k3+k1-k2", "k3+k1-k2=210=k3+k1-k2", "k3+k1-k2=21/=k3+k1-k2"
                                                        };
    
    std::vector<std::string> vector_diagram_template_1s2s = {"200+200-100=k1=200+200-100",

                                                        "200+200-k1=100=200+200-k1",

                                                        "k1+k2-k3=100=k1+k2-k3", "k1+k2-k3=200=k1+k2-k3",

                                                        "100+k1-200=200=100+k1-200",

                                                        "200+k1-100=k2=200+k1-100",

                                                        "k1+k2-100=k3=k1+k2-100", "k1+k2-200=k3=k1+k2-200",

                                                        "100+k1-k2=200=100+k1-k2", "200+k1-k2=100=200+k1-k2"
                                                        };

    std::vector<std::string> vector_diagram_template_forgotten = {"100+k1-k2=200=100+k1-k2",

                                                        "100+k1-200=k2=100+k1-200",
                                                        };

    std::vector<std::string> vector_diagram_template_more = {"200+k1-k2=k3=200+k1-k2",
    
                                                        "100+k1-k2=k3=100+k1-k2"
                                                        };

    std::vector<std::string> vector_diagram_template_extra = {"100+k1-211=k2=100+k1-211", "100+k1-210=k2=100+k1-210", "100+k1-21/=k2=100+k1-21/"};
    
    std::vector<std::string> vector_diagram_template_emergency_fix = {"100+k1-200=200=100+k1-200",
								"200+200-k1=100=200+200-k1",
								"200+200-100=k1=200+200-100",
								"200+k1-100=k2=200+k1-100",
								"100+k1-k2=200=100+k1-k2",
								"200+k1-k2=100=200+k1-k2",
								"100+k1-k2=200=100+k1-k2",
								"100+k1-k2=21/=100+k1-k2",
								"100+k1-k2=210=100+k1-k2",
								"200+k1-k2=210=200+k1-k2",
								"200+k1-k2=211=200+k1-k2",
								"200+k1-k2=21/=200+k1-k2",
								"k1+k2-100=k3=k1+k2-100",
								"k1+k2-k3=100=k1+k2-k3",
								"100+k1-200=k2=100+k1-200",
								"100+k1-210=k2=100+k1-210",
								"100+k1-211=k2=100+k1-211",
								"100+k1-21/=k2=100+k1-21/",
								"k1+k2-200=k3=k1+k2-200",
								"k1+k2-k3=200=k1+k2-k3",
								"100+k1-k2=k3=100+k1-k2",
								"k3+k1-k2=210=k3+k1-k2",
								"k3+k1-k2=211=k3+k1-k2",
								"k3+k1-k2=21/=k3+k1-k2",
								"200+k1-k2=k3=200+k1-k2"};

    std::vector<std::string> vector_diagram_template_high_intermediate = { 
	"211+k1-k2=300=211+k1-k2",
	"210+k1-k2=300=210+k1-k2",

	"211+k1-k2=310=211+k1-k2",
	"210+k1-k2=310=210+k1-k2",

	"211+k1-k2=320=211+k1-k2",
	"210+k1-k2=320=210+k1-k2",

	"211+k1-k2=311=211+k1-k2",
	"210+k1-k2=311=210+k1-k2",
	"21/+k1-k2=311=21/+k1-k2",

	"211+k1-k2=321=211+k1-k2",
	"210+k1-k2=321=210+k1-k2",
	"21/+k1-k2=321=21/+k1-k2",

	"211+k1-k2=322=211+k1-k2",
	"210+k1-k2=322=210+k1-k2",
	"21/+k1-k2=322=21/+k1-k2",


	"100+k1-211=300=100+k1-211",
	"100+k1-210=300=100+k1-210",

	"100+k1-211=310=100+k1-211",
	"100+k1-210=310=100+k1-210",

	"100+k1-211=320=100+k1-211",
	"100+k1-210=320=100+k1-210",

	"100+k1-211=311=100+k1-211",
	"100+k1-210=311=100+k1-210",
	"100+k1-21/=311=100+k1-21/",

	"100+k1-211=321=100+k1-211",
	"100+k1-210=321=100+k1-210",
	"100+k1-21/=321=100+k1-21/",

	"100+k1-211=322=100+k1-211",
	"100+k1-210=322=100+k1-210",
	"100+k1-21/=322=100+k1-21/",


	"100+k1-200=300=100+k1-200",
	"100+k1-200=310=100+k1-200",
	"100+k1-200=320=100+k1-200",
	"100+k1-200=311=100+k1-200",
	"100+k1-200=321=100+k1-200",
	"100+k1-200=322=100+k1-200",

	"100+k1-k2=300=100+k1-k2",
	"100+k1-k2=310=100+k1-k2",
	"100+k1-k2=320=100+k1-k2",
	"100+k1-k2=311=100+k1-k2",
	"100+k1-k2=321=100+k1-k2",
	"100+k1-k2=322=100+k1-k2",

	"200+k1-k2=300=200+k1-k2",
	"200+k1-k2=310=200+k1-k2",
	"200+k1-k2=320=200+k1-k2",
	"200+k1-k2=311=200+k1-k2",
	"200+k1-k2=321=200+k1-k2",
	"200+k1-k2=322=200+k1-k2"};

    std::vector<std::string> vector_diagram_template_fine_tuned = {
	"200+k1-210=<00=200+k1-210"};

    std::vector<std::string> vector_diagram_template = {};
    vector_diagram_template.insert(vector_diagram_template.end(), vector_diagram_template_2p.begin(), vector_diagram_template_2p.end());
    vector_diagram_template.insert(vector_diagram_template.end(), vector_diagram_template_1s2s.begin(), vector_diagram_template_1s2s.end());
    vector_diagram_template.insert(vector_diagram_template.end(), vector_diagram_template_forgotten.begin(), vector_diagram_template_forgotten.end());
    vector_diagram_template.insert(vector_diagram_template.end(), vector_diagram_template_more.begin(), vector_diagram_template_more.end());
    vector_diagram_template.insert(vector_diagram_template.end(), vector_diagram_template_extra.begin(), vector_diagram_template_extra.end());
    vector_diagram_template.insert(vector_diagram_template.end(), vector_diagram_template_high_intermediate.begin(), vector_diagram_template_high_intermediate.end());

    //vector_diagram_template = vector_diagram_template_fine_tuned;//COMMENT THIS OUT
    //std::cout << "Total of  " << vector_diagram_template.size() << " templates" << std::endl;

    std::vector<double> vector_xi = {0.1, 0.5, 0.63, 0.8, 1, 1.25, 1.6, 2, 2.5, 3.15, 4, 5, 6.3, 8, 10, 12.5, 16, 20, 25, 31.5, 40, 50, 63, 80, 100,
    					5.6, 6.7, 7.1, 7.5, 9, 45, 53, 56, 60, 71,
    					57, 58, 59,
    					0.85, 0.9, 0.925, 0.95, 0.975,
    					1.03, 1.06, 1.09, 1.12, 1.18,
    					1.4, 1.7, 1.75, 1.8, 1.85, 1.9, 1.95, 2.06, 2.12, 2.18, 2.24, 2.8, 54.5, 56.25, 56.5, 56.75, 57.25, 57.5, 57.75, 58.5};
    //vector_xi = {103, 104.5, 104.75, 105, 105.25, 105.5, 107}; //COMMENT THIS OUT
    //vector_xi = {3.35, 3.55, 3.75, 4.12, 4.25, 4.37, 4.50, 4.62, 4.75, 4.87, 5.3, 6};
    double dm_velocity = 240.0/(3e5);

    //Thread set up
    int N_threads = 8;

    //std::cout << "Threads used: " << N_threads << std::endl;
    log_file << "Threads used: " << N_threads << std::endl;

    //std::cout << "Points to calculate: ";
    for (int i=0; i<vector_critical_xis.size(); i++){
        //std::cout << 1.0*vector_N_matrix_elements[i]*vector_N_points_per_matrix_element[i] << "@" << 1.0*vector_critical_xis[i] << "/";
    }
    //std::cout << " x " << vector_xi.size()*vector_diagram_template.size() << std::endl;

    log_file << "Points to calculate: ";
    for (int i=0; i<vector_critical_xis.size(); i++){
        log_file << 1.0*vector_N_matrix_elements[i]*vector_N_points_per_matrix_element[i] << "@" << 1.0*vector_critical_xis[i] << "/";
    }
    log_file << " x " << vector_xi.size()*vector_diagram_template.size() << std::endl;
    //std::cout << std::endl;
    log_file << std::endl;

    std::vector<uint64_t> seeds(N_threads);
    std::random_device rd;
    for (auto& s : seeds) s = rd();

    //Start of the actual calculation
    std::vector<double> partial_real_part(N_threads);
    std::vector<double> partial_imag_part(N_threads);
    std::vector<double> partial_square_real_part(N_threads);
    std::vector<double> partial_square_imag_part(N_threads);
    int current_percent = 0;
    int percent_counter = 0;

    for (int xi_counter=0; xi_counter<vector_xi.size(); xi_counter++){
        int template_counter = LSF_task - 1;

        if (template_counter < 0){
            std::cerr << "Template counter is " << template_counter << ", below 0" << std::endl;
            return 2;
        }else if (template_counter >= vector_diagram_template.size()){
            std::cerr << "Template counter is " << template_counter << ", above size " << vector_diagram_template.size() << "." << std::endl;
            return 3;
        }
        if (TYPE_OF_UNBOUND_FUNCTION == 0){
            //std::cout << "This is the code with full unbound wavefunctions." << std::endl;
            log_file << "This is the code with full unbound wavefunctions." << std::endl;
        }else if(TYPE_OF_UNBOUND_FUNCTION == 1){
            //std::cout << "This is the code with full+approx unbound wavefunctions." << std::endl;
            log_file << "This is the code with full+approx unbound wavefunctions." << std::endl;
        }else if(TYPE_OF_UNBOUND_FUNCTION == 2){
            //std::cout << "This is the code with approx (Bessel with phases) unbound wavefunctions." << std::endl;
            log_file << "This is the code with approx (Bessel with phases) unbound wavefunctions." << std::endl;
        }
        CONST_XI = vector_xi[xi_counter];
        CONST_ALPHA = vector_xi[xi_counter]*dm_velocity;
        std::string diagram_template = vector_diagram_template[template_counter];

        std::vector<std::thread> pool;
        std::fill(partial_real_part.begin(), partial_real_part.end(), 0.0);
        std::fill(partial_imag_part.begin(), partial_imag_part.end(), 0.0);
        std::fill(partial_square_real_part.begin(), partial_square_real_part.end(), 0.0);
        std::fill(partial_square_imag_part.begin(), partial_square_imag_part.end(), 0.0);
        current_percent = 0;
        percent_counter = 0;
        //bar.set_option(indicators::option::SavedStartTime{false});
        //bar.set_progress(0.0);
            
        for (int i=0; i<vector_critical_xis.size(); i++){
            if (CONST_XI < vector_critical_xis[i]){
                N_matrix_elements = vector_N_matrix_elements[i];
                N_points_per_matrix_element = vector_N_points_per_matrix_element[i];
                break;
            }
        }

        //std::cout << CONST_XI << " " << 1.0*N_matrix_elements << " " << 1.0*N_points_per_matrix_element << std::endl;
        if (easy_mode){
            N_matrix_elements = N_matrix_elements/((int) 3e2);
            N_points_per_matrix_element = N_points_per_matrix_element/((int) 1.2e1);
        }

        for (int thread_number=0; thread_number<N_threads; thread_number++){
            pool.emplace_back([&, thread_number]{
                RNG.seed(seeds[thread_number]);

                for (int i=thread_number; i<N_matrix_elements; i+=N_threads){
                    auto [first_matrix_element, second_matrix_element] = generate_diagram(diagram_template);                                        
                    auto [value1, error1] = first_matrix_element.evaluate(N_points_per_matrix_element, range);
                    auto [value2, error2] = second_matrix_element.evaluate(N_points_per_matrix_element, range);
                    
                    std::complex<double> weighted_contribution = value1 * value2 * first_matrix_element.weight;
                    partial_real_part[thread_number] += weighted_contribution.real();
                    partial_imag_part[thread_number] += weighted_contribution.imag();
                    partial_square_real_part[thread_number] += std::pow(weighted_contribution.real(), 2);
                    partial_square_imag_part[thread_number] += std::pow(weighted_contribution.imag(), 2);
                
                    if (thread_number == 0){
                        percent_counter++;
                        if (percent_counter >= N_matrix_elements/(100*N_threads)){
                            percent_counter = 0;
                            current_percent++;
                            //bar.set_progress(current_percent-1);
                        }   
                    }
                }
            });
        }
        for (auto& thread : pool) thread.join();

        double total_real_part = std::accumulate(partial_real_part.begin(), partial_real_part.end(), 0.0);
        double total_imag_part = std::accumulate(partial_imag_part.begin(), partial_imag_part.end(), 0.0);
        double total_square_real_part = std::accumulate(partial_square_real_part.begin(), partial_square_real_part.end(), 0.0);
        double total_square_imag_part = std::accumulate(partial_square_imag_part.begin(), partial_square_imag_part.end(), 0.0);

        std::complex<double> mean_value(total_real_part/N_matrix_elements, total_imag_part/N_matrix_elements);
        double real_variance = total_square_real_part/N_matrix_elements - mean_value.real()*mean_value.real();
        double imag_variance = total_square_imag_part/N_matrix_elements - mean_value.imag()*mean_value.imag();
        std::complex<double> error_value(std::sqrt(std::max(0.0, real_variance))/std::sqrt(N_matrix_elements), 
                                         std::sqrt(std::max(0.0, imag_variance))/std::sqrt(N_matrix_elements));

        //std::cout << xi_counter*vector_diagram_template.size()+template_counter+1 << "/" << vector_xi.size()*vector_diagram_template.size() << " " << std::endl;
        log_file << xi_counter*vector_diagram_template.size()+template_counter+1 << "/" << vector_xi.size()*vector_diagram_template.size() << " " << std::endl;
        //std::cout << "Total matrix element for " << diagram_template << ":" << std::endl;
        log_file << "Total matrix element for " << diagram_template << ":" << std::endl;
        //std::cout << "At alpha=" << CONST_ALPHA << " and xi=" << CONST_XI << std::endl;
        log_file << "At alpha=" << CONST_ALPHA << " and xi=" << CONST_XI << std::endl;
        //std::cout << mean_value << std::endl;
        log_file << mean_value << std::endl;
        //std::cout << error_value << std::endl;
        log_file << error_value << std::endl; 
        //std::cout << std::endl;
        log_file << std::endl;
    }

    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_sec = end - start;
    log_file << "Program completed after " <<  std::chrono::duration_cast<std::chrono::seconds>(end - start).count() << " s" << std::endl;
    log_file.close();
    return 0;
}
