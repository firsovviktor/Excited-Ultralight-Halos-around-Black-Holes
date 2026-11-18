#include <iostream>
#include <fstream>
#include <string>
#include <filesystem>
#include <vector>
#include <complex>
#include <map>
#include <set>
#include <tuple>
#include <functional>
#include <cmath>
#include <numbers>
#include <array>
#include <matplot/matplot.h>
#include <indicators/progress_bar.hpp>

indicators::ProgressBar bar{
    indicators::option::BarWidth{50},
    indicators::option::Start{"["}, indicators::option::End{"]"},
    indicators::option::Fill{"="}, indicators::option::Lead{">"},
    indicators::option::ShowPercentage{true},
    indicators::option::ShowElapsedTime{true},
    indicators::option::ShowRemainingTime{true},
};

//NFW density profile in GeV/cm^3. r is measured in pc
double NFW(double r){
    double rho_0=0.38; //NFW density in GeV/cm^3
    double r_0=18600.0; //NFW scale in pc
    return rho_0/(r/r_0 * std::pow(1+r/r_0, 2));
}

//luminous velocity profile in km/s. r is measured in pc
double V_NFW(double r){
    return 240.0;
}

// Like std::stod, but returns denormals (e.g. 1e-320) instead of throwing out_of_range, as Python's float() does
double to_double(const std::string& s) {
    return std::strtod(s.c_str(), nullptr);
}

// Replace the character following each 'k' by '0', so templates differing only in k indices compare equal
std::string strip_template(std::string t) {
    for (size_t i = 1; i < t.length(); ++i) {
        if (t[i - 1] == 'k') {
            t[i] = '0';
        }
    }
    return t;
}

// Number of non-overlapping occurrences of sub in s, like Python's str.count
int count_substr(const std::string& s, const std::string& sub) {
    int count = 0;
    for (size_t pos = s.find(sub); pos != std::string::npos; pos = s.find(sub, pos + sub.length())) {
        count++;
    }
    return count;
}

class Matrix_Element {
public:
    std::string approx;
    std::string templ;
    std::string stripped_templ;
    double alpha;
    double xi;
    std::complex<double> value;
    std::complex<double> error;

    Matrix_Element(std::string approx, std::string templ, double alpha, double xi, std::complex<double> value, std::complex<double> error)
        : approx(approx), templ(templ), stripped_templ(strip_template(templ)), alpha(alpha), xi(xi), value(value), error(error) {}

    void print() const {
        std::cout << "Approx: " << approx << "\n"
                  << "Template: " << templ << "\n"
                  << "Alpha: " << alpha << "\n"
                  << "Xi: " << xi << "\n"
                  << "Value: " << value << "\n"
                  << "Error: " << error << "\n";
    }

    bool is_the_same(const Matrix_Element& other) const {
        return (approx == other.approx) && 
               (stripped_templ == other.stripped_templ) && 
               (alpha == other.alpha) && 
               (xi == other.xi);
    }
};

const bool USE_HIGH_INTERMEDIATE = false;
const double ASSUMED_DECAY_PHASE = 1.0; // 1.0 for maximal decay, 0.0 for average and -1 for minimal

// Replace all non-overlapping occurrences of from by to, like Python's str.replace
std::string replace_all(std::string s, const std::string& from, const std::string& to) {
    for (size_t pos = s.find(from); pos != std::string::npos; pos = s.find(from, pos + to.length())) {
        s.replace(pos, from.length(), to);
    }
    return s;
}

bool contains(const std::string& s, const std::string& sub) {
    return s.find(sub) != std::string::npos;
}

// Replace the n=3 intermediate states by "nlm"
std::string replace_high_intermediate(std::string s) {
    for (const char* state : {"300", "310", "320", "311", "321", "322"}) {
        s = replace_all(s, state, "nlm");
    }
    return s;
}

class coefficients_N {
public:
    double Nkkk = 0.;
    double Nkk2p = 0.;
    double Nkk2s = 0.;
    double Nkk1s = 0.;
    double Nk2p2p = 0.;
    double Nk2p2s = 0.;
    double Nk2p1s = 0.;
    double Nk2s2s = 0.;
    double Nk2s1s = 0.;
    double N2p2p1s = 0.;
    double N2p2s1s = 0.;
    double N2s2s1s = 0.;
};

class coefficients {
public:
    double xi;
    double alpha;
    coefficients_N N1s;
    coefficients_N N2s;
    coefficients_N N2p;

    coefficients(double xi, double alpha) : xi(xi), alpha(alpha) {}

    void populate(const Matrix_Element& matrix_element) {
        if (xi != matrix_element.xi || alpha != matrix_element.alpha) return;

        const std::string& templ = matrix_element.templ;
        std::string tem = replace_high_intermediate(templ);
        int counts_1s = count_substr(tem, "100") - count_substr(templ, "=100=");
        int counts_2s = count_substr(tem, "200") - count_substr(templ, "=200=");
        int counts_2p = count_substr(tem, "21") - count_substr(templ, "=211=") - count_substr(templ, "=210=") - count_substr(templ, "=21/=");
        int counts_k = count_substr(tem, "k") - count_substr(templ, "=k1=") - count_substr(templ, "=k2=") - count_substr(templ, "=k3=");

        std::tuple<int, int, int, int> key = {counts_k, counts_2p, counts_2s, counts_1s};
        std::string key_str = "(" + std::to_string(counts_k) + ", " + std::to_string(counts_2p) + ", "
                            + std::to_string(counts_2s) + ", " + std::to_string(counts_1s) + ")";

        double n_f = counts_k / 2.0;
        int d_f = count_substr(templ, "=k1=") + count_substr(templ, "=k2=") + count_substr(templ, "=k3=");
        const double pi = std::numbers::pi;
        double prefactor = 2 * pi * std::pow(pi, -1.5 * n_f) * std::pow(8 * pi * pi * pi * xi * xi, -d_f);

        double value;
        if (key == std::tuple{2, 4, 0, 0}) {
            value = prefactor * std::abs(matrix_element.value);
        } else {
            value = prefactor * matrix_element.value.real();
        }
        [[maybe_unused]] double error = prefactor * std::abs(matrix_element.error);

        tem = replace_all(replace_all(replace_all(templ, "211", "21x"), "210", "21x"), "21/", "21x");
        tem = replace_high_intermediate(tem);
        tem = replace_all(replace_all(replace_all(tem, "k1", "k"), "k2", "k"), "k3", "k");
        tem = tem.substr(0, tem.rfind('='));

        bool high_m0 = contains(templ, "300") || contains(templ, "310") || contains(templ, "320");
        auto exceptional_template = [&]() {
            std::cout << "Exceptional template: " << key_str << ", " << templ << ", " << tem << std::endl;
        };

        if (key == std::tuple{6, 0, 0, 0}) {
            if (tem == "k+k-k=100") {
                N1s.Nkkk += 2 * value;
            } else if (tem == "k+k-k=200") {
                N2s.Nkkk += 2 * value;
            } else if (tem == "k+k-k=21x") {
                N2p.Nkkk += 2 * value;
            } else {
                exceptional_template();
            }

        } else if (key == std::tuple{4, 2, 0, 0}) {
            if (tem == "21x+k-k=100") {
                N1s.Nkk2p += 4 * value;
                N2p.Nkk2p -= 4 * value;
            } else if (tem == "21x+k-k=200") {
                N2s.Nkk2p += 4 * value;
                N2p.Nkk2p -= 4 * value;
            } else if (tem == "21x+k-k=k") {
                N2p.Nkk2p -= 4 * value;
            } else if (tem == "k+k-21x=k") {
                N2p.Nkk2p += 2 * value;
            } else if (tem == "21x+k-k=nlm") {
                if (USE_HIGH_INTERMEDIATE) {
                    if (high_m0) {
                        if (contains(templ, "211")) {
                            N2p.Nkk2p -= 2 * 4 * value;
                        } else {
                            N2p.Nkk2p -= 1 * 4 * value;
                        }
                    } else {
                        N2p.Nkk2p -= 2 * 4 * value;
                    }
                }
            } else {
                exceptional_template();
            }

        } else if (key == std::tuple{4, 0, 2, 0}) {
            if (tem == "k+k-200=k") {
                N2s.Nkk2s += 2 * value;
            } else if (tem == "200+k-k=k") {
                N2s.Nkk2s -= 4 * value;
            } else if (tem == "200+k-k=100") {
                N1s.Nkk2s += 4 * value;
                N2s.Nkk2s -= 4 * value;
            } else if (tem == "200+k-k=21x") {
                N2s.Nkk2s -= 4 * value;
                N2p.Nkk2s += 4 * value;
            } else if (tem == "200+k-k=nlm") {
                if (USE_HIGH_INTERMEDIATE) {
                    if (high_m0) {
                        N2s.Nkk2s -= 1 * 4 * value;
                    } else {
                        N2s.Nkk2s -= 2 * 4 * value;
                    }
                }
            } else {
                exceptional_template();
            }

        } else if (key == std::tuple{4, 0, 0, 2}) {
            if (tem == "k+k-100=k") {
                N1s.Nkk1s += 2 * value;
            } else if (tem == "100+k-k=k") {
                N1s.Nkk1s -= 4 * value;
            } else if (tem == "100+k-k=200") {
                N1s.Nkk1s -= 4 * value;
                N2s.Nkk1s += 4 * value;
            } else if (tem == "100+k-k=21x") {
                N1s.Nkk1s -= 4 * value;
                N2p.Nkk1s += 4 * value;
            } else if (tem == "100+k-k=nlm") {
                if (USE_HIGH_INTERMEDIATE) {
                    if (high_m0) {
                        N1s.Nkk1s -= 1 * 4 * value;
                    } else {
                        N1s.Nkk1s -= 2 * 4 * value;
                    }
                }
            } else {
                exceptional_template();
            }

        } else if (key == std::tuple{2, 4, 0, 0}) {
            if (templ == "211+211-k1=100=211+211-k1") {
                N1s.Nk2p2p += value;
                N2p.Nk2p2p -= 2 * value;
            } else if (templ == "211+210-k1=100=211+210-k1") {
                N1s.Nk2p2p += 4 * value;
                N2p.Nk2p2p -= 8 * value;
            } else if (templ == "211+21/-k1=100=211+21/-k1") {
                N1s.Nk2p2p += 4 * value;
                N2p.Nk2p2p -= 8 * value;
            } else if (templ == "210+210-k1=100=210+210-k1") {
                N1s.Nk2p2p += value;
                N2p.Nk2p2p -= 2 * value;
            } else if (templ == "210+21/-k1=100=210+21/-k1") {
                N1s.Nk2p2p += 4 * value;
                N2p.Nk2p2p -= 8 * value;
            } else if (templ == "21/+21/-k1=100=21/+21/-k1") {
                N1s.Nk2p2p += value;
                N2p.Nk2p2p -= 2 * value;
            } else if (templ == "211+21/-k1=100=210+210-k1") {
                // assuming phase of the decay
                N1s.Nk2p2p += 4 * value * ASSUMED_DECAY_PHASE;
                N2p.Nk2p2p -= 8 * value * ASSUMED_DECAY_PHASE;
            } else {
                exceptional_template();
            }

        } else if (key == std::tuple{2, 2, 2, 0}) {
            if (tem == "21x+200-k=100") {
                N1s.Nk2p2s += 4 * value;
                N2s.Nk2p2s -= 4 * value;
                N2p.Nk2p2s -= 4 * value;
            } else if (tem == "21x+k-200=k") {
                N2s.Nk2p2s += 4 * value;
                N2p.Nk2p2s -= 4 * value;
            } else if (tem == "200+k-21x=k") {
                N2s.Nk2p2s -= 4 * value;
                N2p.Nk2p2s += 4 * value;
            } else {
                exceptional_template();
            }

        } else if (key == std::tuple{2, 2, 0, 2}) {
            if (tem == "100+k-21x=200") {
                N1s.Nk2p1s -= 4 * value;
                N2s.Nk2p1s += 4 * value;
                N2p.Nk2p1s += 4 * value;
            } else if (tem == "100+k-21x=k") {
                N1s.Nk2p1s -= 4 * value;
                N2p.Nk2p1s += 4 * value;
            } else if (tem == "21x+k-100=k") {
                N1s.Nk2p1s += 4 * value;
                N2p.Nk2p1s -= 4 * value;
            } else if (templ == "100+k1-211=211=100+k1-211" ||
                       templ == "100+k1-211=210=100+k1-211" ||
                       templ == "100+k1-211=21/=100+k1-211" ||
                       templ == "100+k1-210=211=100+k1-210" ||
                       templ == "100+k1-210=210=100+k1-210" ||
                       templ == "100+k1-210=21/=100+k1-210" ||
                       templ == "100+k1-21/=211=100+k1-21/" ||
                       templ == "100+k1-21/=210=100+k1-21/" ||
                       templ == "100+k1-21/=21/=100+k1-21/") {
                N1s.Nk2p1s -= 4 * value;
                N2p.Nk2p1s += 8 * value;
            } else if (tem == "100+k-21x=nlm") {
                if (USE_HIGH_INTERMEDIATE) {
                    if (high_m0) {
                        if (contains(templ, "211")) {
                            N1s.Nk2p1s -= 2 * 4 * value;
                            N2p.Nk2p1s += 2 * 4 * value;
                        } else {
                            N1s.Nk2p1s -= 1 * 4 * value;
                            N2p.Nk2p1s += 1 * 4 * value;
                        }
                    } else {
                        N1s.Nk2p1s -= 2 * 4 * value;
                        N2p.Nk2p1s += 2 * 4 * value;
                    }
                }
            } else {
                exceptional_template();
            }

        } else if (key == std::tuple{2, 0, 4, 0}) {
            if (tem == "200+200-k=100") {
                N1s.Nk2s2s += value;
                N2s.Nk2s2s -= 2 * value;
            } else {
                exceptional_template();
            }

        } else if (key == std::tuple{2, 0, 2, 2}) {
            if (tem == "100+k-200=200") {
                N1s.Nk2s1s -= 4 * value;
                N2s.Nk2s1s += 8 * value;
            } else if (tem == "100+k-200=k") {
                N1s.Nk2s1s -= 4 * value;
                N2s.Nk2s1s += 4 * value;
            } else if (tem == "200+k-100=k") {
                N1s.Nk2s1s += 4 * value;
                N2s.Nk2s1s -= 4 * value;
            } else if (tem == "100+k-200=21x") {
                N1s.Nk2s1s -= 4 * value;
                N2s.Nk2s1s += 4 * value;
                N2p.Nk2s1s += 4 * value;
            } else if (tem == "100+k-200=nlm") {
                if (USE_HIGH_INTERMEDIATE) {
                    if (high_m0) {
                        N1s.Nk2s1s -= 1 * 4 * value;
                        N2s.Nk2s1s += 1 * 4 * value;
                    } else {
                        N1s.Nk2s1s -= 2 * 4 * value;
                        N2s.Nk2s1s += 2 * 4 * value;
                    }
                }
            } else {
                exceptional_template();
            }

        } else if (key == std::tuple{0, 4, 0, 2}) {
            if (templ == "211+211-100=k1=211+211-100") {
                N1s.N2p2p1s += value;
                N2p.N2p2p1s -= 2 * value;
            } else if (templ == "211+210-100=k1=211+210-100") {
                N1s.N2p2p1s += 4 * value;
                N2p.N2p2p1s -= 8 * value;
            } else if (templ == "211+21/-100=k1=211+21/-100") {
                N1s.N2p2p1s += 4 * value;
                N2p.N2p2p1s -= 8 * value;
            } else if (templ == "210+210-100=k1=210+210-100") {
                N1s.N2p2p1s += value;
                N2p.N2p2p1s -= 2 * value;
            } else if (templ == "210+21/-100=k1=210+21/-100") {
                N1s.N2p2p1s += 4 * value;
                N2p.N2p2p1s -= 8 * value;
            } else if (templ == "21/+21/-100=k1=21/+21/-100") {
                N1s.N2p2p1s += value;
                N2p.N2p2p1s -= 2 * value;
            } else if (templ == "211+21/-100=k1=210+210-100") {
                // assuming phase of the decay
                N1s.N2p2p1s += 4 * value * ASSUMED_DECAY_PHASE;
                N2p.N2p2p1s -= 8 * value * ASSUMED_DECAY_PHASE;
            } else {
                exceptional_template();
            }

        } else if (key == std::tuple{0, 2, 2, 2}) {
            if (tem == "21x+200-100=k") {
                N1s.N2p2s1s += 4 * value;
                N2s.N2p2s1s -= 4 * value;
                N2p.N2p2s1s -= 4 * value;
            } else {
                exceptional_template();
            }

        } else if (key == std::tuple{0, 0, 4, 2}) {
            if (tem == "200+200-100=k") {
                N1s.N2s2s1s += value;
                N2s.N2s2s1s -= 2 * value;
            } else {
                exceptional_template();
            }

        } else {
            std::cout << "Exceptional key: " << key_str << ", " << templ << ", " << tem << std::endl;
        }
    }
};

const bool DIAGNOSTIC_MODE = false;
std::vector<std::map<std::string, double>> N_1s_ins_and_outs;
std::vector<std::map<std::string, double>> N_2s_ins_and_outs;
std::vector<std::map<std::string, double>> N_2p_ins_and_outs;

const std::array<std::string, 12> COLLISION_TERM_NAMES = {
    "Nkkk", "Nkk2p", "Nkk2s", "Nkk1s",
    "Nk2p2p", "Nk2p2s", "Nk2p1s", "Nk2s2s",
    "Nk2s1s", "N2p2p1s", "N2p2s1s", "N2s2s1s"
};

// The 12 three-body terms of dN/dt for a state with coefficients C, in the order of COLLISION_TERM_NAMES
std::array<double, 12> collision_terms(const coefficients_N& C, double N_1s, double N_2s, double N_2p, double N_k) {
    return {
        N_k * N_k * N_k * C.Nkkk,
        N_k * N_k * N_2p * C.Nkk2p,
        N_k * N_k * N_2s * C.Nkk2s,
        N_k * N_k * N_1s * C.Nkk1s,

        N_k * N_2p * N_2p * C.Nk2p2p,
        N_k * N_2p * N_2s * C.Nk2p2s,
        N_k * N_2p * N_1s * C.Nk2p1s,
        N_k * N_2s * N_2s * C.Nk2s2s,

        N_k * N_2s * N_1s * C.Nk2s1s,
        N_2p * N_2p * N_1s * C.N2p2p1s,
        N_2p * N_2s * N_1s * C.N2p2s1s,
        N_2s * N_2s * N_1s * C.N2s2s1s
    };
}

// Record the individual terms of one dN/dt evaluation
void record_ins_and_outs(std::vector<std::map<std::string, double>>& ins_and_outs, const std::array<double, 12>& terms,
                         const std::string& BH_name, double BH_term, double total) {
    ins_and_outs.emplace_back();
    for (size_t i = 0; i < terms.size(); ++i) {
        ins_and_outs.back()[COLLISION_TERM_NAMES[i]] = terms[i];
    }
    ins_and_outs.back()[BH_name] = BH_term;
    ins_and_outs.back()["total"] = total;
}

double dot_N_1s(double N_1s, double N_2s, double N_2p, double N_k, bool DIAGNOSE_SWITCH,
                const coefficients& coeff, double Gamma_BH_1s) {
    std::array<double, 12> terms = collision_terms(coeff.N1s, N_1s, N_2s, N_2p, N_k);
    double dot = 0.0;
    for (double term : terms) dot += term;

    dot -= N_1s * Gamma_BH_1s;

    if (DIAGNOSTIC_MODE && DIAGNOSE_SWITCH) {
        record_ins_and_outs(N_1s_ins_and_outs, terms, "1sBH", -N_1s * Gamma_BH_1s, dot);
    }

    return dot;
}

double dot_N_2s(double N_1s, double N_2s, double N_2p, double N_k, bool DIAGNOSE_SWITCH,
                const coefficients& coeff, double Gamma_BH_2s) {
    std::array<double, 12> terms = collision_terms(coeff.N2s, N_1s, N_2s, N_2p, N_k);
    double dot = 0.0;
    for (double term : terms) dot += term;

    dot -= N_2s * Gamma_BH_2s;

    if (DIAGNOSTIC_MODE && DIAGNOSE_SWITCH) {
        record_ins_and_outs(N_2s_ins_and_outs, terms, "2sBH", -N_2s * Gamma_BH_2s, dot);
    }

    return dot;
}

double dot_N_2p(double N_1s, double N_2s, double N_2p, double N_k, bool DIAGNOSE_SWITCH,
                const coefficients& coeff, double Gamma_BH_2p) {
    std::array<double, 12> terms = collision_terms(coeff.N2p, N_1s, N_2s, N_2p, N_k);
    double dot = 0.0;
    for (double term : terms) dot += term;

    dot -= 3 * N_2p * Gamma_BH_2p; //factor of 3 as decay acts on all multiplets

    if (DIAGNOSTIC_MODE && DIAGNOSE_SWITCH) {
        record_ins_and_outs(N_2p_ins_and_outs, terms, "2pBH", -3 * N_2p * Gamma_BH_2p, dot / 3);
    }

    return dot / 3; //N_2p is population of 1 state, so it changes slower
}

bool in_tolerance(double delta_N_1s, double delta_N_2s, double delta_N_2p, double N_1s, double N_2s, double N_2p,
                  double tol_rel, double tol_abs) {
    if (std::abs(delta_N_1s) > tol_rel * std::abs(N_1s) + tol_abs) return false;
    if (std::abs(delta_N_2s) > tol_rel * std::abs(N_2s) + tol_abs) return false;
    if (std::abs(delta_N_2p) > tol_rel * std::abs(N_2p) + tol_abs) return false;
    return true;
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
    std::string base_directory_path = "./logs";

    int file_counter = 0;
    std::vector<Matrix_Element> ALL_MATRIX_ELEMENTS;

    // Reading the files
    for (const auto& inner_directory : std::filesystem::directory_iterator(base_directory_path)) {
        if (inner_directory.path().string().starts_with("./logs/logs")){
            for (const auto& file : std::filesystem::directory_iterator(inner_directory)) {
                //std::cout << file.path() << std::endl;
                file_counter ++;
                std::ifstream read_file(file.path());
                if (!read_file.is_open()) {
                    std::cerr << "Error: Could not open file " << file.path() << std::endl;
                }

                // Parse "(re,im)" into a complex number
                auto parse_complex = [](const std::string& s) {
                    std::string inner = s.substr(1, s.length() - 2);
                    size_t comma = inner.find(',');
                    return std::complex<double>(to_double(inner.substr(0, comma)), to_double(inner.substr(comma + 1)));
                };

                std::string line;
                int trigger_counter = -1;
                std::string approx, templ;
                double alpha = 0, xi = 0;
                std::complex<double> value;
                while (std::getline(read_file, line)) {
                    if (line.find("This is the code") != std::string::npos) {
                        trigger_counter = 0;
                        approx = "";
                        if (line.find("full unbound") != std::string::npos) {
                            approx = "full";
                        } else if (line.find("full+approx") != std::string::npos) {
                            approx = "full+approx";
                        } else if (line.find("approx") != std::string::npos) {
                            approx = "approx";
                        }
                        continue;
                    }
                    if (trigger_counter < 0) continue;
                    trigger_counter++;

                    if (trigger_counter == 2) {
                        // "Total matrix element for <template>:" -> last token without the trailing ':'
                        templ = line.substr(line.rfind(' ') + 1);
                        templ.pop_back();
                    } else if (trigger_counter == 3) {
                        // "At alpha=<alpha> and xi=<xi>"
                        size_t index = line.find('=');
                        alpha = to_double(line.substr(index + 1, line.find(' ', index) - index - 1));
                        index = line.find('=', index + 1);
                        xi = to_double(line.substr(index + 1, line.find(' ', index) - index - 1));
                    } else if (trigger_counter == 4) {
                        value = parse_complex(line);
                    } else if (trigger_counter == 5) {
                        ALL_MATRIX_ELEMENTS.emplace_back(approx, templ, alpha, xi, value, parse_complex(line));
                        trigger_counter = -1;
                    }
                }
            }
        }
    }

    std::cout << "Total files " << file_counter << std::endl;
    std::cout << "Total matrix elements " << ALL_MATRIX_ELEMENTS.size() << std::endl;

    // Remove duplicate entries, keeping the one with the smallest error
    std::vector<bool> to_be_deleted(ALL_MATRIX_ELEMENTS.size(), false);
    for (size_t i = 0; i < ALL_MATRIX_ELEMENTS.size(); ++i) {
        for (size_t j = 0; j < i; ++j) {
            const Matrix_Element& me_i = ALL_MATRIX_ELEMENTS[i];
            const Matrix_Element& me_j = ALL_MATRIX_ELEMENTS[j];
            if (me_i.is_the_same(me_j)) {
                std::complex<double> dif = me_i.value - me_j.value;
                std::complex<double> err = me_i.error + me_j.error;
                if ((std::abs(dif) > 2 * std::abs(err)) && (std::abs(dif) > 1e-60)) {
                    std::cout << "Possible inconsistency:" << std::endl;
                    me_i.print();
                    me_j.print();
                    std::cout << std::endl;
                }

                if (std::abs(me_i.error) < std::abs(me_j.error)) {
                    to_be_deleted[j] = true;
                } else {
                    to_be_deleted[i] = true;
                }
            }
        }
    }

    std::vector<Matrix_Element> unique_elements;
    for (size_t i = 0; i < ALL_MATRIX_ELEMENTS.size(); ++i) {
        if (!to_be_deleted[i]) unique_elements.push_back(ALL_MATRIX_ELEMENTS[i]);
    }
    ALL_MATRIX_ELEMENTS = std::move(unique_elements);

    std::cout << "Unique matrix elements " << ALL_MATRIX_ELEMENTS.size() << std::endl;

    // Count matrix elements by (k, 2p, 2s, 1s) content of the template, per (xi, alpha) point
    std::map<std::tuple<int, int, int, int>, int, std::greater<>> ME_counts;
    std::set<std::pair<double, double>> xi_counts;
    int count_high_intermediate = 0;
    for (const Matrix_Element& matrix_element : ALL_MATRIX_ELEMENTS) {
        const std::string& t = matrix_element.templ;
        if (t.find("32") != std::string::npos || t.find("31") != std::string::npos || t.find("30") != std::string::npos) {
            count_high_intermediate++;
        } else {
            int counts_1s = count_substr(t, "100") - count_substr(t, "=100=");
            int counts_2s = count_substr(t, "200") - count_substr(t, "=200=");
            int counts_2p = count_substr(t, "21") - count_substr(t, "=211=") - count_substr(t, "=210=") - count_substr(t, "=21/=");
            int counts_k = count_substr(t, "k") - count_substr(t, "=k1=") - count_substr(t, "=k2=") - count_substr(t, "=k3=");

            ME_counts[{counts_k, counts_2p, counts_2s, counts_1s}]++;
            xi_counts.insert({matrix_element.xi, matrix_element.alpha});
        }
    }

    double total = 0;
    for (const auto& [key, count] : ME_counts) {
        auto [k, p2, s2, s1] = key;
        double per_point = static_cast<double>(count) / xi_counts.size();
        std::cout << "(" << k << ", " << p2 << ", " << s2 << ", " << s1 << "): " << per_point << std::endl;
        total += per_point;
    }
    std::cout << "Total: " << total << " out of 81" << std::endl;
    std::cout << "High intermediate: " << static_cast<double>(count_high_intermediate) / xi_counts.size() << " out of 48" << std::endl;

    std::vector<coefficients> ALL_COEFFICIENTS;
    for (const auto& [xi, alpha] : xi_counts) {
        coefficients coeff(xi, alpha);

        int n_current = 0;
        for (const Matrix_Element& matrix_element : ALL_MATRIX_ELEMENTS) {
            if (matrix_element.xi == xi && matrix_element.alpha == alpha) {
                coeff.populate(matrix_element);
                n_current++;
            }
        }

        ALL_COEFFICIENTS.push_back(coeff);
        if (n_current == 81 + 48) {
            std::cout << xi << " " << n_current << "✓ ";
        } else {
            std::cout << std::endl << xi << " " << n_current << "✗ " << std::endl;
        }
    }
    std::cout << std::endl;
    std::cout << std::endl;


    std::vector<double> vector_xi = {0.5, 0.63, 0.8, 
                                    1, 1.25, 1.6, 2, 2.5, 3.15, 4, 5, 6.3, 8, 
                                    10, 12.5, 16, 20, 25, 31.5, 40, 50, 63, 80, 
                                    100};

    std::vector<double> vector_f = {1e4, 1.25e4, 1.6e4, 2e4, 2.5e4, 3.15e4, 4e4, 5e4, 6.3e4, 8e4, 
                                    1e5, 1.25e5, 1.6e5, 2e5, 2.5e5, 3.15e5, 4e5, 5e5, 6.3e5, 8e5, 
                                    1e6, 1.25e6, 1.6e6, 2e6, 2.5e6, 3.15e6, 4e6, 5e6, 6.3e6, 8e6,
                                    1e7, 1.25e7, 1.6e7, 2e7, 2.5e7, 3.15e7, 4e7, 5e7, 6.3e7, 8e7, 
                                    1e8, 1.25e8, 1.6e8, 2e8, 2.5e8, 3.15e8, 4e8, 5e8, 6.3e8, 8e8, 
                                    1e9, 1.25e9, 1.6e9, 2e9, 2.5e9, 3.15e9, 4e9, 5e9, 6.3e9, 8e9,
                                    1e10, 1.25e10, 1.6e10, 2e10, 2.5e10, 3.15e10, 4e10, 5e10, 6.3e10, 8e10,
                                    1e11, 1.25e11, 1.6e11, 2e11, 2.5e11, 3.15e11, 4e11, 5e11, 6.3e11, 8e11,
                                    1e12, 1.25e12, 1.6e12, 2e12, 2.5e12, 3.15e12, 4e12, 5e12, 6.3e12, 8e12,
                                    1e13};

    //vector_xi = {12.5, 16, 20, 25, 31.5, 40, 50, 63, 80, 100};
    /*vector_f = {1e15, 1.25e15, 1.6e15, 2e15, 2.5e15, 3.15e15, 4e15, 5e15, 6.3e15, 8e15,
                1e16, 1.25e16, 1.6e16, 2e16, 2.5e16, 3.15e16, 4e16, 5e16, 6.3e16, 8e16,
                1e17, 1.25e17, 1.6e17, 2e17, 2.5e17, 3.15e17, 4e17, 5e17, 6.3e17, 8e17,
                1e18, 1.25e18, 1.6e18, 2e18, 2.5e18, 3.15e18, 4e18, 5e18, 6.3e18, 8e18,
                1e19};*/

    //vector_xi = {80};
    //vector_f = {1e5};
    bool superradiance_test = false;
    std::string astro_regime = "BH";
    //Sun has M=1, r=8e3, Gamma = 0
    //BH has M=10, r=1e3, Gamma=Gamma_BH
    //IBH has M=1e3, r=1e3, Gamma=Gamma_BH
    //SBH has M=1e5, r=1e3, Gamma=Gamma_BH
    //star has M=10, r=1e3, Gamma=0

    std::string plotting_regime = "true";
    //false has no plot
    //true shows a plot
    //log shows a plot and logs sampled data

    std::string filename;
    if (astro_regime == "Sun"){
        filename = generate_timestamped_filename("Sun_solve_ODE_logs", ".txt"); 
    }else if (astro_regime == "BH"){
        filename = generate_timestamped_filename("BH_solve_ODE_logs", ".txt"); 
    }else if (astro_regime == "IBH"){
        filename = generate_timestamped_filename("IBH_solve_ODE_logs", ".txt");
    }else if (astro_regime == "SBH"){
        filename = generate_timestamped_filename("SBH_solve_ODE_logs", ".txt");
    }else if (astro_regime == "star"){
        filename = generate_timestamped_filename("star_solve_ODE_logs", ".txt");
    }else{
        std::cerr << "Unknown astro_regime" << std::endl;
        return 1;
    }
    std::ofstream log_file(filename, std::ios::app);

    //bar.set_option(indicators::option::SavedStartTime{false});
    //bar.set_progress(0.0);
    for (int xi_index=0; xi_index<vector_xi.size(); xi_index++){
        for (int f_index=0; f_index<vector_f.size(); f_index++){
            bar.set_progress(100*(xi_index*vector_f.size()+f_index)/(vector_xi.size()*vector_f.size()));
            double xi = vector_xi[xi_index]; //focusing parameter
            double f = vector_f[f_index]; //Axion decay constant in GeV
            int number_of_points_in_ODE = int(1e8); //Number of points in the ODE integration
            //1.25, 1e8 4e4 or 100, 5e4, 1.5e6

            double M, r;
            if (astro_regime == "Sun"){
                M = 1;
                r = 8e3;
            }else if (astro_regime == "BH"){
                M = 10;
                r = 1e3;
            }else if (astro_regime == "IBH"){
                M = 1e3;
                r = 1e3;
            }else if (astro_regime == "SBH"){
                M = 1e5;
                r = 1e3;
            }else if (astro_regime == "star"){
                M = 10;
                r = 1e3;
            }else{
                std::cerr << "Unknown astro_regime" << std::endl;
                return 1;
            }

            double alpha = xi*V_NFW(r)/3e5; //Gravitational fine structure constant

            if (xi_counts.contains({xi, alpha})) {
                int n_elements = 0;
                for (const Matrix_Element& me : ALL_MATRIX_ELEMENTS) {
                    if (me.xi == xi && me.alpha == alpha) n_elements++;
                }
                log_file << "Found a match for (" << xi << ", " << alpha << ") with " << n_elements << " elements" << std::endl;
            } else {
                std::cout << "No match found for (xi, alpha) = (" << xi << ", " << alpha << ")" << std::endl;
                std::cout << "Available (xi, alpha):" << std::endl;
                for (const auto& [key_xi, key_alpha] : xi_counts) {
                    std::cout << "(" << key_xi << ", " << key_alpha << ")" << std::endl;
                }
                return 1;
            }

            coefficients coeff(xi, alpha);
            for (const coefficients& coefficient : ALL_COEFFICIENTS) {
                if (coefficient.xi == xi && coefficient.alpha == alpha) {
                    coeff = coefficient;
                    break;
                }
            }

            double alpha_0 = 7.5e-5;
            double v_tilde = V_NFW(r)/240.0;
            double rho_tilde = NFW(r)/0.4;
            double alpha_tilde = alpha/7.5e-5;
            double f_tilde = f/1e8;
            double M_tilde = M;
            double m_tilde = alpha/alpha_0/M; //particle mass in 10**-14 eV

            double Gamma_BH = 0;
            double Gamma_BH_1s = 0; 
            double Gamma_BH_2s = 0;
            double Gamma_BH_2p = 0;

            if ((astro_regime == "BH")||(astro_regime == "IBH")||(astro_regime == "SBH")){
                Gamma_BH = 2*9.32268 * std::pow(rho_tilde, -2) * std::pow(M_tilde, -4) * std::pow(v_tilde, 4) * std::pow(f_tilde, 4) * std::pow(alpha_tilde, 7); //2*\Gamma'_1s in the overleaf
                Gamma_BH_1s = Gamma_BH;
                Gamma_BH_2s = Gamma_BH/8;
                Gamma_BH_2p = 2*3.04507e-18 * std::pow(rho_tilde, -2) * std::pow(M_tilde, -4) * std::pow(v_tilde, 4) * std::pow(f_tilde, 4) * std::pow(alpha_tilde, 11); //2*\Gamma'_1s in the overleaf
                if (superradiance_test){
                    Gamma_BH_2p = 0.0;
                }
            }

            double time_conversion = 1.03658e12 * std::pow(rho_tilde, -2) * std::pow(M_tilde, -3) * std::pow(v_tilde, 4) * std::pow(f_tilde, 4) * alpha_tilde; //in years, t/\tau in the overleaf
            double actual_N_k = 6.85591e61 * rho_tilde * std::pow(m_tilde, -6) * std::pow(M_tilde, -2) * std::pow(v_tilde, -1);
            double max_N_1s = 3e67 * std::pow(f_tilde, 2) * std::pow(m_tilde, -3) * std::pow(M_tilde, -1); //Gilad's value, maybe should have x2
            double max_N_2s = 1.6e69 * std::pow(f_tilde, 2) * std::pow(m_tilde, -3) * std::pow(M_tilde, -1); //never reached, doesn't matter
            double max_N_2p = 1.6e69 * std::pow(f_tilde, 2) * std::pow(m_tilde, -3) * std::pow(M_tilde, -1);

            double tau_rel = 9.05699e9 * std::pow(f_tilde, 4) * std::pow(m_tilde, 3) * std::pow(rho_tilde, -2) * std::pow(v_tilde, 2);
            double tau_1s = 1.11189e11 * std::pow(m_tilde, -6) * std::pow(M_tilde, -5);
            double tau_coherent = 0.0065 * std::pow(m_tilde, -1) * std::pow(v_tilde, -2);
            double tau_Compton = 2.09e-9 * std::pow(m_tilde, -1);
            double tau_mixing = 21.0 * std::pow(m_tilde, -3) * std::pow(M_tilde, -2);

            double rho_crit_to_BG_1s = 2.9e4 * std::pow(f_tilde, 2) * std::pow(m_tilde, 4) * std::pow(M_tilde, 2) * std::pow(rho_tilde, -1);

            double N_k = 1.0;

            std::vector<double> time = {0};
            double dt = 1e-1;

            int tolerance_streak = 0;
            double tol_rel = 1e-3;
            double tol_abs = 1e-6;

            std::vector<double> N_1s = {1e-4};
            std::vector<double> N_2s = {1e-4};
            std::vector<double> N_2p = {1e-4};

            if (superradiance_test){
                N_1s = {1e-50};
                N_2s = {1e-50};
                N_2p = {M_tilde/m_tilde * 1.1e78/actual_N_k};
            }

            int N = number_of_points_in_ODE;
            int saturation_check = -1;

            for (int i=0; i<N; i++) {
                bool DIAGNOSE_SWITCH = true;//(i % 100 == 0);
                double delta_N_1s = dt * dot_N_1s(N_1s.back(), N_2s.back(), N_2p.back(), N_k, DIAGNOSE_SWITCH, coeff, Gamma_BH_1s);
                double delta_N_2s = dt * dot_N_2s(N_1s.back(), N_2s.back(), N_2p.back(), N_k, DIAGNOSE_SWITCH, coeff, Gamma_BH_2s);
                double delta_N_2p = dt * dot_N_2p(N_1s.back(), N_2s.back(), N_2p.back(), N_k, DIAGNOSE_SWITCH, coeff, Gamma_BH_2p);

                if (in_tolerance(delta_N_1s, delta_N_2s, delta_N_2p, N_1s.back(), N_2s.back(), N_2p.back(), tol_rel, tol_abs)) {
                    tolerance_streak++;
                    time.push_back(time.back() + dt);
                    N_1s.push_back(std::min(std::max(1 / actual_N_k, N_1s.back() + delta_N_1s), max_N_1s / actual_N_k));
                    N_2s.push_back(std::min(std::max(1 / actual_N_k, N_2s.back() + delta_N_2s), max_N_2s / actual_N_k));
                    N_2p.push_back(std::min(std::max(1 / actual_N_k, N_2p.back() + delta_N_2p), max_N_2p / actual_N_k));

                    if (time.back() * time_conversion > 1.5e10) break;//CHANGE BACK TO 1.5e10!!!
                } else {
                    if (DIAGNOSTIC_MODE && DIAGNOSE_SWITCH) {
                        N_1s_ins_and_outs.pop_back();
                        N_2s_ins_and_outs.pop_back();
                        N_2p_ins_and_outs.pop_back();
                    }
                    tolerance_streak = 0;
                    dt = dt / 1.1;
                }

                if (tolerance_streak > 100) {
                    dt = dt * 1.1;
                    tolerance_streak = 0;
                }

                if (saturation_check > 0){
                    saturation_check++;
                }else{
                    if (i%int(1e5)==0){
                        saturation_check = 1;
                        int last_index = time.size()-1;
                        for (int j=-4999; j<0; j++){
                            if (std::abs(N_1s[last_index+j]-N_1s[last_index]) > 1e-6){
                                saturation_check = -1; 
                            }
                            if (std::abs(N_2s[last_index+j]-N_2s[last_index]) > 1e-6){
                                saturation_check = -1; 
                            }
                            if (std::abs(N_2p[last_index+j]-N_2p[last_index]) > 1e-6){
                                saturation_check = -1; 
                            }
                        }
                        /*if (saturation_check == 1){
                            for (int j=0; j>-1e4; j--){
                                if (N_2p[i+j] != 0){
                                    std::cout << i+j << " " << N_1s[i+j] << " " << N_2s[i+j] << " " << N_2p[i+j] << std::endl;
                                    std::cout << i+j+1 << " " << N_1s[i+j+1] << " " << N_2s[i+j+1] << " " << N_2p[i+j+1] << std::endl;
                                    break;
                                }
                            }
                        }*/
                    }
                }
                if (saturation_check > 1e5){
                    if ((plotting_regime == "true")||(plotting_regime == "log")) std::cout << "stopping due to saturation" << std::endl; 
                    break;
                }
            }
            std::vector<double> converted_time = {};
            for (int i=0; i<time.size(); i++){
                converted_time.push_back(time[i]*time_conversion);
            }
            log_file << std::scientific << std::setprecision(3);
            std::cout << std::scientific << std::setprecision(3);
            if (DIAGNOSTIC_MODE){
                std::cout << std::showpos;
            }
            //std::cout << "ODE: " << time.size() - 1 << " accepted steps, t = " << time.back() * time_conversion << " years" << std::endl;

            log_file << "xi=" << xi << ", f=" << f << "GeV, T=" << converted_time.back() << "years" << std::endl;
            log_file << "N1s=" << N_1s.back() << ", N2s=" << N_2s.back() << ", N2p=" << N_2p.back() << std::endl;
            
            
            if ((plotting_regime == "true")||(plotting_regime == "log")){
                int index = 0;
                std::cout << index << " " << converted_time[index] << " " << N_1s[index] << " " << N_2s[index] << " " << N_2p[index] << std::endl;

                if (DIAGNOSTIC_MODE){
                    std::cout << "N_1s: ";
                    for (const auto& [key, value] : N_1s_ins_and_outs[index]){
                        std::cout << key << ": " << value << ",   ";
                    }
                    std::cout << std::endl;
                    std::cout << "N_2s: ";
                    for (const auto& [key, value] : N_2s_ins_and_outs[index]){
                        std::cout << key << ": " << value << ",   ";
                    }
                    std::cout << std::endl;
                    std::cout << "N_2p: ";
                    for (const auto& [key, value] : N_2p_ins_and_outs[index]){
                        std::cout << key << ": " << value << ",   ";
                    }
                    std::cout << std::endl;
                    std::cout << std::endl;
                }

                for (int i=0; i<100; i++){
                    if ( int(std::exp(std::log(time.size())/100*i)) > index ){
                        index = int(std::exp(std::log(time.size())/100*i));
                        std::cout << index << " " << converted_time[index] << " " << N_1s[index] << " " << N_2s[index] << " " << N_2p[index] << std::endl;

                        if (DIAGNOSTIC_MODE){
                            std::cout << "N_1s: ";
                            for (const auto& [key, value] : N_1s_ins_and_outs[index]){
                                std::cout << key << ": " << value << ",   ";
                            }
                            std::cout << std::endl;
                            std::cout << "N_2s: ";
                            for (const auto& [key, value] : N_2s_ins_and_outs[index]){
                                std::cout << key << ": " << value << ",   ";
                            }
                            std::cout << std::endl;
                            std::cout << "N_2p: ";
                            for (const auto& [key, value] : N_2p_ins_and_outs[index]){
                                std::cout << key << ": " << value << ",   ";
                            }
                            std::cout << std::endl;
                            std::cout << std::endl;
                        }
                    }
                }
                index = time.size()-1;
                std::cout << index << " " << converted_time[index] << " " << N_1s[index] << " " << N_2s[index] << " " << N_2p[index] << std::endl;

                if (DIAGNOSTIC_MODE){
                    std::cout << "N_1s: ";
                    for (const auto& [key, value] : N_1s_ins_and_outs[index]){
                        std::cout << key << ": " << value << ",   ";
                    }
                    std::cout << std::endl;
                    std::cout << "N_2s: ";
                    for (const auto& [key, value] : N_2s_ins_and_outs[index]){
                        std::cout << key << ": " << value << ",   ";
                    }
                    std::cout << std::endl;
                    std::cout << "N_2p: ";
                    for (const auto& [key, value] : N_2p_ins_and_outs[index]){
                        std::cout << key << ": " << value << ",   ";
                    }
                    std::cout << std::endl;
                    std::cout << std::endl;
                }

                std::cout << std::endl;

                std::vector<double> sampled_converted_time = {};
                std::vector<double> sampled_N_1s = {};
                std::vector<double> sampled_N_2s = {};
                std::vector<double> sampled_N_2p = {};
                int NUMBER_OF_POINTS_IN_PLOT = int(1e4);
                for (int i=0; i<NUMBER_OF_POINTS_IN_PLOT+1; i++){
                    index = (time.size()-1)*i/NUMBER_OF_POINTS_IN_PLOT;
                    sampled_converted_time.push_back(converted_time[index]);
                    sampled_N_1s.push_back(N_1s[index]);
                    sampled_N_2s.push_back(N_2s[index]);
                    sampled_N_2p.push_back(N_2p[index]);
                }

                if (plotting_regime == "log"){
                    std::string filename_plot;
                    if (astro_regime == "Sun"){
                        filename_plot = generate_timestamped_filename("Sun_solve_ODE_plot", ".txt"); 
                    }else if (astro_regime == "BH"){
                        filename_plot = generate_timestamped_filename("BH_solve_ODE_plot", ".txt"); 
                    }else if (astro_regime == "IBH"){
                        filename_plot = generate_timestamped_filename("IBH_solve_ODE_plot", ".txt");
                    }else if (astro_regime == "SBH"){
                        filename_plot = generate_timestamped_filename("SBH_solve_ODE_plot", ".txt"); 
                    }else if (astro_regime == "star"){
                        filename_plot = generate_timestamped_filename("star_solve_ODE_plot", ".txt");
                    }else{
                        std::cerr << "Unknown astro_regime" << std::endl;
                        return 1;
                    }

                    std::ofstream log_file_plot(filename_plot, std::ios::app);
                    log_file_plot << std::setprecision(16) << "xi= " << xi << " f= " << f << std::endl;  
                    for (int i=0; i<sampled_converted_time.size(); i++){
                        log_file_plot << sampled_converted_time[i] << " ";
                    }
                    log_file_plot << std::endl;
                    for (int i=0; i<sampled_converted_time.size(); i++){
                        log_file_plot << sampled_N_1s[i] << " ";
                    }
                    log_file_plot << std::endl;
                    for (int i=0; i<sampled_converted_time.size(); i++){
                        log_file_plot << sampled_N_2s[i] << " ";
                    }
                    log_file_plot << std::endl;
                    for (int i=0; i<sampled_converted_time.size(); i++){
                        log_file_plot << sampled_N_2p[i] << " ";
                    }
                }
                
                auto fig = matplot::gcf();
                fig->backend()->run_command("unset warnings");
                auto lines = matplot::plot(sampled_converted_time, sampled_N_1s,
                                           sampled_converted_time, sampled_N_2s,
                                           sampled_converted_time, sampled_N_2p);
                
                lines[0]->display_name("1s");
                lines[1]->display_name("2s");
                lines[2]->display_name("2p");

                lines[0]->color("blue");
                lines[1]->color("red");
                lines[2]->color("green");

                for (auto& l : lines){
                    l->line_width(2.5);
                }

                matplot::xlim({0, converted_time.back()});
                matplot::ylim({0, 1.05*std::max(*std::max_element(N_2p.begin(), N_2p.end()), *std::max_element(N_1s.begin(), N_1s.end()))});
                matplot::legend();
                matplot::show();
            }
        }
    }

    bar.set_progress(100);
    return 0;
}