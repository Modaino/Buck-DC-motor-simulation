#include <iostream>
#include <vector>
#include <boost/numeric/odeint.hpp>
#include <fstream>

using namespace boost::numeric::odeint;
using namespace std;

typedef std::vector< double > state_type;

double Vi = 24;
double wref = 0.8;
double A = 15;
//konstansok
double L = 0.0077;
double R = 4;
double psi = 0.115;
double theta = 3.5/(10*10*10*10*10);
double B = 4.26/(10*10*10*10*10);
double frq = 346;
double Tload = 0.01;
double Amplitudo = 8;

const int RESOLUTION = 1000;
const double LOWERBOUND = 0.0;
const double UPPERBOUND =1;

double sawtooth(double t) {
    double result = t;
    while (result > 1) {
        result = result-1;
    }
    return result;
}

double sgn(double x){
    if(x > 0){ return 1; }
    if(x < 0){ return -1; }
    if(x == 0){ return 0; }
}

double q(double t, double w){ //a kapcsoló függvénye w - szögsebesség, t - idő
    return (sgn(Amplitudo*sawtooth(t*frq)-A*((w/Amplitudo)-wref))+1)/2;
}

void buck_motor(const state_type &x, state_type &dxdt, const double t)
{
    //egyenlet
    dxdt[0] = -x[0]*(R/L) - x[1]*(psi/L) + (Vi/L)*q(t, x[1]); //x[0] == I valamint x[1] == w
    dxdt[1] = x[0]*(psi/theta) - x[1]*(B/theta) - Tload/theta;
}

struct push_back_state_and_time
{
    std::vector< state_type >& m_states;
    std::vector< double >& m_times;

    push_back_state_and_time( std::vector< state_type > &states , std::vector< double > &times )
            : m_states( states ) , m_times( times ) { }

    void operator()( const state_type &x , double t )
    {
        m_states.push_back( x );
        m_times.push_back( t );
    }
};

int main() {

    state_type x(2);
    x[0] = 0.0; // start at i=0.0, w=0.0
    x[1] = 0.0;

    int arraySize = (int)(UPPERBOUND-LOWERBOUND)*RESOLUTION;

    size_t steps[arraySize] = {0};

    vector<state_type> x_vec;
    vector<double> times;

    for (double j = 0.0; j <arraySize; j++) {
        A = (double)LOWERBOUND + j/RESOLUTION;
        steps[(int)j] = integrate(buck_motor,
                                 x, 0.0, 1.0, 0.00001,
                                 push_back_state_and_time(x_vec, times));
    }




/* output */
    ofstream myfile1;
    myfile1.open("results.csv");
    for (double k = 0; k < arraySize; ++k) {
        cout<< steps[(int)k] << endl;
        A = (double)LOWERBOUND + k/RESOLUTION;
        bool flag = true; int counter = 0;
        for( size_t i=1; i<=steps[(int)k]; i++ ) {
            if(fmod(times[i], 1/frq)*frq > 0.999999  && flag){
                myfile1 << A << "," << times[i] << "," << x_vec[i][0] << "," << x_vec[i][1] <<'\n';// A, time, i, W
                flag = false;
            }else{
                counter++;
            }
            if(counter >= 100){
                flag = true;
                counter = 0;
            }
        }
    }

    myfile1.close();
    return 0;
}
