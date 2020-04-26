#include <iostream>
#include <vector>
#include <stdio.h>
#include <fstream>
#include <algorithm>
#include <valarray>
#include <boost/numeric/odeint.hpp>

using namespace std;
using namespace boost::numeric::odeint;

typedef std::vector< double > state_type;
typedef vector<vector<double> > data_type;

void write_simulation(const data_type &y_vec, const data_type &betas);
void write_phase_matrices(const data_type &beta_ph0,const data_type &beta_ph1,const data_type &beta_ph2,const data_type &beta_ph3);
void write_fast_results(const data_type &fast_results);
void write_optimisations(data_type optim_results);

//Write a vector<vector<class T> > to an output stream. Output is formatted as a matrix.
template < class T >
std::ostream& operator << (std::ostream& os, const std::vector<std::vector<T> >& v)
{
    for (typename std::vector<std::vector<T> >::const_iterator i = v.begin(); i != v.end(); ++i){
        for (typename std::vector<T>::const_iterator ii = i->begin(); ii != i->end(); ++ii)
        {
            if (ii == i->begin()){
                os << *ii;
            } else {
                os << "," << *ii;
            }
        }
        os << "\n";
    }
    return os;
}

template < class T >
std::ostream& operator << (std::ostream& os, const std::vector<T> & v)
{
    for (typename std::vector<T>::const_iterator i = v.begin(); i != v.end(); ++i){
        if (i == v.begin()){
            os << *i;
        } else {
            os << "," << *i;
        }
    }
    //os << "\n";
    return os;
}

double GenTime(double T2, double R0){
  double G = T2 * ((R0-1)/log(2));
  return(G);
}

data_type parse2DCsvFile(string inputFileName) {

    data_type data;
    ifstream inputFile(inputFileName);
    int l = 0;

    while (inputFile) {
        l++;
        string s;
        if (!getline(inputFile, s)) break;
        if (s[0] != '#') {
            istringstream ss(s);
            vector<double> record;

            while (ss) {
                string line;
                if (!getline(ss, line, ','))
                    break;
                try {
                    record.push_back(stof(line));
                }
                catch (const std::invalid_argument e) {
                    cout << "NaN found in file " << inputFileName << " line " << l
                         << endl;
                    e.what();
                }
            }

            data.push_back(record);
        }
    }

    if (!inputFile.eof()) {
        cerr << "Could not read file " << inputFileName << "\n";
        __throw_invalid_argument("File not found.");
    }

    return data;
}


double linearBeta(double time, double betaBegin, double betaEnd, double timeBegin, double duration){
  double inc = (betaEnd - betaBegin) / duration;
  double betaOut = betaBegin + inc * (time - timeBegin);
  return(betaOut);
}

data_type calc_beta_3Comp_taper(double t, int num_comps, vector<double> sd, data_type beta_ph0,data_type beta_ph1,data_type beta_ph2,data_type beta_ph3 ){
    data_type result;
    if (t<=sd[0]){
        result = beta_ph0;
    } else if ((t>sd[0]) && (t<=sd[1])){
        result = beta_ph1;
    } else if((t>sd[1]) && (t<=sd[2])) {
        for (int i=0; i<num_comps; ++i){
            vector<double> temp;
            for (int j=0; j<num_comps; ++j){
                temp.push_back(linearBeta(t, beta_ph1[i][j], beta_ph2[i][j], sd[1],sd[2]-sd[1]));
            }
            result.push_back(temp);
        }
    } else if(t>sd[2] && (t<=sd[3])) {
        result = beta_ph3;
    } else if(t>sd[3]){
        result = beta_ph3;
    }
    return result;
}


class SIRS_3Comp
{
    int m_num_comps;
    double m_gamma, m_omega,m_sd1, m_sd2, m_sd3, m_sd4;
    data_type m_beta_ph0,m_beta_ph1,m_beta_ph2,m_beta_ph3;
    vector<double> m_fj;

    public:
        SIRS_3Comp(int num_comps, double gamma, double omega, double sd1, double sd2, double sd3,double sd4, data_type beta_ph0,data_type beta_ph1,data_type beta_ph2,data_type beta_ph3, vector<double> fj):
            m_num_comps(num_comps), m_gamma(gamma), m_omega(omega), m_sd1(sd1), m_sd2(sd2),m_sd3(sd3),m_sd4(sd4),m_beta_ph0(beta_ph0),m_beta_ph1(beta_ph1),m_beta_ph2(beta_ph2),m_beta_ph3(beta_ph3),m_fj(fj) {}
        void operator() (const state_type &y, state_type &f, const double t)
        {
            data_type m_beta = calc_beta_3Comp_taper(t, m_num_comps, {m_sd1,m_sd2,m_sd3,m_sd4}, m_beta_ph0, m_beta_ph1,m_beta_ph2,m_beta_ph3);
            for (int i=0;i<m_num_comps;++i){
                double sum_beta=0;
                for (int j=0;j<m_num_comps; ++j){
                    sum_beta += m_beta[i][j]*y[i]*y[j+m_num_comps];
                }
                f[i] = -1.0*sum_beta + m_omega*y[i+(2*m_num_comps)];
                f[i+m_num_comps] = sum_beta - m_gamma*y[i+m_num_comps];
                f[i+(2*m_num_comps)] = m_gamma*y[i+m_num_comps] - m_omega*y[i+(2*m_num_comps)];
                f[i+(3*m_num_comps)] = sum_beta;
             }
        }
};


struct push_back_state_and_time
{
    vector< state_type >& m_states;
    vector< double >& m_times;
    vector< double > m_fj;
    int m_num_comp;

    push_back_state_and_time(vector< state_type > &states , vector< double > &times, vector<double> fj, int num_comp)
    : m_states( states ) , m_times( times ), m_fj(fj), m_num_comp(num_comp) { }

    void operator()( const state_type &x , double t )
    {
        state_type x1 = x;
        x1.insert(x1.begin(),t);
        //0 1  2  3   4   5   6  7  8   9   10  11 12 13  14  15  16    17    18     19     20     21    22    23    24  25  26
        //t,Sv,Sh,Sr1,Sr2,Sr3,Iv,Ih,Ir1,Ir2,Ir3,Rv,Rh,Rr1,Rr2,Rr3,cumIv,cumIh,cumIr1,cumIr2,cumIr3,sumSr,sumIr,sumRr,pIv,pIh,pIr
        //x1.insert(x1.end(),{x1[3]+x1[4]+x1[5],  x1[8]+x1[9]+x1[10], x1[13]+x1[14]+x1[15], x1[6]/m_fj[0], x1[7]/m_fj[1]});
        //x1.insert(x1.end(), x1[22]/m_fj[2]);
        m_states.push_back( x1 );
        m_times.push_back( t );
    }
};

struct push_back_time_Iv_prop
{
    state_type & m_states;
    std::vector< double >& m_times;

    push_back_time_Iv_prop(state_type  &states , vector< double > &times )
    : m_states( states ) , m_times( times ) { }

    void operator()( const state_type &x , double t )
    {

            m_states.push_back( (x[2]) );
            m_times.push_back( t );


    }
};

int main(int /* argc */ , char** /* argv */ )
{
    double omega;
    double  sd1 = 0, sd2 = 0, sd3 = 0,sd4=0;
    int simtime;
    string output_dir, filename;

    double R0 = 2.8;
    double T2 = 3.3;

    double gamma = 1.0/GenTime(T2,R0);

    double beta1_ph0 = 1.7*gamma;
    double beta2_ph0 = 1.7*gamma;
    double beta3_ph0 = 1.7*gamma;
    double beta4_ph0 = 1.7*gamma;


    double beta1_ph1 = 0.8*gamma;
    double beta2_ph1 = 0.9*gamma;
    double beta3_ph1 = 0.9*gamma;
    double beta4_ph1 = 0.8*gamma;

    double beta1_ph2 = 0.4*gamma;
    double beta2_ph2 = 1.85*gamma;
    double beta3_ph2 = 2.25*gamma;
    double beta4_ph2 = 0.4*gamma;

    double beta1_ph3 = 0.4*gamma;
    double beta2_ph3 = 1.85*gamma;
    double beta3_ph3 = 2.25*gamma;
    double beta4_ph3 = 0.4*gamma;


    int n_comp=50;
    data_type beta_ph0(n_comp, vector<double>(n_comp,0)),beta_ph1(n_comp, vector<double>(n_comp,0)),beta_ph2(n_comp, vector<double>(n_comp,0)),beta_ph3(n_comp, vector<double>(n_comp,0));

    for (int i=0; i<n_comp;++i){
        for (int j=0; j<n_comp;++j){
            if (i==0){
                if (j<2){
                    beta_ph0[i][j] = beta1_ph0;
                    beta_ph1[i][j] = beta1_ph1;
                    beta_ph2[i][j] = beta1_ph2;
                    beta_ph3[i][j] = beta1_ph3;
                } else {
                    beta_ph0[i][j] = beta4_ph0;
                    beta_ph1[i][j] = beta4_ph1;
                    beta_ph2[i][j] = beta4_ph2;
                    beta_ph3[i][j] = beta4_ph3;
                }
            } else if (i==1){
                if (j==0){
                    beta_ph0[i][j] = beta1_ph0;
                    beta_ph1[i][j] = beta1_ph1;
                    beta_ph2[i][j] = beta1_ph2;
                    beta_ph3[i][j] = beta1_ph3;
                }
                if (j==1){
                    beta_ph0[i][j] = beta1_ph0;
                    beta_ph1[i][j] = beta1_ph1;
                    beta_ph2[i][j] = beta1_ph2;
                    beta_ph3[i][j] = beta1_ph3;
                }
                if(j>1){
                    beta_ph0[i][j] = beta2_ph0;
                    beta_ph1[i][j] = beta2_ph1;
                    beta_ph2[i][j] = beta2_ph2;
                    beta_ph3[i][j] = beta2_ph3;
                }
            } else if (i>1){
                if (j==0){
                    beta_ph0[i][j] = beta4_ph0;
                    beta_ph1[i][j] = beta4_ph1;
                    beta_ph2[i][j] = beta4_ph2;
                    beta_ph3[i][j] = beta4_ph3;
                } else if (j==1){
                    beta_ph0[i][j] = beta2_ph0;
                    beta_ph1[i][j] = beta2_ph1;
                    beta_ph2[i][j] = beta2_ph2;
                    beta_ph3[i][j] = beta2_ph3;
                } else {
                    beta_ph0[i][j] = beta3_ph0;
                    beta_ph1[i][j] = beta3_ph1;
                    beta_ph2[i][j] = beta3_ph2;
                    beta_ph3[i][j] = beta3_ph3;
                }
            }
        }
    }

    data_type beta_ph0_1 = {
        {beta1_ph0, beta1_ph0, beta4_ph0,beta4_ph0,beta4_ph0},
        {beta1_ph0, beta1_ph0/2.0, beta2_ph0,beta2_ph0,beta2_ph0},
        {beta4_ph0, beta2_ph0, beta3_ph0,beta3_ph0,beta3_ph0},
        {beta4_ph0, beta2_ph0, beta3_ph0,beta3_ph0,beta3_ph0},
        {beta4_ph0, beta2_ph0, beta3_ph0,beta3_ph0,beta3_ph0}
    };

    data_type beta_ph1_1 = {
        {beta1_ph1, beta1_ph1, beta4_ph1,beta4_ph1,beta4_ph1},
        {beta1_ph1, beta1_ph1/2.0, beta2_ph1,beta2_ph1,beta2_ph1},
        {beta4_ph1, beta2_ph1, beta3_ph1,beta3_ph1,beta3_ph1},
        {beta4_ph1, beta2_ph1, beta3_ph1,beta3_ph1,beta3_ph1},
        {beta4_ph1, beta2_ph1, beta3_ph1,beta3_ph1,beta3_ph1}
    };

    data_type beta_ph2_1 = {
        {beta1_ph2, beta1_ph2, beta4_ph2,beta4_ph2,beta4_ph2},
        {beta1_ph2, beta1_ph2/2.0, beta2_ph2,beta2_ph2,beta2_ph2},
        {beta4_ph2, beta2_ph2, beta3_ph2,beta3_ph2,beta3_ph2},
        {beta4_ph2, beta2_ph2, beta3_ph2,beta3_ph2,beta3_ph2},
        {beta4_ph2, beta2_ph2, beta3_ph2,beta3_ph2,beta3_ph2}
    };

    data_type beta_ph3_1 = {
        {beta1_ph3, beta1_ph3, beta4_ph3,beta4_ph3,beta4_ph3},
        {beta1_ph3, beta1_ph3/2.0, beta2_ph3,beta2_ph3,beta2_ph3},
        {beta4_ph3, beta2_ph3, beta3_ph3,beta3_ph3,beta3_ph3},
        {beta4_ph3, beta2_ph3, beta3_ph3,beta3_ph3,beta3_ph3},
        {beta4_ph3, beta2_ph3, beta3_ph3,beta3_ph3,beta3_ph3}
    };
    omega = 1.0/365.0;

    simtime = 365*2;

    sd1 = 71.0;
    sd2 = sd1+(6*7);
    sd3 = sd2+(12*7);
    sd4 = simtime;
    vector<double> fj = {0.2,0.2,0.6};
    int tot_comps = n_comp;
    int r_comps = tot_comps-2;
    double fr_per_comp = fj[2]/r_comps;
    double i0=0.0001;
    bool simulation=false, fast_analysis=false, optimisation = false;
    simulation = true;
    //fast_analysis = true;
    //optimisation = true;

    state_type y(4*tot_comps);
    if(simulation){
        for(int i=0;i<tot_comps-r_comps;++i){
            y[i]= fj[i]-(fj[i]*i0);
            y[i+tot_comps] = fj[i]*i0;
            y[i+2*tot_comps] = 0;
            y[i+3*tot_comps] = 0;
        }
        for(int i=tot_comps-r_comps; i<tot_comps; ++i){
            y[i]= fr_per_comp - (fr_per_comp*i0);
            y[i+tot_comps] = fr_per_comp*i0;
            y[i+2*tot_comps] = 0;
            y[i+3*tot_comps] = 0;
        }
        runge_kutta4< state_type > stepper;
        vector<double> times;
        SIRS_3Comp vuln(tot_comps,gamma,omega,sd1, sd2, sd3, sd4,beta_ph0, beta_ph1,beta_ph2,beta_ph3,fj);
        vector<state_type> y_vec;
        integrate_n_steps(stepper,vuln, y, 0.0 , 1.0, simtime, push_back_state_and_time(y_vec, times,fj,tot_comps));
        data_type betas;
        for (size_t t = 0; t<y_vec.size(); ++t){
            data_type tempbetas = calc_beta_3Comp_taper(t,tot_comps,{sd1,sd2,sd3,sd4}, beta_ph0, beta_ph1,beta_ph2,beta_ph3);
            vector<double> temp;
            for (auto i : tempbetas){
                for (auto j : i){
                    temp.push_back(j);
                }
            }
            betas.push_back(temp);
        }
        write_simulation(y_vec, betas);
        write_phase_matrices(beta_ph0,beta_ph1, beta_ph2, beta_ph3);
    }

    if(fast_analysis){
        data_type fast_results;
        string file("c:/temp/output/FAST_paras_9betas_all_phases.csv");
        vector<vector<double> > paras = parse2DCsvFile(file);
        for (size_t i=0; i<paras.size();++i){
            for(int i=0;i<tot_comps-r_comps;++i){
                y[i]= fj[i]-(fj[i]*i0);
                y[i+tot_comps] = fj[i]*i0;
                y[i+2*tot_comps] = 0;
                y[i+3*tot_comps] = 0;
            }
            for(int i=tot_comps-r_comps; i<tot_comps; ++i){
                y[i]= fr_per_comp - (fr_per_comp*i0);
                y[i+tot_comps] = fr_per_comp*i0;
                y[i+2*tot_comps] = 0;
                y[i+3*tot_comps] = 0;
            }

            beta_ph0 = {
                {beta1_ph0*paras[i][0], beta1_ph0*paras[i][1], beta4_ph0*paras[i][2],beta4_ph0*paras[i][2],beta4_ph0*paras[i][2]},
                {beta1_ph0*paras[i][3], beta1_ph0*paras[i][4], beta2_ph0*paras[i][5],beta2_ph0*paras[i][5],beta2_ph0*paras[i][5]},
                {beta4_ph0*paras[i][6], beta2_ph0*paras[i][7], beta3_ph0*paras[i][8],beta3_ph0*paras[i][8],beta3_ph0*paras[i][8]},
                {beta4_ph0*paras[i][6], beta2_ph0*paras[i][7], beta3_ph0*paras[i][8],beta3_ph0*paras[i][8],beta3_ph0*paras[i][8]},
                {beta4_ph0*paras[i][6], beta2_ph0*paras[i][7], beta3_ph0*paras[i][8],beta3_ph0*paras[i][8],beta3_ph0*paras[i][8]}
            };
            beta_ph1 = {
                {beta1_ph1*paras[i][0], beta1_ph1*paras[i][1], beta4_ph1*paras[i][2],beta4_ph1*paras[i][2],beta4_ph1*paras[i][2]},
                {beta1_ph1*paras[i][3], beta1_ph1*paras[i][4], beta2_ph1*paras[i][5],beta2_ph1*paras[i][5],beta2_ph1*paras[i][5]},
                {beta4_ph1*paras[i][6], beta2_ph1*paras[i][7], beta3_ph1*paras[i][8],beta3_ph1*paras[i][8],beta3_ph1*paras[i][8]},
                {beta4_ph1*paras[i][6], beta2_ph1*paras[i][7], beta3_ph1*paras[i][8],beta3_ph1*paras[i][8],beta3_ph1*paras[i][8]},
                {beta4_ph1*paras[i][6], beta2_ph1*paras[i][7], beta3_ph1*paras[i][8],beta3_ph1*paras[i][8],beta3_ph1*paras[i][8]}
            };
            beta_ph2 = {
                {beta1_ph2*paras[i][0], beta1_ph2*paras[i][1], beta4_ph2*paras[i][2],beta4_ph2*paras[i][2],beta4_ph2*paras[i][2]},
                {beta1_ph2*paras[i][3], beta1_ph2*paras[i][4], beta2_ph2*paras[i][5],beta2_ph2*paras[i][5],beta2_ph2*paras[i][5]},
                {beta4_ph2*paras[i][6], beta2_ph2*paras[i][7], beta3_ph2*paras[i][8],beta3_ph2*paras[i][8],beta3_ph2*paras[i][8]},
                {beta4_ph2*paras[i][6], beta2_ph2*paras[i][7], beta3_ph2*paras[i][8],beta3_ph2*paras[i][8],beta3_ph2*paras[i][8]},
                {beta4_ph2*paras[i][6], beta2_ph2*paras[i][7], beta3_ph2*paras[i][8],beta3_ph2*paras[i][8],beta3_ph2*paras[i][8]}
           };
            beta_ph3 = beta_ph2;
            runge_kutta4< state_type > stepper;
            vector<double> times;

            SIRS_3Comp vuln(tot_comps,gamma,omega,sd1, sd2, sd3, sd4,beta_ph0, beta_ph1,beta_ph2,beta_ph3,fj);
            vector<state_type> y_vec;
            integrate_n_steps(stepper,vuln, y, 0.0 , 1.0, simtime, push_back_state_and_time(y_vec, times,fj,tot_comps));
            std::vector<double> temp;
            //0 1  2  3   4   5   6  7  8   9   10  11 12 13  14  15  16    17    18     19     20     21    22    23    24  25  26
            //t,Sv,Sh,Sr1,Sr2,Sr3,Iv,Ih,Ir1,Ir2,Ir3,Rv,Rh,Rr1,Rr2,Rr3,cumIv,cumIh,cumIr1,cumIr2,cumIr3,sumSr,sumIr,sumRr,pIv,pIh,pIr
            for (size_t i=0; i<y_vec.size();++i){
                temp.push_back(y_vec[i][24]);
            }
            double  incrIv = 0;
            int sd_2nd_peak_Iv = static_cast<int>(sd1)+300;
            for (size_t i=72; i<sd1+400;++i){
                if (temp[i]>=temp[i-1]){
                    incrIv=1;
                    sd_2nd_peak_Iv = i;
                    break;
                }
            }
            auto it = *max_element(std::begin(temp)+sd_2nd_peak_Iv, std::begin(temp)+(sd1+365));
            double  higher;
            if (it>temp[71]){
                higher = 1.0;
            } else {
                higher = 0.0;
            }
            vector<double> temp1 = {it,y_vec[71+365][16]/fj[0], higher};
            fast_results.push_back(temp1);
        }
        write_fast_results(fast_results);
    }
    if (optimisation){
        data_type optim_results;
        //for 3 values of beta1 (0.0, 0.4, 0.8) and 3 values of beta2 (0.9, 1.85, 2.8) vary beta3 from 0.9 to 2.8 and beta2 from 0.0 to 0.8

        for (double beta1=0.0; beta1<0.81;beta1+=0.4){
            for (double beta2=0.9; beta2<2.81;beta2+=0.95){
                for (double beta3=0.9; beta3<2.801;beta3+=0.01){
                    for (double beta4=0.0; beta4<0.81;beta4+=0.01){
                        //create phase2+3 matrix
//                        double beta1 = 0.8, beta2=2.8,beta3=2.8,beta4=0.8;
                        beta_ph2_1 = {
                                {beta1*gamma, beta1*gamma, beta4*gamma,beta4*gamma,beta4*gamma},
                                {beta1*gamma, beta1*gamma, beta2*gamma,beta2*gamma,beta2*gamma},
                                {beta4*gamma, beta2*gamma, beta3*gamma,beta3*gamma,beta3*gamma},
                                {beta4*gamma, beta2*gamma, beta3*gamma,beta3*gamma,beta3*gamma},
                                {beta4*gamma, beta2*gamma, beta3*gamma,beta3*gamma,beta3*gamma}
                            };
                        beta_ph3_1 = beta_ph2_1;
                        //Initialise starting values
                        for(int i=0;i<tot_comps-r_comps;++i){
                            y[i]= fj[i]-(fj[i]*i0);
                            y[i+tot_comps] = fj[i]*i0;
                            y[i+2*tot_comps] = 0;
                            y[i+3*tot_comps] = 0;
                        }
                        for(int i=tot_comps-r_comps; i<tot_comps; ++i){
                            y[i]= fr_per_comp - (fr_per_comp*i0);
                            y[i+tot_comps] = fr_per_comp*i0;
                            y[i+2*tot_comps] = 0;
                            y[i+3*tot_comps] = 0;
                        }
                        runge_kutta4< state_type > stepper;
                        vector<double> times;

                        SIRS_3Comp vuln(tot_comps,gamma,omega,sd1, sd2, sd3, sd4,beta_ph0_1, beta_ph1_1,beta_ph2_1,beta_ph3_1,fj);
                        vector<state_type> y_vec;
                        integrate_n_steps(stepper,vuln, y, 0.0 , 1.0, simtime, push_back_state_and_time(y_vec, times,fj,tot_comps));
                        std::vector<double> tempIv, tempIh,tempIr;
                        //0 1  2  3   4   5   6  7  8   9   10  11 12 13  14  15  16    17    18     19     20     21    22    23    24  25  26
                        //t,Sv,Sh,Sr1,Sr2,Sr3,Iv,Ih,Ir1,Ir2,Ir3,Rv,Rh,Rr1,Rr2,Rr3,cumIv,cumIh,cumIr1,cumIr2,cumIr3,sumSr,sumIr,sumRr,pIv,pIh,pIr
                        //Outputs of interest:
                        //1) peak 2 Iv <= peak 1 Iv;
                        //2) peak 2 Iv <= peak 1 Iv & peak 2 Is <= peak 1 Is & peak 2 Ir <= peak 1 Ir;
                        //3) No increase at all in Iv, Is or Ir after peak 1
                        for (size_t i=0; i<y_vec.size();++i){
                            tempIv.push_back(y_vec[i][24]);
                            tempIh.push_back(y_vec[i][25]);
                            tempIr.push_back(y_vec[i][26]);
                        }
//                        auto maxIv = *max_element(std::begin(tempIv)+sd1, std::begin(tempIv)+(sd1+365));
//                        auto maxIh = *max_element(std::begin(tempIh)+sd1, std::begin(tempIh)+(sd1+365));
//                        auto maxIr = *max_element(std::begin(tempIr)+sd1, std::begin(tempIr)+(sd1+365));
                        //No increase at all in Iv,Is,Ir means difference in Ix[t+1] - Ix[t] must be 0 or negative from day 72 onwards.
                        double  incrIv = 0,incrIh = 0, incrIr =0;
                        int sd_2nd_peak_Iv = static_cast<int>(sd1)+300,sd_2nd_peak_Ih = static_cast<int>(sd1)+300,sd_2nd_peak_Ir = static_cast<int>(sd1)+300;
                        for (size_t i=72; i<sd1+365;++i){
                            if (tempIv[i]>=tempIv[i-1]){
                                incrIv=1;
                                sd_2nd_peak_Iv = i;
                                break;
                            }
                        }
                        for (size_t i=72; i<sd1+365;++i){
                            if (tempIh[i]>=tempIh[i-1]){
                                incrIh=1;
                                 sd_2nd_peak_Ih = i;
                                break;
                            }
                        }
                        for (size_t i=72; i<sd1+365;++i){
                            if (tempIr[i]>=tempIr[i-1]){
                                incrIr=1;
                                sd_2nd_peak_Ir = i;
                                break;
                            }
                        }
                        auto maxIv = *max_element(std::begin(tempIv)+sd_2nd_peak_Iv, std::begin(tempIv)+(sd1+365));
                        auto maxIh = *max_element(std::begin(tempIh)+sd_2nd_peak_Ih, std::begin(tempIh)+(sd1+365));
                        auto maxIr = *max_element(std::begin(tempIr)+sd_2nd_peak_Ir, std::begin(tempIr)+(sd1+365));

                        //output consists of: beta3,beta4,beta1,beta2,cumIv,1stpeakIv, 2ndpeakIv,1stpeakIh, 2ndpeakIh,1stpeakIr,2ndpeakIr,no_incrIv,no_incrIh,no_incrIr
                        vector<double> temp = {beta1, beta2, beta3, beta4,y_vec[71+365][16]/fj[0], tempIv[71], maxIv,tempIh[71], maxIh,tempIr[71], maxIr, incrIv,incrIh,incrIr};
                        optim_results.push_back(temp);
                        //cout << ".";
                    }
                }
            }
        }
        write_optimisations(optim_results);
    }
    return 0;
}


void write_simulation(const data_type &y_vec, const data_type &betas){
    string filename = "c:/temp/output/SIRS_MComp_96-2-2_Mark.csv";
    ofstream output(filename);
    //0 1  2  3   4   5   6  7  8   9   10  11 12 13  14  15  16    17    18     19     20     21    22    23    24  25  26
    //t,Sv,Sh,Sr1,Sr2,Sr3,Iv,Ih,Ir1,Ir2,Ir3,Rv,Rh,Rr1,Rr2,Rr3,cumIv,cumIh,cumIr1,cumIr2,cumIr3,sumSr,sumIr,sumRr,pIv,pIh,pIr
    //output << "t,Sv,Sh,Sr1,Sr2,Sr3,Iv,Ih,Ir1,Ir2,Ir3,Rv,Rh,Rr1,Rr2,Rr3,cumIv,cumIh,cumIr1,cumIr2,cumIr3,sumSr,sumIr,sumRr,pIv,pIh,pIr,beta11,beta12,beta13,beta14,beta15,beta21,beta22,beta23,beta24,beta25,beta31,beta32,beta33,beta34,beta35,beta41,beta42,beta43,beta44,beta45,beta51,beta52,beta53,beta54,beta55\n";
    for (size_t i=0; i<y_vec.size();++i){
        output << y_vec[i] << "," << betas[i] << "\n";
    }
//    output << y_vec;
    output.close();
}

void write_phase_matrices(const data_type &beta_ph0,const data_type &beta_ph1,const data_type &beta_ph2,const data_type &beta_ph3){
    string filename = "c:/temp/output/SIRS_MComp_96-2-2_betas_Mark.csv";
    ofstream output(filename);
    output << "v,h,r\n";
    output << beta_ph0 << "\n";
    output << beta_ph1 << "\n";
    output << beta_ph2 << "\n";
    output << beta_ph3 << "\n";
    output.close();
}


void write_fast_results(const data_type &fast_results){
    string filename = "c:/temp/output/SIRS_mult_Comp_FAST_results_9betas_all_phases.csv";
    ofstream output(filename);
    output << "peak_Iv, cum_Iv, higher\n";
    output << fast_results;
    output.close();
}

void write_optimisations(data_type optim_results){
    string filename = "c:/temp/output/output/SIRS_mult_Comp_optim_results2.csv";
    ofstream output(filename);
    output << "beta1,beta2,beta3,beta4,cumIv_1y,first_peakIv,sec_peakIv,first_peakIh,sec_peakIh,first_peakIr,sec_peakIr,incrIv,incrIh,incrIr\n";
    output << optim_results;
}
