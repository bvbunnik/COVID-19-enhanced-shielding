#include <stdio.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>
#include <boost/numeric/odeint.hpp>

using namespace std;
using namespace boost::numeric::odeint;

typedef std::vector< double > state_type;
typedef vector<vector<double> > data_type;

void simulation();
void fast_analysis();
void optimisation();
void P3P4_betas();
void compliance();
void vary_Re_P1();
void vary_Re_P2();
void immunity();
void write_simulation(string filename, const data_type &y_vec, const data_type &betas, int v_comps, int h_comps, int r_comps);
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
        m_states.push_back( x1 );
        m_times.push_back( t );
    }
};

struct push_back_state_time_pIv_pIh_pIr
{
    vector< state_type >& m_states;
    vector< double >& m_times;
    vector< double > m_fj;
    int m_v_comps, m_h_comps, m_r_comps, m_n_comp;

    push_back_state_time_pIv_pIh_pIr(vector< state_type > &states , vector< double > &times, vector<double> fj, int v_comps,int h_comps, int r_comps)
    : m_states( states ) , m_times( times ), m_fj(fj), m_v_comps(v_comps), m_h_comps(h_comps), m_r_comps(r_comps), m_n_comp(v_comps+h_comps+r_comps) { }

    void operator()( const state_type &x , double t )
    {
        state_type x1 = x;
        x1.insert(x1.begin(),t);
        //insert sum of comps
        //I_comps start at i=n_comp, the first v_comp are Iv, then h_comp Ih and lastly r_comp Ir.
        double cumIv,Iv, Ih, Ir;
        cumIv = accumulate(x.begin()+3*m_n_comp, x.begin()+3*m_n_comp+m_v_comps, 0.0);
        Iv = accumulate(x.begin()+m_n_comp, x.begin()+m_n_comp+m_v_comps, 0.0);
        Ih = accumulate(x.begin()+m_n_comp+m_v_comps, x.begin()+m_n_comp+m_v_comps+m_h_comps, 0.0);
        Ir = accumulate(x.begin()+m_n_comp+m_v_comps+m_h_comps, x.begin()+m_n_comp+m_n_comp, 0.0);
        x1.insert(x1.end(), {cumIv, Iv/m_fj[0],Ih/m_fj[1],Ir/m_fj[2]});
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
    simulation();
    //fast_analysis();
    //optimisation();
    //P3P4_betas();
    //compliance();
    //vary_Re_P1();
    //vary_Re_P2();
    //immunity();
    return 0;
}

void simulation(){
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

    int n_comp=5;
    //2-2-96 -> 1v-1v-48r
    //8-8-84 -> 4v-4h-42r
    //14-14-72 -> 7v-7h-36r
    //20-20-60 -> 10v-10h-30r or 1v-1h-3r
    //20-40-40 -> 1v-2h-2r
    //20-10-70 -> 2v-1h-7r
    vector<double> fj = {0.2,0.2,0.6};
    int v_comps = 1;
    int h_comps = 1;
    int r_comps = n_comp - v_comps - h_comps;
    double fv_per_comp = fj[0]/v_comps;
    double fh_per_comp = fj[1]/h_comps;
    double fr_per_comp = fj[2]/r_comps;
    data_type beta_ph0(n_comp, vector<double>(n_comp,0)),beta_ph1(n_comp, vector<double>(n_comp,0)),beta_ph2(n_comp, vector<double>(n_comp,0)),beta_ph3(n_comp, vector<double>(n_comp,0));
    omega = 1.0/365.0;

    simtime = 365*2;

    sd1 = 71.0;
    sd2 = sd1+(6*7);
    sd3 = sd2+(12*7);
    sd4 = simtime;

    double i0=0.0001;
    filename = "SIRS_MComp_simulation_20-20-60.csv";
    for (int i=0; i<n_comp;++i){
        for (int j=0; j<n_comp;++j){
            if (i<v_comps){
                if (j<v_comps+h_comps){
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
            } else if (i>=v_comps & i<v_comps+h_comps){
                if (j<v_comps){
                    beta_ph0[i][j] = beta1_ph0;
                    beta_ph1[i][j] = beta1_ph1;
                    beta_ph2[i][j] = beta1_ph2;
                    beta_ph3[i][j] = beta1_ph3;
                } else if (j>=v_comps & j<v_comps + h_comps){
                    beta_ph0[i][j] = beta1_ph0;
                    beta_ph1[i][j] = beta1_ph1;
                    beta_ph2[i][j] = beta1_ph2;
                    beta_ph3[i][j] = beta1_ph3;
                } else {
                    beta_ph0[i][j] = beta2_ph0;
                    beta_ph1[i][j] = beta2_ph1;
                    beta_ph2[i][j] = beta2_ph2;
                    beta_ph3[i][j] = beta2_ph3;
                }
            } else if (i>=v_comps + h_comps){
                if (j<v_comps){
                    beta_ph0[i][j] = beta4_ph0;
                    beta_ph1[i][j] = beta4_ph1;
                    beta_ph2[i][j] = beta4_ph2;
                    beta_ph3[i][j] = beta4_ph3;
                } else if (j<v_comps + h_comps){
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

    state_type y(4*n_comp);
    for(int i=0;i<v_comps;++i){
        y[i]= fv_per_comp - (fv_per_comp*i0);
        y[i+n_comp] = fv_per_comp*i0;
        y[i+2*n_comp] = 0;
        y[i+3*n_comp] = 0;
    }
    for(int i=v_comps;i<v_comps+h_comps;++i){
        y[i]= fh_per_comp - (fh_per_comp*i0);
        y[i+n_comp] = fh_per_comp*i0;
        y[i+2*n_comp] = 0;
        y[i+3*n_comp] = 0;
    }
    for(int i=v_comps+h_comps; i<n_comp; ++i){
        y[i]= fr_per_comp - (fr_per_comp*i0);
        y[i+n_comp] = fr_per_comp*i0;
        y[i+2*n_comp] = 0;
        y[i+3*n_comp] = 0;
    }
    runge_kutta4< state_type > stepper;
    vector<double> times;
    SIRS_3Comp vuln(n_comp,gamma,omega,sd1, sd2, sd3, sd4,beta_ph0, beta_ph1,beta_ph2,beta_ph3,fj);
    vector<state_type> y_vec;
    integrate_n_steps(stepper,vuln, y, 0.0 , 1.0, simtime, push_back_state_time_pIv_pIh_pIr(y_vec, times,fj,v_comps, h_comps,r_comps));
    data_type betas;
    for (size_t t = 0; t<y_vec.size(); ++t){
        data_type tempbetas = calc_beta_3Comp_taper(t,n_comp,{sd1,sd2,sd3,sd4}, beta_ph0, beta_ph1,beta_ph2,beta_ph3);
        vector<double> temp;
        for (auto i : tempbetas){
            for (auto j : i){
                temp.push_back(j);
            }
        }
        betas.push_back(temp);
    }
    write_simulation(filename, y_vec, betas,v_comps, h_comps, r_comps);
    write_phase_matrices(beta_ph0,beta_ph1, beta_ph2, beta_ph3);
}


void fast_analysis(){
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

    int n_comp=10;
    //2-2-96 -> 1v-1v-48r
    //8-8-84 -> 4v-4h-42r
    //14-14-72 -> 7v-7h-36r
    //20-20-60 -> 10v-10h-30r or 1v-1h-3r
    //20-40-40 -> 1v-2h-2r
    //20-10-70 -> 2v-1h-7r
    vector<double> fj = {0.2,0.1,0.7};
    int v_comps = 2;
    int h_comps = 1;
    int r_comps = n_comp - v_comps - h_comps;
    double fv_per_comp = fj[0]/v_comps;
    double fh_per_comp = fj[1]/h_comps;
    double fr_per_comp = fj[2]/r_comps;
    data_type beta_ph0(n_comp, vector<double>(n_comp,0)),beta_ph1(n_comp, vector<double>(n_comp,0)),beta_ph2(n_comp, vector<double>(n_comp,0)),beta_ph3(n_comp, vector<double>(n_comp,0));
    omega = 1.0/365.0;

    simtime = 365*2;

    sd1 = 71.0;
    sd2 = sd1+(6*7);
    sd3 = sd2+(12*7);
    sd4 = simtime;

    double i0=0.0001;

    for (int i=0; i<n_comp;++i){
        for (int j=0; j<n_comp;++j){
            if (i<v_comps){
                if (j<v_comps+h_comps){
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
            } else if (i>=v_comps & i<v_comps+h_comps){
                if (j<v_comps){
                    beta_ph0[i][j] = beta1_ph0;
                    beta_ph1[i][j] = beta1_ph1;
                    beta_ph2[i][j] = beta1_ph2;
                    beta_ph3[i][j] = beta1_ph3;
                } else if (j>=v_comps & j<v_comps + h_comps){
                    beta_ph0[i][j] = beta1_ph0;
                    beta_ph1[i][j] = beta1_ph1;
                    beta_ph2[i][j] = beta1_ph2;
                    beta_ph3[i][j] = beta1_ph3;
                } else {
                    beta_ph0[i][j] = beta2_ph0;
                    beta_ph1[i][j] = beta2_ph1;
                    beta_ph2[i][j] = beta2_ph2;
                    beta_ph3[i][j] = beta2_ph3;
                }
            } else if (i>=v_comps + h_comps){
                if (j<v_comps){
                    beta_ph0[i][j] = beta4_ph0;
                    beta_ph1[i][j] = beta4_ph1;
                    beta_ph2[i][j] = beta4_ph2;
                    beta_ph3[i][j] = beta4_ph3;
                } else if (j<v_comps + h_comps){
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

    state_type y(4*n_comp);
    data_type fast_results;
    string file("FAST_paras.csv");
    vector<vector<double> > paras = parse2DCsvFile(file);
    for (size_t i=0; i<paras.size();++i){
        for(int i=0;i<v_comps;++i){
            y[i]= fv_per_comp - (fv_per_comp*i0);
            y[i+n_comp] = fv_per_comp*i0;
            y[i+2*n_comp] = 0;
            y[i+3*n_comp] = 0;
        }
        for(int i=v_comps;i<v_comps+h_comps;++i){
            y[i]= fh_per_comp - (fh_per_comp*i0);
            y[i+n_comp] = fh_per_comp*i0;
            y[i+2*n_comp] = 0;
            y[i+3*n_comp] = 0;
        }
        for(int i=v_comps+h_comps; i<n_comp; ++i){
            y[i]= fr_per_comp - (fr_per_comp*i0);
            y[i+n_comp] = fr_per_comp*i0;
            y[i+2*n_comp] = 0;
            y[i+3*n_comp] = 0;
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

        SIRS_3Comp vuln(n_comp,gamma,omega,sd1, sd2, sd3, sd4,beta_ph0, beta_ph1,beta_ph2,beta_ph3,fj);
        vector<state_type> y_vec;
        integrate_n_steps(stepper,vuln, y, 0.0 , 1.0, simtime, push_back_state_and_time(y_vec, times,fj,n_comp));
        std::vector<double> temp;
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

void optimisation(){
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

    int n_comp=10;
    //2-2-96 -> 1v-1v-48r
    //8-8-84 -> 4v-4h-42r
    //14-14-72 -> 7v-7h-36r
    //20-20-60 -> 10v-10h-30r or 1v-1h-3r
    //20-40-40 -> 1v-2h-2r
    //20-10-70 -> 2v-1h-7r
    vector<double> fj = {0.2,0.1,0.7};
    int v_comps = 2;
    int h_comps = 1;
    int r_comps = n_comp - v_comps - h_comps;
    double fv_per_comp = fj[0]/v_comps;
    double fh_per_comp = fj[1]/h_comps;
    double fr_per_comp = fj[2]/r_comps;
    data_type beta_ph0(n_comp, vector<double>(n_comp,0)),beta_ph1(n_comp, vector<double>(n_comp,0)),beta_ph2(n_comp, vector<double>(n_comp,0)),beta_ph3(n_comp, vector<double>(n_comp,0));
    omega = 1.0/365.0;

    simtime = 365*2;

    sd1 = 71.0;
    sd2 = sd1+(6*7);
    sd3 = sd2+(12*7);
    sd4 = simtime;

    double i0=0.0001;

    for (int i=0; i<n_comp;++i){
        for (int j=0; j<n_comp;++j){
            if (i<v_comps){
                if (j<v_comps+h_comps){
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
            } else if (i>=v_comps & i<v_comps+h_comps){
                if (j<v_comps){
                    beta_ph0[i][j] = beta1_ph0;
                    beta_ph1[i][j] = beta1_ph1;
                    beta_ph2[i][j] = beta1_ph2;
                    beta_ph3[i][j] = beta1_ph3;
                } else if (j>=v_comps & j<v_comps + h_comps){
                    beta_ph0[i][j] = beta1_ph0;
                    beta_ph1[i][j] = beta1_ph1;
                    beta_ph2[i][j] = beta1_ph2;
                    beta_ph3[i][j] = beta1_ph3;
                } else {
                    beta_ph0[i][j] = beta2_ph0;
                    beta_ph1[i][j] = beta2_ph1;
                    beta_ph2[i][j] = beta2_ph2;
                    beta_ph3[i][j] = beta2_ph3;
                }
            } else if (i>=v_comps + h_comps){
                if (j<v_comps){
                    beta_ph0[i][j] = beta4_ph0;
                    beta_ph1[i][j] = beta4_ph1;
                    beta_ph2[i][j] = beta4_ph2;
                    beta_ph3[i][j] = beta4_ph3;
                } else if (j<v_comps + h_comps){
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

    state_type y(4*n_comp);
    data_type optim_results;
    //for 3 values of beta1 (0.0, 0.4, 0.8) and 3 values of beta2 (0.9, 1.85, 2.8) vary beta3 from 0.9 to 2.8 and beta2 from 0.0 to 0.8

    for (double beta1=0.0; beta1<0.81;beta1+=0.4){
        for (double beta2=0.9; beta2<2.81;beta2+=0.95){
            for (double beta3=0.9; beta3<2.801;beta3+=0.01){
                for (double beta4=0.0; beta4<0.81;beta4+=0.01){

                    for (int i=0; i<n_comp;++i){
                        for (int j=0; j<n_comp;++j){
                            if (i<v_comps){
                                if (j<v_comps+h_comps){
                                    beta_ph2[i][j] = beta1*gamma;
                                    beta_ph3[i][j] = beta1*gamma;
                                } else {
                                    beta_ph2[i][j] = beta4*gamma;
                                    beta_ph3[i][j] = beta4*gamma;
                                }
                            } else if (i>=v_comps & i<v_comps+h_comps){
                                if (j<v_comps){
                                    beta_ph2[i][j] = beta1*gamma;
                                    beta_ph3[i][j] = beta1*gamma;
                                } else if (j>=v_comps & j<v_comps + h_comps){
                                    beta_ph2[i][j] = beta1*gamma;
                                    beta_ph3[i][j] = beta1*gamma;
                                } else {
                                    beta_ph2[i][j] = beta2*gamma;
                                    beta_ph3[i][j] = beta2*gamma;
                                }
                            } else if (i>=v_comps + h_comps){
                                if (j<v_comps){
                                    beta_ph2[i][j] = beta4*gamma;
                                    beta_ph3[i][j] = beta4*gamma;
                                } else if (j<v_comps + h_comps){
                                    beta_ph2[i][j] = beta2*gamma;
                                    beta_ph3[i][j] = beta2*gamma;
                                } else {
                                    beta_ph2[i][j] = beta3*gamma;
                                    beta_ph3[i][j] = beta3*gamma;
                                }
                            }
                        }
                    }
                    //Initialise starting values
                    for(int i=0;i<v_comps;++i){
                        y[i]= fv_per_comp - (fv_per_comp*i0);
                        y[i+n_comp] = fv_per_comp*i0;
                        y[i+2*n_comp] = 0;
                        y[i+3*n_comp] = 0;
                    }
                    for(int i=v_comps;i<v_comps+h_comps;++i){
                        y[i]= fh_per_comp - (fh_per_comp*i0);
                        y[i+n_comp] = fh_per_comp*i0;
                        y[i+2*n_comp] = 0;
                        y[i+3*n_comp] = 0;
                    }
                    for(int i=v_comps+h_comps; i<n_comp; ++i){
                        y[i]= fr_per_comp - (fr_per_comp*i0);
                        y[i+n_comp] = fr_per_comp*i0;
                        y[i+2*n_comp] = 0;
                        y[i+3*n_comp] = 0;
                    }
                    runge_kutta4< state_type > stepper;
                    vector<double> times;

                    SIRS_3Comp vuln(n_comp,gamma,omega,sd1, sd2, sd3, sd4,beta_ph0, beta_ph1,beta_ph2,beta_ph3,fj);
                    vector<state_type> y_vec;
                    integrate_n_steps(stepper,vuln, y, 0.0 , 1.0, simtime, push_back_state_time_pIv_pIh_pIr(y_vec, times,fj,v_comps, h_comps,r_comps));
                    std::vector<double> tempIv, tempIh,tempIr;
                    //last three cells of y_vec are pIv, pIh & pIr
                    //Outputs of interest:
                    //1) peak 2 Iv <= peak 1 Iv;
                    //2) peak 2 Iv <= peak 1 Iv & peak 2 Is <= peak 1 Is & peak 2 Ir <= peak 1 Ir;
                    //3) No increase at all in Iv, Is or Ir after peak 1
                    int end = y_vec[0].size();
                    for (size_t i=0; i<y_vec.size();++i){
                        tempIv.push_back(y_vec[i][end-3]);
                        tempIh.push_back(y_vec[i][end-2]);
                        tempIr.push_back(y_vec[i][end-1]);
                    }
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

                    vector<double> temp = {beta1, beta2, beta3, beta4,y_vec[71+365][end-4]/fj[0], tempIv[71], maxIv,tempIh[71], maxIh,tempIr[71], maxIr, incrIv,incrIh,incrIr};
                    optim_results.push_back(temp);
                    //cout << ".";
                }
            }
        }
    }
    write_optimisations(optim_results);
}

void P3P4_betas(){
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

    int n_comp=5;
    //2-2-96 -> 1v-1v-48r
    //8-8-84 -> 4v-4h-42r
    //14-14-72 -> 7v-7h-36r
    //20-20-60 -> 10v-10h-30r or 1v-1h-3r
    //20-40-40 -> 1v-2h-2r
    //20-10-70 -> 2v-1h-7r
    vector<double> fj = {0.2,0.2,0.6};
    int v_comps = 1;
    int h_comps = 1;
    int r_comps = n_comp - v_comps - h_comps;
    double fv_per_comp = fj[0]/v_comps;
    double fh_per_comp = fj[1]/h_comps;
    double fr_per_comp = fj[2]/r_comps;
    data_type beta_ph0(n_comp, vector<double>(n_comp,0)),beta_ph1(n_comp, vector<double>(n_comp,0)),beta_ph2(n_comp, vector<double>(n_comp,0)),beta_ph3(n_comp, vector<double>(n_comp,0));
    omega = 1.0/365.0;

    simtime = 365*2;

    sd1 = 71.0;
    sd2 = sd1+(6*7);
    sd3 = sd2+(12*7);
    sd4 = simtime;
    data_type results;
    double i0=0.0001;
    for (double perc=0.75; perc<1.251; perc+=0.01){
        for (int i=0; i<n_comp;++i){
            for (int j=0; j<n_comp;++j){
                if (i<v_comps){
                    if (j<v_comps+h_comps){
                        beta_ph0[i][j] = beta1_ph0;
                        beta_ph1[i][j] = beta1_ph1;
                        beta_ph2[i][j] = perc*beta1_ph2;
                        beta_ph3[i][j] = perc*beta1_ph3;
                    } else {
                        beta_ph0[i][j] = beta4_ph0;
                        beta_ph1[i][j] = beta4_ph1;
                        beta_ph2[i][j] = perc*beta4_ph2;
                        beta_ph3[i][j] = perc*beta4_ph3;
                    }
                } else if (i>=v_comps & i<v_comps+h_comps){
                    if (j<v_comps){
                        beta_ph0[i][j] = beta1_ph0;
                        beta_ph1[i][j] = beta1_ph1;
                        beta_ph2[i][j] = perc*beta1_ph2;
                        beta_ph3[i][j] = perc*beta1_ph3;
                    } else if (j>=v_comps & j<v_comps + h_comps){
                        beta_ph0[i][j] = beta1_ph0;
                        beta_ph1[i][j] = beta1_ph1;
                        beta_ph2[i][j] = perc*beta1_ph2;
                        beta_ph3[i][j] = perc*beta1_ph3;
                    } else {
                        beta_ph0[i][j] = beta2_ph0;
                        beta_ph1[i][j] = beta2_ph1;
                        beta_ph2[i][j] = perc*beta2_ph2;
                        beta_ph3[i][j] = perc*beta2_ph3;
                    }
                } else if (i>=v_comps + h_comps){
                    if (j<v_comps){
                        beta_ph0[i][j] = beta4_ph0;
                        beta_ph1[i][j] = beta4_ph1;
                        beta_ph2[i][j] = perc*beta4_ph2;
                        beta_ph3[i][j] = perc*beta4_ph3;
                    } else if (j<v_comps + h_comps){
                        beta_ph0[i][j] = beta2_ph0;
                        beta_ph1[i][j] = beta2_ph1;
                        beta_ph2[i][j] = perc*beta2_ph2;
                        beta_ph3[i][j] = perc*beta2_ph3;
                    } else {
                        beta_ph0[i][j] = beta3_ph0;
                        beta_ph1[i][j] = beta3_ph1;
                        beta_ph2[i][j] = perc*beta3_ph2;
                        beta_ph3[i][j] = perc*beta3_ph3;
                    }
                }
            }
        }

        state_type y(4*n_comp);

        //Initialise starting values
        for(int i=0;i<v_comps;++i){
            y[i]= fv_per_comp - (fv_per_comp*i0);
            y[i+n_comp] = fv_per_comp*i0;
            y[i+2*n_comp] = 0;
            y[i+3*n_comp] = 0;
        }
        for(int i=v_comps;i<v_comps+h_comps;++i){
            y[i]= fh_per_comp - (fh_per_comp*i0);
            y[i+n_comp] = fh_per_comp*i0;
            y[i+2*n_comp] = 0;
            y[i+3*n_comp] = 0;
        }
        for(int i=v_comps+h_comps; i<n_comp; ++i){
            y[i]= fr_per_comp - (fr_per_comp*i0);
            y[i+n_comp] = fr_per_comp*i0;
            y[i+2*n_comp] = 0;
            y[i+3*n_comp] = 0;
        }
        runge_kutta4< state_type > stepper;
        vector<double> times;

        SIRS_3Comp vuln(n_comp,gamma,omega,sd1, sd2, sd3, sd4,beta_ph0, beta_ph1,beta_ph2,beta_ph3,fj);
        vector<state_type> y_vec;
        integrate_n_steps(stepper,vuln, y, 0.0 , 1.0, simtime, push_back_state_time_pIv_pIh_pIr(y_vec, times,fj,v_comps, h_comps,r_comps));
        //Calc height first peak and second peak
        vector<double> tempIv;
        int end = y_vec[0].size();
        for (size_t i=0; i<y_vec.size();++i){
            tempIv.push_back(y_vec[i][end-3]);
        }
        double  incrIv = 0;
        int sd_2nd_peak_Iv = static_cast<int>(sd1)+300;
        for (size_t i=72; i<sd1+365;++i){
            if (tempIv[i]>=tempIv[i-1]){
                incrIv=1;
                sd_2nd_peak_Iv = i;
                break;
            }
        }
        auto maxIv = *max_element(std::begin(tempIv)+sd_2nd_peak_Iv, std::begin(tempIv)+(sd1+365));
        auto rel_height = maxIv/tempIv[71];
        results.push_back({perc,rel_height});
    }
    cout << results;
}

void compliance(){
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

    int n_comp=5;
    //2-2-96 -> 1v-1v-48r
    //8-8-84 -> 4v-4h-42r
    //14-14-72 -> 7v-7h-36r
    //20-20-60 -> 10v-10h-30r or 1v-1h-3r
    //20-40-40 -> 1v-2h-2r
    //20-10-70 -> 2v-1h-7r
    vector<double> fj = {0.2,0.2,0.6};
    int v_comps = 1;
    int h_comps = 1;
    int r_comps = n_comp - v_comps - h_comps;
    double fv_per_comp = fj[0]/v_comps;
    double fh_per_comp = fj[1]/h_comps;
    double fr_per_comp = fj[2]/r_comps;
    data_type beta_ph0(n_comp, vector<double>(n_comp,0)),beta_ph1(n_comp, vector<double>(n_comp,0)),beta_ph2(n_comp, vector<double>(n_comp,0)),beta_ph3(n_comp, vector<double>(n_comp,0));
    omega = 1.0/365.0;

    simtime = 365*2;

    sd1 = 71.0;
    sd2 = sd1+(6*7);
    sd3 = sd2+(12*7);
    sd4 = simtime;
    data_type results;
    double i0=0.0001;
    for (double comp=0.0; comp<1.01; comp+=0.01){
        for (int i=0; i<n_comp;++i){
            for (int j=0; j<n_comp;++j){
                if (i<v_comps){
                    if (j<v_comps+h_comps){
                        beta_ph0[i][j] = beta1_ph0;
                        beta_ph1[i][j] = beta1_ph1;
                        beta_ph2[i][j] = beta1_ph2+(1.0-comp)*(beta1_ph0-beta1_ph2);
                        beta_ph3[i][j] = beta1_ph3+(1.0-comp)*(beta1_ph0-beta1_ph3);
                    } else {
                        beta_ph0[i][j] = beta4_ph0;
                        beta_ph1[i][j] = beta4_ph1;
                        beta_ph2[i][j] = beta4_ph2+(1.0-comp)*(beta4_ph0-beta4_ph3);
                        beta_ph3[i][j] = beta4_ph3+(1.0-comp)*(beta4_ph0-beta4_ph3);
                    }
                } else if (i>=v_comps & i<v_comps+h_comps){
                    if (j<v_comps){
                        beta_ph0[i][j] = beta1_ph0;
                        beta_ph1[i][j] = beta1_ph1;
                        beta_ph2[i][j] = beta1_ph2+(1.0-comp)*(beta1_ph0-beta1_ph2);
                        beta_ph3[i][j] = beta1_ph3+(1.0-comp)*(beta1_ph0-beta1_ph3);
                    } else if (j>=v_comps & j<v_comps + h_comps){
                        beta_ph0[i][j] = beta1_ph0;
                        beta_ph1[i][j] = beta1_ph1;
                        beta_ph2[i][j] = beta1_ph2+(1.0-comp)*(beta1_ph0-beta1_ph2);
                        beta_ph3[i][j] = beta1_ph3+(1.0-comp)*(beta1_ph0-beta1_ph3);
                    } else {
                        beta_ph0[i][j] = beta2_ph0;
                        beta_ph1[i][j] = beta2_ph1;
                        beta_ph2[i][j] = beta2_ph2;
                        beta_ph3[i][j] = beta2_ph3;
                    }
                } else if (i>=v_comps + h_comps){
                    if (j<v_comps){
                        beta_ph0[i][j] = beta4_ph0;
                        beta_ph1[i][j] = beta4_ph1;
                        beta_ph2[i][j] = beta4_ph2+(1.0-comp)*(beta4_ph0-beta4_ph3);
                        beta_ph3[i][j] = beta4_ph3+(1.0-comp)*(beta4_ph0-beta4_ph3);
                    } else if (j<v_comps + h_comps){
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

        state_type y(4*n_comp);

        //Initialise starting values
        for(int i=0;i<v_comps;++i){
            y[i]= fv_per_comp - (fv_per_comp*i0);
            y[i+n_comp] = fv_per_comp*i0;
            y[i+2*n_comp] = 0;
            y[i+3*n_comp] = 0;
        }
        for(int i=v_comps;i<v_comps+h_comps;++i){
            y[i]= fh_per_comp - (fh_per_comp*i0);
            y[i+n_comp] = fh_per_comp*i0;
            y[i+2*n_comp] = 0;
            y[i+3*n_comp] = 0;
        }
        for(int i=v_comps+h_comps; i<n_comp; ++i){
            y[i]= fr_per_comp - (fr_per_comp*i0);
            y[i+n_comp] = fr_per_comp*i0;
            y[i+2*n_comp] = 0;
            y[i+3*n_comp] = 0;
        }
        runge_kutta4< state_type > stepper;
        vector<double> times;

        SIRS_3Comp vuln(n_comp,gamma,omega,sd1, sd2, sd3, sd4,beta_ph0, beta_ph1,beta_ph2,beta_ph3,fj);
        vector<state_type> y_vec;
        integrate_n_steps(stepper,vuln, y, 0.0 , 1.0, simtime, push_back_state_time_pIv_pIh_pIr(y_vec, times,fj,v_comps, h_comps,r_comps));
        //Calc height first peak and second peak
        vector<double> tempIv;
        int end = y_vec[0].size();
        for (size_t i=0; i<y_vec.size();++i){
            tempIv.push_back(y_vec[i][end-3]);
        }
        double  incrIv = 0;
        int sd_2nd_peak_Iv = static_cast<int>(sd1)+300;
        for (size_t i=72; i<sd1+365;++i){
            if (tempIv[i]>=tempIv[i-1]){
                incrIv=1;
                sd_2nd_peak_Iv = i;
                break;
            }
        }
        auto maxIv = *max_element(std::begin(tempIv)+sd_2nd_peak_Iv, std::begin(tempIv)+(sd1+365));
        auto rel_height = maxIv/tempIv[71];
        results.push_back({comp,rel_height});
    }
    cout << results;
}

void vary_Re_P1(){
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

    int n_comp=5;
    //2-2-96 -> 1v-1v-48r
    //8-8-84 -> 4v-4h-42r
    //14-14-72 -> 7v-7h-36r
    //20-20-60 -> 10v-10h-30r or 1v-1h-3r
    //20-40-40 -> 1v-2h-2r
    //20-10-70 -> 2v-1h-7r
    vector<double> fj = {0.2,0.2,0.6};
    int v_comps = 1;
    int h_comps = 1;
    int r_comps = n_comp - v_comps - h_comps;
    double fv_per_comp = fj[0]/v_comps;
    double fh_per_comp = fj[1]/h_comps;
    double fr_per_comp = fj[2]/r_comps;
    data_type beta_ph0(n_comp, vector<double>(n_comp,0)),beta_ph1(n_comp, vector<double>(n_comp,0)),beta_ph2(n_comp, vector<double>(n_comp,0)),beta_ph3(n_comp, vector<double>(n_comp,0));
    omega = 1.0/365.0;

    simtime = 365*2;

    sd1 = 71.0;
    sd2 = sd1+(6*7);
    sd3 = sd2+(12*7);
    sd4 = simtime;
    data_type results;
    double i0=0.0001;
    for (double ph0=1.4; ph0<2.01; ph0+=0.01){
        beta1_ph0 = beta2_ph0 = beta3_ph0 = beta4_ph0 = ph0*gamma;
        for (int i=0; i<n_comp;++i){
            for (int j=0; j<n_comp;++j){
                if (i<v_comps){
                    if (j<v_comps+h_comps){
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
                } else if (i>=v_comps & i<v_comps+h_comps){
                    if (j<v_comps){
                        beta_ph0[i][j] = beta1_ph0;
                        beta_ph1[i][j] = beta1_ph1;
                        beta_ph2[i][j] = beta1_ph2;
                        beta_ph3[i][j] = beta1_ph3;
                    } else if (j>=v_comps & j<v_comps + h_comps){
                        beta_ph0[i][j] = beta1_ph0;
                        beta_ph1[i][j] = beta1_ph1;
                        beta_ph2[i][j] = beta1_ph2;
                        beta_ph3[i][j] = beta1_ph3;
                    } else {
                        beta_ph0[i][j] = beta2_ph0;
                        beta_ph1[i][j] = beta2_ph1;
                        beta_ph2[i][j] = beta2_ph2;
                        beta_ph3[i][j] = beta2_ph3;
                    }
                } else if (i>=v_comps + h_comps){
                    if (j<v_comps){
                        beta_ph0[i][j] = beta4_ph0;
                        beta_ph1[i][j] = beta4_ph1;
                        beta_ph2[i][j] = beta4_ph2;
                        beta_ph3[i][j] = beta4_ph3;
                    } else if (j<v_comps + h_comps){
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

        state_type y(4*n_comp);

        //Initialise starting values
        for(int i=0;i<v_comps;++i){
            y[i]= fv_per_comp - (fv_per_comp*i0);
            y[i+n_comp] = fv_per_comp*i0;
            y[i+2*n_comp] = 0;
            y[i+3*n_comp] = 0;
        }
        for(int i=v_comps;i<v_comps+h_comps;++i){
            y[i]= fh_per_comp - (fh_per_comp*i0);
            y[i+n_comp] = fh_per_comp*i0;
            y[i+2*n_comp] = 0;
            y[i+3*n_comp] = 0;
        }
        for(int i=v_comps+h_comps; i<n_comp; ++i){
            y[i]= fr_per_comp - (fr_per_comp*i0);
            y[i+n_comp] = fr_per_comp*i0;
            y[i+2*n_comp] = 0;
            y[i+3*n_comp] = 0;
        }
        runge_kutta4< state_type > stepper;
        vector<double> times;

        SIRS_3Comp vuln(n_comp,gamma,omega,sd1, sd2, sd3, sd4,beta_ph0, beta_ph1,beta_ph2,beta_ph3,fj);
        vector<state_type> y_vec;
        integrate_n_steps(stepper,vuln, y, 0.0 , 1.0, simtime, push_back_state_time_pIv_pIh_pIr(y_vec, times,fj,v_comps, h_comps,r_comps));
        //Calc height first peak and second peak
        vector<double> tempIv;
        int end = y_vec[0].size();
        for (size_t i=0; i<y_vec.size();++i){
            tempIv.push_back(y_vec[i][end-3]);
        }
        double  incrIv = 0;
        int sd_2nd_peak_Iv = static_cast<int>(sd1)+300;
        for (size_t i=72; i<sd1+365;++i){
            if (tempIv[i]>=tempIv[i-1]){
                incrIv=1;
                sd_2nd_peak_Iv = i;
                break;
            }
        }
        auto maxIv = *max_element(std::begin(tempIv)+sd_2nd_peak_Iv, std::begin(tempIv)+(sd1+365));
        auto rel_height = maxIv/tempIv[71];
        results.push_back({ph0,rel_height});
    }
    cout << results;
}

void vary_Re_P2(){
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

    int n_comp=5;
    //2-2-96 -> 1v-1v-48r
    //8-8-84 -> 4v-4h-42r
    //14-14-72 -> 7v-7h-36r
    //20-20-60 -> 10v-10h-30r or 1v-1h-3r
    //40-20-40 -> 2v-1h-2r
    //10-20-70 -> 1v-2h-7r
    vector<double> fj = {0.2,0.2,0.6};
    int v_comps = 1;
    int h_comps = 1;
    int r_comps = n_comp - v_comps - h_comps;
    double fv_per_comp = fj[0]/v_comps;
    double fh_per_comp = fj[1]/h_comps;
    double fr_per_comp = fj[2]/r_comps;
    data_type beta_ph0(n_comp, vector<double>(n_comp,0)),beta_ph1(n_comp, vector<double>(n_comp,0)),beta_ph2(n_comp, vector<double>(n_comp,0)),beta_ph3(n_comp, vector<double>(n_comp,0));
    omega = 1.0/365.0;

    simtime = 365*2;

    sd1 = 71.0;
    sd2 = sd1+(6*7);
    sd3 = sd2+(12*7);
    sd4 = simtime;
    data_type results;
    double i0=0.0001;
    //Phase2 beta1&4=0.6, beta2&3=0.7 or beta1&4=1.0, beta2&3=1.1
    beta1_ph1 = beta4_ph1=1.0*gamma;
    beta2_ph1 = beta3_ph1=1.1*gamma;
    for (int i=0; i<n_comp;++i){
        for (int j=0; j<n_comp;++j){
            if (i<v_comps){
                if (j<v_comps+h_comps){
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
            } else if (i>=v_comps & i<v_comps+h_comps){
                if (j<v_comps){
                    beta_ph0[i][j] = beta1_ph0;
                    beta_ph1[i][j] = beta1_ph1;
                    beta_ph2[i][j] = beta1_ph2;
                    beta_ph3[i][j] = beta1_ph3;
                } else if (j>=v_comps & j<v_comps + h_comps){
                    beta_ph0[i][j] = beta1_ph0;
                    beta_ph1[i][j] = beta1_ph1;
                    beta_ph2[i][j] = beta1_ph2;
                    beta_ph3[i][j] = beta1_ph3;
                } else {
                    beta_ph0[i][j] = beta2_ph0;
                    beta_ph1[i][j] = beta2_ph1;
                    beta_ph2[i][j] = beta2_ph2;
                    beta_ph3[i][j] = beta2_ph3;
                }
            } else if (i>=v_comps + h_comps){
                if (j<v_comps){
                    beta_ph0[i][j] = beta4_ph0;
                    beta_ph1[i][j] = beta4_ph1;
                    beta_ph2[i][j] = beta4_ph2;
                    beta_ph3[i][j] = beta4_ph3;
                } else if (j<v_comps + h_comps){
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

    state_type y(4*n_comp);

    //Initialise starting values
    for(int i=0;i<v_comps;++i){
        y[i]= fv_per_comp - (fv_per_comp*i0);
        y[i+n_comp] = fv_per_comp*i0;
        y[i+2*n_comp] = 0;
        y[i+3*n_comp] = 0;
    }
    for(int i=v_comps;i<v_comps+h_comps;++i){
        y[i]= fh_per_comp - (fh_per_comp*i0);
        y[i+n_comp] = fh_per_comp*i0;
        y[i+2*n_comp] = 0;
        y[i+3*n_comp] = 0;
    }
    for(int i=v_comps+h_comps; i<n_comp; ++i){
        y[i]= fr_per_comp - (fr_per_comp*i0);
        y[i+n_comp] = fr_per_comp*i0;
        y[i+2*n_comp] = 0;
        y[i+3*n_comp] = 0;
    }
    runge_kutta4< state_type > stepper;
    vector<double> times;

    SIRS_3Comp vuln(n_comp,gamma,omega,sd1, sd2, sd3, sd4,beta_ph0, beta_ph1,beta_ph2,beta_ph3,fj);
    vector<state_type> y_vec;
    integrate_n_steps(stepper,vuln, y, 0.0 , 1.0, simtime, push_back_state_time_pIv_pIh_pIr(y_vec, times,fj,v_comps, h_comps,r_comps));
    vector<vector<double>> betas;

    write_simulation("SIRS_MComp_20-20-60-P2_10.csv", y_vec, betas, v_comps,h_comps, r_comps);
}


void immunity(){
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

    int n_comp=5;
    //2-2-96 -> 1v-1v-48r
    //8-8-84 -> 4v-4h-42r
    //14-14-72 -> 7v-7h-36r
    //20-20-60 -> 10v-10h-30r or 1v-1h-3r
    //20-40-40 -> 1v-2h-2r
    //20-10-70 -> 2v-1h-7r
    vector<double> fj = {0.2,0.2,0.6};
    int v_comps = 1;
    int h_comps = 1;
    int r_comps = n_comp - v_comps - h_comps;
    double fv_per_comp = fj[0]/v_comps;
    double fh_per_comp = fj[1]/h_comps;
    double fr_per_comp = fj[2]/r_comps;
    data_type beta_ph0(n_comp, vector<double>(n_comp,0)),beta_ph1(n_comp, vector<double>(n_comp,0)),beta_ph2(n_comp, vector<double>(n_comp,0)),beta_ph3(n_comp, vector<double>(n_comp,0));

    omega = 1.0/365.0;

    simtime = 365*2;

    sd1 = 71.0;
    sd2 = sd1+(6*7);
    sd3 = sd2+(12*7);
    sd4 = simtime;
    data_type results;
    double i0=0.0001;
    for (int i=0; i<n_comp;++i){
        for (int j=0; j<n_comp;++j){
            if (i<v_comps){
                if (j<v_comps+h_comps){
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
            } else if (i>=v_comps & i<v_comps+h_comps){
                if (j<v_comps){
                    beta_ph0[i][j] = beta1_ph0;
                    beta_ph1[i][j] = beta1_ph1;
                    beta_ph2[i][j] = beta1_ph2;
                    beta_ph3[i][j] = beta1_ph3;
                } else if (j>=v_comps & j<v_comps + h_comps){
                    beta_ph0[i][j] = beta1_ph0;
                    beta_ph1[i][j] = beta1_ph1;
                    beta_ph2[i][j] = beta1_ph2;
                    beta_ph3[i][j] = beta1_ph3;
                } else {
                    beta_ph0[i][j] = beta2_ph0;
                    beta_ph1[i][j] = beta2_ph1;
                    beta_ph2[i][j] = beta2_ph2;
                    beta_ph3[i][j] = beta2_ph3;
                }
            } else if (i>=v_comps + h_comps){
                if (j<v_comps){
                    beta_ph0[i][j] = beta4_ph0;
                    beta_ph1[i][j] = beta4_ph1;
                    beta_ph2[i][j] = beta4_ph2;
                    beta_ph3[i][j] = beta4_ph3;
                } else if (j<v_comps + h_comps){
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


    for (double dur = 1.0; dur<365.1; dur+=1.0){
        state_type y(4*n_comp);

        //Initialise starting values
        for(int i=0;i<v_comps;++i){
            y[i]= fv_per_comp - (fv_per_comp*i0);
            y[i+n_comp] = fv_per_comp*i0;
            y[i+2*n_comp] = 0;
            y[i+3*n_comp] = 0;
        }
        for(int i=v_comps;i<v_comps+h_comps;++i){
            y[i]= fh_per_comp - (fh_per_comp*i0);
            y[i+n_comp] = fh_per_comp*i0;
            y[i+2*n_comp] = 0;
            y[i+3*n_comp] = 0;
        }
        for(int i=v_comps+h_comps; i<n_comp; ++i){
            y[i]= fr_per_comp - (fr_per_comp*i0);
            y[i+n_comp] = fr_per_comp*i0;
            y[i+2*n_comp] = 0;
            y[i+3*n_comp] = 0;
        }
        runge_kutta4< state_type > stepper;

        vector<double> times;
        double omega1=1.0/dur;

        SIRS_3Comp vuln(n_comp,gamma,omega1,sd1, sd2, sd3, sd4,beta_ph0, beta_ph1,beta_ph2,beta_ph3,fj);
        vector<state_type> y_vec;
        integrate_n_steps(stepper,vuln, y, 0.0 , 1.0, simtime, push_back_state_time_pIv_pIh_pIr(y_vec, times,fj,v_comps, h_comps,r_comps));
        //Calc height first peak and second peak
        vector<double> tempIv;
        int end = y_vec[0].size();
        for (size_t i=0; i<y_vec.size();++i){
            tempIv.push_back(y_vec[i][end-3]);
        }
        double  incrIv = 0;
        int sd_2nd_peak_Iv = static_cast<int>(sd1)+360;
        for (size_t i=72; i<sd1+365;++i){
            if (tempIv[i]>=tempIv[i-1]){
                incrIv=1;
                sd_2nd_peak_Iv = i;
                break;
            }
        }
        auto maxIv = *max_element(std::begin(tempIv)+sd_2nd_peak_Iv, std::begin(tempIv)+(sd1+365));
        auto rel_height = maxIv/tempIv[71];
        results.push_back({dur,rel_height});
    }
    cout << results;
}

void test_shielders(){
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

    int n_comp=5;
    //2-2-96 -> 1v-1v-48r
    //8-8-84 -> 4v-4h-42r
    //14-14-72 -> 7v-7h-36r
    //20-20-60 -> 10v-10h-30r or 1v-1h-3r
    //20-40-40 -> 1v-2h-2r
    //20-10-70 -> 2v-1h-7r
    vector<double> fj = {0.2,0.2,0.6};
    int v_comps = 1;
    int h_comps = 1;
    int r_comps = n_comp - v_comps - h_comps;
    double fv_per_comp = fj[0]/v_comps;
    double fh_per_comp = fj[1]/h_comps;
    double fr_per_comp = fj[2]/r_comps;
    data_type beta_ph0(n_comp, vector<double>(n_comp,0)),beta_ph1(n_comp, vector<double>(n_comp,0)),beta_ph2(n_comp, vector<double>(n_comp,0)),beta_ph3(n_comp, vector<double>(n_comp,0));
    omega = 1.0/365.0;

    simtime = 365*2;

    sd1 = 71.0;
    sd2 = sd1+(6*7);
    sd3 = sd2+(12*7);
    sd4 = simtime;
    data_type results;
    double i0=0.0001;
    for (double perc=0.0; perc<0.51; perc+=0.5){
        //beta_BB, EE & HH are affected
        beta_ph0 = {
            {beta1_ph0,beta1_ph0,beta4_ph0,beta4_ph0,beta4_ph0},
            {beta1_ph0,beta1_ph0,beta2_ph0,beta2_ph0,beta2_ph0},
            {beta4_ph0,beta2_ph0,beta3_ph0,beta3_ph0,beta3_ph0},
            {beta4_ph0,beta2_ph0,beta3_ph0,beta3_ph0,beta3_ph0},
            {beta4_ph0,beta2_ph0,beta3_ph0,beta3_ph0,beta3_ph0}
        };
        beta_ph1 = {
            {beta1_ph1,beta1_ph1,beta4_ph1,beta4_ph1,beta4_ph1},
            {beta1_ph1,beta1_ph1,beta2_ph1,beta2_ph1,beta2_ph1},
            {beta4_ph1,beta2_ph1,beta3_ph1,beta3_ph1,beta3_ph1},
            {beta4_ph1,beta2_ph1,beta3_ph1,beta3_ph1,beta3_ph1},
            {beta4_ph1,beta2_ph1,beta3_ph1,beta3_ph1,beta3_ph1}
        };
        beta_ph2 = {
            {beta1_ph2,beta1_ph2,beta4_ph2,beta4_ph2,beta4_ph2},
            {beta1_ph2,beta1_ph2,beta2_ph2,beta2_ph2,beta2_ph2},
            {beta4_ph2,beta2_ph2,beta3_ph2,beta3_ph2,beta3_ph2},
            {beta4_ph2,beta2_ph2,beta3_ph2,beta3_ph2,beta3_ph2},
            {beta4_ph2,beta2_ph2,beta3_ph2,beta3_ph2,beta3_ph2}
       };
        beta_ph3 = beta_ph2;
        state_type y(4*n_comp);

        //Initialise starting values
        for(int i=0;i<v_comps;++i){
            y[i]= fv_per_comp - (fv_per_comp*i0);
            y[i+n_comp] = fv_per_comp*i0;
            y[i+2*n_comp] = 0;
            y[i+3*n_comp] = 0;
        }
        for(int i=v_comps;i<v_comps+h_comps;++i){
            y[i]= fh_per_comp - (fh_per_comp*i0);
            y[i+n_comp] = fh_per_comp*i0;
            y[i+2*n_comp] = 0;
            y[i+3*n_comp] = 0;
        }
        for(int i=v_comps+h_comps; i<n_comp; ++i){
            y[i]= fr_per_comp - (fr_per_comp*i0);
            y[i+n_comp] = fr_per_comp*i0;
            y[i+2*n_comp] = 0;
            y[i+3*n_comp] = 0;
        }
        runge_kutta4< state_type > stepper;
        vector<double> times;

        SIRS_3Comp vuln(n_comp,gamma,omega,sd1, sd2, sd3, sd4,beta_ph0, beta_ph1,beta_ph2,beta_ph3,fj);
        vector<state_type> y_vec;
        integrate_n_steps(stepper,vuln, y, 0.0 , 1.0, simtime, push_back_state_time_pIv_pIh_pIr(y_vec, times,fj,v_comps, h_comps,r_comps));
        //Calc height first peak and second peak
        vector<double> tempIv;
        int end = y_vec[0].size();
        for (size_t i=0; i<y_vec.size();++i){
            tempIv.push_back(y_vec[i][end-3]);
        }
        double  incrIv = 0;
        int sd_2nd_peak_Iv = static_cast<int>(sd1)+300;
        for (size_t i=72; i<sd1+365;++i){
            if (tempIv[i]>=tempIv[i-1]){
                incrIv=1;
                sd_2nd_peak_Iv = i;
                break;
            }
        }
        auto maxIv = *max_element(std::begin(tempIv)+sd_2nd_peak_Iv, std::begin(tempIv)+(sd1+365));
        auto rel_height = maxIv/tempIv[71];
        results.push_back({perc,rel_height});
    }
    cout << results;
}

void write_simulation(string filename, const data_type &y_vec, const data_type &betas, int v_comps, int h_comps, int r_comps){

    ofstream output(filename);
    output << "t,";
    for (int i=0; i < v_comps;++i){
        output << "Sv" + to_string(i) + ",";
    }
    for (int i=0; i < h_comps;++i){
        output << "Ss" + to_string(i) + ",";
    }
    for (int i=0; i < r_comps;++i){
        output << "Sg" + to_string(i) + ",";
    }
    for (int i=0; i < v_comps;++i){
        output << "Iv" + to_string(i) + ",";
    }
    for (int i=0; i < h_comps;++i){
        output << "Is" + to_string(i) + ",";
    }
    for (int i=0; i < r_comps;++i){
        output << "Ig" + to_string(i) + ",";
    }
    for (int i=0; i < v_comps;++i){
        output << "Rv" + to_string(i) + ",";
    }
    for (int i=0; i < h_comps;++i){
        output << "Rs" + to_string(i) + ",";
    }
    for (int i=0; i < r_comps;++i){
        output << "Rg" + to_string(i) + ",";
    }
    for (int i=0; i < v_comps;++i){
        output << "cumIv" + to_string(i) + ",";
    }
    for (int i=0; i < h_comps;++i){
        output << "cumIs" + to_string(i) + ",";
    }
    for (int i=0; i < r_comps;++i){
        output << "cumIg" + to_string(i) + ",";
    }
    output << "cumIv, pIv,pIs,pIg\n";
//    for (size_t i=0; i<y_vec.size();++i){
//        output << y_vec[i] << "," << betas[i] << "\n";
//    }
    output << y_vec;
    output.close();
}

void write_phase_matrices(const data_type &beta_ph0,const data_type &beta_ph1,const data_type &beta_ph2,const data_type &beta_ph3){
    string filename = "SIRS_MComp_simulation_betas.csv";
    ofstream output(filename);
    output << "v,h,r\n";
    output << beta_ph0 << "\n";
    output << beta_ph1 << "\n";
    output << beta_ph2 << "\n";
    output << beta_ph3 << "\n";
    output.close();
}


void write_fast_results(const data_type &fast_results){
    string filename = "SIRS_mult_Comp_FAST_results.csv";
    ofstream output(filename);
    output << "peak_Iv, cum_Iv, higher\n";
    output << fast_results;
    output.close();
}

void write_optimisations(data_type optim_results){
    string filename = "SIRS_mult_Comp_optim_results.csv";
    ofstream output(filename);
    output << "beta1,beta2,beta3,beta4,cumIv_1y,first_peakIv,sec_peakIv,first_peakIh,sec_peakIh,first_peakIr,sec_peakIr,incrIv,incrIh,incrIr\n";
    output << optim_results;
}
