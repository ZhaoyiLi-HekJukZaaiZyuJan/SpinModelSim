#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <random>
#include <gmp.h>
#include <stdio.h>

#include <cxxopts/cxxopts.hpp>

using namespace std;

clock_t start_t = clock();
auto start = chrono::steady_clock::now();

struct lattice {
        vector<vector<int> > spins;
        vector<vector<vector<int> > > bonds;
        int S[2];
        double T;
        double J;
        int energy;

        lattice(const int S[2], const  double T, const double J) {
            this->S[0] = S[0];
            this->S[1] = S[1];
            this->T = T;
            this->J = J;
            vector<vector<int> > spins(S[0], vector<int>(S[1], 1));
            this->spins = spins;
            vector<vector<vector<int> > > bonds(2, vector<vector<int> >(S[0], vector<int>(S[1], J)));
            this->bonds = bonds;

            this->energy = update_energy();
        }

        int update_energy() {
            int energy = 0;
            for (int i = 0; i < S[0]; i++) {
                for (int j = 0; j < S[1]; j++) {
                    energy += bonds[0][i][j] * spins[i][j] * spins[(1 + i) % S[0]][j];
                    energy += bonds[1][i][j] * spins[i][j] * spins[i][(1 + j) % S[1]];
                }
            }
            this->energy = energy;
            return energy;
        }

        int get_state_number() {
            int conversionMat[S[0] * S[1]];
            int powerOf2 = 1;
            for (int i = S[0] * S[1] - 1; i >= 0; i--) {
                conversionMat[i] = powerOf2;
                powerOf2 *= 2;
            }

            int stateNumber = 0;
            for (int i = 0; i < S[0]; i++) {
                for (int j = 0; j < S[1]; j++) {
                    stateNumber += conversionMat[i * S[0] + j] * ((spins[i][j] + 1) / 2);
                }
            }

            return stateNumber;
        }

        void random_initial_condition(int seed = 0) {
            mt19937 rng(seed);
            uniform_real_distribution<double> distribution(0.0, 1.0);

            for (int i = 0; i < S[0]; i++) {
                for (int j = 0; j < S[1]; j++) {
                    if (distribution(rng) < 0.5) {
                        spins[i][j] = -1;
                    } else {
                        spins[i][j] = 1;
                    }
                }
            }
        }

        void print() {
            for (int i = 0; i < S[0]; i++) {
                for (int j = 0; j < S[1]; j++) {
                    cout << spins[i][j] << " ";
                }
                cout << endl;
            }
        }

        void make_step(mt19937& rng) {
            
            uniform_int_distribution<int> distribution(0, S[0] - 1);

            int x = distribution(rng);
            int y = distribution(rng);

            int d_energy = 0;
            d_energy += bonds[0][x][y] * spins[x][y] * spins[(x + 1) % S[0]][y];
            d_energy += bonds[0][(x - 1 + S[0]) % S[0]][y] * spins[x][y] * spins[(x - 1 + S[0]) % S[0]][y];
            d_energy += bonds[1][x][y] * spins[x][y] * spins[x][(y + 1) % S[1]];
            d_energy += bonds[1][x][(y - 1 + S[1]) % S[1]] * spins[x][y] * spins[x][(y - 1 + S[1]) % S[1]];
            d_energy *= -2;

            double r = generate_canonical<double, 10>(rng);
            if (r < exp(-d_energy / T)) {
                spins[x][y] *= -1;
                energy += d_energy;
            }
        }
        
        void make_sweep(int nsteps, mt19937& rng) {
            for(int i=0; i <  nsteps; i++){
                this->make_step(rng);
            }
        }
};

void simulate(const int S[2], const double T, const double J, const int waitSweeps, const int rptSweeps, const bool out, const bool track_state, const string fname){
    random_device rd;
    mt19937 rng(rd());
    ofstream outfile_m;
    ofstream outfile_e;
    ofstream outfile_s;
    ofstream outfile_stats_m;
    ofstream outfile_stats_e;

    long long int sum_e = 0, sum_m = 0;
    unsigned long long int sum_e2 = 0, sum_m2 = 0;


	if(out != 0){
		outfile_m.open(fname +"_m.out", fstream::app);
        outfile_e.open(fname +"_e.out", fstream::app);
        outfile_stats_m.open(fname +"_stats_m.out", fstream::app);
        outfile_stats_e.open(fname +"_stats_e.out", fstream::app);
		
        outfile_m << T << "," << flush;
        outfile_e << T << "," << flush;
        outfile_stats_m << T << "," << flush;
        outfile_stats_e << T << "," << flush;
        if(track_state != 0){
             outfile_s.open(fname +"_s.out", fstream::app);
             outfile_s << T << "," << flush;
        }
	}
    
    lattice testlattice(S, T, J);
    testlattice.random_initial_condition();
    testlattice.energy = testlattice.update_energy();


    // Wait for waitSweeps
    for (int i = 0; i < waitSweeps; i++) {
        testlattice.make_sweep(S[0]*S[1], rng);
    }

    // Repeat for rptSweeps and sample for each sweep
    for (int i = 0; i < rptSweeps; i++) {
        int energy = 0, state = 0, magnetization = 0; 
        if (i % (rptSweeps / 100) == 0) {
            cout << ".";
            cout.flush();
        }
        testlattice.make_sweep(S[0]*S[1], rng);

        // Recording state, energy, and magnetization
        energy = testlattice.energy;
        magnetization = 0;
        for (int j = 0; j < S[0]; j++) {
            for (int k = 0; k < S[1]; k++) {
                magnetization += testlattice.spins[j][k];
            }  
        }
        if (track_state){
              state = testlattice.get_state_number();
        }
        sum_e += energy;
        sum_e2 += energy*energy;
        sum_m += magnetization;
        sum_m2 += magnetization*magnetization;

        if(out != 0){
            outfile_m << magnetization << "," << flush;
            outfile_e << energy << ","  << flush;
            if(track_state != 0){
                outfile_s << state << "," << flush;
            }
        }
    }
    if(out != 0){
        outfile_m << endl;
        outfile_e << endl;
        if(track_state != 0){
            outfile_s << endl;
        }
        outfile_stats_m << sum_m << "," << sum_m2 << endl;
        outfile_stats_e << sum_e << "," << sum_e2 << endl;
    }
    outfile_m.close();
    outfile_s.close();
    outfile_e.close();
}

int main(int argc, const char *argv[]) {
    start_t = clock();
    string fname;
    int n, waitSweep, rptSweep, nt;
    double tmin, tmax, J;
    bool out, track_state;

    //test
    cout << "test" << endl;
    int test_int = 1000;
    cout << test_int << endl;
    mpz_t long_int;
    mpz_set_si(long_int,test_int);
    mpz_pow_ui(long_int, long_int,10) ;
    cout << long_int << endl;


    cxxopts::Options options(*argv,
							 "Simulator for 2D ising model"
							 );

    options.add_options()
	("fname", "filename", cxxopts::value(fname)->default_value(""))
    ("track_state", "keep track of the state as an int", cxxopts::value(track_state)->default_value("0"))
    ("out", "output in this directory as an .out file", cxxopts::value(out)->default_value("0"))
    ("waitSweep", "number of sweeps to wit before sampling", cxxopts::value(waitSweep)->default_value("10"))
    ("rptSweep", "number of sweeps sampled", cxxopts::value(rptSweep)->default_value("10000"))
	("J", "interaction strength", cxxopts::value(J)->default_value("-1"))
    ("tmin", "temperature", cxxopts::value(tmin)->default_value("0.1"))
    ("tmax", "temperature", cxxopts::value(tmax)->default_value("2.1"))
    ("nt", "temperature", cxxopts::value(nt)->default_value("10"))
	("n", "size of lattice", cxxopts::value(n)->default_value("3"));
    options.parse(argc, argv);

    //outputing options
	cout << "tmin:" << tmin << ";tmax:" << tmax << endl;
	cout << "n:" << n << endl;
    cout << "waitSweep:" << waitSweep << ", rptSweep:" << rptSweep <<  endl;

    int binst = 0 == nt ? 1: nt;

    if (fname == "") {
		fname = "../data/n=" + to_string(n) + ",T=(" + to_string(tmin).substr(0,3)+ "," +to_string(tmax).substr(0,3) + "),J=" +to_string(J).substr(0,4) + ",rptSweep=" + to_string(rptSweep);
	}
    cout << "fname:" << fname << endl;
    int S[2] = {n,n};
    
	for (int j = 0; j <= nt; j++) {
        double t = tmin + (tmax-tmin)/binst * j;
        simulate(S, t, J, waitSweep, rptSweep, out, track_state,  fname);
     }   

    cout << "t(Generation):"<< double(clock()-start_t)/CLOCKS_PER_SEC << endl; //timing
    return 0;
}

//simple testing commands
//./simulate -n 20 --out --tmin 0.1./simulate -n 20 --out --tmin 0.1 --tmax 3.1 --nt 10