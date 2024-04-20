#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <time.h>
#include <string>
#include <iomanip>
#include <map>
#include <algorithm>
#include <mpi.h>
#include <vector>
#include "calc_c.h"
#include "calc_energy.h"

const double et0{0.3 * pow(2.25 * M_PI, 2. / 3.)};
const double ex0{-0.75 * pow(2.25 * M_PI, 1. / 3.) / M_PI};
double ls;
double le;
int Nl;
double dl;
double zs;
double ze;
int Nz;
double dz;
int Njob;
const std::string result_filename{"results.csv"};
bool out_err;

void myhandler(const char *reason, const char *file, int line, int gsl_errno)
{
    if (out_err)
    {
        std::cout << "===================================================" << std::endl;
        std::cout << reason << std::endl;
    }
}

void parser(std::string infile)
{
    std::map<std::string, std::string> indata;

    std::ifstream ifs(infile);
    if (ifs.fail())
        throw std::runtime_error("Failed to open file.");
    std::string str;
    while (std::getline(ifs, str))
    {
        std::string name, val;
        std::istringstream iss(str);
        iss >> name >> val;
        indata[name] = val;
    }
    ifs.close();

    auto indata_check = [&](std::string name)
    {
        try
        {
            return indata.at(name);
        }
        catch (...)
        {
            throw std::runtime_error(name + " is not defined.");
        }
    };

    out_err = (std::stoi(indata_check("out_err")) == 1);

    ls = std::stod(indata_check("ls"));
    le = std::stod(indata_check("le"));
    Nl = std::stoi(indata_check("Nl"));
    dl = (le - ls) / (double)(Nl + (int)(Nl == 0));

    zs = std::stod(indata_check("zs"));
    ze = std::stod(indata_check("ze"));
    Nz = std::stoi(indata_check("Nz"));
    dz = (ze - zs) / (double)(Nz + (int)(Nz == 0));

    Njob = (Nl + 1) * (Nz + 1);
}

void write_out_results(std::vector<std::pair<int, double>> result_vec,
                       std::vector<std::pair<int, double>> result_gb_vec)
{
    auto calc_Etot{[](double rs, double zeta, double Ecorr)
                   {
                       double fu{pow(1. + zeta, 1. / 3.)};
                       double fd{pow(1. - zeta, 1. / 3.)};
                       double Et{0.5 * (gsl_pow_5(fu) + gsl_pow_5(fd)) * et0 / (rs * rs)};
                       double Ex{0.5 * (gsl_pow_4(fu) + gsl_pow_4(fd)) * ex0 / rs};
                       return Et + Ex + Ecorr;
                   }};

    std::ofstream ofs(result_filename, std::ofstream::app);

    for (int i{0}; i < Njob; i++)
    {
        const int iz{i % (Nz + 1)};
        const int il{i / (Nz + 1)};
        double l{ls + il * dl};
        double rs{pow(10., l)};
        const double zeta{zs + (double)iz * dz};

        const double result{result_vec[i].second};
        const double result_gb{result_gb_vec[i].second};

        ofs << std::showpos << std::scientific << std::setprecision(18) << rs;
        ofs << ",";
        ofs << std::showpos << std::scientific << std::setprecision(18) << zeta;
        ofs << ",";
        ofs << std::showpos << std::scientific << std::setprecision(18) << result;
        ofs << ",";
        ofs << std::showpos << std::scientific << std::setprecision(18) << result_gb;
        ofs << ",";
        ofs << std::showpos << std::scientific << std::setprecision(18) << calc_Etot(rs, zeta, result);
        ofs << ",";
        ofs << std::showpos << std::scientific << std::setprecision(18) << calc_Etot(rs, zeta, result_gb);
        ofs << "\n";
    }
    ofs.close();
}

void main_MPI(int argc, char *argv[])
{
    int rank, proc;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &proc);

    const int Nloop{(Njob + proc - 1) / proc};
    int *ijob_local{new int[Nloop]};
    for (int i{0}; i < Nloop; i++)
        ijob_local[i] = -1;
    double *result_local{new double[Nloop]};
    double *result_gb_local{new double[Nloop]};

    for (int i{0}; i < Nloop || rank + i * proc < Njob; i++)
    {
        const int ijob{rank + i * proc};
        const int iz{ijob % (Nz + 1)};
        const int il{ijob / (Nz + 1)};
        double l{ls + il * dl};
        double rs{pow(10., l)};
        const double zeta{zs + (double)iz * dz};

        ijob_local[i] = ijob;
        result_local[i] = calc_energy(rs, zeta);
        result_gb_local[i] = calc_energy_gb(rs, zeta);

        std::cout << "Job i=" << ijob << " (rs=" << rs
                  << ", zeta=" << zeta << ") done" << std::endl;
    }

    int *ijob_global{new int[Nloop * proc]};
    double *result_global{new double[Nloop * proc]};
    double *result_gb_global{new double[Nloop * proc]};
    MPI_Allgather(ijob_local, Nloop, MPI_INT, ijob_global, Nloop, MPI_INT, MPI_COMM_WORLD);
    MPI_Allgather(result_local, Nloop, MPI_DOUBLE, result_global, Nloop, MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Allgather(result_gb_local, Nloop, MPI_DOUBLE, result_gb_global, Nloop, MPI_DOUBLE, MPI_COMM_WORLD);

    if (rank == 0)
    {
        std::vector<std::pair<int, double>> result_vec;
        std::vector<std::pair<int, double>> result_gb_vec;

        for (int i{0}; i < Nloop * proc; i++)
        {
            const int ijob{ijob_global[i]};
            if (ijob >= 0)
            {
                result_vec.push_back({ijob, result_global[i]});
                result_gb_vec.push_back({ijob, result_gb_global[i]});
            }
        }

        std::sort(result_vec.begin(), result_vec.end());
        std::sort(result_gb_vec.begin(), result_gb_vec.end());
        write_out_results(result_vec, result_gb_vec);
    }

    delete[] ijob_local;
    delete[] ijob_global;
    delete[] result_local;
    delete[] result_global;
    delete[] result_gb_local;
    delete[] result_gb_global;

    MPI_Finalize();
}

int main(int argc, char *argv[])
{
    gsl_set_error_handler(&myhandler);

    try
    {
        if (argc == 1)
            throw std::runtime_error("input file not given.");
        parser(argv[1]);
        main_MPI(argc, argv);
    }
    catch (std::runtime_error &e)
    {
        std::cerr << "runtime_error: " << e.what() << std::endl;
        return -1;
    }
}
