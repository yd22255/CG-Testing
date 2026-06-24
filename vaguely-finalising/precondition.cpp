#include <cstdio>
#include <string>
#include <thread>
#include <array>
#include <vector>
#include <numeric> 
#include <functional>
#include <iostream>
#include <algorithm>
#include <random>
#include <mdspan>
#include <ranges>
#include <ctime>
#include <execution>
//#include <exec/static_thread_pool.hpp>
#include "libraryfuncs.hpp"
#include "matrix_gen.hpp"
#include "precondition.hpp"

namespace stdexec=std::execution;
using namespace std::this_thread;
using namespace mylibs;
using namespace mymats;
using namespace precon;

//just going to throw ideas out here and see where i get
    //i think i need to pass number of rows/number of columns in each of these because it would make sqcoo_spgemm just coo_spgemm

std::pair<std::pair<std::vector<double>,std::vector<double>>,std::vector<double>> precon::ilu_preconditionMatrix(std::vector<double> Arow,std::vector<double> Acolumn,std::vector<double> Adata,int Anumrows){
    auto [LUrowcols,LUdata] = coo_ilu(Arow,Acolumn,Adata);
    auto[LUArowcols,LUAdata] = coo_spgemm(LUrowcols.first,LUrowcols.second,LUdata,Anumrows,Anumrows,Arow,Acolumn,Adata,Anumrows,Anumrows);
    std::pair<std::pair<std::vector<double>,std::vector<double>>,std::vector<double>> toreturn = std::make_pair(LUArowcols,LUAdata);
    return toreturn;
}
std::vector<double> precon::ilu_preconditionVector(std::vector<double> Arow,std::vector<double> Acolumn,std::vector<double> Adata,std::vector<double>b){
    auto [LUrowcols,LUdata] = coo_ilu(Arow,Acolumn,Adata);
    std::vector<double> LU_b = sparsemvdot_product(LUrowcols.first,LUrowcols.second,LUdata,b);
    return LU_b;
}

std::pair<std::pair<std::vector<double>,std::vector<double>>,std::vector<double>> precon::jacobi_preconditionMatrix(std::vector<double> Arow,std::vector<double> Acolumn,std::vector<double> Adata, int Anumrows, int Anumcols){
    auto [Jrowcols,Jdata] = coo_jacobi(Arow,Acolumn,Adata);
    auto[JArowcols,JAdata] = coo_spgemm(Jrowcols.first,Jrowcols.second,Jdata,Anumrows,Anumcols,Arow,Acolumn,Adata,Anumrows,Anumcols);
    std::pair<std::pair<std::vector<double>,std::vector<double>>,std::vector<double>> toreturn = std::make_pair(JArowcols,JAdata);
    return toreturn;
}
std::vector<double> precon::jacobi_preconditionVector(std::vector<double> Arow,std::vector<double> Acolumn,std::vector<double> Adata,std::vector<double>b){
    auto [Jrowcols,Jdata] = coo_jacobi(Arow,Acolumn,Adata);
    std::vector<double> J_b = sparsemvdot_product(Jrowcols.first,Jrowcols.second,Jdata,b);
    return J_b;
}