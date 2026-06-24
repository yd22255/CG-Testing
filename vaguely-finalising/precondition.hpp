#ifndef PRECON
#define PRECON
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

namespace precon{
    std::pair<std::pair<std::vector<double>,std::vector<double>>,std::vector<double>> ilu_preconditionMatrix(std::vector<double> Arow,std::vector<double> Acolumn,std::vector<double> Adata,int Anumrows);
    std::vector<double> ilu_preconditionVector(std::vector<double> Arow,std::vector<double> Acolumn,std::vector<double> Adata,std::vector<double>b);
    std::pair<std::pair<std::vector<double>,std::vector<double>>,std::vector<double>> jacobi_preconditionMatrix(std::vector<double> Arow,std::vector<double> Acolumn,std::vector<double> Adata, int Anumrows, int Anumcols);
    std::vector<double> jacobi_preconditionVector(std::vector<double> Arow,std::vector<double> Acolumn,std::vector<double> Adata,std::vector<double>b);
}   

#endif