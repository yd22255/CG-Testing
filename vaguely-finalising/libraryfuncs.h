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


std::vector<double> sparsemvdot_product(std::vector<double> Arow, std::vector<double> Acolumn,std::vector<double>Adata, const std::span<double> vector){}

std::vector<double> vectortranspose(const std::vector<double> vector){} 

double vvdot_product(const std::vector<double> v1, const std::vector<double> v2){}

std::pair<std::pair<std::vector<double>,std::vector<double>>,std::vector<double>> matrixtranspose(std::vector<double> Mrow,std::vector<double> Mcolumn,std::vector<double> Mdata){}

std::vector<double> densematrixtranspose(std::vector<double> matrix){}

std::vector<double> vectorize(std::vector<double>Arow,std::vector<double>Acolumn,std::vector<double>Adata){}

double norm(std::vector<double> vector){}

std::pair<std::pair<std::vector<double>,std::vector<double>>,std::vector<double>> cootocsr(std::vector<double> Arow,std::vector<double> Acolumn,std::vector<double> Adata,int numrows){}

int csrvalueaccess(std::vector<double> Ainds, std::vector<double>Adata,int x, int y){}

std::pair<std::pair<std::vector<double>,std::vector<double>>,std::vector<double>> coo_spgemm(std::vector<double> Arow,std::vector<double> Acolumn,std::vector<double> Adata,int Anumrows,int Anumcols,std::vector<double> Brow,std::vector<double> Bcolumn,std::vector<double> Bdata,int Bnumrows,int Bnumcols){}

std::pair<std::pair<std::vector<double>,std::vector<double>>,std::vector<double>> coo_jacobi(std::vector<double> &Arow, std::vector<double>&Acolumn,std::vector<double>&Adata){}

std::pair<std::pair<std::vector<double>,std::vector<double>>,std::vector<double>> coo_ilu(std::vector<double> &Arow, std::vector<double>&Acolumn,std::vector<double>&Adata){}

int coo(std::vector<double> A,std::vector<double> &Arow, std::vector<double>&Acolumn,std::vector<double>&Adata,int rows,int cols){} 

bool roundtoepsilon(const double x, double epsilon){}

