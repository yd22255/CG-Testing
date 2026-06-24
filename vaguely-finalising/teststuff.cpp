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
#include "libraryfuncs.hpp"
#include "matrix_gen.hpp"

//namespace stdexec=std::execution;
// using namespace std::this_thread;

using namespace mylibs;
using namespace mymats;

int main(){
    std::vector<double>A={0.5,0.5,0.5,1.0};
    std::vector<double> Arow(0);
    std::vector<double> Acolumn(0);
    std::vector<double> Adata(0);
    mylibs::coo(A,Arow,Acolumn,Adata,sqrt(A.size()),sqrt(A.size()));
    std::vector<double>B={0.5,0.5,0.5,1.0};
    std::vector<double> Brow(0);
    std::vector<double> Bcolumn(0);
    std::vector<double> Bdata(0);
    mylibs::coo(B,Brow,Bcolumn,Bdata,sqrt(B.size()),sqrt(B.size()));
    auto [Cindcols,Cdata]=mylibs::coo_spgemm(Arow,Acolumn,Adata,sqrt(A.size()),sqrt(A.size()),Brow,Bcolumn,Bdata,sqrt(B.size()),sqrt(B.size()));
    int cdatsize = Cdata.size();
    std::cout<<cdatsize;
    //auto v = std::views::iota(0, cdatsize);
    // std::for_each(v.begin(),v.end(),[&](auto iv) {
    //     std::cout<<Cdata[iv]<<std::endl;
    // });
    return 0;
}