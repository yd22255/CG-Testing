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
#include "matrix_gen.hpp"


double mymats::random_element(){
    int min = 1;
    int max = 5;

    // Initialize a random number generator
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> range(min, max);

    // Generate random number in the range [min, max]
    double randomValue = range(gen);
    return randomValue; 
}

std::vector<double> mymats::diagonally_dominant(std::vector<double> matrix, int vectorsize){ //must be square
    std::mdspan<double,std::dextents<size_t,2>> matrixspan {matrix.data(),vectorsize,vectorsize};
                                                        //{ 
    for (int i = 0; i < matrixspan.extent(0); i++){     //|
        for (int j = 0; j < matrixspan.extent(1); j++){ //|
            std::cout<<matrixspan[i,j]<<" ";            //|
        }                                               //|
        std::cout<<"\n";                                //|
    }                                                   //|
    std::cout<<"\n";                                    //} purely for dev checking
    // exec::static_thread_pool pool(8);
    // auto sched = pool.get_scheduler();
    auto v=std::views::iota(0, vectorsize);
    // auto sumspansnd = stdexec::schedule(sched)
    // | stdexec::bulk(std::execution::par,matrixspan.extent(0), [&](int i) { 
    //     std::vector<double> sumtotal(v.size(),0);
    //     std::span sumspan{sumtotal.data(),sumtotal.size()};
    //     std::for_each(v.begin(),v.end(),[sumspan,sumtotal,matrixspan,i](auto iv) { 
    //         auto j = iv; 
    //         sumspan[i] += (matrixspan[i,j]);});});
    // stdexec::sync_wait(std::move(sumspansnd));
    //attempt at stdexecing the below:
    for (int i=0; i<matrixspan.extent(0); i++) {
        std::vector<double> sumtotal(v.size(),0);
        std::span sumspan{sumtotal.data(),sumtotal.size()};
        std::for_each(v.begin(),v.end(),[sumspan,sumtotal,matrixspan,i](auto iv) { 
            auto j = iv; 
            sumspan[i] += (matrixspan[i,j]); 
        });
        double addupsumtotal = std::accumulate(sumtotal.begin(),sumtotal.end(),0); //PROBLEM - learn transform_reduce
        matrixspan[i,i]=addupsumtotal;   
    }
    //PROBLEM:i do not know how i convert this to stdexec. do i just bulk the first for loop and set up a bunch of schedulers inside? 
    //that seems awkward though and potentially system-breaking

                                                        //{
    for (int i = 0; i < matrixspan.extent(0); i++){     //|
        for (int j = 0; j < matrixspan.extent(1); j++){ //|
            std::cout<<matrixspan[i,j]<<" ";            //|
        }                                               //|
        std::cout<<"\n";                                //|
    }                                                   //|
    std::cout<<"\n";                                    //} purely for dev checking
    return matrix;
}

std::vector<double> mymats::symmetric(std::vector<double> matrix, int vectorsize){
    std::mdspan<double,std::dextents<size_t,2>> matrixspan {matrix.data(),vectorsize,vectorsize};
    int vs2 = vectorsize*vectorsize;
    auto v = std::views::iota(0, vs2);
    std::for_each(v.begin(),v.end(),[matrixspan,vectorsize](auto iv) {
            int i = iv / vectorsize;
            int j = iv % vectorsize; 
            matrixspan[j,i] = matrixspan[i,j];
        });
    return matrix;
}

std::vector<double> mymats::generate_sparse_matrix(const double size, const double nonzeros){ 
    //this is currently generating the wrong number of nonzeros, because of symmetric/diag. 
    //need to figure out a better generation algorithm that guarantees the nonzero count at the end?
    //going to be really hard because of DD, nonzero count becomes minimum sqrt(size) plus the symmetric nonzeros
    int vectorsize = sqrt(size);
    std::vector<double> a(size,0.0);
    //generate the 0 matrix
    //PROBLEM:not sure how to stdexec this
    std::for_each_n(a.begin(), nonzeros+1, [](auto& n) { n = random_element(); });
    //fill the requisite number of random digits
    std::mt19937 rng(std::time(nullptr)); //generate the rng and pass it to shuffle
	std::shuffle(a.begin(), a.end(), rng); //shuffle
    //a should contain (nonzeros) random elements and (size-nonzeros) zeros, in a random order
    std::vector<double> aS=symmetric(a,vectorsize);
    std::vector<double> aDD=diagonally_dominant(aS,vectorsize);
    return aDD; 
}

std::vector<double> mymats::generate_vector(const double size){ //needs fixing to be able to generate 0s, but not hard
    int vectorsize = sqrt(size);
    std::vector<double> a(size,0.0);
    //generate 0 matrix
    //PROBLEM:not sure how to stdexec this
    std::for_each_n(a.begin(), size, [](auto& n) { n = random_element(); });
    std::mt19937 rng(std::time(nullptr)); //generate the rng and pass it to shuffle
	std::shuffle(a.begin(), a.end(), rng);//shuffle
    // a should contain (nonzeros) random elements and (size-nonzeros) zeros, in a random order
    return a; 
    
}