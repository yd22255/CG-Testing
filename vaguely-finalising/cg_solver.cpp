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
#include <exec/static_thread_pool.hpp>
#include "libraryfuncs.hpp"
#include "matrix_gen.hpp"
#include "precondition.hpp"

namespace stdexec=std::execution;
using namespace std::this_thread;
using namespace mylibs;
using namespace mymats;
using namespace precon;
//WHEN YOU EVENTUALLY COME BACK HERE: RUN THROUGH AND ADD NUMROW/NUMCOLUMN KNOWLEDGE ON EACH MATRIX
//IN THEORY, WE SET A NUMBER OF ROWS AND COLS WHEN WE MAKE A SPARSE MATRIX. Reuse those numbers.

// int coo(std::vector<double> A,std::vector<double> &Arow, std::vector<double>&Acolumn,std::vector<double>&Adata,int rows,int cols); 
// std::pair<std::pair<std::vector<double>,std::vector<double>>,std::vector<double>> coo_spgemm(std::vector<double> Arow,std::vector<double> Acolumn,std::vector<double> Adata,int Anumrows,int Anumcols,std::vector<double> Brow,std::vector<double> Bcolumn,std::vector<double> Bdata,int Bnumrows,int Bnumcols);
// std::vector<double> vectortranspose(const std::vector<double> vector); 
// bool roundtoepsilon(const double x, double epsilon);


std::vector<double> testmvdot_product(std::vector<double> Arow, std::vector<double> Acolumn,std::vector<double>Adata, const std::span<double> vector){   
    
    int vectorsize = vector.size();
    int vs2 = vectorsize*vectorsize;
    auto v = std::views::iota(0, vs2);
    std::vector<double> result(vector.size(),0.0);
    std::span resultspan{result.data(),result.size()};
    
    std::for_each(v.begin(),v.end(),[resultspan,Arow,Acolumn,Adata,vector,vectorsize](auto iv) {
            int i = iv / vectorsize; 
            int j = iv % vectorsize; 
            int arowsize = Arow.size();
            auto v3 = std::views::iota(0, arowsize);
            std::for_each(v3.begin(),v3.end(),[resultspan,vector,Arow,Acolumn,Adata,i,j](auto z){
                if(Arow[z]==i && Acolumn[z]==j){ 
                    resultspan[i] += (Adata[z] * vector[j]);
                }
            });
            });
    
    return result;
}

std::pair<std::pair<std::vector<double>,int>,bool>  sparsecgsolver(std::vector<double> Arow, std::vector<double>Acolumn,std::vector<double>Adata, std::span<double>b, std::span<double>x0,double epsilon,const int max_iters){
    //PROBLEM: 
    //i feel like ive violated some principles here but s&r aren't even out yet so theres no principles to stop me :)
    exec::static_thread_pool pool(8);
    auto sched = pool.get_scheduler();
    int vectorsize = x0.size();
    std::vector<double> solutions(x0.size());
    auto v1 = std::views::iota(0, vectorsize);
    // auto begin = stdexec::schedule(sched);
    auto solsnd = stdexec::schedule(sched)
    | stdexec::bulk(std::execution::par,vectorsize, [&](int i) { 
        solutions[i]=x0[i];
      });
    stdexec::sync_wait(std::move(solsnd));
    //this doesnt really want to be & but i dont have a good workaround right now
    // std::for_each(v1.begin(),v1.end(),[&solutions,x0](auto i){
    //     solutions[i]=x0[i];
    // }); //OI:1V/16V = 1/16
    
    std::span solspan{solutions.data(),solutions.size()};
    std::vector<double> r(x0.size(),0.0);
    std::span rspan{r.data(),r.size()};
    std::cout<<"test";
    std::vector<double>Ax0=testmvdot_product(Arow,Acolumn,Adata,x0);
    std::cout<<"test";
    auto v = std::views::iota(0, vectorsize);
    auto bax0snd = stdexec::schedule(sched)
    | stdexec::bulk(std::execution::par,vectorsize, [&](int i) { 
        rspan[i]=(b[i]-Ax0[i]);});
    stdexec::sync_wait(std::move(bax0snd));
    // std::for_each(v.begin(),v.end(),[rspan,b,Ax0,vectorsize](auto i) {
    //         rspan[i]=(b[i]-Ax0[i]);}); 
        //OI: 2V/24V = 1/12

    std::vector<double> d(r.size(),0.0);
    std::span dspan{d.data(),d.size()};
    auto dspansnd = stdexec::schedule(sched)
    | stdexec::bulk(std::execution::par,vectorsize, [&](int i) { 
        dspan[i]=rspan[i];});
    stdexec::sync_wait(std::move(dspansnd));
    // std::for_each(v.begin(),v.end(),[dspan,rspan,vectorsize](auto i) {
    //         dspan[i]=rspan[i];});
        //OI: 1V/16V = 1/16

    //std::vector<double> rT = vectortranspose(r); //OI:0 (i think) 
    // double deltanew=vvdot_product(rT,r); //OI: 1/3x8 = 1/24?
    // double delta0=deltanew;
    // double deltaold=0;

    int iters=0;
//     while (iters<max_iters && deltanew>(epsilon*epsilon)*delta0){
//         std::vector<double> q = sparsemvdot_product(Arow,Acolumn,Adata,dspan); //OI:?
//         std::span qspan{q.data(),q.size()};
//         std::vector<double> dT = vectortranspose(d); //OI:?
//         double alpha = deltanew/vvdot_product(dT,q); //OI:?

//         std::vector<double>alphad(d.size());
//         std::span alphadspan{alphad.data(),alphad.size()}; 
//         auto alphadspansnd = stdexec::schedule(sched)
//         | stdexec::bulk(std::execution::par,vectorsize, [&](int i) { 
//             alphadspan[i]=(dspan[i]*alpha);});
//         stdexec::sync_wait(std::move(alphadspansnd));
//         // std::for_each(v.begin(),v.end(),[alphadspan,dspan,alpha,vectorsize](auto i) {
//         //         alphadspan[i]=(dspan[i]*alpha);});
        
        
// //this is messing with things as for_each doesn't seem to like changing of the data you're using inline ---------
// //i have found a solution but it is abhorrent. I want to change it but i have no idea how to.
// //------------------------------------------------------
//         std::vector<double> sol2(solutions.size(),0.0);
//         std::span sol2span{sol2.data(),sol2.size()};
//         auto solupdatesnd = stdexec::schedule(sched)
//         | stdexec::bulk(std::execution::par,vectorsize, [&](int i) { 
//             solspan[i]+=alphad[i];
//         });
//         stdexec::sync_wait(std::move(solupdatesnd));
//         // std::for_each(v.begin(),v.end(),[solspan,sol2span,vectorsize](auto i) {
//         //         sol2span[i]=solspan[i];});
//         // std::for_each(v.begin(),v.end(),[solspan,sol2span,alphad,vectorsize](auto i) {
//         //         solspan[i]=sol2span[i]+alphad[i];});

//         Ax0=sparsemvdot_product(Arow,Acolumn,Adata,solutions);
//         if(iters % 50==0){
//             auto rupdatesnd = stdexec::schedule(sched)
//             | stdexec::bulk(std::execution::par,vectorsize, [&](int i) { 
//                 rspan[i]=(b[i]-Ax0[i]);  });
//             stdexec::sync_wait(std::move(rupdatesnd));
//             // std::for_each(v.begin(),v.end(),[rspan,b,Ax0,vectorsize](auto i) {
//             //     rspan[i]=(b[i]-Ax0[i]);  });
//         }else{
// //this is also messing with things as for_each doesn't seem to like changing of the data you're using inline ---------
//             std::vector<double> r2(r.size(),0.0);
//             std::span r2span{r2.data(),r2.size()};
//             auto rupdatesnd = stdexec::schedule(sched)
//             | stdexec::bulk(std::execution::par,vectorsize, [&](int i) { 
//                 rspan[i]-=alpha*qspan[i];
//             });
//             stdexec::sync_wait(std::move(rupdatesnd));
//             // std::for_each(v.begin(),v.end(),[rspan,r2span,vectorsize](auto i) {
//             //     r2span[i]=rspan[i];});
//             // std::for_each(v.begin(),v.end(),[rspan,r2span,alpha,qspan,vectorsize](auto i) {
//             //     rspan[i]=(r2span[i]-alpha*qspan[i]);  });
//         }
//         deltaold=deltanew;
//         std::vector<double> rT = vectortranspose(r);
//         deltanew=vvdot_product(rT,r);
//         double beta = deltanew/deltaold;
// //got a feeling this will also mess for the same reason...
// //there's a super disgusting solution to just copy d before performing the iteration? but that seems vile. 
//         std::vector<double> d2(d.size(),0.0);
//         std::span d2span{d2.data(),d2.size()};
//         auto dupdatesnd = stdexec::schedule(sched)
//         | stdexec::bulk(std::execution::par,vectorsize, [&](int i) { 
//             dspan[i]=rspan[i]+beta*dspan[i];
//         });
//         stdexec::sync_wait(std::move(dupdatesnd));
//         // std::for_each(v.begin(),v.end(),[dspan,d2span,vectorsize](auto i) {
//         // d2span[i]=dspan[i];});
//         // std::for_each(v.begin(),v.end(),[rspan,d2span,dspan,beta,vectorsize](auto i) {
//         // dspan[i]=(rspan[i]+beta*d2span[i]);  });
//         iters+=1;
//     }
    std::pair<std::vector<double>,int> tuple = std::make_pair(solutions,iters);
    std::pair<std::pair<std::vector<double>,int>,bool> toreturn = std::make_pair(tuple,true);
    return toreturn;
}

std::vector<double> sparseverify(std::vector<double> &Arow, std::vector<double>&Acolumn,std::vector<double>&Adata,std::span<double>b, std::vector<double>xfinal){
    int vectorsize = xfinal.size();
    int vs2 = vectorsize*vectorsize;
    auto v = std::views::iota(0, vs2);
    std::vector<double> result(xfinal.size(),0.0);
    std::span resultspan{result.data(),result.size()};
    std::for_each(v.begin(),v.end(),[resultspan,Arow,Acolumn,Adata,xfinal,vectorsize](auto iv) {
            int i = iv / vectorsize; 
            int j = iv % vectorsize; 
            int arowsize = Arow.size();
            auto v3 = std::views::iota(0, arowsize);
            std::for_each(v3.begin(),v3.end(),[resultspan,xfinal,Arow,Acolumn,Adata,i,j](auto z){
                if(Arow[z]==i && Acolumn[z]==j){ 
                    resultspan[i] += (Adata[z] * xfinal[j]);
                }
            });
            });
    std::vector<double> error(vectorsize,0);
    for (size_t i = 0; i < vectorsize; i++){
        if(result[i]!=b[i]){
            error[i]=(abs(result[i]-b[i]));
        }else{error[i]=0;}
    } 
    return error;
}


int sparsemode(std::vector<double>A,std::vector<double>b,std::vector<double>x0,double epsilon,int max_iters){
    // int sumtotal = 0;
    // std::for_each(A.begin(),A.end(),[&sumtotal](auto x){
    //     if(x!=0){
    //         sumtotal+=1;
    //     }
    // });
    std::vector<double> Arow(0);
    std::vector<double> Acolumn(0);
    std::vector<double> Adata(0);
    coo(A,Arow,Acolumn,Adata,sqrt(A.size()),sqrt(A.size()));
    // std::for_each(Arow.begin(),Arow.end(),[](auto i) {
    //         std::cout<<i<<" ";});
    // std::cout<<std::endl;
    // std::for_each(Acolumn.begin(),Acolumn.end(),[](auto i) {
    //         std::cout<<i<<" ";});
    // std::cout<<std::endl;
    // std::for_each(Adata.begin(),Adata.end(),[](auto i) {
    //         std::cout<<i<<" ";});
    // std::cout<<std::endl;
    std::span<double> bspan{b.data(),b.size()};
    std::span<double> x0span{x0.data(),x0.size()};
    std::vector<double>xfinal={};
    int vectorsize = x0.size();
    // std::pair<std::pair<std::vector<double>,int>,bool> returned=cgsolver(Aspan,bspan,x0span,epsilon);
    auto [solsiters,complete]=sparsecgsolver(Arow,Acolumn,Adata,bspan,x0span,epsilon,max_iters);
    
    for (int i = 1 ; i <=vectorsize ; i++){
        xfinal.push_back(solsiters.first[(solsiters.first.size())-i]);
        //iterate back through matrix to get the last vector so i can verify the answer ax-b=0 (and print the error) 
        }
    std::reverse(xfinal.begin(),xfinal.end());
    std::cout<<"Final Solution: ";
    for (int i = 1 ; i <=xfinal.size() ; i++){
        std::cout<<xfinal[i-1]<<" ";
    }
    std::cout<<std::endl;
    std::vector<double> error = sparseverify(Arow,Acolumn,Adata,bspan,xfinal);
    std::cout<<"Total error: ";
    double errormag = 0;
    for (int i = 0 ; i <=error.size()-1 ; i++){
        if ((round(error[i])==0)&&(roundtoepsilon(x0[i],epsilon))){ 
                error[i]=0;
        } 
        errormag += abs(error[i]);
        std::cout<<error[i]<<" ";
    }
    
    std::cout<<"\nError Magnitude: "<<errormag<<std::endl;
    
    //std::cout<<"end"<<std::endl;
    return 0;
}

int runtest(std::vector<double>A,std::vector<double>b,std::vector<double>x0){
    std::vector<double> Arow(0);
    std::vector<double> Acolumn(0);
    std::vector<double> Adata(0);
    coo(A,Arow,Acolumn,Adata,sqrt(A.size()),sqrt(A.size()));
    std::vector<double> Ax0 = testmvdot_product(Arow,Acolumn,Adata,x0);
    std::vector<double> sparsetest = sparsemvdot_product(Arow,Acolumn,Adata,x0);
    return 0;
}

int main(){
    //std::vector<double> testA = generate_sparse_matrix(36,5); //has to be square size
    // std::vector<double> testb = generate_vector(85);
    // std::vector<double> testx0 = generate_vector(85);
    std::vector<double>A={{0.5,0.5,0.5,1.0}};
    std::vector<double>Asparse={0,0,0.5,0,1,0.5,1,0,0.5,1,1,1.0};
    std::mdspan<double,std::dextents<size_t,2>> Asparsespan {Asparse.data(),3,Asparse.size()/3};

    std::vector<double>b={{0.0,2.0}};
    std::vector<double>x0={{2.3,-2.2}};
    std::vector<double>xfinal(x0.size(),0);
    int vectorsize = x0.size();
    double epsilon =1e-6;
    int max_iters = 100;
    sparsemode(A,b,x0,epsilon,max_iters);
    runtest(A,b,x0);
    
    return 0;
}