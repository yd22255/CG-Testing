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
//#include <exec/static_thread_pool.hpp>
#define assert
//all this needs to be parallelised.

using namespace mylibs;

std::vector<double> mylibs::sparsemvdot_product(std::vector<double> Arow, std::vector<double> Acolumn,std::vector<double>Adata, const std::span<double> vector){   
    
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

std::vector<double> mylibs::vectortranspose(const std::vector<double> vector){   
    return vector; //returns 0 on compilerexplorer, here mostly for the sake of developer mathematical understanding
} 

double mylibs::vvdot_product(const std::vector<double> v1, const std::vector<double> v2){   
    double result = inner_product(v1.begin(), v1.end(), v2.begin(), 0.0l);
    return result;
}

std::pair<std::pair<std::vector<double>,std::vector<double>>,std::vector<double>> mylibs::matrixtranspose(std::vector<double> Mrow,std::vector<double> Mcolumn,std::vector<double> Mdata){  
    std::vector<double> resultRow(Mrow.size(),0);
    std::vector<double> resultColumn(Mcolumn.size(),0);
    //in theory
    //transposing a sparse matrix is just flipping Row/Column right. 
    //so i just loop through length(Adata) and swap them around?
    int vectorsize=Mdata.size();
    int vs2 = vectorsize*vectorsize;
    auto v = std::views::iota(0, vs2);
    std::for_each(v.begin(),v.end(),[&](auto iv) {
            int i = iv / vectorsize;
            int j = iv % vectorsize; 
            resultRow[i] = Mcolumn[i];
            resultColumn[i]=Mrow[i];});

    std::pair<std::vector<double>,std::vector<double>> tuple = std::make_pair(resultRow,resultColumn);
    std::pair<std::pair<std::vector<double>,std::vector<double>>,std::vector<double>> toreturn = std::make_pair(tuple,Mdata);
    return toreturn;
}

std::vector<double> mylibs::densematrixtranspose(std::vector<double> matrix){
    int vectorsize=sqrt(matrix.size());  
    std::vector<double> result(matrix.size(),0);
    std::mdspan<double,std::dextents<size_t,2>> matrixspan {matrix.data(),vectorsize,vectorsize};
    std::mdspan<double,std::dextents<size_t,2>> resultspan {result.data(),vectorsize,vectorsize};
    int vs2 = vectorsize*vectorsize;
    auto v = std::views::iota(0, vs2);
    std::for_each(v.begin(),v.end(),[=](auto iv) {
            int i = iv / vectorsize;
            int j = iv % vectorsize; 
            resultspan[j,i] = matrixspan[i,j];});
    return result;
}

std::vector<double> mylibs::vectorize(std::vector<double>Arow,std::vector<double>Acolumn,std::vector<double>Adata){
    std::vector<double> v(Arow.size());
    for (int i=1; i <= Arow.size(); i++){ // just copy the values from the first column
        v[i]=Arow[i];
    }
    return v;
}

double mylibs::norm(std::vector<double> vector){
    int vectorsize = vector.size();
    auto v = std::views::iota(0, vectorsize);
    std::vector<double> sumtotal(v.size(),0);
    std::span sumspan{sumtotal.data(),sumtotal.size()};
    std::for_each(v.begin(),v.end(),[sumspan,sumtotal,vector](auto i) {
        sumspan[i] += (vector[i]); 
    });
    double addupsumtotal = std::accumulate(sumtotal.begin(),sumtotal.end(),0);
    return sqrt(addupsumtotal);
}

std::pair<std::pair<std::vector<double>,std::vector<double>>,std::vector<double>> mylibs::cootocsr(std::vector<double> Arow,std::vector<double> Acolumn,std::vector<double> Adata,int numrows,int numcols){
    int nonzeros=Adata.size();
    std::vector<double>csrRows(0,numrows+1);
    for (int i = 0; i < nonzeros; i++){
        csrRows[Arow[i]+1]++;
    }
    for (int i = 0; i < numrows; i++){
        csrRows[i + 1] += csrRows[i];
    }
    std::pair<std::vector<double>,std::vector<double>> tuple = std::make_pair(csrRows,Acolumn);
    std::pair<std::pair<std::vector<double>,std::vector<double>>,std::vector<double>> toreturn = std::make_pair(tuple,Adata);
    return toreturn;  
}

int mylibs::csrvalueaccess(std::vector<double> Ainds, std::vector<double>Adata,int x, int y){
    double remainder=Ainds[y+1]-Ainds[y];
    return Adata[remainder+(x-1)];
}

int mylibs::csrvalueupdate(std::vector<double> Ainds,std::vector<double> Acol, std::vector<double>Adata,int x, int y,double total){
    double remainder=Ainds[y+1]-Ainds[y];
    Adata[remainder+(x-1)]=total;
    //the data just goes into the row of A we're operating on. 
    
    //update column and Cinds in here
}

std::pair<std::pair<std::vector<double>,std::vector<double>>,std::vector<double>> mylibs::coo_spgemm(std::vector<double> Arow,std::vector<double> Acolumn,std::vector<double> Adata,int Anumrows,int Anumcols,std::vector<double> Brow,std::vector<double> Bcolumn,std::vector<double> Bdata,int Bnumrows,int Bnumcols){
    //this is just spmm
    //for some reason people call it SPGEMM - sparse generalised matrix multiplication
    //coo is the requirements - takes in a COO-stored matrix
    assert(Acolumn.size()==Brow.size());

    std::vector<double> Cinds(0,Adata.size());
    std::vector<double> Ccolumn(0,Adata.size()); //Bcolumn.size()?
    std::vector<double> Cdata(0,Anumrows-1); //HOW BIG IS THIS? 
    auto[Aindcols,Adata]=cootocsr(Arow,Acolumn,Adata,Anumrows,Anumcols);
    auto[Bindcols,Bdata]=cootocsr(Brow,Bcolumn,Bdata,Bnumrows,Bnumcols);
    //figure how to do mm with csr
    //it is SUPPOSED to be easier

    //Gustavson's Algorithm.
    auto Ainds = Aindcols.first;
    auto Acols = Aindcols.second;
    auto Binds = Bindcols.first;
    auto Bcols = Bindcols.second;
    
    auto v = std::views::iota(0, Anumrows);
    auto v1 = std::views::iota(0, Anumcols); //could also be Bnumrows
    auto v2 = std::views::iota(0, Bnumcols);
    std::for_each(v.begin(),v.end(),[=](auto m){
        std::for_each(v1.begin(),v1.end(),[=](auto k){
            std::for_each(v2.begin(),v2.end(),[=](auto n){
                //figure out how to access C[m,n] += A[m,k] * B[k,n] in csr format
                //going to do this the way I've done it before, it's incredibly messy but I just want something that Works.
                //cant do that bc it's CSR not COO.
                
                auto Cmn=csrvalueaccess(Cinds,Cdata,m,n);
                auto Amk = csrvalueaccess(Ainds,Adata,m,k);
                auto Bkn=csrvalueaccess(Binds,Bdata,k,n);
                auto total=Cmn+(Amk*Bkn);
                csrvalueupdate(Cinds,Ccolumn,Cdata,m,n,total);

                //during iteration through indices

                std::cout<<"printloop"<<std::endl;
            });
        });
    });
    //CURRENTLY THIS WON'T UPDATE CColumn!!! VERY RELEVANT PROBLEM
    std::pair<std::vector<double>,std::vector<double>> tuple = std::make_pair(Cinds,Ccolumn);
    std::pair<std::pair<std::vector<double>,std::vector<double>>,std::vector<double>> toreturn = std::make_pair(tuple,Cdata);
    return toreturn;
}

// int nonsquarecoo(std::vector<double> A,std::vector<double> &Arow, std::vector<double>&Acolumn,std::vector<double>&Adata,int rows,int cols){ 
//     std::mdspan<double,std::dextents<size_t,2>> Aspan {A.data(),rows,cols};
//     Arow.clear();
//     Acolumn.clear();
//     Adata.clear();
//     int vs2 = rows*cols;
//     auto v = std::views::iota(0, vs2);
//     std::for_each(v.begin(),v.end(),[Aspan,&Arow,&Acolumn,&Adata,rows,cols](auto iv) {
//             int i = iv / rows; 
//             int j = iv % cols; //i do not know if this works.
//             if(Aspan[i,j]!=0){
//                 Arow.push_back(i);
//                 Acolumn.push_back(j);
//                 Adata.push_back(Aspan[i,j]);};
//             });
//             //this only works for square matrices ffs.
//     return 0;
// }

std::pair<std::pair<std::vector<double>,std::vector<double>>,std::vector<double>> mylibs::coo_jacobi(std::vector<double> &Arow, std::vector<double>&Acolumn,std::vector<double>&Adata){
    //take the matrix and create a diagonal matrix
    std::vector<double> retRow (Arow.size(),0);
    std::vector<double> retCol (Acolumn.size(),0);
    std::vector<double> retData(Adata.size(),0);
    int rowlen=Arow.size();
    int collen=Acolumn.size();
    auto vrow=std::views::iota(0,rowlen);
    auto vcol=std::views::iota(0,collen);
    std::for_each(vrow.begin(),vrow.end(),[&](auto i){
        std::for_each(vcol.begin(),vcol.end(),[&](auto j){
            if (Arow[i]==Acolumn[j]){
                retRow[i]=Arow[i];
                retCol[j]=Acolumn[j];
                retData[i]=Adata[i];
            }
        });
    });
    std::pair<std::vector<double>,std::vector<double>> tuple = std::make_pair(retRow,retCol);
    std::pair<std::pair<std::vector<double>,std::vector<double>>,std::vector<double>> toreturn = std::make_pair(tuple,retData);
}

std::pair<std::pair<std::vector<double>,std::vector<double>>,std::vector<double>> mylibs::coo_ilu(std::vector<double> &Arow, std::vector<double>&Acolumn,std::vector<double>&Adata){
    int size = Arow.size(); //think i should make this nonsquare
    std::vector<double> LRow (Arow.size(),0);
    std::vector<double> LCol (Acolumn.size(),0);
    std::vector<double> LData(Adata.size(),0);
    std::vector<double> URow (Arow.size(),0);
    std::vector<double> UCol (Acolumn.size(),0);
    std::vector<double> UData(Adata.size(),0);
    auto vsize=std::views::iota(0,size);
    std::for_each(vsize.begin(),vsize.end(),[&](auto i){
        //U
        std::for_each(vsize.begin(),vsize.end(),[&](auto j){
            int sum = 0;
            std::for_each(vsize.begin(),vsize.end(),[&](auto k){
                int datapts=Adata.size();
                auto vdsize=std::views::iota(0,datapts);
                double lowerik = 0;
                double upperkj=0;
                std::for_each(vdsize.begin(),vdsize.end(),[&](auto z){
                    if(LRow[z]==i && LCol[z]==k){
                        lowerik=LData[z];
                    }
                    if(URow[z]==k && UCol[z]==j){
                        upperkj=UData[z];
                    }
                });
                sum+=(lowerik*upperkj);
            });
            int datapts=Adata.size();
            auto vdsize=std::views::iota(0,datapts);
            std::for_each(vdsize.begin(),vdsize.end(),[&](auto z){
                    if(URow[z]==i && UCol[z]==j){
                        UData[z]=Adata[z]-sum;
                    }
                });
        });
        
        //L
        std::for_each(vsize.begin(),vsize.end(),[&](auto j){
            if(i==j){
                int datapts=Adata.size();
                auto vdsize=std::views::iota(0,datapts);
                std::for_each(vdsize.begin(),vdsize.end(),[&](auto z){
                    if(URow[z]==i && UCol[z]==i){
                        UData[z]=1;
                    }
                });
            }
            int sum = 0;
            std::for_each(vsize.begin(),vsize.end(),[&](auto k){
                int datapts=Adata.size();
                auto vdsize=std::views::iota(0,datapts);
                double lowerjk = 0;
                double upperki=0;
                std::for_each(vdsize.begin(),vdsize.end(),[&](auto z){
                    if(LRow[z]==j && LCol[z]==k){
                        lowerjk=LData[z];
                    }
                    if(URow[z]==k && UCol[z]==i){
                        upperki=UData[z];
                    }
                });
                sum+=(lowerjk*upperki);
            });
            int datapts=Adata.size();
            auto vdsize=std::views::iota(0,datapts);
            //this is a bit disgusting but i cant find a better way of accessing the different variables
            std::for_each(vdsize.begin(),vdsize.end(),[&](auto z){
                    if(LRow[z]==i && LCol[z]==j){
                        std::for_each(vdsize.begin(),vdsize.end(),[&](auto y){
                            if(Arow[y]==i && Acolumn[y]==j){
                                std::for_each(vdsize.begin(),vdsize.end(),[&](auto x){
                                    if(URow[x]==i && UCol[x]==i){
                                        LData[z]=(Adata[y]-sum)/UData[x];
                                    } 
                                });
                            }
                        });
                    }
                });
            });
        });
    auto [LUrowcols,LUdata]=coo_spgemm(LRow,LCol,LData,size,size,URow,UCol,UData,size,size);
    //convert this to a format that works :)
    std::pair<std::vector<double>,std::vector<double>> tuple = std::make_pair(LUrowcols.first,LUrowcols.second);
    std::pair<std::pair<std::vector<double>,std::vector<double>>,std::vector<double>> toreturn = std::make_pair(tuple,LUdata);
}



int mylibs::coo(std::vector<double> A,std::vector<double> &Arow, std::vector<double>&Acolumn,std::vector<double>&Adata,int rows,int cols){ 
    std::mdspan<double,std::dextents<size_t,2>> Aspan {A.data(),rows,cols};
    Arow.clear();
    Acolumn.clear();
    Adata.clear();
    int vs2 = rows*cols;
    auto v = std::views::iota(0, vs2);
    std::for_each(v.begin(),v.end(),[Aspan,&Arow,&Acolumn,&Adata,rows,cols](auto iv) {
            int i = iv / rows; 
            int j = iv % cols; 
            if(Aspan[i,j]!=0){
                Arow.push_back(i);
                Acolumn.push_back(j);
                Adata.push_back(Aspan[i,j]);};
            });
            //can i fix these being push_backs?
    return 0;
} 

int mylibs::csr(std::vector<double> A,std::vector<double> &Arow, std::vector<double>&Acolumn,std::vector<double>&Adata,int rows,int cols){
    return 0;
}

bool mylibs::roundtoepsilon(const double x, double epsilon){
    if (x<round(x)&&x>(round(x)-abs(epsilon))){
        return true;
    } else if (x>round(x)&&x<(round(x)+abs(epsilon))){
        return true;
    } else{
        return false;
    }
}

