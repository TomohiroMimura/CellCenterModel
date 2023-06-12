#include<iostream>
#include<cmath>
#include<fstream>
#include<vector>
#include<time.h>

#include"vector3.hpp"
#include"class.hpp"
#include"force.hpp"
#include"ODE_solver.hpp"
#include"parameters.hpp"
#include"ShapeSampling.hpp"
#include"Random.hpp"
#include"restructure.hpp"
#include<omp.h>
int main(){
	#ifdef _OPENMP
	double start,end;
	start=omp_get_wtime();//clock関数だと並列計算で全てのスレッドの総和になるため
	#endif

	Random* p_rnd = new Random;
	std::cout<<"Random seed:"<<p_rnd->GetRandomSeed()<<std::endl;

	ShapeSampling* p_ss=new ShapeSampling(p_rnd);
	p_ss->CalcMDstep();
	//p_ss->Lloyd_Algorithm();
	std::cout << "Simulation finished!." << std::endl;

	delete p_rnd;
	delete p_ss;

	#ifdef _OPENMP
	end=omp_get_wtime();
	std::cout<<"duration="<<end-start<<"sec("<<(end-start)/60.0<<"min)"<<"\n";
	std::ofstream fout("d_Result.dat",std::ios::app);
	fout<<"duration="<<end-start<<"sec("<<(end-start)/60.0<<"min)"<<"\n";
	fout.close();
	#endif
	
	return 0;
}
