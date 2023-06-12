#include"ODE_solver.hpp"
#include"parameters.hpp"
#include"force.hpp"
#include<omp.h>
#include<vector>
namespace ODE_solver{
	//一次精度で常微分方程式を解く(オイラー法)
	void motionVertexFirst(Shape* p_s){
		force::CalcTotalForce(p_s,0);
		force::OMP_Reduction_Frc(p_s,0);
	#ifdef _OPENMP
	#pragma omp parallel for num_threads(NUM_THREAD)
	#endif
		for(int i=0;i<p_s->GetVertex().size();i++){
			if(p_s->GetVertex(i)->GetFixFlag()==0){
				double viscosity=VISCOSITY;
				if(p_s->GetVertex(i)->GetBoundaryFlag()==1)viscosity=BOUNDARY_VISCOSITY;
				Vector3<double> v;
				v=p_s->GetVertex(i)->GetLoc(0)+p_s->GetVertex(i)->GetForce(0)*(DELTA_TIME)/(viscosity);
				p_s->GetVertex(i)->SetLoc(p_s->periodize_vector(v),0);//0に直接代入
			}
			p_s->GetVertex(i)->resetForce(0);
		}
	}
	//2次精度で微分方程式を解く(中点法).ただし、中点の座標と力は一次精度で求める
	void motionVertexSecond(Shape* p_s){
		force::CalcTotalForce(p_s,0);
		force::OMP_Reduction_Frc(p_s,0);
		//Δt/2後の位置を求める
		#ifdef _OPENMP
		#pragma omp parallel for num_threads(NUM_THREAD)
		#endif	
		for(int i=0;i<p_s->GetVertex().size();i++){
			if(p_s->GetVertex(i)->GetFixFlag()==0){
				double viscosity=VISCOSITY;
				if(p_s->GetVertex(i)->GetBoundaryFlag()==1)viscosity=BOUNDARY_VISCOSITY;
				Vector3<double>v;
				v=p_s->GetVertex(i)->GetLoc(0)+p_s->GetVertex(i)->GetForce(0)*(DELTA_TIME)/(2.0*viscosity);
				p_s->GetVertex(i)->SetLoc(p_s->periodize_vector(v),1);//1に代入
			}
			p_s->GetVertex(i)->resetForce(0);
		}
		force::CalcTotalForce(p_s,1);//座標1を基準に
		force::OMP_Reduction_Frc(p_s,1);
#ifdef _OPENMP
#pragma omp parallel for num_threads(NUM_THREAD)
#endif
		for(int i=0;i<p_s->GetVertex().size();i++){
			if(p_s->GetVertex(i)->GetFixFlag()==0){
				double viscosity=VISCOSITY;
				if(p_s->GetVertex(i)->GetBoundaryFlag()==1)viscosity=BOUNDARY_VISCOSITY;
				Vector3<double>v;
				v=p_s->GetVertex(i)->GetLoc(0)+p_s->GetVertex(i)->GetForce(1)*(DELTA_TIME)/(viscosity);
				p_s->GetVertex(i)->SetLoc(p_s->periodize_vector(v),0);//0
			}
			p_s->GetVertex(i)->resetForce(1);
		}
	}
	//3次精度(3次のルンゲクッタ)
	void motionVertexThird(Shape* p_s){
		force::CalcTotalForce(p_s,0);
		force::OMP_Reduction_Frc(p_s,0);
		//Δt/2後の位置を求める
		#ifdef _OPENMP
		#pragma omp parallel for num_threads(NUM_THREAD)
		#endif	
		for(int i=0;i<p_s->GetVertex().size();i++){
			if(p_s->GetVertex(i)->GetFixFlag()==0){
				double viscosity=VISCOSITY;
				if(p_s->GetVertex(i)->GetBoundaryFlag()==1){
					viscosity=BOUNDARY_VISCOSITY;
				}
				Vector3<double>v;
				v=p_s->GetVertex(i)->GetLoc(0)+p_s->GetVertex(i)->GetForce(0)*(DELTA_TIME)/(2*viscosity);
				p_s->GetVertex(i)->SetLoc(p_s->periodize_vector(v),1);
			}
			//p_s->GetVertex(i)->resetForce(0);
		}

		force::CalcTotalForce(p_s,1);
		force::OMP_Reduction_Frc(p_s,1);
		#ifdef _OPENMP
		#pragma omp parallel for num_threads(NUM_THREAD)
		#endif
		for(int i=0;i<p_s->GetVertex().size();i++){
			if(p_s->GetVertex(i)->GetFixFlag()==0){
				double viscosity=VISCOSITY;
				if(p_s->GetVertex(i)->GetBoundaryFlag()==1){
					viscosity=BOUNDARY_VISCOSITY;
				}
				Vector3<double>v;
				v=p_s->GetVertex(i)->GetLoc(0)+(2.0*p_s->GetVertex(i)->GetForce(1)-p_s->GetVertex(i)->GetForce(0))*(DELTA_TIME)/(2*viscosity);
				p_s->GetVertex(i)->SetLoc(p_s->periodize_vector(v),2);
			}
			//p_s->GetVertex(i)->resetForce(1);
		}
		
		force::CalcTotalForce(p_s,2);
		force::OMP_Reduction_Frc(p_s,2);
		#ifdef _OPENMP
		#pragma omp parallel for num_threads(NUM_THREAD)
		#endif
		for(int i=0;i<p_s->GetVertex().size();i++){
			if(p_s->GetVertex(i)->GetFixFlag()==0){
				double viscosity=VISCOSITY;
				if(p_s->GetVertex(i)->GetBoundaryFlag()==1){
					viscosity=BOUNDARY_VISCOSITY;
				}
				Vector3<double>v;
				v=p_s->GetVertex(i)->GetLoc(0)+(p_s->GetVertex(i)->GetForce(2)/6.0+p_s->GetVertex(i)->GetForce(1)/1.5+p_s->GetVertex(i)->GetForce(0)/6.0)*(DELTA_TIME)/(viscosity);
				p_s->GetVertex(i)->SetLoc(p_s->periodize_vector(v),0);
			}
			p_s->GetVertex(i)->resetForce(1);
			p_s->GetVertex(i)->resetForce(0);
			p_s->GetVertex(i)->resetForce(2);
		}
	}
	//4次精度(4次のルンゲクッタ)
	void motionVertexFourth(Shape* p_s){
		force::CalcTotalForce(p_s,0);
		force::OMP_Reduction_Frc(p_s,0);
		//Δt/2後の位置を求める
		#ifdef _OPENMP
		#pragma omp parallel for num_threads(NUM_THREAD)
		#endif	
		for(int i=0;i<p_s->GetVertex().size();i++){
			if(p_s->GetVertex(i)->GetFixFlag()==0){
				double viscosity=VISCOSITY;
				if(p_s->GetVertex(i)->GetBoundaryFlag()==1){
					viscosity=BOUNDARY_VISCOSITY;
				}
				Vector3<double>v;
				v=p_s->GetVertex(i)->GetLoc(0)+p_s->GetVertex(i)->GetForce(0)*(DELTA_TIME)/(2*viscosity);
				p_s->GetVertex(i)->SetLoc(p_s->periodize_vector(v),1);
			}
			//p_s->GetVertex(i)->resetForce(0);
		}
		//dt/2後の位置を求める
		force::CalcTotalForce(p_s,1);
		force::OMP_Reduction_Frc(p_s,1);
		#ifdef _OPENMP
		#pragma omp parallel for num_threads(NUM_THREAD)
		#endif
		for(int i=0;i<p_s->GetVertex().size();i++){
			if(p_s->GetVertex(i)->GetFixFlag()==0){
				double viscosity=VISCOSITY;
				if(p_s->GetVertex(i)->GetBoundaryFlag()==1){
					viscosity=BOUNDARY_VISCOSITY;
				}
				Vector3<double>v;
				v=p_s->GetVertex(i)->GetLoc(0)+(p_s->GetVertex(i)->GetForce(1))*(DELTA_TIME)/(2*viscosity);
				p_s->GetVertex(i)->SetLoc(p_s->periodize_vector(v),2);
			}
			//p_s->GetVertex(i)->resetForce(1);
		}
		
		force::CalcTotalForce(p_s,2);
		force::OMP_Reduction_Frc(p_s,2);
		#ifdef _OPENMP
		#pragma omp parallel for num_threads(NUM_THREAD)
		#endif
		for(int i=0;i<p_s->GetVertex().size();i++){
			if(p_s->GetVertex(i)->GetFixFlag()==0){
				double viscosity=VISCOSITY;
				if(p_s->GetVertex(i)->GetBoundaryFlag()==1){
					viscosity=BOUNDARY_VISCOSITY;
				}
				Vector3<double>v;
				v=p_s->GetVertex(i)->GetLoc(0)+(p_s->GetVertex(i)->GetForce(2))*(DELTA_TIME)/(viscosity);
				p_s->GetVertex(i)->SetLoc(p_s->periodize_vector(v),3);
			}
		}
		force::CalcTotalForce(p_s,3);
		force::OMP_Reduction_Frc(p_s,3);
		#ifdef _OPENMP
		#pragma omp parallel for num_threads(NUM_THREAD)
		#endif
		for(int i=0;i<p_s->GetVertex().size();i++){
			if(p_s->GetVertex(i)->GetFixFlag()==0){
				double viscosity=VISCOSITY;
				if(p_s->GetVertex(i)->GetBoundaryFlag()==1){
					viscosity=BOUNDARY_VISCOSITY;
				}
				Vector3<double>v;
				v=p_s->GetVertex(i)->GetLoc(0)+(p_s->GetVertex(i)->GetForce(0)/6.0+p_s->GetVertex(i)->GetForce(1)/3.0+p_s->GetVertex(i)->GetForce(2)/3.0+p_s->GetVertex(i)->GetForce(3)/6.0)*(DELTA_TIME)/(viscosity);
				p_s->GetVertex(i)->SetLoc(p_s->periodize_vector(v),0);
			}
			p_s->GetVertex(i)->resetForce(1);
			p_s->GetVertex(i)->resetForce(0);
			p_s->GetVertex(i)->resetForce(2);
			p_s->GetVertex(i)->resetForce(3);
		}
	}
	

}//namespace ODE_solver


