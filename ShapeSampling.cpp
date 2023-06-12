#include<omp.h>
#include<fstream>
#include"ShapeSampling.hpp"
#include"ODE_solver.hpp"
#include"class.hpp"
#include"Random.hpp"
#include"parameters.hpp"
ShapeSampling::ShapeSampling(Random*p_r){
    p_rnd=p_r;
    p_s=new Shape(p_rnd);
}

ShapeSampling::~ShapeSampling(){
    delete p_s;
}

void ShapeSampling::CalcMDstep(){
    p_s->OutputParameters();

    if(SurfaceApicalConstriction_circleMode==1)SurfaceApicalConstriction_circle();
    if(SurfaceApicalConstriction_circumferenceMode==1)SurfaceApicalConstriction_circumference();
    if(LineApicalConstrictionMode==1)LineApicalConstriction();
    if(SurfaceBeltApicalConstrictionMode==1)SurfaceBeltApicalConstriction();
    if(LineBeltApicalConstrictionMode==1)LineBeltApicalConstriction();


    for(int i=0;i<STEP_END+1;i++){
        double start,end;
        //p_s->CalcProximityPoint();//隣接点のcosの統計平均
        p_s->isConsistent();
        if(i%OUTPUT_VTK_STEP==0){
            p_s->CalcCrossBoundaryFlag();
            p_s->OutputVTK(i);
            p_s->OutputOFF(i);
            //p_s->OutputOFF_Mapping(i);//周期境界だとバグる
            p_s->OutputVTK_Voronoi(i);
        }
        /*----------------------------*/
        //細胞分裂
        #ifdef _OPENMP
        start=omp_get_wtime();
        #endif
        if(CELL_DIVISION_MODE!=0){
            //if(i%CellDivisionStep==0){
                restructure::cellDivision(p_rnd,p_s,i);//細胞分裂
                //p_s->OutputOFF(i);
                
            //}
        }
        #ifdef _OPENMP
        end=omp_get_wtime();
        std::ofstream fout2("d_TIME_CellDiv.dat",std::ios::app);
        fout2<<i<<" "<<end-start<<"\n";
        fout2.close();
        #endif
        /*------------------------------*/
        //p_s->isSuitable(i);
        //Fluctuation();
        
        std::ofstream ofs("d_Force.dat",std::ios::app);
        ofs<<"step-"<<i<<"\n";
        ofs.close();

        /*力学計算-------------------*/
        #ifdef _OPENMP
        start=omp_get_wtime();
        #endif
        if(DEGREE_ACCURACY==1){
            ODE_solver::motionVertexFirst(p_s);//力学計算
        }
        else if(DEGREE_ACCURACY==2){
            ODE_solver::motionVertexSecond(p_s);//力学計算
        }
        else if(DEGREE_ACCURACY==3){
            ODE_solver::motionVertexThird(p_s);
        }
        else if(DEGREE_ACCURACY==4){
            ODE_solver::motionVertexFourth(p_s);
        }
        #ifdef _OPENMP
        end=omp_get_wtime();
        std::ofstream fout3("d_TIME_Dynamics.dat",std::ios::app);
        fout3<<i<<" "<<end-start<<"\n";
        fout3.close();
        #endif
        /*---------------------------*/

        if(BOUNDARY_FIX_MODE==0){
            p_s->EraseBoundaryFaceAndEdge();//自由端での境界の処理(長くなった境界線を消す)
        }

        std::cout<<"step-"<<i<<std::endl;
        
        #ifdef _OPENMP
        start=omp_get_wtime();
        #endif
        
        if(FLIP_OPERATION_MODE==1){//フリップの判定と実行
            char fname[100]="d_Flip.dat";
            std::ofstream fout(fname,std::ios::app);
            fout<<"step-"<<i<<"\n";
            fout.close();
            p_s->FlipOperation();
        }
        #ifdef _OPENMP
        end=omp_get_wtime();
        std::ofstream fout4("d_TIME_Flip.dat",std::ios::app);
        fout4<<i<<" "<<end-start<<"\n";
        fout4.close();
        #endif


        p_s->updateCellTime();//細胞時間の更新
        if(i<TurnOffStep)p_s->ForcedDisplacement();//強制変位

        /*-----------出力ファイル系----------------*/
        if(i%Cell_Time_Output_Step==0){
            p_s->OutputCellTime(i);
            p_s->OutputOFF_periodic(i);
        }

        if(i%(CELL_DIVISION_STEP_AVE/100)==0){
            p_s->OutputFlipCount(i);//flipのカウントを出力
            for(int j=0;j<9;j++){
                p_s->OutputDivisionCount(i,j);
            }
            /*
            for(int j=0;j<5;j++){
                p_s->OutputInternalCount(i,j);
            }
            */
            p_s->OutputRange(i);
            p_s->OutputEnergy(i);
            p_s->OutputDistortion(i);
        
            p_s->VoronoiCenterEvaluation(i);
            p_s->CellCenterEvaluation(i);
        }
    }
    p_s->OutputResult();//実行結果のカウントを出力
}
/*
void ShapeSampling::Fluctuation(){//z方向の揺らぎを与える
    std::cout<<"Fluctuation"<<std::endl;
	for(int i=0;i<p_s->GetFace().size();i++){
		p_s->GetFace(i)->CalcArea(0);
        p_s->GetFace(i)->CalcNormalVector(0);
        Vector3<double> n=p_s->GetFace(i)->GetNormalVector();
		if(p_s->GetFace(i)->GetArea()<p_s->GetFace(i)->GetBaseArea()/2.0){
			for(int j=0;j<3;j++){
                if(p_s->GetVertex(p_s->GetFace(i)->GetVertex(j))->GetFixFlag()==0){
                    p_s->GetVertex(p_s->GetFace(i)->GetVertex(j))->SetLoc(p_s->periodize_vector(p_s->GetVertex(p_s->GetFace(i)->GetVertex(j))->GetLoc(0)+n*p_rnd->RAND_R(- BASE_LENGTH/100.0, BASE_LENGTH/100.0)),0);
                }
			}
		}
	}
}
*/
void ShapeSampling::SurfaceApicalConstriction_circumference(){
    for(int i=0;i<p_s->GetVertex().size();i++){
            if(p_s->GetVertex(i)->GetLoc(0).CalcNorm()<ApicalRadius_Max+0.1&&p_s->GetVertex(i)->GetLoc(0).CalcNorm()>ApicalRadius_Min+0.1-1.0){
                p_s->GetVertex(i)->SetApicalFlag(1);
            }
            if(CellDivisionRegionMode==2){
                if(p_s->GetVertex(i)->GetLoc(0).CalcNorm()<ApicalRadius_Max+0.1){
                    p_s->GetVertex(i)->SetFlagRec();
                }
            }
            else if(CellDivisionRegionMode==3){
                if(p_s->GetVertex(i)->GetLoc(0).CalcNorm()>ApicalRadius_Max-0.9){
                    p_s->GetVertex(i)->SetFlagRec();
                }
            }
            else if(CellDivisionRegionMode==0){
                p_s->GetVertex(i)->SetFlagRec();
            }
    }
}
void ShapeSampling::SurfaceApicalConstriction_circle(){
    for(int i=0;i<p_s->GetVertex().size();i++){
        //if(p_s->GetVertex(i)->GetBoundaryFlag()==0){
            //if(p_s->GetVertex(i)->GetLoc(0).CalcNorm()<ApicalRadius+0.1&&p_s->GetVertex(i)->GetLoc(0).CalcNorm()>ApicalRadius+0.1-1.0){
            if(p_s->GetVertex(i)->GetLoc(0).CalcNorm()<ApicalRadius+0.1){
                p_s->GetVertex(i)->SetApicalFlag(1);
            }
            if(CellDivisionRegionMode==2){
                if(p_s->GetVertex(i)->GetLoc(0).CalcNorm()<ApicalRadius+0.1){
                    p_s->GetVertex(i)->SetFlagRec();
                }
            }
            else if(CellDivisionRegionMode==3){
                if(p_s->GetVertex(i)->GetLoc(0).CalcNorm()>ApicalRadius+0.1-1.0){
                    p_s->GetVertex(i)->SetFlagRec();
                }
            }
            else if(CellDivisionRegionMode==0){
                p_s->GetVertex(i)->SetFlagRec();
            }
    }
}

void ShapeSampling::LineApicalConstriction(){//円形に指定
    std::vector<int> v1_list;//内側
    std::vector<int> v2_list;//外側
    for(int i=0;i<p_s->GetVertex().size();i++){
        if(p_s->GetVertex(i)->GetLoc(0).CalcNorm()>ApicalRadius+0.1-1.0&&p_s->GetVertex(i)->GetLoc(0).CalcNorm()<ApicalRadius+0.1){
            v2_list.push_back(i);
        }
        else if(p_s->GetVertex(i)->GetLoc(0).CalcNorm()>ApicalRadius+0.1-2.0&&p_s->GetVertex(i)->GetLoc(0).CalcNorm()<ApicalRadius+0.1-1.0){
            v1_list.push_back(i);
        }
    }
    for(int i=0;i<v1_list.size();i++){
        for(int j=0;j<p_s->GetVertex(v1_list[i])->GetEdge().size();j++){
            int n1=p_s->GetEdge(p_s->GetVertex(v1_list[i])->GetEdge(j))->GetVertex(0);
            if(n1==v1_list[i]){
                n1=p_s->GetEdge(p_s->GetVertex(v1_list[i])->GetEdge(j))->GetVertex(1);
            }
            auto result=std::find(v2_list.begin(),v2_list.end(),n1);
            if(result!=v2_list.end()){//v1_listとv2_list間を結ぶ線が見つかったとき

                Edge* p_e1= p_s->GetEdge(p_s->GetVertex(v1_list[i])->GetEdge(j));
                p_e1->SetApicalFlag(2);
            }
        }
    }
    for(int i=0;i<p_s->GetVertex().size();i++){
        if(CellDivisionRegionMode==2){
            if(p_s->GetVertex(i)->GetLoc(0).CalcNorm()<ApicalRadius+0.1-1.0){
                p_s->GetVertex(i)->SetFlagRec();
            }
        }
        else if(CellDivisionRegionMode==3){
            if(p_s->GetVertex(i)->GetLoc(0).CalcNorm()>ApicalRadius+0.1-1.0){
                p_s->GetVertex(i)->SetFlagRec();
            }
        }
        else if(CellDivisionRegionMode==0){
            p_s->GetVertex(i)->SetFlagRec();
        }
    }
}
void ShapeSampling::SurfaceBeltApicalConstriction(){
    for(int i=0;i<p_s->GetVertex().size();i++){
        //if(p_s->GetVertex(i)->GetBoundaryFlag()==0){
            //if(fabs(p_s->GetVertex(i)->GetLoc(0).GetX())<ApicalRadius+0.1&&fabs(p_s->GetVertex(i)->GetLoc(0).GetX())>ApicalRadius+0.1-1.0){
            if(fabs(p_s->GetVertex(i)->GetLoc(0).GetX())<ApicalRadius+0.1){
                p_s->GetVertex(i)->SetApicalFlag(1);
            }
            if(CellDivisionRegionMode==2){
                if(fabs(p_s->GetVertex(i)->GetLoc(0).GetX())<ApicalRadius+0.1){
                    p_s->GetVertex(i)->SetFlagRec();
                }
            }
            else if(CellDivisionRegionMode==3){
                if(fabs(p_s->GetVertex(i)->GetLoc(0).GetX())>ApicalRadius-0.9){
                    p_s->GetVertex(i)->SetFlagRec();
                }
            }
            else if(CellDivisionRegionMode==0){
                p_s->GetVertex(i)->SetFlagRec();
            }
    }
}
void ShapeSampling::LineBeltApicalConstriction(){
    std::vector<int> v1_list;//内側
    std::vector<int> v2_list;//外側
    for(int i=0;i<p_s->GetVertex().size();i++){
        if(fabs(p_s->GetVertex(i)->GetLoc(0).GetX())>ApicalRadius+0.1-1.0&&fabs(p_s->GetVertex(i)->GetLoc(0).GetX())<ApicalRadius+0.1){
            v2_list.push_back(i);
        }
        else if(fabs(p_s->GetVertex(i)->GetLoc(0).GetX())>ApicalRadius+0.1-2.0&&fabs(p_s->GetVertex(i)->GetLoc(0).GetX())<ApicalRadius+0.1-1.0){
            v1_list.push_back(i);
        }
    }
    for(int i=0;i<v1_list.size();i++){
        for(int j=0;j<p_s->GetVertex(v1_list[i])->GetEdge().size();j++){
            int n1=p_s->GetEdge(p_s->GetVertex(v1_list[i])->GetEdge(j))->GetVertex(0);
            if(n1==v1_list[i]){
                n1=p_s->GetEdge(p_s->GetVertex(v1_list[i])->GetEdge(j))->GetVertex(1);
            }
            auto result=std::find(v2_list.begin(),v2_list.end(),n1);
            if(result!=v2_list.end()){//v1_listとv2_list間を結ぶ線が見つかったとき

                Edge* p_e1= p_s->GetEdge(p_s->GetVertex(v1_list[i])->GetEdge(j));
                p_e1->SetApicalFlag(2);
            }
        }
    }
    for(int i=0;i<p_s->GetVertex().size();i++){
        if(CellDivisionRegionMode==2){
            if(fabs(p_s->GetVertex(i)->GetLoc(0).GetX())<ApicalRadius+0.1-1.0){
                p_s->GetVertex(i)->SetFlagRec();
            }
        }
        else if(CellDivisionRegionMode==3){
            if(fabs(p_s->GetVertex(i)->GetLoc(0).GetX())>ApicalRadius+0.1-1.0){
                p_s->GetVertex(i)->SetFlagRec();
            }
        }
        else if(CellDivisionRegionMode==0){
            p_s->GetVertex(i)->SetFlagRec();
        }
    }
}
void ShapeSampling::Lloyd_Algorithm(){
    for(int i=0;i<STEP_END+1;i++){
        p_s->OutputOFF(i);
        p_s->Lloyd_Algorithm();
        p_s->FlipOperation();
    }
}
Shape* ShapeSampling::GetShape()const{
    return p_s;
}
Random* ShapeSampling::GetRandom()const{
    return p_rnd;
}