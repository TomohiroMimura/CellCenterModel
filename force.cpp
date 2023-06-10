#include<vector>
#include<cmath>
#include"force.hpp"
#include"parameters.hpp"
#include"vector3.hpp"
#include"class.hpp"
#include<omp.h>
#include<fstream>

namespace force{
    void CalcTotalForce(Shape*p_s,int deg){
        //重複して計算しないように前もって計算
        #ifdef _OPENMP
        #pragma omp parallel for schedule(guided) num_threads(NUM_THREAD)
        #endif
        for(int i=0;i<NUM_THREAD;i++){
            p_s->LineEnergy[i]=0.0;
            p_s->AreaEnergy[i]=0.0;
            p_s->DihedralAngleEnergy[i]=0.0;
            p_s->Z_Energy[i]=0.0;
            p_s->V_Center_Energy[i]=0.0;
            p_s->RepulsiveEnergy[i]=0.0;
            p_s->ApicalConstrictionEnergy[i]=0.0;

            p_s->VoronoiAreaEnergy[i]=0.0;
            p_s->VoronoiPerimeterEnergy[i]=0.0;
        }

        #ifdef _OPENMP
        #pragma omp parallel for schedule(guided) num_threads(NUM_THREAD)
        #endif
        for(int i=0;i<p_s->GetFace().size();i++){
            p_s->GetFace(i)->CalcNormalVector(deg);
            p_s->GetFace(i)->CalcArea(deg);
        }
        #ifdef _OPENMP
        #pragma omp parallel for schedule(guided) num_threads(NUM_THREAD)
        #endif
        for(int i=0;i<p_s->GetVertex().size();i++){
            p_s->GetVertex(i)->CalcNormalVector2(deg);
        }
        #ifdef _OPENMP
        #pragma omp parallel for schedule(guided) num_threads(NUM_THREAD)
        #endif
        for(int i=0;i<p_s->GetEdge().size();i++){
            p_s->GetEdge(i)->CalcLength(deg);
            if(p_s->GetEdge(i)->GetFace().size()==2)p_s->GetEdge(i)->CalcDihedralAngle(deg);
        }
        #ifdef _OPENMP
        #pragma omp parallel for schedule(guided) num_threads(NUM_THREAD)
        #endif
        for(int i=0;i<p_s->GetFace().size();i++){
            p_s->GetFace(i)->CalcCenter(deg);//重心
        }

            if(K1_T_LENGTH!=0||K2_T_LENGTH!=0){//lineの長さ保存
                CalcLineForce(p_s,deg);
            }
            if(K_T_AREA!=0){//三角形の面積保存
                CalcAreaForce(p_s,deg);
            }
            if(K_T_INTERIOR_ANGLE!=0){//三角形の内角保存(没)
                CalcInteriorAngleForce(p_s,deg);
            }
            
            if(K_T_THETA!=0){//三角形同士の二面角保存
                CalcDihedralAngleForce(p_s,deg);
            }

            if(K_Z!=0){//z方向変形束縛力
                CalcZ_axisRestraintForce(p_s,deg);
            }
            if(K_Y!=0){
                CalcY_axisRestraintForce(p_s,deg);
            }

            if(K_V_CENTER!=0||K_V_OUT_CENTER!=0){//重心の重心を母点に近づける
                CalcVoronoiCenterForce(p_s,deg);
            }

            if(K_V_AREA!=0||K_V_PERIMETER!=0){
                #ifdef _OPENMP
                #pragma omp parallel for schedule(guided) num_threads(NUM_THREAD)
                #endif
                for(int i=0;i<p_s->GetVertex().size();i++){
                    p_s->GetVertex(i)->CalcVoronoiArea(deg);
                    p_s->GetVertex(i)->CalcVoronoiPerimeter(deg);
                }
                CalcVoronoiAreaForce(p_s,deg);
                CalcVoronoiPerimeterForce(p_s,deg);
            }

            CalcApicalConstrictionForce(p_s,deg);
            //apical収縮によるedgeの収縮

            if(K_REP!=0){
                CalcRepulsiveForce(p_s,deg);
            }
            if(K_lateral!=0&&Layer_NUM>=2)CalcLateralForce(p_s,deg);
            if(Fluct!=0){//(没)作用反作用を満たすのかが不明瞭
                Fluctuation(p_s,deg);
            }
            if(pressure!=0)CalcPressure(p_s,deg);

            if(Bending_Elastic_Constant!=0){//後に置いといた方が確実(voronoi_areaに代入するから)
                if(CurvatureType==1){//julicher
            #ifdef _OPENMP
            #pragma omp parallel for schedule(guided) num_threads(NUM_THREAD)
            #endif
                    for(int i=0;i<p_s->GetVertex().size();i++){
                        p_s->GetVertex(i)->CalcVoronoiArea(deg);
                        p_s->GetVertex(i)->CalcMeanCurvature();//julicher
                    }
                    CalcCurvatureForce(p_s,deg);
                }
                else if(CurvatureType==2){//laplace-beltrami
            #ifdef _OPENMP
            #pragma omp parallel for schedule(guided) num_threads(NUM_THREAD)
            #endif
                    for(int i=0;i<p_s->GetFace().size();i++){
                        p_s->GetFace(i)->CalcCircumCenter(deg);
                    }
            #ifdef _OPENMP
            #pragma omp parallel for schedule(guided) num_threads(NUM_THREAD)
            #endif
                    for(int i=0;i<p_s->GetVertex().size();i++){
                        p_s->GetVertex(i)->CalcMeanCurvature(deg);//laplace-beltrami
                    }                    
                    CalcCurvatureForce2(p_s,deg);
                }
            }
    }

    void CalcLineForce(Shape*p_s,int deg){
        //deg=degree accuracy
        
#ifdef _OPENMP
#pragma omp parallel for num_threads(NUM_THREAD)
#endif
        for(int i=0;i<p_s->GetEdge().size();i++){
            CalcLineForce(p_s->GetEdge(i),deg);
        }
    }

    void CalcAreaForce(Shape*p_s,int deg){
#ifdef _OPENMP
#pragma omp parallel for num_threads(NUM_THREAD)
#endif
        for(int i=0;i<p_s->GetFace().size();i++){
            CalcAreaForce(p_s->GetFace(i),deg);
        }
    }
    void CalcDihedralAngleForce(Shape*p_s,int deg){
#ifdef _OPENMP
#pragma omp parallel for num_threads(NUM_THREAD)
#endif
        for(int i=0;i<p_s->GetEdge().size();i++){
            CalcDihedralAngleForce(p_s->GetEdge(i),deg);
        }
    }

    void CalcLineForce(Edge* p_e,int deg){
        int id_tnum = omp_get_thread_num();//並列処理用
        //p_e->CalcLength(deg);//必要(9/24追記:エラーの原因になる)
        Vector3<double> line_vec=p_e->p_s->periodize_vector(p_e->p_s->GetVertex(p_e->GetVertex(0))->GetLoc(deg)-p_e->p_s->GetVertex(p_e->GetVertex(1))->GetLoc(deg));
        Vector3<double>l_grad=line_vec/(line_vec.CalcNorm()+std::pow(10.0,-5));
        Vector3<double>frc_tmp=l_grad*(p_e->Get_k1_length()+p_e->Get_k2_length()*(p_e->GetLength()-p_e->GetBaseLength()));
        //1次ポテンシャルと2次ポテンシャル
        p_e->p_s->GetVertex(p_e->GetVertex(0))->addFrcThread(-frc_tmp,id_tnum);
        p_e->p_s->GetVertex(p_e->GetVertex(1))->addFrcThread(frc_tmp,id_tnum);
        p_e->p_s->LineEnergy[id_tnum]+=(p_e->Get_k1_length()*+0.5*p_e->Get_k2_length()*(p_e->GetLength()-p_e->GetBaseLength()))*(p_e->GetLength()-p_e->GetBaseLength());
        //double r=(DELTA_TIME/VISCOSITY);
        /*
        if(frc_tmp.CalcNorm()*r>1){
            std::ofstream fout("d_Force.dat",std::ios::app);
            fout<<"LineForce\n";
            fout<<"Edge-"<<p_e->GetVertex(0)<<" and "<<p_e->GetVertex(1)<<"\n";
            fout<<"frc0 = "<<-frc_tmp<<"\n";
            fout<<"frc1 = "<<frc_tmp<<"\n";
            fout.close();
        }
        */
    }
    void CalcAreaForce(Face* p_f,int deg){
        int id_tnum = (int)omp_get_thread_num();//並列処理用
        Vector3<double> a=p_f->p_s->periodize_vector(p_f->p_s->GetVertex(p_f->GetVertex(1))->GetLoc(deg)-p_f->p_s->GetVertex(p_f->GetVertex(0))->GetLoc(deg));
        Vector3<double> b=p_f->p_s->periodize_vector(p_f->p_s->GetVertex(p_f->GetVertex(2))->GetLoc(deg)-p_f->p_s->GetVertex(p_f->GetVertex(0))->GetLoc(deg));
        //p_f->CalcNormalVector(deg);
        //p_f->CalcArea(deg);

        //x0について
        Vector3<double>s_gradx0=0.5*p_f->GetNormalVector()%p_f->p_s->periodize_vector(b-a);
        Vector3<double> frc_tmpx0=s_gradx0*p_f->Get_k_area()*(p_f->GetArea()-p_f->GetBaseArea());
        p_f->p_s->GetVertex(p_f->GetVertex(0))->addFrcThread(-frc_tmpx0,id_tnum);
        //x1について
        Vector3<double>s_gradx1=-0.5*p_f->GetNormalVector()%b;
        Vector3<double>frc_tmpx1=s_gradx1*p_f->Get_k_area()*(p_f->GetArea()-p_f->GetBaseArea());
        p_f->p_s->GetVertex(p_f->GetVertex(1))->addFrcThread(-frc_tmpx1,id_tnum);
        //x2について
        Vector3<double>s_gradx2=0.5*p_f->GetNormalVector()%a;
        Vector3<double>frc_tmpx2=s_gradx2*p_f->Get_k_area()*(p_f->GetArea()-p_f->GetBaseArea());
        p_f->p_s->GetVertex(p_f->GetVertex(2))->addFrcThread(-frc_tmpx2,id_tnum);

        p_f->p_s->AreaEnergy[id_tnum]+=0.5*p_f->Get_k_area()*(p_f->GetArea()-p_f->GetBaseArea())*(p_f->GetArea()-p_f->GetBaseArea());
        //面積弾性
        /*
        double r=(DELTA_TIME/VISCOSITY);
        if(frc_tmpx0.CalcNorm()*r>1||frc_tmpx1.CalcNorm()*r>1||frc_tmpx2.CalcNorm()*r>1){
            std::ofstream fout("d_Force.dat",std::ios::app);
            fout<<"AreaForce\n";
            fout<<"Triangle-"<<p_f->GetVertex(0)<<" and "<<p_f->GetVertex(1)<<" and "<<p_f->GetVertex(2)<<"\n";
            fout<<"frc0 = "<<-frc_tmpx0<<"\n";
            fout<<"frc1 = "<<-frc_tmpx1<<"\n";
            fout<<"frc2 = "<<-frc_tmpx2<<"\n";
            fout.close();
        }
        */
    }
    void CalcDihedralAngleForce(Edge* p_e,int deg){
        int id_tnum = (int)omp_get_thread_num();//並列処理用
        if(p_e->GetFace().size()==2){
            Vector3<double> x0=p_e->p_s->GetVertex(p_e->GetVertex(0))->GetLoc(deg);
            Vector3<double> x1=p_e->p_s->GetVertex(p_e->GetVertex(1))->GetLoc(deg);
            Vertex* v2;
            Vertex* v3;
            Vector3<double>x2;
            Vector3<double>x3;
            for(int i=0;i<p_e->GetFace().size();i++){
                for(int k=0;k<3;k++){
                    if(p_e->p_s->GetFace(p_e->GetFace(i))->GetVertex(k)==p_e->GetVertex(0) && p_e->p_s->GetFace(p_e->GetFace(i))->GetVertex((k+1)%3)==p_e->GetVertex(1)){
                        v2=p_e->p_s->GetVertex(p_e->p_s->GetFace(p_e->GetFace(i))->GetVertex((k+2)%3));
                        x2=v2->GetLoc(deg);
                        break;
                    }
                    else if(p_e->p_s->GetFace(p_e->GetFace(i))->GetVertex(k)==p_e->GetVertex(0) && p_e->p_s->GetFace(p_e->GetFace(i))->GetVertex((k+2)%3)==p_e->GetVertex(1)){
                        v3=p_e->p_s->GetVertex(p_e->p_s->GetFace(p_e->GetFace(i))->GetVertex((k+1)%3));
                        x3=v3->GetLoc(deg);
                        break;
                    }
                }
            }

            Vector3<double> e0=p_e->p_s->periodize_vector(x1-x0);
            Vector3<double> e1=p_e->p_s->periodize_vector(x2-x0);
            Vector3<double> e2=p_e->p_s->periodize_vector(x3-x0);
            Vector3<double> n1=e0%e1;
            Vector3<double> n2=e2%e0;

            //p_e->CalcDihedralAngle(deg);

            Vector3<double> theta_gradx0=-(p_e->p_s->periodize_vector(x2-x1)*e0)*n1/(e0.CalcNorm()*(n1*n1)+std::pow(10.0,-7))-(p_e->p_s->periodize_vector(x3-x1)*e0)*n2/(e0.CalcNorm()*(n2*n2)+std::pow(10.0,-7));
            Vector3<double> theta_gradx1=(e1*e0)*n1/(e0.CalcNorm()*(n1*n1)+std::pow(10.0,-7))+(e2*e0)*n2/(e0.CalcNorm()*(n2*n2)+std::pow(10.0,-7));
            Vector3<double> theta_gradx2=-(e0.CalcNorm())*n1/((n1*n1)+std::pow(10.0,-7));
            Vector3<double> theta_gradx3=-(e0.CalcNorm())*n2/((n2*n2)+std::pow(10.0,-7));


           double k=p_e->Get_k_theta();
           //double k=p_e->Get_k_theta()*3*(e0.CalcSqr())/(0.5*n1.CalcNorm()+0.5*n2.CalcNorm());
           //double k=p_e->Get_k_theta()*e0.CalcNorm()*(p_e->p_s->GetVertex(p_e->GetVertex(0))->GetMeanCurvature()+p_e->p_s->GetVertex(p_e->GetVertex(1))->GetMeanCurvature());
           double BaseAngle=p_e->GetBaseDihedralAngle()+p_e->p_s->p_rnd->RAND_Normal(0,fluctSD);
           if(p_e->GetApicalFlag()==2){
               BaseAngle-=APICAL_CONSTRICTION_ANGLE_L;
               k*=M_APICAL_THETA_L;
           }
           if(p_e->p_s->GetVertex(p_e->GetVertex(0))->GetApicalFlag()==1||p_e->p_s->GetVertex(p_e->GetVertex(0))->GetApicalFlag()==1){
               BaseAngle-=APICAL_CONSTRICTION_ANGLE_S;
               k*=M_APICAL_THETA_S;
           }
           /*
            Vector3<double> frc0=-k*theta_gradx0;
           Vector3<double> frc1=-k*theta_gradx1;
           Vector3<double> frc2=-k*theta_gradx2;
           Vector3<double> frc3=-k*theta_gradx3;
           */
           Vector3<double> frc0=-k*(p_e->GetDihedralAngle()-BaseAngle)*theta_gradx0;
           Vector3<double> frc1=-k*(p_e->GetDihedralAngle()-BaseAngle)*theta_gradx1;
           Vector3<double> frc2=-k*(p_e->GetDihedralAngle()-BaseAngle)*theta_gradx2;
           Vector3<double> frc3=-k*(p_e->GetDihedralAngle()-BaseAngle)*theta_gradx3;

           p_e->p_s->DihedralAngleEnergy[id_tnum]+=0.5*k*(p_e->GetDihedralAngle()-BaseAngle)*(p_e->GetDihedralAngle()-BaseAngle);
           
            /*
           Vector3<double> frc0=-(8*k*(tan((p_e->GetDihedralAngle()-M_PI)/2.0)-tan((BaseAngle-M_PI)/2.0))/(cos((p_e->GetDihedralAngle()-M_PI)/2.0)*cos((p_e->GetDihedralAngle()-M_PI)/2.0)))*theta_gradx0;
           Vector3<double> frc1=-(8*k*(tan((p_e->GetDihedralAngle()-M_PI)/2.0)-tan((BaseAngle-M_PI)/2.0))/(cos((p_e->GetDihedralAngle()-M_PI)/2.0)*cos((p_e->GetDihedralAngle()-M_PI)/2.0)))*theta_gradx1;
           Vector3<double> frc2=-(8*k*(tan((p_e->GetDihedralAngle()-M_PI)/2.0)-tan((BaseAngle-M_PI)/2.0))/(cos((p_e->GetDihedralAngle()-M_PI)/2.0)*cos((p_e->GetDihedralAngle()-M_PI)/2.0)))*theta_gradx2;
           Vector3<double> frc3=-(8*k*(tan((p_e->GetDihedralAngle()-M_PI)/2.0)-tan((BaseAngle-M_PI)/2.0))/(cos((p_e->GetDihedralAngle()-M_PI)/2.0)*cos((p_e->GetDihedralAngle()-M_PI)/2.0)))*theta_gradx3;
            */
            //if(e0.CalcNorm()*(n1*n1)>1e-4&&e0.CalcNorm()*(n2*n2)>1e-4){
                p_e->p_s->GetVertex(p_e->GetVertex(0))->addFrcThread(frc0,id_tnum);
                p_e->p_s->GetVertex(p_e->GetVertex(1))->addFrcThread(frc1,id_tnum);
                v2->addFrcThread(frc2,id_tnum);
                v3->addFrcThread(frc3,id_tnum);
            //}
            /*
            double r=(DELTA_TIME/VISCOSITY);
            if(frc0.CalcNorm()*r>1||frc1.CalcNorm()*r>1||frc2.CalcNorm()*r>1||frc3.CalcNorm()*r>1){
                std::ofstream fout("d_Force.dat",std::ios::app);
                fout<<"DihedralAngleForce\n";
                fout<<"Edge-"<<p_e->GetVertex(0)<<" and "<<p_e->GetVertex(1)<<"\n";
                fout<<"frc0 = "<<frc0<<"\n";
                fout<<"frc1 = "<<frc1<<"\n";
                fout<<"frc2 = "<<frc2<<"\n";
                fout<<"frc3 = "<<frc3<<"\n";
                fout.close();
               //char lname[100];
               //sprintf(lname,"%d-%d",p_e->GetVertex(0),p_e->GetVertex(1));
               //p_e->p_s->OutputOFF(0,lname); 
            }
            */
        }
    }
    void CalcCurvatureForce(Shape* p_s,int deg){//作用反作用が成り立つ
#ifdef _OPENMP
#pragma omp parallel for num_threads(NUM_THREAD)
#endif
        for(int m=0;m<p_s->GetEdge().size();m++){
            Edge* p_e=p_s->GetEdge(m);
            if(p_e->GetFace().size()==1)continue;
            if(p_s->GetVertex(p_e->GetVertex(0))->GetBoundaryFlag()==1||p_s->GetVertex(p_e->GetVertex(1))->GetBoundaryFlag()==1)continue;
            int id_tnum = omp_get_thread_num();//並列処理用
            Vector3<double> line_vec=p_s->periodize_vector(p_s->GetVertex(p_e->GetVertex(0))->GetLoc(deg)-p_s->GetVertex(p_e->GetVertex(1))->GetLoc(deg));
            //v0-v1
            Vector3<double>l_grad=line_vec/(line_vec.CalcNorm()+std::pow(10.0,-5));
            Vector3<double> v=-Bending_Elastic_Constant*(p_e->GetDihedralAngle()-M_PI)*(p_e->p_s->GetVertex(p_e->GetVertex(0))->GetMeanCurvature()*p_s->GetVertex(p_e->GetVertex(0))->GetVoronoiArea()+p_s->GetVertex(p_e->GetVertex(1))->GetMeanCurvature()*p_s->GetVertex(p_e->GetVertex(1))->GetVoronoiArea())*l_grad;
            p_s->GetVertex(p_e->GetVertex(0))->addFrcThread(v,id_tnum);
            p_s->GetVertex(p_e->GetVertex(1))->addFrcThread(-v,id_tnum);


            Vector3<double> x0=p_e->p_s->GetVertex(p_e->GetVertex(0))->GetLoc(deg);
            Vector3<double> x1=p_e->p_s->GetVertex(p_e->GetVertex(1))->GetLoc(deg);
            Vertex* v2;
            Vertex* v3;
            Vector3<double>x2;
            Vector3<double>x3;
            for(int i=0;i<p_e->GetFace().size();i++){
                for(int k=0;k<3;k++){
                    if(p_e->p_s->GetFace(p_e->GetFace(i))->GetVertex(k)==p_e->GetVertex(0) && p_e->p_s->GetFace(p_e->GetFace(i))->GetVertex((k+1)%3)==p_e->GetVertex(1)){
                        v2=p_e->p_s->GetVertex(p_e->p_s->GetFace(p_e->GetFace(i))->GetVertex((k+2)%3));
                        x2=v2->GetLoc(deg);
                        break;
                    }
                    else if(p_e->p_s->GetFace(p_e->GetFace(i))->GetVertex(k)==p_e->GetVertex(0) && p_e->p_s->GetFace(p_e->GetFace(i))->GetVertex((k+2)%3)==p_e->GetVertex(1)){
                        v3=p_e->p_s->GetVertex(p_e->p_s->GetFace(p_e->GetFace(i))->GetVertex((k+1)%3));
                        x3=v3->GetLoc(deg);
                        break;
                    }
                }
            }
            Vector3<double> e0=p_e->p_s->periodize_vector(x1-x0);
            Vector3<double> e1=p_e->p_s->periodize_vector(x2-x0);
            Vector3<double> e2=p_e->p_s->periodize_vector(x3-x0);
            Vector3<double> n1=e0%e1;
            Vector3<double> n2=e2%e0;

            Vector3<double> theta_gradx[4];

            theta_gradx[0]=-(p_e->p_s->periodize_vector(x2-x1)*e0)*n1/(e0.CalcNorm()*(n1*n1)+std::pow(10.0,-7))-(p_e->p_s->periodize_vector(x3-x1)*e0)*n2/(e0.CalcNorm()*(n2*n2)+std::pow(10.0,-7));
            theta_gradx[1]=(e1*e0)*n1/(e0.CalcNorm()*(n1*n1)+std::pow(10.0,-7))+(e2*e0)*n2/(e0.CalcNorm()*(n2*n2)+std::pow(10.0,-7));
            theta_gradx[2]=-(e0.CalcNorm())*n1/((n1*n1)+std::pow(10.0,-7));
            theta_gradx[3]=-(e0.CalcNorm())*n2/((n2*n2)+std::pow(10.0,-7));
            Vector3<double> frc[4];
            
            for(int j=0;j<4;j++){  
                frc[j]=-Bending_Elastic_Constant*theta_gradx[j]*p_e->GetLength()*(p_e->p_s->GetVertex(p_e->GetVertex(0))->GetMeanCurvature()*p_s->GetVertex(p_e->GetVertex(0))->GetVoronoiArea()+p_s->GetVertex(p_e->GetVertex(1))->GetMeanCurvature()*p_s->GetVertex(p_e->GetVertex(1))->GetVoronoiArea());
            }
            
            p_e->p_s->GetVertex(p_e->GetVertex(0))->addFrcThread(frc[0],id_tnum);
            p_e->p_s->GetVertex(p_e->GetVertex(1))->addFrcThread(frc[1],id_tnum);
            v2->addFrcThread(frc[2],id_tnum);
            v3->addFrcThread(frc[3],id_tnum);
        }
#ifdef _OPENMP
#pragma omp parallel for num_threads(NUM_THREAD)
#endif
        for(int i=0;i<p_s->GetFace().size();i++){
            Face* p_f=p_s->GetFace(i);
            if(p_s->GetVertex(p_f->GetVertex(0))->GetBoundaryFlag()==1||p_s->GetVertex(p_f->GetVertex(1))->GetBoundaryFlag()==1||p_s->GetVertex(p_f->GetVertex(2))->GetBoundaryFlag()==1)continue;
            int id_tnum = (int)omp_get_thread_num();//並列処理用
            Vector3<double> a=p_f->p_s->periodize_vector(p_f->p_s->GetVertex(p_f->GetVertex(1))->GetLoc(deg)-p_f->p_s->GetVertex(p_f->GetVertex(0))->GetLoc(deg));
            Vector3<double> b=p_f->p_s->periodize_vector(p_f->p_s->GetVertex(p_f->GetVertex(2))->GetLoc(deg)-p_f->p_s->GetVertex(p_f->GetVertex(0))->GetLoc(deg));

            //x0について
            Vector3<double>s_gradx[3];
            Vector3<double>frc_tmpx[3];
            s_gradx[0]=0.5*p_f->GetNormalVector()%p_f->p_s->periodize_vector(b-a);
            s_gradx[1]=-0.5*p_f->GetNormalVector()%b;
            s_gradx[2]=0.5*p_f->GetNormalVector()%a;
            for(int j=0;j<3;j++){
                frc_tmpx[j]=2.0*Bending_Elastic_Constant*(p_s->GetVertex(p_f->GetVertex(0))->GetMeanCurvature()*p_s->GetVertex(p_f->GetVertex(0))->GetMeanCurvature()+p_s->GetVertex(p_f->GetVertex(1))->GetMeanCurvature()*p_s->GetVertex(p_f->GetVertex(1))->GetMeanCurvature()+p_s->GetVertex(p_f->GetVertex(2))->GetMeanCurvature()*p_s->GetVertex(p_f->GetVertex(2))->GetMeanCurvature())*s_gradx[j]/3.0;
                p_s->GetVertex(p_f->GetVertex(j))->addFrcThread(frc_tmpx[j],id_tnum);
            }
        }
    }

    void CalcCurvatureForce2(Shape* p_s,int deg){
#ifdef _OPENMP
#pragma omp parallel for num_threads(NUM_THREAD)
#endif
        for(int m=0;m<p_s->GetEdge().size();m++){
            if(p_s->GetEdge(m)->GetFace().size()==1)continue;
            Edge* p_e=p_s->GetEdge(m);
            if(p_s->GetVertex(p_e->GetVertex(0))->GetBoundaryFlag()==1||p_s->GetVertex(p_e->GetVertex(1))->GetBoundaryFlag()==1)continue;
            int id_tnum = omp_get_thread_num();//並列処理用
            int id_i=p_e->GetVertex(0);
            int id_j=p_e->GetVertex(1);
            int id_k,id_l;
            
            for(int i=0;i<p_e->GetFace().size();i++){
                for(int k=0;k<3;k++){
                    if(p_e->p_s->GetFace(p_e->GetFace(i))->GetVertex(k)==id_i && p_e->p_s->GetFace(p_e->GetFace(i))->GetVertex((k+1)%3)==id_j){
                        id_l=p_e->p_s->GetFace(p_e->GetFace(i))->GetVertex((k+2)%3);
                        break;
                    }
                    else if(p_e->p_s->GetFace(p_e->GetFace(i))->GetVertex(k)==id_i && p_e->p_s->GetFace(p_e->GetFace(i))->GetVertex((k+2)%3)==id_j){
                        id_k=p_e->p_s->GetFace(p_e->GetFace(i))->GetVertex((k+1)%3);
                        break;
                    }
                }
            }
            Vector3<double>r_ik=p_s->periodize_vector(p_s->GetVertex(id_i)->GetLoc(deg)-p_s->GetVertex(id_k)->GetLoc(deg));
            Vector3<double>r_jk=p_s->periodize_vector(p_s->GetVertex(id_j)->GetLoc(deg)-p_s->GetVertex(id_k)->GetLoc(deg));
            Vector3<double>r_il=p_s->periodize_vector(p_s->GetVertex(id_i)->GetLoc(deg)-p_s->GetVertex(id_l)->GetLoc(deg));
            Vector3<double>r_jl=p_s->periodize_vector(p_s->GetVertex(id_j)->GetLoc(deg)-p_s->GetVertex(id_l)->GetLoc(deg));
            Vector3<double> r_ij=p_s->periodize_vector(p_s->GetVertex(id_i)->GetLoc(deg)-p_s->GetVertex(id_j)->GetLoc(deg));

            double cos_k=(r_ik*r_jk)/(r_ik.CalcNorm()*r_jk.CalcNorm());
            if(abs(cos_k>1))cos_k/=abs(cos_k);
            double sin_k=sqrt(1-cos_k*cos_k);
            double cot_k=cos_k/(sin_k+1e-4);
            double cos_l=(r_jl*r_il)/(r_il.CalcNorm()*r_jl.CalcNorm());
            if(abs(cos_l>1))cos_l/=abs(cos_l);
            double sin_l=sqrt(1-cos_l*cos_l);
            double cot_l=cos_l/(sin_l+1e-4);

            double k=Bending_Elastic_Constant;
            Vector3<double> fi2=k*(p_s->GetVertex(id_i)->GetMeanCurvature()*p_s->GetVertex(id_i)->GetNormalVector()-p_s->GetVertex(id_j)->GetMeanCurvature()*p_s->GetVertex(id_j)->GetNormalVector())*(cot_k+cot_l)+(k/2.0)*(cot_k+cot_l)*(r_ij)*(p_s->GetVertex(id_i)->GetMeanCurvature()*p_s->GetVertex(id_i)->GetMeanCurvature()+p_s->GetVertex(id_j)->GetMeanCurvature()*p_s->GetVertex(id_j)->GetMeanCurvature());

            p_s->GetVertex(id_i)->addFrcThread(fi2,id_tnum);
            p_s->GetVertex(id_j)->addFrcThread(-fi2,id_tnum);

            Vector3<double> grad_k_i=r_ik*cot_k/(r_ik*r_ik)-r_jk/(r_jk.CalcNorm()*r_ik.CalcNorm()*sin_k);
            Vector3<double> grad_k_j=r_jk*cot_k/(r_jk*r_jk)-r_ik/(r_ik.CalcNorm()*r_jk.CalcNorm()*sin_k);
            Vector3<double> grad_k_k=cot_k*((-r_ik)/(r_ik*r_ik)+(-r_jk)/(r_jk*r_jk))+(r_jk+r_ik)/(r_ik.CalcNorm()*r_jk.CalcNorm()*sin_k);

            Vector3<double> grad_l_i=r_il*cot_l/(r_il*r_il)-r_jl/(r_jl.CalcNorm()*r_il.CalcNorm()*sin_l);
            Vector3<double> grad_l_j=r_jl*cot_l/(r_jl*r_jl)-r_il/(r_il.CalcNorm()*r_jl.CalcNorm()*sin_l);
            Vector3<double> grad_l_l=cot_l*((-r_il)/(r_il*r_il)+(-r_jl)/(r_jl*r_jl))+(r_jl+r_il)/(r_il.CalcNorm()*r_jl.CalcNorm()*sin_l);

            double f1=k*(p_s->GetVertex(id_i)->GetMeanCurvature()*p_s->GetVertex(id_i)->GetNormalVector()-p_s->GetVertex(id_j)->GetMeanCurvature()*p_s->GetVertex(id_j)->GetNormalVector())*(r_ij)+(k/4.0)*(p_s->GetVertex(id_i)->GetMeanCurvature()*p_s->GetVertex(id_i)->GetMeanCurvature()+p_s->GetVertex(id_j)->GetMeanCurvature()*p_s->GetVertex(id_j)->GetMeanCurvature())*(r_ij*r_ij);

            p_s->GetVertex(id_i)->addFrcThread(f1*(-grad_k_i/(sin_k*sin_k)-grad_l_i/(sin_l*sin_l)),id_tnum);
            p_s->GetVertex(id_j)->addFrcThread(f1*(-grad_k_j/(sin_k*sin_k)-grad_l_j/(sin_l*sin_l)),id_tnum);
            p_s->GetVertex(id_k)->addFrcThread(f1*(-grad_k_k/(sin_k*sin_k)),id_tnum);
            p_s->GetVertex(id_l)->addFrcThread(f1*(-grad_l_l/(sin_l*sin_l)),id_tnum);


        }
    }

    void CalcInteriorAngleForce(Shape* p_s,int deg){
        int k_interiorAngle=K_T_INTERIOR_ANGLE;
#ifdef _OPENMP
#pragma omp parallel for num_threads(NUM_THREAD)
#endif
        for(int i=0;i<p_s->GetFace().size();i++){
            int id_tnum = (int)omp_get_thread_num();//並列処理用
            //p_s->GetFace(i)->CalcArea(deg);
            for(int j=0;j<3;j++){
                Vector3<double> a=p_s->periodize_vector(p_s->GetVertex(p_s->GetFace(i)->GetVertex((j+1)%3))->GetLoc(deg)-p_s->GetVertex(p_s->GetFace(i)->GetVertex(j))->GetLoc(deg));
                Vector3<double> b=p_s->periodize_vector(p_s->GetVertex(p_s->GetFace(i)->GetVertex((j+2)%3))->GetLoc(deg)-p_s->GetVertex(p_s->GetFace(i)->GetVertex(j))->GetLoc(deg));

                double cos_interiorAngle=(a*b)/(a.CalcNorm()*b.CalcNorm());
                //勾配
                Vector3<double> grad_j0=2*p_s->GetFace(i)->GetArea()*a/((a*b)*(a*a))+2*p_s->GetFace(i)->GetArea()*b/((b*b)*(a*b))+((a%b)%(p_s->periodize_vector(b-a)))/(2*(a*b)*p_s->GetFace(i)->GetArea());
                Vector3<double> grad_j1=-2*p_s->GetFace(i)->GetArea()*a/((a*b)*(a*a))+((a%b)%(-b))/(2*(a*b)*p_s->GetFace(i)->GetArea());
                Vector3<double> grad_j2=-2*p_s->GetFace(i)->GetArea()*b/((b*b)*(a*b))+((a%b)%a)/(2*(a*b)*p_s->GetFace(i)->GetArea());
                //力
                p_s->GetVertex(p_s->GetFace(i)->GetVertex(j))->addFrcThread( k_interiorAngle*(cos_interiorAngle-0.5)*sqrt(1-cos_interiorAngle*cos_interiorAngle)*grad_j0,id_tnum);
                p_s->GetVertex(p_s->GetFace(i)->GetVertex((j+1)%3))->addFrcThread( k_interiorAngle*(cos_interiorAngle-0.5)*sqrt(1-cos_interiorAngle*cos_interiorAngle)*grad_j1,id_tnum);
                p_s->GetVertex(p_s->GetFace(i)->GetVertex((j+2)%3))->addFrcThread( k_interiorAngle*(cos_interiorAngle-0.5)*sqrt(1-cos_interiorAngle*cos_interiorAngle)*grad_j2,id_tnum);
            }
        }
    }
    void CalcZ_axisRestraintForce(Shape* p_s,int deg){
#ifdef _OPENMP
#pragma omp parallel for num_threads(NUM_THREAD)
#endif
        for(int i=0;i<p_s->GetVertex().size();i++){
            int id_tnum = (int)omp_get_thread_num();//並列処理用
            if(KzType==0){
                p_s->GetVertex(i)->CalcVoronoiArea(deg);
                p_s->GetVertex(i)->addFrcThread(Vector3<double>(0.0,0.0,-K_Z*p_s->GetVertex(i)->GetVoronoiArea()*p_s->GetVertex(i)->GetLoc(deg).GetZ()),id_tnum);
                p_s->Z_Energy[id_tnum]+=0.5*K_Z*p_s->GetVertex(i)->GetVoronoiArea()*(p_s->GetVertex(i)->GetLoc(deg).GetZ())*(p_s->GetVertex(i)->GetLoc(deg).GetZ());
            }
            else if(KzType==1){
                if(p_s->GetVertex(i)->GetLoc(deg).GetZ()>0){
                    p_s->GetVertex(i)->CalcVoronoiArea(deg);
                    p_s->GetVertex(i)->addFrcThread(Vector3<double>(0.0,0.0,-K_Z*p_s->GetVertex(i)->GetVoronoiArea()*p_s->GetVertex(i)->GetLoc(deg).GetZ()),id_tnum);
                    p_s->Z_Energy[id_tnum]+=0.5*K_Z*p_s->GetVertex(i)->GetVoronoiArea()*(p_s->GetVertex(i)->GetLoc(deg).GetZ())*(p_s->GetVertex(i)->GetLoc(deg).GetZ());
                }
            }
            else if(KzType==-1){
                if(p_s->GetVertex(i)->GetLoc(deg).GetZ()<0){
                    p_s->GetVertex(i)->CalcVoronoiArea(deg);
                    p_s->GetVertex(i)->addFrcThread(Vector3<double>(0.0,0.0,-K_Z*p_s->GetVertex(i)->GetVoronoiArea()*p_s->GetVertex(i)->GetLoc(deg).GetZ()),id_tnum);
                    p_s->Z_Energy[id_tnum]+=0.5*K_Z*p_s->GetVertex(i)->GetVoronoiArea()*(p_s->GetVertex(i)->GetLoc(deg).GetZ())*(p_s->GetVertex(i)->GetLoc(deg).GetZ());
                }
            }
        }
    }
    void CalcY_axisRestraintForce(Shape* p_s,int deg){
#ifdef _OPENMP
#pragma omp parallel for num_threads(NUM_THREAD)
#endif
        for(int i=0;i<p_s->GetVertex().size();i++){
            int id_tnum = (int)omp_get_thread_num();//並列処理用
                //p_s->GetVertex(i)->CalcVoronoiArea(deg);
                p_s->GetVertex(i)->addFrcThread(Vector3<double>(0.0,-K_Y*p_s->GetVertex(i)->GetLoc(deg).GetY(),0.0),id_tnum);
                p_s->Y_Energy[id_tnum]+=0.5*K_Y*(p_s->GetVertex(i)->GetLoc(deg).GetY())*(p_s->GetVertex(i)->GetLoc(deg).GetY());
        }
    }

    void CalcApicalConstrictionForce(Shape* p_s,int deg){//辺の収縮力
#ifdef _OPENMP
#pragma omp parallel for num_threads(NUM_THREAD)
#endif
        for(int i=0;i<p_s->GetEdge().size();i++){
            int id_tnum = (int)omp_get_thread_num();//並列処理用
            if(p_s->GetEdge(i)->GetFace().size()==2){
                //Line型
                if(p_s->GetEdge(i)->GetApicalFlag()==2){//Line型
                Edge* p_e1=p_s->GetEdge(i);
                int num1;//辺をはさんだ対角同士の点のインデックス
                int num2;
                for(int fidx:p_e1->GetFace()){
                    for(int j=0;j<3;j++){
                        if(p_s->GetFace(fidx)->GetVertex(j)==p_e1->GetVertex(0)&&p_s->GetFace(fidx)->GetVertex((j+1)%3)==p_e1->GetVertex(1)){
                            num1=p_s->GetFace(fidx)->GetVertex((j+2)%3);
                            break;
                        }
                        else if(p_s->GetFace(fidx)->GetVertex(j)==p_e1->GetVertex(1)&&p_s->GetFace(fidx)->GetVertex((j+1)%3)==p_e1->GetVertex(0)){
                            num2=p_s->GetFace(fidx)->GetVertex((j+2)%3);
                            break;
                        }
                    }
                }

                //voronoi点同士を近づける力
                double length=(p_s->periodize_vector(p_s->GetVertex(num1)->GetLoc(deg)-p_s->GetVertex(num2)->GetLoc(deg))).CalcNorm()/3.0;
                Vector3<double> l_grad=(p_s->periodize_vector(p_s->GetVertex(num1)->GetLoc(deg)-p_s->GetVertex(num2)->GetLoc(deg)))/(3.0*(p_s->periodize_vector(p_s->GetVertex(num1)->GetLoc(deg)-p_s->GetVertex(num2)->GetLoc(deg))).CalcNorm());
                p_s->GetVertex(num1)->addFrcThread(-K_APICAL_DIS_L*length*l_grad,id_tnum);
                p_s->GetVertex(num2)->addFrcThread(K_APICAL_DIS_L*length*l_grad,id_tnum);
                
                p_s->ApicalConstrictionEnergy[id_tnum]+=0.5*K_APICAL_DIS_L*length*length;
                }



                //surface型
                if(p_s->GetVertex(p_s->GetEdge(i)->GetVertex(0))->GetApicalFlag()==1||p_s->GetVertex(p_s->GetEdge(i)->GetVertex(1))->GetApicalFlag()==1){//Surface型
                Edge* p_e1=p_s->GetEdge(i);
                int num1;//辺をはさんだ対角同士の点のインデックス
                int num2;
                for(int fidx:p_e1->GetFace()){
                    for(int j=0;j<3;j++){
                        if(p_s->GetFace(fidx)->GetVertex(j)==p_e1->GetVertex(0)&&p_s->GetFace(fidx)->GetVertex((j+1)%3)==p_e1->GetVertex(1)){
                            num1=p_s->GetFace(fidx)->GetVertex((j+2)%3);
                            break;
                        }
                        else if(p_s->GetFace(fidx)->GetVertex(j)==p_e1->GetVertex(1)&&p_s->GetFace(fidx)->GetVertex((j+1)%3)==p_e1->GetVertex(0)){
                            num2=p_s->GetFace(fidx)->GetVertex((j+2)%3);
                            break;
                        }
                    }
                }
                //voronoi点同士を近づける力
                double length=(p_s->periodize_vector(p_s->GetVertex(num1)->GetLoc(deg)-p_s->GetVertex(num2)->GetLoc(deg))).CalcNorm()/3.0;
                Vector3<double> l_grad=(p_s->periodize_vector(p_s->GetVertex(num1)->GetLoc(deg)-p_s->GetVertex(num2)->GetLoc(deg)))/(3.0*(p_s->periodize_vector(p_s->GetVertex(num1)->GetLoc(deg)-p_s->GetVertex(num2)->GetLoc(deg))).CalcNorm());
                p_s->GetVertex(num1)->addFrcThread(-K_APICAL_DIS_S*length*l_grad,id_tnum);
                p_s->GetVertex(num2)->addFrcThread(K_APICAL_DIS_S*length*l_grad,id_tnum);

                p_s->ApicalConstrictionEnergy[id_tnum]+=0.5*K_APICAL_DIS_S*length*length;
                }
            }
        }
    }
    void CalcRepulsiveForce(Shape* p_s,int deg){//同じレイヤーのみ
        #ifdef _OPENMP
        double start,end;
	    start=omp_get_wtime();
        #endif
        p_s->DivideMeshList(deg);

        #ifdef _OPENMP
        #pragma omp parallel for schedule(guided) num_threads(NUM_THREAD)
        #endif
        for(int i=0;i<p_s->GetEdge().size();i++){//隣り合う点のみ先に計算
            int id_tnum = (int)omp_get_thread_num();//並列処理用
            Vector3<double> r=p_s->periodize_vector(p_s->GetVertex(p_s->GetEdge(i)->GetVertex(0))->GetLoc(deg)-p_s->GetVertex(p_s->GetEdge(i)->GetVertex(1))->GetLoc(deg));
            if(RC>r.CalcNorm()){//rcよりも距離が近いときの判定
                Vector3<double> v=K_REP*((RC-r.CalcNorm())/(RC*RC))*r/r.CalcNorm();
                if(p_s->GetVertex(p_s->GetEdge(i)->GetVertex(0))->GetFixFlag()!=1){//固定するときは要らない
                    p_s->GetVertex(p_s->GetEdge(i)->GetVertex(0))->addFrcThread(v,id_tnum);//足し合わせる
                }
                if(p_s->GetVertex(p_s->GetEdge(i)->GetVertex(1))->GetFixFlag()!=1){//固定するときはいらない
                    p_s->GetVertex(p_s->GetEdge(i)->GetVertex(1))->addFrcThread(-v,id_tnum);
                }
                double r=(DELTA_TIME/VISCOSITY);
                if(v.CalcNorm()*r>1){
                    std::ofstream fout("d_Force.dat",std::ios::app);
                    fout<<"RepulsiveForce\n";
                    fout<<"v-v"<<p_s->GetEdge(i)->GetVertex(0)<<" and "<<p_s->GetEdge(i)->GetVertex(1)<<"\n";
                    fout<<"frc0 = "<<v<<"\n";
                    fout<<"frc1 = "<<-v<<"\n";
                    fout.close();
                }
            }
        }
        //p_s->CalcVertexNormalVector(deg);

        #ifdef _OPENMP
        #pragma omp parallel for schedule(guided) num_threads(NUM_THREAD)
        #endif
        //for(int i:v_list){//隣り合う点を除く
        for(int i=0;i<p_s->GetVertex().size();i++){
            int id_tnum = (int)omp_get_thread_num();//並列処理用
            
            int numx=p_s->GetVertex(i)->GetOffsetLoc().GetX();
            int numy=p_s->GetVertex(i)->GetOffsetLoc().GetY();
            int numz=p_s->GetVertex(i)->GetOffsetLoc().GetZ();
            for(int x=numx-1;x<=numx+1;x++){//27格子のみを考える
                if(x==-1 ||x==p_s->GetMAXX()+1){
                        continue;
                    }
                for(int y=numy-1;y<=numy+1;y++){
                    if(y==-1 ||y==p_s->GetMAXY()+1){
                        continue;
                    }                
                    for(int z=numz-1;z<=numz+1;z++){
                        if(z==-1 ||z==p_s->GetMAXZ()+1){
                            continue;
                        }
                        for(int j:p_s->GetNodeIndex(p_s->GetVertex(i)->GetLayerFlag(),x+y*(p_s->GetMAXX()+1)+z*(p_s->GetMAXX()+1)*(p_s->GetMAXY()+1))){
                            if(j>i){
                                bool found=false;
                                for(int vidx:p_s->GetVertex(i)->GetNextVertex()){//隣にある粒子のときは
                                    if(vidx==j){//隣同士は含めない
                                        found=true;
                                        break;
                                    }
                                }
                                if(found){
                                    continue;
                                }
                                Vector3<double>r=p_s->periodize_vector(p_s->GetVertex(i)->GetLoc(deg)-p_s->GetVertex(j)->GetLoc(deg));
                                if(RC>r.CalcNorm()){//rcよりも距離が近いときの判定
                                    Vector3<double>v=K_REP*((RC-r.CalcNorm())/(RC*RC))*r/r.CalcNorm();
                                    if(p_s->GetVertex(i)->GetFixFlag()!=1){//固定するときは要らない
                                        p_s->GetVertex(i)->addFrcThread(v,id_tnum);//足し合わせる
                                    }
                                    if(p_s->GetVertex(j)->GetFixFlag()!=1){//固定するときはいらない
                                        p_s->GetVertex(j)->addFrcThread(-v,id_tnum);
                                    }
                                    p_s->RepulsiveEnergy[id_tnum]+=0.5*K_REP*((r.CalcNorm()/RC)-1)*((r.CalcNorm()/RC)-1);
                                    /*
                                    double r=(DELTA_TIME/VISCOSITY);
                                    if(v.CalcNorm()*r>1){
                                        std::ofstream fout("d_Force.dat",std::ios::app);
                                        fout<<"RepulsiveForce\n";
                                        fout<<"v-v"<<i<<" and "<<j<<"\n";
                                        fout<<"frc0 = "<<v<<"\n";
                                        fout<<"frc1 = "<<-v<<"\n";
                                        fout.close();
                                    }
                                    */
                                }
                            }
                        }
                    }
                } 
            }
        }

        #ifdef _OPENMP
        end=omp_get_wtime();
        std::ofstream fout("d_Collision_Time.dat",std::ios::app);
	    fout<<end-start<<"\n";
        #endif
    }
    void CalcLateralForce(Shape* p_s,int deg){
        if(Layer_NUM<=1)return;
        int thres=std::ceil(R_UP/RC);//天井関数
        #ifdef _OPENMP
        #pragma omp parallel for schedule(guided) num_threads(NUM_THREAD)
        #endif
        for(int i=0;i<p_s->GetVertex().size();i++){
            int id_tnum = (int)omp_get_thread_num();//並列処理用
        
            int numx=p_s->GetVertex(i)->GetOffsetLoc().GetX();
            int numy=p_s->GetVertex(i)->GetOffsetLoc().GetY();
            int numz=p_s->GetVertex(i)->GetOffsetLoc().GetZ();
            for(int x=numx-thres;x<=numx+thres;x++){//(2thres+1)^3格子のみを考える
                if(x<=-1 ||x>=p_s->GetMAXX()+1){
                        continue;
                    }
                for(int y=numy-thres;y<=numy+thres;y++){
                    if(y<=-1 ||y>=p_s->GetMAXY()+1){
                        continue;
                    }                
                    for(int z=numz-thres;z<=numz+thres;z++){
                        if(z<=-1 ||z>=p_s->GetMAXZ()+1){
                            continue;
                        }
                        for(int j:p_s->GetNodeIndex((p_s->GetVertex(i)->GetLayerFlag()+1)%Layer_NUM,x+y*(p_s->GetMAXX()+1)+z*(p_s->GetMAXX()+1)*(p_s->GetMAXY()+1))){//自分とは別のレイヤー
                            if(j>i){
                                Vector3<double>r=p_s->periodize_vector(p_s->GetVertex(i)->GetLoc(deg)-p_s->GetVertex(j)->GetLoc(deg));
                                if((r*p_s->GetVertex(i)->GetNormalVector()/r.CalcNorm())*(r*p_s->GetVertex(j)->GetNormalVector()/r.CalcNorm())<0.75)continue;
                                //力を加える点を選択する
                                if(R_UP>r.CalcNorm()&&R_LOW<=r.CalcNorm()){
                                    Vector3<double>v=K_lateral*((R_EQ-r.CalcNorm())/(R_EQ*R_EQ))*r/r.CalcNorm();
                                    if(p_s->GetVertex(i)->GetFixFlag()!=1){//固定するときは要らない
                                        p_s->GetVertex(i)->addFrcThread(v,id_tnum);//足し合わせる
                                    }
                                    if(p_s->GetVertex(j)->GetFixFlag()!=1){//固定するときはいらない
                                        p_s->GetVertex(j)->addFrcThread(-v,id_tnum);
                                    }
                                    p_s->LateralEnergy[id_tnum]+=0.5*K_lateral*((r.CalcNorm()/R_EQ)-1)*((r.CalcNorm()/R_EQ)-1);
                                }
                            }
                        }
                    }
                } 
            }        
        }
    }
    void CalcPressure(Shape* p_s,int deg){
#ifdef _OPENMP
#pragma omp parallel for num_threads(NUM_THREAD)
#endif
        for(int i=0;i<p_s->GetVertex().size();i++){
            int id_tnum = (int)omp_get_thread_num();//並列処理用
                p_s->GetVertex(i)->CalcVoronoiArea(deg);
                p_s->GetVertex(i)->addFrcThread(Vector3<double>(p_s->GetVertex(i)->GetVoronoiArea()*p_s->GetVertex(i)->GetNormalVector()*pressure),id_tnum);
        }
    }

    void CalcVoronoiCenterForce(Shape*p_s,int deg){//重心に母点を近づける力
#ifdef _OPENMP
#pragma omp parallel for num_threads(NUM_THREAD)
#endif
        for(int i=0;i<p_s->GetVertex().size();i++){
            int id_tnum = (int)omp_get_thread_num();//並列処理用
            if(p_s->GetVertex(i)->GetBoundaryFlag()==0){
                //p_s->GetVertex(i)->CalcNormalVector(deg);
                Vector3<double> g(0.0,0.0,0.0);
                for(int j=1;j<p_s->GetVertex(i)->GetNextVertex().size();j++){
                    g+=p_s->periodize_vector(p_s->GetVertex(p_s->GetVertex(i)->GetNextVertex(j))->GetLoc(deg)-p_s->GetVertex(p_s->GetVertex(i)->GetNextVertex(0))->GetLoc(deg))/(double)p_s->GetVertex(i)->GetNextVertex().size();
                }
                g=p_s->periodize_vector(g+p_s->GetVertex(p_s->GetVertex(i)->GetNextVertex(0))->GetLoc(deg));
                Vector3<double> n=p_s->GetVertex(i)->GetNormalVector();
                Vector3<double> x=p_s->GetVertex(i)->GetLoc(deg);
                Vector3<double> x_g=p_s->periodize_vector(x-g);
                Vector3<double> x_g_n=(x_g*n)*n;
                Vector3<double> f=-K_V_CENTER*(x_g-x_g_n)-K_V_OUT_CENTER*(x_g_n);
                p_s->GetVertex(i)->addFrcThread(f,id_tnum);
                for(int j:p_s->GetVertex(i)->GetNextVertex()){
                    p_s->GetVertex(j)->addFrcThread(-f/(double)p_s->GetVertex(i)->GetNextVertex().size(),id_tnum);
                }
                p_s->V_Center_Energy[id_tnum]+=0.5*K_V_CENTER*(x_g-x_g_n).CalcSqr()+0.5*K_V_OUT_CENTER*(x_g_n).CalcSqr();
            }
        }
    }
    
    void CalcVoronoiAreaForce(Shape* p_s,int deg){//calccenter不用
        double k=K_V_AREA;
        #ifdef _OPENMP
        #pragma omp parallel for num_threads(NUM_THREAD)
        #endif
        for(int i=0;i<p_s->GetVertex().size();i++){
            int id_tnum = (int)omp_get_thread_num();//並列処理用
            if(p_s->GetVertex(i)->GetBoundaryFlag()==1)continue;
            Vertex* p_v=p_s->GetVertex(i);
            Vector3<double>n=p_v->GetNormalVector();
            //p_v->CalcVoronoiArea(deg);
            
            if(voronoi_area_type==2){//微分項はnの射影ではない
                for(int j=0;j<p_v->GetNextVertex().size();j++){
                    Vector3<double> a=p_s->periodize_vector(p_s->GetFace(p_v->GetFace(j))->GetCenter()-p_s->GetFace(p_v->GetFace((j-2+p_v->GetNextVertex().size())%p_v->GetNextVertex().size()))->GetCenter());
                    Vector3<double> b=p_s->periodize_vector(p_s->GetFace(p_v->GetFace((j+1)%p_v->GetNextVertex().size()))->GetCenter()-p_s->GetFace(p_v->GetFace((j-1+p_v->GetNextVertex().size())%p_v->GetNextVertex().size()))->GetCenter());
                    
                    p_s->GetVertex(p_v->GetNextVertex(j))->addFrcThread((-k*(p_v->GetVoronoiArea()-BaseVoronoiArea)*(a+b)%n)/6.0,id_tnum);
                }
            }
            else if(voronoi_area_type==1){//nの射影だがラフな計算
                for(int j=0;j<p_v->GetNextVertex().size();j++){
                    Vector3<double> c=p_s->periodize_vector(p_s->GetVertex(p_v->GetNextVertex(j))->GetLoc(deg)-p_s->GetVertex(p_v->GetNextVertex((j+1)%p_v->GetNextVertex().size()))->GetLoc(deg));
                    Vector3<double> d=p_s->periodize_vector(p_s->GetVertex(p_v->GetNextVertex((j+1)%p_v->GetNextVertex().size()))->GetLoc(deg)-p_s->GetVertex(i)->GetLoc(deg));
                    Vector3<double> e=p_s->periodize_vector(p_s->GetVertex(i)->GetLoc(deg)-p_s->GetVertex(p_v->GetNextVertex(j))->GetLoc(deg));

                    p_s->GetVertex(i)->addFrcThread(-k*(p_v->GetVoronoiArea()-BaseVoronoiArea)*0.5*(c%n),id_tnum);
                    p_s->GetVertex(p_v->GetNextVertex(j))->addFrcThread(-k*(p_v->GetVoronoiArea()-BaseVoronoiArea)*0.5*(d%n),id_tnum);
                    p_s->GetVertex(p_v->GetNextVertex((j+1)%p_v->GetNextVertex().size()))->addFrcThread(-k*(p_v->GetVoronoiArea()-BaseVoronoiArea)*0.5*(e%n),id_tnum);
                    
                }
            }
            p_s->VoronoiAreaEnergy[id_tnum]+=0.5*k*(p_v->GetVoronoiArea()-BaseVoronoiArea)*(p_v->GetVoronoiArea()-BaseVoronoiArea);
        }
    }
    void CalcVoronoiPerimeterForce(Shape* p_s,int deg){
        double k=K_V_PERIMETER;
        #ifdef _OPENMP
        #pragma omp parallel for num_threads(NUM_THREAD)
        #endif
        for(int i=0;i<p_s->GetVertex().size();i++){//edgeで回すと接していない点同士を出すのがめんどそう
            if(p_s->GetVertex(i)->GetBoundaryFlag()==1)continue;
            int id_tnum = (int)omp_get_thread_num();//並列処理用
            //p_s->GetVertex(i)->CalcVoronoiPerimeter(deg);
            double voronoi_perimeter=p_s->GetVertex(i)->GetVoronoiPerimeter();
            for(int j=0;j<p_s->GetVertex(i)->GetNextVertex().size();j++){
                Vector3<double> lj=p_s->periodize_vector(p_s->GetVertex((j+2)%p_s->GetVertex(i)->GetNextVertex().size())->GetLoc(deg)-p_s->GetVertex(j)->GetLoc(deg));
                Vector3<double> ljj=p_s->periodize_vector(p_s->GetVertex(j)->GetLoc(deg)-p_s->GetVertex((j-2+p_s->GetVertex(i)->GetNextVertex().size())%p_s->GetVertex(i)->GetNextVertex().size())->GetLoc(deg));
                
                p_s->GetVertex(p_s->GetVertex(i)->GetNextVertex(j))->addFrcThread(-k*(voronoi_perimeter-BaseVoronoiPerimeter)*((-lj/(3.0*lj.CalcNorm())+ljj/(3.0*ljj.CalcNorm()))),id_tnum);
            }
            p_s->VoronoiPerimeterEnergy[id_tnum]+=0.5*k*(voronoi_perimeter-BaseVoronoiPerimeter)*(voronoi_perimeter-BaseVoronoiPerimeter);
        }
    }
    
    void Fluctuation(Shape* p_s,int deg){//没
        #ifdef _OPENMP
        #pragma omp parallel for num_threads(NUM_THREAD)
        #endif
        for(int i=0;i<p_s->GetFace().size();i++){
            int id_tnum = (int)omp_get_thread_num();//並列処理用
            for(int j=0;j<3;j++){
                p_s->GetVertex(p_s->GetFace(i)->GetVertex(j))->addFrcThread(p_s->GetFace(i)->GetNormalVector()*p_s->p_rnd->RAND_R(-Fluct, Fluct)*VISCOSITY/(DELTA_TIME*p_s->GetFace(i)->GetArea()),id_tnum);
            }
        }    
    }

    void OMP_Reduction_Frc(Shape* p_s,int deg){
        #ifdef _OPENMP
        #pragma omp parallel for num_threads(NUM_THREAD)
        #endif
        for(int i=0;i<p_s->GetVertex().size();i++){
            p_s->GetVertex(i)->resetForce(deg);
            for(int tnum=0;tnum<NUM_THREAD;tnum++){
                p_s->GetVertex(i)->addForce(p_s->GetVertex(i)->GetFrcThread(tnum),deg);
                //std::cout<<"frc_thread"<<tnum<<"="<<p_s->GetVertex(i)->GetForce(deg)<<std::endl;;
                p_s->GetVertex(i)->resetFrcThread(tnum);
            }
        }      
    }
    Vector3<double> CalcDihedralGlad(Shape* p_s,Vector3<double> x0,Vector3<double> x1,Vector3<double>x2,Vector3<double> x3,int i){//0と1がe,2が左,3が右
        Vector3<double> e0=p_s->periodize_vector(x1-x0);
        Vector3<double> e1=p_s->periodize_vector(x2-x0);
        Vector3<double> e2=p_s->periodize_vector(x3-x0);
        Vector3<double> n1=e0%e1;
        Vector3<double> n2=e2%e0;
        if(i==0){
            return -(p_s->periodize_vector(x2-x1)*e0)*n1/(e0.CalcNorm()*(n1*n1)+std::pow(10.0,-7))-(p_s->periodize_vector(x3-x1)*e0)*n2/(e0.CalcNorm()*(n2*n2)+std::pow(10.0,-7));
        }
        else if(i==1){
            return (e1*e0)*n1/(e0.CalcNorm()*(n1*n1)+std::pow(10.0,-7))+(e2*e0)*n2/(e0.CalcNorm()*(n2*n2)+std::pow(10.0,-7));
        }
        else if(i==2){
            return -(e0.CalcNorm())*n1/((n1*n1)+std::pow(10.0,-7));
        }
        else if(i==3){
            return -(e0.CalcNorm())*n2/((n2*n2)+std::pow(10.0,-7));
        }
        else{
            std::cout<<"int i(CalcDihedralGlad) must be 0~3."<<std::endl;
            exit(1);
        }
    }

    void CalcCurvatureEnergyDifference(Shape* p_s,int deg){

        #ifdef _OPENMP
        #pragma omp parallel for schedule(guided) num_threads(NUM_THREAD)
        #endif
        for(int i=0;i<p_s->GetFace().size();i++){
            p_s->GetFace(i)->CalcNormalVector(deg);
            p_s->GetFace(i)->CalcArea(deg);
        }
        #ifdef _OPENMP
        #pragma omp parallel for schedule(guided) num_threads(NUM_THREAD)
        #endif
        for(int i=0;i<p_s->GetVertex().size();i++){
            p_s->GetVertex(i)->CalcNormalVector2(deg);
        }
        #ifdef _OPENMP
        #pragma omp parallel for schedule(guided) num_threads(NUM_THREAD)
        #endif
        for(int i=0;i<p_s->GetEdge().size();i++){
            p_s->GetEdge(i)->CalcLength(deg);
            if(p_s->GetEdge(i)->GetFace().size()==2)p_s->GetEdge(i)->CalcDihedralAngle(deg);
        }
        #ifdef _OPENMP
        #pragma omp parallel for schedule(guided) num_threads(NUM_THREAD)
        #endif
        for(int i=0;i<p_s->GetFace().size();i++){
            p_s->GetFace(i)->CalcCenter(deg);//重心
        }
        #ifdef _OPENMP
        #pragma omp parallel for schedule(guided) num_threads(NUM_THREAD)
        #endif
        for(int i=0;i<p_s->GetVertex().size();i++){
                p_s->GetVertex(i)->CalcVoronoiArea(deg);
                p_s->GetVertex(i)->CalcMeanCurvature();//julicher
        }
        CalcCurvatureForce(p_s,deg);
        OMP_Reduction_Frc(p_s,deg);

        double dh=1e-5;
        for(int m=0;m<p_s->GetVertex().size();m++){

            if(p_s->GetVertex(m)->GetBoundaryFlag()==1)continue;
            double frc_x=CalcCurvatureEnergyDifference(p_s,m,0,dh);
            double frc_y=CalcCurvatureEnergyDifference(p_s,m,1,dh);
            double frc_z=CalcCurvatureEnergyDifference(p_s,m,2,dh);

            Vector3<double> force1=p_s->GetVertex(m)->GetForce(0);

            std::cout<<force1.GetX()-frc_x<<" "<<force1.GetY()-frc_y<<" "<<force1.GetZ()-frc_z<<std::endl;
        }
    }

    double CalcCurvatureEnergyDifference(Shape* p_s,int m,int n,double dh){
        if(n==0)p_s->GetVertex(m)->SetLoc(Vector3<double>(p_s->GetVertex(m)->GetLoc(0).GetX()+dh,p_s->GetVertex(m)->GetLoc(0).GetY(),p_s->GetVertex(m)->GetLoc(0).GetZ()),0);
        else if(n==1)p_s->GetVertex(m)->SetLoc(Vector3<double>(p_s->GetVertex(m)->GetLoc(0).GetX(),p_s->GetVertex(m)->GetLoc(0).GetY()+dh,p_s->GetVertex(m)->GetLoc(0).GetZ()),0);
        else if(n==2)p_s->GetVertex(m)->SetLoc(Vector3<double>(p_s->GetVertex(m)->GetLoc(0).GetX(),p_s->GetVertex(m)->GetLoc(0).GetY(),p_s->GetVertex(m)->GetLoc(0).GetZ()+dh),0);
            for(int i:p_s->GetVertex(m)->GetFace()){
                p_s->GetFace(i)->CalcArea(0);
                p_s->GetFace(i)->CalcCenter(0);
                p_s->GetFace(i)->CalcNormalVector(0);
            }
            for(int i:p_s->GetVertex(m)->GetEdge()){
                p_s->GetEdge(i)->CalcDihedralAngle(0);
                p_s->GetEdge(i)->CalcLength(0);
            }
            p_s->GetVertex(m)->CalcNormalVector2(0);
            p_s->GetVertex(m)->CalcVoronoiArea(0);
            p_s->GetVertex(m)->CalcMeanCurvature();
            double ene1=2*Bending_Elastic_Constant*p_s->GetVertex(m)->GetMeanCurvature()*p_s->GetVertex(m)->GetMeanCurvature()*p_s->GetVertex(m)->GetVoronoiArea();

        if(n==0)p_s->GetVertex(m)->SetLoc(Vector3<double>(p_s->GetVertex(m)->GetLoc(0).GetX()-2*dh,p_s->GetVertex(m)->GetLoc(0).GetY(),p_s->GetVertex(m)->GetLoc(0).GetZ()),0);
        else if(n==1)p_s->GetVertex(m)->SetLoc(Vector3<double>(p_s->GetVertex(m)->GetLoc(0).GetX(),p_s->GetVertex(m)->GetLoc(0).GetY()-2*dh,p_s->GetVertex(m)->GetLoc(0).GetZ()),0);
        else if(n==2)p_s->GetVertex(m)->SetLoc(Vector3<double>(p_s->GetVertex(m)->GetLoc(0).GetX(),p_s->GetVertex(m)->GetLoc(0).GetY(),p_s->GetVertex(m)->GetLoc(0).GetZ()-2*dh),0);
            for(int i:p_s->GetVertex(m)->GetFace()){
                p_s->GetFace(i)->CalcArea(0);
                p_s->GetFace(i)->CalcCenter(0);
                p_s->GetFace(i)->CalcNormalVector(0);
            }
            for(int i:p_s->GetVertex(m)->GetEdge()){
                p_s->GetEdge(i)->CalcDihedralAngle(0);
                p_s->GetEdge(i)->CalcLength(0);
            }
            p_s->GetVertex(m)->CalcNormalVector2(0);
            p_s->GetVertex(m)->CalcVoronoiArea(0);
            p_s->GetVertex(m)->CalcMeanCurvature();
            double ene2=2*Bending_Elastic_Constant*p_s->GetVertex(m)->GetMeanCurvature()*p_s->GetVertex(m)->GetMeanCurvature()*p_s->GetVertex(m)->GetVoronoiArea();

        if(n==0)p_s->GetVertex(m)->SetLoc(Vector3<double>(p_s->GetVertex(m)->GetLoc(0).GetX()+dh,p_s->GetVertex(m)->GetLoc(0).GetY(),p_s->GetVertex(m)->GetLoc(0).GetZ()),0);
        else if(n==1)p_s->GetVertex(m)->SetLoc(Vector3<double>(p_s->GetVertex(m)->GetLoc(0).GetX(),p_s->GetVertex(m)->GetLoc(0).GetY()+dh,p_s->GetVertex(m)->GetLoc(0).GetZ()),0);
        else if(n==2)p_s->GetVertex(m)->SetLoc(Vector3<double>(p_s->GetVertex(m)->GetLoc(0).GetX(),p_s->GetVertex(m)->GetLoc(0).GetY(),p_s->GetVertex(m)->GetLoc(0).GetZ()+dh),0);

            return (ene2-ene1)/(2.0*dh);
    }
   
}//namespace force
