#include"class.hpp"
#include<cmath>
#include<vector>
#include"vector3.hpp"
#include<iostream>
#include<fstream>
#include<string>
#include<omp.h>
#include"Random.hpp"
#include<algorithm>
#include"restructure.hpp"
#include<set>
#include<utility>
Shape::Shape(Random* p_r){
	p_rnd=p_r;
	std::string OffFileName=OFF_FILE_NAME;
	ReadOffFile(OffFileName);//OFFファイルの読み込み
	if(PERIODIC_MODE==1){
		assert(p_v.size()*2==p_f.size());
		std::cout<<"This is a periodic Mode."<<std::endl;
		//size_system=Vector3<double>(sqrt(p_v.size()),sqrt(p_v.size())*sqrt(3)/2.0,0);
		//size_system=Vector3<double>(400,20*sqrt(3)/2.0,0);
		size_system=Vector3<double>(size_x,size_y,0.0);
		std::cout<<"system_size="<<size_system<<std::endl;
	}
	else if(PERIODIC_MODE==2){
		std::cout<<"This is a Y-axis periodic Mode."<<std::endl;
		size_system=Vector3<double>(0.0,size_y,0.0);
		std::cout<<"system_size="<<size_system<<std::endl;
	}
	else{
		size_system=Vector3<double>(0.0,0.0,0.0);
		std::cout<<"This is not a periodic Mode."<<std::endl;
	}
	EdgeConstruction();//頂点と面の情報から辺を構成
	AssignNextVertex();//頂点に隣り合う頂点の配列を作成
	if(!isConsistent()){
		std::cout<<"Invalid assingment"<<std::endl;
		exit(1);
	}
	SetBaseConstant();//初期値の設定
	OutputVTK(0);//VTKファイルの出力(初期状態)
	SetFixAndBoundary();//境界の設定

	CalcVertexNormalVector(0);//頂点での法線ベクトルを計算
	if(CELL_TIME_INITIALIZE_MODE==0)CellTimeInitialize();//細胞周期を任意に与える
	else if(CELL_TIME_INITIALIZE_MODE==1)CellTimeInitialize(CELL_TIME_FILE_NAME);
	SetLayerFlag();//layer_flagの設定
	sortCounterClockwise();
	//for(int i=0;i<p_v.size();i++){
		//p_v[i]->sortCounterClockwise();
	//}
	flip_count=0;
	for(int i=0;i<9;i++){
		division_count[i]=0;
	}
	dihedral_count=0;
	for(int i=0;i<5;i++){
		internal_angle_error_count[i]=0;
	}
	for(int i=0;i<15;i++){
		TotalEnergy[i]=0.0;
	}

	for(int i=0;i<NUM_THREAD;i++){
		LineEnergy[i]=0.0;
        AreaEnergy[i]=0.0;
        DihedralAngleEnergy[i]=0.0;
        Z_Energy[i]=0.0;
		Y_Energy[i]=0.0;
        V_Center_Energy[i]=0.0;
        RepulsiveEnergy[i]=0.0;
        ApicalConstrictionEnergy[i]=0.0;
		VoronoiAreaEnergy[i]=0.0;
		VoronoiPerimeterEnergy[i]=0.0;
		LateralEnergy[i]=0.0;
	}

}
Shape::~Shape(){
	//メモリの解放
	int n_v=p_v.size();
	int n_e=p_e.size();
	int n_f=p_f.size();
	for(int i=0;i<n_v;i++){
		delete p_v[i];
	}
	for(int i=0;i<n_e;i++){
		delete p_e[i];
	}
	for(int i=0;i<n_f;i++){
		delete p_f[i];
	}
	p_v.clear();
	p_e.clear();
	p_f.clear();
}
void Shape::ReadOffFile(std::string fname){
	std::ifstream fin(fname);
	if(!fin){
		std::cout<<"Error: cannot open the off file."<<std::endl;
		exit(1);
	}
	std::string str;
	int numVertices,numFaces,numEdges;
	fin>>str;
	fin>>numVertices;
	fin>>numFaces;
	fin>>numEdges;//辺は定義されていないのでダミー
	for(int i=0;i<numVertices;i++){
		double tmp_x,tmp_y,tmp_z;
		fin>>tmp_x;fin>>tmp_y;fin>>tmp_z;
		Vertex* tmp_v=new Vertex(tmp_x,tmp_y,tmp_z,this);//そのまま代入される
		p_v.push_back(tmp_v);//shapeクラスのp_vに入れる
		p_v[i]->SetColorFlag(i);
	}
	for(int i=0;i<numFaces;i++){
		int tmp_n;
		int tmp0,tmp1,tmp2;
		fin>>tmp_n;
		fin>>tmp0;
		fin>>tmp1;
		fin>>tmp2;
		Face* tmp_f=new Face(tmp0,tmp1,tmp2,this);//faceのコンストラクタにより、fに対応するviは指定される⑤
		p_f.push_back(tmp_f);//shapeクラスのp_fに入れる
		//頂点を含む面の登録(vに対応するfを指定②)
		p_v[tmp0]->PushFace(i);
		p_v[tmp1]->PushFace(i);
		p_v[tmp2]->PushFace(i);
	}
	fin.close();
	std::cout<<"ReadOffFile finished.\n";
}
void Shape::EdgeConstruction(){
	for(int i=0;i<p_f.size();i++){
		//面についてのループ
		for(int j=0;j<3;j++){//3辺について
			int tmp_vi0=p_f[i]->GetVertex(j);
			int tmp_vi1=p_f[i]->GetVertex((j+1)%3);
			//Vertex* tmp_v0=p_f[i]->GetVertex(j);
			//Vertex* tmp_v0=p_v[p_f[i]->GetVertex(j)];
			//Vertex* tmp_v1=p_v[p_f[i]->GetVertex((j+1)%3)];
			//Vertex* tmp_v1=p_f[i]->GetVertex((j+1)%3);
			
			if(tmp_vi0>tmp_vi1){
				//辺を構成する2頂点はidの若い順になるようにする
				//後で判定するときに楽
				int tmp=tmp_vi0;
				tmp_vi0=tmp_vi1;
				tmp_vi1=tmp;
				/*Vertex* tmp=tmp_v0;
				tmp_v0=tmp_v1;
				tmp_v1=tmp;
				*/
			}

			//tmp_eは一時的に保持しているedge型インスタンス
			Edge* tmp_e=new Edge(tmp_vi0,tmp_vi1,this);
			//edgeのコンストラクタ内でeに対応するviは指定されている③
			
			tmp_e->PushFace(i);//辺が含まれる面を登録④

			//同じものがあるかどうかの判定
			int flag=0;
			for(int k=0;k<p_e.size();k++){
				
				//同じものがあるとき
				if(p_e[k]->GetVertex(0)==tmp_e->GetVertex(0)&&p_e[k]->GetVertex(1)==tmp_e->GetVertex(1)){
					flag=1;
					p_e[k]->PushFace(i);//同じものがすでにあった場合は、辺に接している面を登録④
					p_f[i]->PushEdge(k);
					delete tmp_e;
					break;
				}
				/*
				if(p_e[k]->GetVertex(0)==tmp_e->GetVertex(0)&&p_e[k]->GetVertex(1)==tmp_e->GetVertex(1)){
					flag=1;
					p_e[k]->PushFace(p_f[i]);//同じものがすでにあった場合は、辺に接している面を登録④
					delete tmp_e;
					break;
				}
				*/
			}

			//同じものがないとき
			if(flag!=1){
				p_e.push_back(tmp_e);//shapeクラスのp_eに登録
				//頂点に接している辺の登録①
				p_v[tmp_e->GetVertex(0)]->PushEdge(p_e.size()-1);
				p_v[tmp_e->GetVertex(1)]->PushEdge(p_e.size()-1);
				//面に接している辺の登録
				p_f[i]->PushEdge(p_e.size()-1);					
				/*
				tmp_e->GetVertex(0)->PushEdge(tmp_e);
				tmp_e->GetVertex(1)->PushEdge(tmp_e);
				*/
			}
		}
	}
	for(int i=0;i<p_e.size();i++){
		//p_e[i]->SetId(i);
		if(p_e[i]->GetFace().size()!=1)p_e[i]->CalcDihedralAngle(0);
	}
	std::cout<<"EdgeConstruction finished.\n";
}
void Shape::AssignNextVertex(){
#ifdef _OPENMP
#pragma omp parallel for num_threads(NUM_THREAD)
#endif
	for(int vidx=0;vidx<p_v.size();vidx++){
		//p_v[vidx].clear();//後で使うとき用
		for(int fidx:p_v[vidx]->GetFace()){
			for(int vidy:p_f[fidx]->GetVertex()){
				if(vidy!=vidx){
					p_v[vidx]->PushNextVertex(vidy);
				}
			}
		}
		p_v[vidx]->SortAndErase();
	}
}

void Shape::SetBaseConstant(){
#ifdef _OPENMP
#pragma omp parallel for num_threads(NUM_THREAD)
#endif
	for(int i=0;i<p_f.size();i++){
		p_f[i]->SetBaseArea(BASE_AREA);
		p_f[i]->SetK_Area(K_T_AREA);
	}
#ifdef _OPENMP
#pragma omp parallel for num_threads(NUM_THREAD)
#endif
	for(int i=0;i<p_e.size();i++){
		p_e[i]->SetBaseLength(BASE_LENGTH);
		p_e[i]->SetBaseDihedralAngle(BASE_DIHEDRAL_ANGLE);
		p_e[i]->SetK1_Length(K1_T_LENGTH);
		p_e[i]->SetK2_Length(K2_T_LENGTH);
		p_e[i]->SetK_Theta(K_T_THETA);
	}
}
//新しく作った時用
void Shape::SetBaseConstant(Edge* _p_e){
	_p_e->SetBaseLength(BASE_LENGTH);
	_p_e->SetBaseDihedralAngle(BASE_DIHEDRAL_ANGLE);
	_p_e->SetK1_Length(K1_T_LENGTH);
	_p_e->SetK2_Length(K2_T_LENGTH);
	_p_e->SetK_Theta(K_T_THETA);
}
void Shape::SetBaseConstant(Face* _p_f){
	_p_f->SetK_Area(K_T_AREA);
	_p_f->SetBaseArea(BASE_AREA);
}

//境界での処理
void Shape::SetFixAndBoundary(){
	for(int i=0;i<p_e.size();i++){
		if(p_e[i]->GetFace().size()==1){
			p_v[p_e[i]->GetVertex(0)]->SetBoundaryFlag(1);
			p_v[p_e[i]->GetVertex(1)]->SetBoundaryFlag(1);
			if(BOUNDARY_FIX_MODE==1){
				p_v[p_e[i]->GetVertex(0)]->SetFixFlag(1);
				p_v[p_e[i]->GetVertex(1)]->SetFixFlag(1);
			}
		}
	}
}
//layer_flagを設定
void Shape::SetLayerFlag(){
	if(Layer_NUM==1)return;
	std::set<int> st;
	for(Vertex* _p_v:p_v){
		st.insert((int)_p_v->GetLoc(0).GetZ());
	}
	if(st.size()!=Layer_NUM){
		std::cout<<"layer_flag is invalid."<<std::endl;
		exit(1);
	}
	for(int i=0;i<p_v.size();i++){
		p_v[i]->SetLayerFlag(std::distance(st.begin(),st.find((int)p_v[i]->GetLoc(0).GetZ())));
		if(LayerDivFlag[p_v[i]->GetLayerFlag()]==1)p_v[i]->SetFlagRec();
	}
}
//細胞周期を任意に与える
void Shape::CellTimeInitialize(){
#ifdef _OPENMP
#pragma omp parallel for num_threads(NUM_THREAD)
#endif
	for(int i=0;i<p_v.size();i++){
		p_v[i]->SetCellTime(p_rnd->RAND_R(0,CELL_DIVISION_STEP_AVE*DELTA_TIME));
	}
}
//細胞周期をファイルから読み込む
void Shape::CellTimeInitialize(const char* fname){
	std::ifstream fin(fname);
	if(!fin){
		std::cout<<"Error: cannot open the Initial Time file."<<std::endl;
		exit(1);
	}
	int n;fin>>n;
	if(n!=p_v.size()){
		std::cout<<"write vertex.size() in first line."<<std::endl;
		exit(1);
	}
	for(int i=0;i<n;i++){
		double t;fin>>t;
		p_v[i]->SetCellTime(t);
	}
	fin.close();
}
//並び替え
void Shape::sortCounterClockwise(){
	std::vector<std::vector<int>> tmp_vi(p_v.size());
	std::vector<std::vector<int>> tmp_fi(p_f.size());
	std::vector<std::vector<int>> tmp_ei(p_e.size());
#ifdef _OPENMP
#pragma omp parallel for num_threads(NUM_THREAD)
#endif
	for(int l=0;l<p_v.size();l++){
		if(p_v[l]->GetBoundaryFlag()==0){   
        int curr_vidx=p_v[l]->GetNextVertex(0);//隣の頂点を一つ指定
        tmp_vi[l].push_back(curr_vidx);//push_backする
        for(int i=0;i<p_v[l]->GetEdge().size();i++){//周囲の辺についてのループ
            bool ans =false;
            for(int fidx:p_v[curr_vidx]->GetFace()){//curr_vidx周りの面について
                for(int k=0;k<3;k++){//三角形のどの頂点か?
                    //三角形の点が反時計回りに登録されているという前提
                    if(p_f[fidx]->GetVertex(k)==curr_vidx&&p_f[fidx]->GetVertex((k+2)%3)==l){
                        curr_vidx=p_f[fidx]->GetVertex((k+1)%3);
                        tmp_fi[l].push_back(fidx);//面をpush_back
                        ans=true;
                        if(curr_vidx!=tmp_vi[l][0]){//一周して
                            tmp_vi[l].push_back(curr_vidx);//被らないため
                        }
                        break;//更新したのでfor文を一つ抜ける
                    }
                }
                if(ans){
                break;
                }
            }   
        }
        assert(tmp_fi[l].size()==p_v[l]->GetFace().size());
        assert(tmp_vi[l].size()==p_v[l]->GetEdge().size());

        //tmp_eiへの代入(ローカルのインデックスとグローバルのインデックスを区別すること)
        for(int i=0;i<tmp_vi[l].size();i++){
            for(int j=0;j<p_v[l]->GetEdge().size();j++){
                if(p_e[p_v[l]->GetEdge(j)]->GetVertex(1)==tmp_vi[l][i]||p_e[p_v[l]->GetEdge(j)]->GetVertex(0)==tmp_vi[l][i]){
                    tmp_ei[l].push_back(p_v[l]->GetEdge(j));
                    //std::cout<<i<<" "<<j<<std::endl;
                }
            }
        }  
        assert(tmp_ei[l].size()==p_v[l]->GetEdge().size());

    	}
    	else{
		std::vector<int> list1=p_v[l]->GetFace();
		std::sort(list1.begin(),list1.end());
		bool ans=false;
		for(int fidx:p_v[l]->GetFace()){
			for(int j=0;j<3;j++){
				if(p_f[fidx]->GetVertex(j)==l){
					std::vector<int>list2=p_v[p_f[fidx]->GetVertex((j+1)%3)]->GetFace();
					std::sort(list2.begin(),list2.end());
					std::vector<int> list_inter;
                    std::set_intersection(list1.begin(), list1.end(), list2.begin(), list2.end(), std::back_inserter(list_inter));
					if(list_inter.size()==1){
						tmp_vi[l].push_back(p_f[fidx]->GetVertex((j+1)%3));
						tmp_vi[l].push_back(p_f[fidx]->GetVertex((j+2)%3));
						tmp_fi[l].push_back(fidx);
						ans=true;
					}
					break;
				}
			}
			if(ans){
				break;
			}
		}
		
        int curr_vidx=tmp_vi[l][1];
        while(tmp_vi[l].size()<p_v[l]->GetEdge().size()){//サイズが等しくなるまで
            bool ans=false;
            for(int fidx:p_v[curr_vidx]->GetFace()){
                //curr_vidx周りの面について
                for(int k=0;k<3;k++){
                    //三角形の点が反時計回りに登録されているという前提
                    if(p_f[fidx]->GetVertex(k)==curr_vidx&&p_f[fidx]->GetVertex((k+2)%3)==l){
                        curr_vidx=p_f[fidx]->GetVertex((k+1)%3);
                        if(tmp_fi[l].size()<p_v[l]->GetFace().size()){
                            tmp_fi[l].push_back(fidx);//面をpush_back
                        }
                        ans=true;
                        if(curr_vidx!=tmp_vi[l][0]){//一周して
                            tmp_vi[l].push_back(curr_vidx);//被らないため
                        }
                        break;//更新したのでfor文を一つ抜ける
                    }
                }
                if(ans){
                   break;//さらに抜ける
                }    
            }
        }
        assert(tmp_vi[l].size()==p_v[l]->GetEdge().size());
        assert(tmp_fi[l].size()==p_v[l]->GetFace().size());
        //tmp_eiへの代入
        for(int i=0;i<tmp_vi[l].size();i++){
            for(int j=0;j<p_v[l]->GetEdge().size();j++){
                if(p_e[p_v[l]->GetEdge(j)]->GetVertex(1)==tmp_vi[l][i]||p_e[p_v[l]->GetEdge(j)]->GetVertex(0)==tmp_vi[l][i]){
                    tmp_ei[l].push_back(p_v[l]->GetEdge(j));
                }
            }
        }
        assert(tmp_ei[l].size()==p_v[l]->GetEdge().size());
    	}
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(NUM_THREAD)
#endif
	for(int l=0;l<p_v.size();l++){
    	for(int i=0;i<tmp_vi[l].size();i++){
        	p_v[l]->SetNextVertex(i,tmp_vi[l][i]);
        	p_v[l]->SetEdge(i,tmp_ei[l][i]);
        	p_v[l]->SetFace(i,tmp_fi[l][i]);
    	}
	}
}
//出力
void Shape::OutputVTK(int step){
	CalcCrossBoundaryFlag();
	int num=0;
	for(int i=0;i<p_f.size();i++){
		if(p_f[i]->GetCrossBoundaryFlag()==0){
			num++;
		}
	}

	char fname[100];
	sprintf(fname,"d%09d.vtk",step);
	std::ofstream fout(fname);
	if(!fout){
		std::cout<<"Error:cannot open output vtk file."<<std::endl;
		exit(1);
	}
	fout<<"# vtk DataFile Version 3.0\n";
	fout<<"Motion\n";
	fout<<"ASCII\n";
	fout<<"DATASET UNSTRUCTURED_GRID\n";
	
	//点
	fout<<"POINTS "<<p_v.size()<<" float\n";
	for(int i=0;i<p_v.size();i++){
		fout<<p_v[i]->GetLoc(0)<<"\n";
	}
	//形状
	fout<<"CELLS "<<num<<" "<<4*num<<"\n";
	for(int i=0;i<p_f.size();i++){
		if(p_f[i]->GetCrossBoundaryFlag()==0){
			fout << "3 ";
			fout << p_f[i]->GetVertex(0) << " ";
			fout << p_f[i]->GetVertex(1) << " ";
			fout << p_f[i]->GetVertex(2) << "\n";
		}
	}
	fout<<"CELL_TYPES "<<num<<"\n";
	for(int i=0;i<p_f.size();i++){
		if(p_f[i]->GetCrossBoundaryFlag()==0){
			fout<<"5\n";//三角形
		}
	}
	/*
	fout<<"POINT_DATA "<<p_v.size()<<"\n";
	fout<<"SCALARS c_apical_flag int\n";
	fout<<"LOOKUP_TABLE default\n";
	for(int i=0;i<p_v.size();i++){
		fout<<(int)((-p_v[i]->GetBaseDihedralAngle()+BaseDihedralAngle)/ApicalConstrictionAngle)<<"\n";
	}
	*/
	fout.close();
}
void Shape::OutputVTK_Voronoi(int step){
	CalcCrossBoundaryFlag();
	int v_num=0;
	for(int i=0;i<p_v.size();i++){
		if(p_v[i]->GetCrossBoundaryFlag()==0){
			v_num++;
		}
	}

	char fname[100];
	sprintf(fname,"e%09d.vtk",step);
	std::ofstream fout(fname);
	if(!fout){
		std::cout<<"Error:cannot open output vtk file."<<std::endl;
		exit(1);
	}
	fout<<"# vtk DataFile Version 3.0\n";
	fout<<"Motion\n";
	fout<<"ASCII\n";
	fout<<"DATASET UNSTRUCTURED_GRID\n";
	
	//点
	fout<<"POINTS "<<p_f.size()<<" float\n";
	for(int i=0;i<p_f.size();i++){
		p_f[i]->CalcCenter(0);
		fout<<p_f[i]->GetCenter()<<"\n";
	}
	//形状
	int num=0;//カウント用
	int edge_num=0;
	for(int i=0;i<p_v.size();i++){
		if(p_v[i]->GetCrossBoundaryFlag()==0){
			num++;
			for(int j=0;j<p_v[i]->GetFace().size();j++){
				num++;
			}
		}
	}
	for(int i=0;i<p_e.size();i++){
		if(p_e[i]->GetApicalFlag()!=0&&p_e[i]->GetFace().size()==2){
			num=num+3;
			edge_num++;
		}
	}
	fout<<"CELLS "<<v_num+edge_num<<" "<<num<<"\n";

	for(int i=0;i<p_v.size();i++){
		if(p_v[i]->GetCrossBoundaryFlag()==0){
			//p_v[i]->sortCounterClockwise();
			fout<<p_v[i]->GetFace().size()<<" ";
			for(int j=0;j<p_v[i]->GetFace().size()-1;j++){
				fout<<p_v[i]->GetFace(j)<<" ";
			}
			fout<<p_v[i]->GetFace(p_v[i]->GetFace().size()-1)<<"\n";
		}
	}
	for(int i=0;i<p_e.size();i++){
		if(p_e[i]->GetApicalFlag()!=0&&p_e[i]->GetFace().size()==2){
			fout<<"2 "<<p_e[i]->GetFace(0)<<" "<<p_e[i]->GetFace(1)<<"\n";
		}
	}

	fout<<"CELL_TYPES "<<v_num+edge_num<<"\n";
	for(int i=0;i<p_v.size();i++){
		if(p_v[i]->GetCrossBoundaryFlag()==0){
			fout<<"7\n";//VTK_POLYGON
		}
	}
	for(int i=0;i<p_e.size();i++){
		if(p_e[i]->GetApicalFlag()!=0&&p_e[i]->GetFace().size()==2){
			fout<<"3\n";//VTK_LINE
		}
	}

	fout<<"CELL_DATA "<<v_num+edge_num<<"\n";
	fout<<"SCALARS APICAL_flag int\n";
	fout<<"LOOKUP_TABLE default\n";
	for(int i=0;i<p_v.size();i++){
		if(p_v[i]->GetCrossBoundaryFlag()==0){
			fout<<p_v[i]->GetApicalFlag()<<"\n";
			//fout<<"0\n";
		}
	}
	for(int i=0;i<p_e.size();i++){
		if(p_e[i]->GetApicalFlag()!=0&&p_e[i]->GetFace().size()==2){
			fout<<p_e[i]->GetApicalFlag()<<"\n";//VTK_LINE
		}
	}
	fout<<"SCALARS Color_Flag int\n";
	fout<<"LOOKUP_TABLE default\n";
	for(int i=0;i<p_v.size();i++){
		if(p_v[i]->GetCrossBoundaryFlag()==0){
			fout<<p_v[i]->GetColorFlag()<<"\n";
		}
	}
	for(int i=0;i<p_e.size();i++){
		if(p_e[i]->GetApicalFlag()!=0&&p_e[i]->GetFace().size()==2){
			fout<<0<<"\n";//VTK_LINE
		}
	}
	//曲率を出すとき用
	/*
#ifdef _OPENMP
#pragma omp parallel for num_threads(NUM_THREAD)
#endif
	for(int i=0;i<p_v.size();i++){
		if(p_v[i]->GetBoundaryFlag()==0){
			p_v[i]->CalcMeanCurvature();//Julicher
		}
	}
	//fout<<"CELL_DATA "<<v_num+edge_num<<"\n";
	fout<<"SCALARS Mean_CURVATURE float\n";
	fout<<"LOOKUP_TABLE default\n";
	for(int i=0;i<p_v.size();i++){
		if(p_v[i]->GetCrossBoundaryFlag()==0){
			fout<<p_v[i]->GetMeanCurvature()<<"\n";
		}
	}
	for(int i=0;i<p_e.size();i++){
		if(p_e[i]->GetApicalFlag()!=0&&p_e[i]->GetFace().size()==2){
			fout<<"0\n";//VTK_LINE
		}
	}

#ifdef _OPENMP
#pragma omp parallel for num_threads(NUM_THREAD)
#endif
	for(int i=0;i<p_v.size();i++){
		if(p_v[i]->GetBoundaryFlag()==0){
			p_v[i]->CalcCurvature();
		}
	}

	//fout<<"CELL_DATA "<<v_num+edge_num<<"\n";
	fout<<"SCALARS GAUSSIAN_CURVATURE float\n";
	fout<<"LOOKUP_TABLE default\n";
	for(int i=0;i<p_v.size();i++){
		if(p_v[i]->GetCrossBoundaryFlag()==0){
			fout<<p_v[i]->GetGaussianCurvature()<<"\n";
		}
	}
	for(int i=0;i<p_e.size();i++){
		if(p_e[i]->GetApicalFlag()!=0&&p_e[i]->GetFace().size()==2){
			fout<<"0\n";//VTK_LINE
		}
	}
	//fout<<"CELL_DATA "<<v_num+edge_num<<"\n";
	fout<<"SCALARS MEAN_CURVATURE2 float\n";
	fout<<"LOOKUP_TABLE default\n";
	for(int i=0;i<p_v.size();i++){
		if(p_v[i]->GetCrossBoundaryFlag()==0){
			fout<<p_v[i]->GetMeanCurvature()<<"\n";
		}
	}
	for(int i=0;i<p_e.size();i++){
		if(p_e[i]->GetApicalFlag()!=0&&p_e[i]->GetFace().size()==2){
			fout<<"0\n";//VTK_LINE
		}
	}
	*/

	fout.close();
}
void Shape::OutputOFF(int step){
	CalcCrossBoundaryFlag();
	int num=0;
	for(int i=0;i<p_f.size();i++){
		if(p_f[i]->GetCrossBoundaryFlag()==0){
			num++;
		}
	}

	char fname[100];
	sprintf(fname, "d%09d.off", step);
	std::ofstream fout(fname);
	if(!fout){
		std::cout << "Error: cannot open output off file." << std::endl;
		exit(1);
	}

	fout << "OFF\n";
	fout << p_v.size() << " " << num << " " << 0 << "\n";

	for(int i = 0; i < p_v.size(); i++){
		fout << p_v[i]->GetLoc(0) << "\n";
	}

	for(int i = 0; i < p_f.size(); i++){
		if(p_f[i]->GetCrossBoundaryFlag()==0){
			fout << "3 ";
			fout << p_f[i]->GetVertex(0) << " ";
			fout << p_f[i]->GetVertex(1) << " ";
			fout << p_f[i]->GetVertex(2) << "\n";
		}
	}
	fout.close();
}
void Shape::OutputOFF_periodic(int step){
	char fname[100];
	sprintf(fname, "d%09d_periodic.off", step);
	std::ofstream fout(fname);
	if(!fout){
		std::cout << "Error: cannot open output off file." << std::endl;
		exit(1);
	}

	fout << "OFF\n";
	fout << p_v.size() << " " << p_f.size() << " " << 0 << "\n";

	for(int i = 0; i < p_v.size(); i++){
		fout << p_v[i]->GetLoc(0) << "\n";
	}

	for(int i = 0; i < p_f.size(); i++){
			fout << "3 ";
			fout << p_f[i]->GetVertex(0) << " ";
			fout << p_f[i]->GetVertex(1) << " ";
			fout << p_f[i]->GetVertex(2) << "\n";
	}
	fout.close();
}
void Shape::OutputOFF_Mapping(int step){
	CalcCrossBoundaryFlag();
	int num=0;
	for(int i=0;i<p_f.size();i++){
		if(p_f[i]->GetCrossBoundaryFlag()==0){
			num++;
		}
	}
	for(int i=0;i<p_v.size();i++){
		p_v[i]->CalcLocMapping();
	}

	char fname[100];
	sprintf(fname, "d%09dm.off", step);
	std::ofstream fout(fname);
	if(!fout){
		std::cout << "Error: cannot open output off file." << std::endl;
		exit(1);
	}

	fout << "OFF\n";
	fout << p_v.size() << " " << num << " " << 0 << "\n";

	for(int i = 0; i < p_v.size(); i++){
		fout << p_v[i]->GetLocMapping() << "\n";
	}

	for(int i = 0; i < p_f.size(); i++){
		if(p_f[i]->GetCrossBoundaryFlag()==0){
			fout << "3 ";
			fout << p_f[i]->GetVertex(0) << " ";
			fout << p_f[i]->GetVertex(1) << " ";
			fout << p_f[i]->GetVertex(2) << "\n";
		}
	}
	fout.close();
}
void Shape::OutputOFF(int step,char* str){
	CalcCrossBoundaryFlag();
	int num=0;
	for(int i=0;i<p_f.size();i++){
		if(p_f[i]->GetCrossBoundaryFlag()==0){
			num++;
		}
	}
	char fname[100];
	sprintf(fname, "d%09d_%s.off",step,str);
	std::ofstream fout(fname);
	if(!fout){
		std::cout << "Error: cannot open output off file." << std::endl;
		exit(1);
	}
	fout << "OFF\n";
	fout << p_v.size() << " " << num << " " << 0 << "\n";
	for(int i = 0; i < p_v.size(); i++){
		fout << p_v[i]->GetLoc(0) << "\n";
	}
	for(int i = 0; i < p_f.size(); i++){
		if(p_f[i]->GetCrossBoundaryFlag()==0){
			fout << "3 ";
			fout << p_f[i]->GetVertex(0) << " ";
			fout << p_f[i]->GetVertex(1) << " ";
			fout << p_f[i]->GetVertex(2) << "\n";
		}
	}
	fout.close();
}

void Shape::Outputarray(const std::vector<int> list){//配列の出力をするだけ
	std::ofstream fout("array.dat");
	for(int i=0;i<list.size();i++){
		fout<<i<<" "<<list[i]<<"\n";
	}
	fout.close();
}
void Shape::OutputFlipCount(int step){//Flipのカウント
	std::ofstream fout("d_FlipCount.dat",std::ios::app);
	fout<<(double)step/(double)CELL_DIVISION_STEP_AVE<<" "<<flip_count<<"\n";
	fout.close();
}
void Shape::OutputDivisionCount(int step,int type){//分裂のカウント
	char fname[100];
	sprintf(fname,"d_DivisionCount_type%d.dat",type);
	std::ofstream fout(fname,std::ios::app);
	fout<<(double)step/(double)CELL_DIVISION_STEP_AVE<<" "<<division_count[type]<<"\n";
	fout.close();
}
void Shape::OutputInternalCount(int step,int type){//内角のエラー数
	char fname[100];
	sprintf(fname,"d_InternalCount_type%d.dat",type);
	std::ofstream fout(fname,std::ios::app);
	fout<<(double)step/(double)CELL_DIVISION_STEP_AVE<<" "<<internal_angle_error_count[type]<<"\n";
	fout.close();
}
void Shape::OutputParameters(){
	std::ofstream fout("d_Parameters.dat");
	fout<<"vertex_size="<<p_v.size()<<"\n";
	fout<<"face_size="<<p_f.size()<<"\n";
	fout<<"edge_size="<<p_e.size()<<"\n";
	fout<<"\n";
	
	fout<<"CELL_DivisionStep_AVE="<<CELL_DIVISION_STEP_AVE<<"\n";
	fout<<"CELL_DivisionStep_SD="<<CELL_DIVISION_STEP_SD<<"\n";
	fout<<"DivisionDistributionMode="<<DivisionDistributionMode<<"\n";
	fout<<"DIV_alpha="<<DIV_alpha<<"\n";
	fout<<"DIV_beta="<<DIV_beta<<"\n";

	fout<<"\n";
	fout<<"DegreeAccuracy="<<DEGREE_ACCURACY<<"\n";//計算精度(1次:1,2次:2)
	fout<<"numThread="<<NUM_THREAD<<"\n";//並列計算の使用スレッド数

	fout<<"Viscosity="<<VISCOSITY<<"\n";//粘性係数
	fout<<"BoundaryViscosity="<<BOUNDARY_VISCOSITY<<"\n";//境界での粘性係数
	fout<<"\n";
	fout<<"k_T_area="<<K_T_AREA<<"\n";//面積弾性係数
	fout<<"k_T_interiorAngle="<<K_T_INTERIOR_ANGLE<<"\n";//内角保存
	fout<<"K1_T_LENGTH="<< K1_T_LENGTH<<"\n";//1次に比例する辺の弾性係数
	fout<<"K2_T_LENGTH="<< K2_T_LENGTH<<"\n";//2次に比例する辺の弾性係数
	fout<<"K_T_THETA="<< K_T_THETA<<"\n";//二面角
	fout<<"K_Z="<< K_Z<<"\n";//z方向の変形を抑制
	fout<<"KzType="<<KzType<<"\n";
	fout<<"voronoi_area_type="<<voronoi_area_type<<"\n";//面積の求め方(1:rough,2:strict)
	fout<<"K_APICAL_DIS_S="<< K_APICAL_DIS_S<<"\n";//apical収縮の収縮力(surface)
	fout<<"K_APICAL_DIS_L="<< K_APICAL_DIS_L<<"\n";//apical収縮の収縮力(line)
	fout<<"M_APICAL_THETA_S="<< M_APICAL_THETA_S<<"\n";//apical収縮の際の二面角エネルギーの比例係数の重みづけ(surface)
	fout<<"M_APICAL_THETA_L="<< M_APICAL_THETA_L<<"\n";//apical収縮の際の二面角エネルギーの比例係数の重みづけ(line)
	fout<<"K_REP="<< K_REP<<"\n";//排除効果
	fout<<"rc="<< RC<<"\n";//排除半径

	fout<<"LayetNumber="<<Layer_NUM<<"\n";//層の数
	fout<<"K_lateral="<<K_lateral<<"\n";//比例係数
	fout<<"R_low="<<R_LOW<<"\n";
	fout<<"R_EQ="<<R_EQ<<"\n";
	fout<<"R_UP="<<R_UP<<"\n";


	fout<<"k_v_center="<< K_V_CENTER<<"\n";//重心に寄せるエネルギー
	fout<<"K_V_OUT_CENTER="<<K_V_OUT_CENTER<<"\n";
	fout<<"k_v_area="<< K_V_AREA<<"\n";//細胞の面積保存
	fout<<"k_v_perimeter="<< K_V_PERIMETER<<"\n";//細胞の周長保存

	fout<<"Fluct="<<Fluct<<"\n";

	fout<<"BaseLength="<< BASE_LENGTH<<"\n";//平衡長さ
	fout<<"BaseArea="<< BASE_AREA<<"\n";//平衡面積
	fout<<"BaseDihedralAngle="<< BASE_DIHEDRAL_ANGLE<<"\n";//平衡二面角
	fout<<"ApicalConstrictionAngle_s="<< APICAL_CONSTRICTION_ANGLE_S<<"\n";//apical収縮の際に変化する二面角量(surface型)
	fout<<"ApicalConstrictionAngle_l="<< APICAL_CONSTRICTION_ANGLE_L<<"\n";//apical収縮の際に変化する二面角量(line型)

	fout<<"FLIP_THRESHOLD="<< FLIP_THRESHOLD_LENGTH<<"\n";//flip判定の閾値
	fout<<"FLIP_COOLTIME_STEP="<<FLIP_COOLTIME_STEP<<"\n";//一度flipした後のcool_timeの長さ

	fout<<"PERIODIC_MODE="<<PERIODIC_MODE<<"\n";//periodic:1,non-periodic:0
	fout<<"boundary_fix_mode="<<BOUNDARY_FIX_MODE<<"\n";//境界を固定するかどうか?
	fout<<"Cell_Division_Mode="<<CELL_DIVISION_MODE<<"\n";//0,1,2,3,4
	fout<<"flip_operation_mode="<<FLIP_OPERATION_MODE<<"\n";//1でflipする,0でflipしない
	fout<<"Cell_Time_Mode="<<CELL_TIME_MODE<<"\n";//細胞分裂する細胞の選び方
	fout<<"\n";
	fout.close();
}
void Shape::OutputCellTime(int step){
	char fname[100];
	sprintf(fname,"d_CELL_TIME%d.dat",step);
	std::ofstream fout(fname);
	fout<<p_v.size()<<"\n";
	for(int i=0;i<p_v.size();i++){
		fout<<p_v[i]->GetCellTime()<<"\n";
	}
	fout.close();
}

void Shape::OutputResult(){
	std::ofstream fout("d_Result.dat");
	fout<<"vertex_size="<<p_v.size()<<"\n";
	fout<<"face_size="<<p_f.size()<<"\n";
	fout<<"edge_size="<<p_e.size()<<"\n";
	fout<<"\n";
	fout<<"flip_count="<<flip_count<<"\n";
	fout<<"dihedral_error_count="<<dihedral_count<<"\n";
	fout<<"division_try_count="<<division_count[0]<<"\n";
	fout<<"division_count="<<division_count[1]<<"\n";
	fout<<"NormalVectorError="<<division_count[2]<<"\n";
	fout<<"LongestAxisError="<<division_count[3]<<"\n";
	fout<<"crosspoint_Error="<<division_count[4]<<"\n";
	fout<<"closeError="<<division_count[5]<<"\n";
	fout<<"DihedralAngleError(new triangle)="<<division_count[6]<<"\n";
	fout<<"connectionError="<<division_count[7]<<"\n";
	fout<<"DihedralAngleError(around)="<<division_count[8]<<"\n";
	fout<<"\n";
	fout<<"InteriorAngleCount"<<"\n";
	fout<<"cos5="<<internal_angle_error_count[0]<<"\n";
	fout<<"cos10="<<internal_angle_error_count[1]<<"\n";
	fout<<"cos15="<<internal_angle_error_count[2]<<"\n";
	fout<<"cos20="<<internal_angle_error_count[3]<<"\n";
	fout<<"cos30="<<internal_angle_error_count[4]<<"\n";
	fout<<"\n";
	fout.close();
}
void Shape::OutputEnergy(int step){
	for(int i=0;i<15;i++){
		TotalEnergy[i]=0.0;
	}
	for(int i=0;i<NUM_THREAD;i++){
		TotalEnergy[1]+=LineEnergy[i];
		TotalEnergy[2]+=AreaEnergy[i];
		TotalEnergy[3]+=DihedralAngleEnergy[i];
		TotalEnergy[4]+=Z_Energy[i];
		TotalEnergy[5]+=V_Center_Energy[i];
		TotalEnergy[6]+=RepulsiveEnergy[i];
		TotalEnergy[7]+=ApicalConstrictionEnergy[i];
		TotalEnergy[8]+=VoronoiAreaEnergy[i];
		TotalEnergy[9]+=VoronoiPerimeterEnergy[i];
		TotalEnergy[10]+=Y_Energy[i];
		TotalEnergy[11]+=LateralEnergy[i];

		TotalEnergy[0]+=LineEnergy[i]+AreaEnergy[i]+DihedralAngleEnergy[i]+Z_Energy[i]+V_Center_Energy[i]+RepulsiveEnergy[i]+ApicalConstrictionEnergy[i]+VoronoiAreaEnergy[i]+VoronoiPerimeterEnergy[i]+
		Y_Energy[i]+LateralEnergy[i];
	}
	std::ofstream fout0("d_totalEnergy.dat",std::ios::app);
    fout0<<(double)step/(double)CELL_DIVISION_STEP_AVE<<" "<<TotalEnergy[0]<<"\n";
	fout0.close();
	std::ofstream fout1("d_LineEnergy.dat",std::ios::app);
    fout1<<(double)step/(double)CELL_DIVISION_STEP_AVE<<" "<<TotalEnergy[1]<<"\n";
	fout1.close();
	std::ofstream fout2("d_AreaEnergy.dat",std::ios::app);
    fout2<<(double)step/(double)CELL_DIVISION_STEP_AVE<<" "<<TotalEnergy[2]<<"\n";
	fout2.close();
	std::ofstream fout3("d_DihedralAngleEnergy.dat",std::ios::app);
    fout3<<(double)step/(double)CELL_DIVISION_STEP_AVE<<" "<<TotalEnergy[3]<<"\n";
	fout3.close();
	std::ofstream fout4("d_Z_Energy.dat",std::ios::app);
    fout4<<(double)step/(double)CELL_DIVISION_STEP_AVE<<" "<<TotalEnergy[4]<<"\n";
	fout4.close();
	std::ofstream fout5("d_V_CenterEnergy.dat",std::ios::app);
    fout5<<(double)step/(double)CELL_DIVISION_STEP_AVE<<" "<<TotalEnergy[5]<<"\n";
	fout5.close();
	std::ofstream fout6("d_RepulsiveEnergy.dat",std::ios::app);
    fout6<<(double)step/(double)CELL_DIVISION_STEP_AVE<<" "<<TotalEnergy[6]<<"\n";
	fout6.close();
	std::ofstream fout7("d_ApicalConstrictionEnergy.dat",std::ios::app);
    fout7<<(double)step/(double)CELL_DIVISION_STEP_AVE<<" "<<TotalEnergy[7]<<"\n";
	fout7.close();
	std::ofstream fout8("d_VoronoiAreaEnergy.dat",std::ios::app);
    fout8<<(double)step/(double)CELL_DIVISION_STEP_AVE<<" "<<TotalEnergy[8]<<"\n";
	fout8.close();
	std::ofstream fout9("d_VoronoiPerimeterEnergy.dat",std::ios::app);
    fout9<<(double)step/(double)CELL_DIVISION_STEP_AVE<<" "<<TotalEnergy[9]<<"\n";
	fout9.close();
	std::ofstream fout10("d_Y_Energy.dat",std::ios::app);
    fout10<<(double)step/(double)CELL_DIVISION_STEP_AVE<<" "<<TotalEnergy[10]<<"\n";
	fout10.close();
	std::ofstream fout11("d_LateralEnergy.dat",std::ios::app);
	fout11<<(double)step/(double)CELL_DIVISION_STEP_AVE<<" "<<TotalEnergy[11]<<"\n";
	fout11.close();

}
void Shape::OutputRange(int step){
	std::vector<double> vx(p_v.size());
	std::vector<double> vy(p_v.size());
	std::vector<double> vz(p_v.size());
	for(int j=0;j<p_v.size();j++){
		vx[j]=(p_v[j]->GetLoc(0).GetX());
		vy[j]=(p_v[j]->GetLoc(0).GetY());
		vz[j]=(p_v[j]->GetLoc(0).GetZ());
	}
	double x_range=*max_element(vx.begin(),vx.end())-*min_element(vx.begin(),vx.end());
	double y_range=*max_element(vy.begin(),vy.end())-*min_element(vy.begin(),vy.end());
	double z_range=*max_element(vz.begin(),vz.end())-*min_element(vz.begin(),vz.end());
	std::ofstream fout("d_XRange.dat",std::ios::app);
	fout<<(double)step/(double)CELL_DIVISION_STEP_AVE<<" "<<x_range<<"\n";
	fout.close();
	std::ofstream fout2("d_YRange.dat",std::ios::app);
	fout2<<(double)step/(double)CELL_DIVISION_STEP_AVE<<" "<<y_range<<"\n";
	fout2.close();
	std::ofstream fout3("d_ZRange.dat",std::ios::app);
	fout3<<(double)step/(double)CELL_DIVISION_STEP_AVE<<" "<<z_range<<"\n";
	fout3.close();
}
void Shape::OutputDistortion(int step){
	int x1=0,x2=0,y1=0,y2=0;
	double avgx1=0.0,avgx2=0.0,avgy1=0.0,avgy2=0.0;
	for(int i=0;i<p_v.size();i++){
		if(p_v[i]->id_x==1){
			avgx1+=p_v[i]->GetLoc(0).GetX();
			x1++;
		}
		else if(p_v[i]->id_x==-1){
			avgx2+=p_v[i]->GetLoc(0).GetX();
			x2++;
		}
		if(p_v[i]->id_y==1){
			avgy1+=p_v[i]->GetLoc(0).GetY();
			y1++;
		}
		else if(p_v[i]->id_y==-1){
			avgy2+=p_v[i]->GetLoc(0).GetY();
			y2++;
		}
	}
	avgx1/=(double)x1;
	avgx2/=(double)x2;
	avgy1/=(double)y1;
	avgy2/=(double)y2;
	std::cout<<avgx1<<" "<<avgx2<<" "<<avgy1<<" "<<avgy2<<std::endl;

	std::ofstream fout("d_XDistortion.dat",std::ios::app);
	fout<<step<<" "<<((avgx1-avgx2)-39.0)/39.0<<"\n";
	fout.close();
	std::ofstream fout2("d_YDistortion.dat",std::ios::app);
	fout2<<step<<" "<<((avgy1-avgy2)-33.7746)/33.7746<<"\n";
	fout2.close();
}

void Shape::FlipCount(){
	flip_count++;
}
void Shape::DivisionCount(int num){
	division_count[num]++;
}
void Shape::DihedralCount(){
	dihedral_count++;
}
void Shape::InternalAngleErrorCount(int num){
	internal_angle_error_count[num]++;
}
//ちゃんと登録されているかどうか
bool Shape::isConsistent(){
	bool flag=false;
#ifdef _OPENMP
#pragma omp parallel for num_threads(NUM_THREAD)
#endif
    for(int vidx=0;vidx<p_v.size();vidx++){
        if(p_v[vidx]->GetFace().size()==0){
            std::cout<<"Vertex-"<<vidx<<" has no Face."<<std::endl;
			flag=true;
        }
        for(int lidx:GetVertex(vidx)->GetEdge()){//Edgeの中に自身があるかを確かめる
	        std::vector<int> ll={GetEdge(lidx)->GetVertex(0),GetEdge(lidx)->GetVertex(1)};
            if(find(ll.begin(),ll.end(),vidx)==ll.end()){
                std::cout<<"Edge-"<<lidx<<" should have vertex-"<<vidx<<std::endl;
				flag=true;
            }
        }
        for(int cidx:GetVertex(vidx)->GetFace()){
            std::vector<int>ff=GetFace(cidx)->GetVertex();
            if(find(ff.begin(),ff.end(),vidx)==ff.end()){
                std::cout<<"Face-"<<cidx<<" should have vertex-"<<vidx<<std::endl;
				flag=true;
            }
        }
        for(int vidy:GetVertex(vidx)->GetNextVertex()){
            std::vector<int>vv=GetVertex(vidy)->GetNextVertex();
            if(find(vv.begin(),vv.end(),vidx)==vv.end()){
                std::cout<<"NextVertex-"<<vidy<<" should have vertex-"<<vidx<<std::endl;
				flag=true;
            }
        }
    }
#ifdef _OPENMP
#pragma omp parallel for num_threads(NUM_THREAD)
#endif
    for(int lidx = 0; lidx <GetEdge().size(); ++lidx) {
        if(GetEdge(lidx)->GetFace().size() == 0) continue;
        for(int i = 0; i < 2; ++i) {
            std::vector<int> vv = GetVertex(GetEdge(lidx)->GetVertex(i))->GetEdge();
			if(find(vv.begin(), vv.end(), lidx) == vv.end()) {
                std::cout << "Vertex-" <<GetEdge(lidx)->GetVertex(i)<< " should have line-" << lidx << std::endl;
				flag=true;
            }
        }
	    for(int cidx: GetEdge(lidx)->GetFace()) {
            std::vector<int> cc = GetFace(cidx)->GetEdge();
            if(find(cc.begin(), cc.end(), lidx) == cc.end()) {
                std::cout << "Face-" <<cidx<< " should have line-" <<lidx<< std::endl;
				flag=true;
            }
        }
    }
#ifdef _OPENMP
#pragma omp parallel for num_threads(NUM_THREAD)
#endif
    for(int cidx = 0; cidx < GetFace().size(); ++cidx) {
        for(int vidx: GetFace(cidx)->GetVertex()) {
            std::vector<int> vv = GetVertex(vidx)->GetFace();
            if(find(vv.begin(), vv.end(), cidx) == vv.end()) {
                std::cout << "Vertex-" <<vidx<< " should have Face-" <<cidx<< std::endl;
				flag=true;
            }
        }
        for(int lidx: GetFace(cidx)->GetEdge()) {
            std::vector<int> ll = GetEdge(lidx)->GetFace();
            if(find(ll.begin(), ll.end(), cidx) == ll.end()) {
                std::cout << "Line-" <<lidx<< " should have Face-" <<cidx<< std::endl;
				flag=true;
            }
        }
    }        
    return !flag;
}
//二面角でおかしいところがないか
bool Shape::isSuitable(int num){//ステップ数を引数にする
	for(int i=0;i<p_e.size();i++){
		if(p_e[i]->GetFace().size()==2){
			p_e[i]->CalcDihedralAngle(0);
			if(fabs(p_e[i]->GetDihedralAngle())<M_PI/2.0){//二面角が90°を超えているのはおかしい
				std::cout<<"invalid dihedral angle"<<std::endl;
				std::ofstream fout("d_Dihedral_Angle_Error.dat",std::ios::app);
				fout<<"Step-"<<num<<"\n";
				fout<<"Edge-"<<i<<"\n";
				fout<<"Vertex-"<<p_e[i]->GetVertex(0)<<"and"<<p_e[i]->GetVertex(1)<<"\n";
				fout.close();
				char fname[100]="dihedral_error";
				OutputOFF(num,fname);
				dihedral_count++;
				return false;
			}
		}
	}
	return true;
}

//計算
Vector3<double> Shape::periodize_vector(Vector3<double> v){//[-s/2,s/2]
	Vector3<double> a=v;
	if(PERIODIC_MODE==1){
		if(size_system.CalcNorm()==0){
			std::cout<<"size_system.calcNorm=0"<<std::endl;
			exit(1);
		}
		Vector3<double> b=a;
		if(a.GetX()>(size_system.GetX()/2.0)){
			a.SetX(a.GetX()-size_system.GetX());
		}
		else if(a.GetX()<-size_system.GetX()/2.0){
			//std::cout<<"w"<<std::endl;
			a.SetX(a.GetX()+size_system.GetX());
		}
		if(a.GetY()>size_system.GetY()/2.0){
			a.SetY(a.GetY()-size_system.GetY());
		}
		else if(a.GetY()<-size_system.GetY()/2.0){
			a.SetY(a.GetY()+size_system.GetY());
		}
		
		if(a==b){
			if(fabs(a.GetX())>size_system.GetX()||fabs(a.GetY())>size_system.GetY()){
				std::cout<<"b="<<b<<std::endl;
				std::cout<<"a="<<a<<std::endl;
				exit(1);
			}
		}
	}
	else if(PERIODIC_MODE==2){
		if(size_system.CalcNorm()==0){
			std::cout<<"size_system.calcNorm=0"<<std::endl;
			exit(1);
		}
		Vector3<double> b=a;
		if(a.GetY()>size_system.GetY()/2.0){
			a.SetY(a.GetY()-size_system.GetY());
		}
		else if(a.GetY()<-size_system.GetY()/2.0){
			a.SetY(a.GetY()+size_system.GetY());
		}
		if(a==b){
			if(fabs(a.GetY())>size_system.GetY()){
				std::cout<<"b="<<b<<std::endl;
				std::cout<<"a="<<a<<std::endl;
				exit(1);
			}
		}
	}
	return a;
}
Vector3<double> Shape::periodized_center(std::vector<Vector3<double>> v){
	Vector3<double> s(0.0,0.0,0.0);
	for(int i=1;i<v.size();i++){
		Vector3<double> b=v[i]-v[0];
		b=periodize_vector(b);
		s+=b;
	}
	Vector3<double> g=v[0]+s/(double)v.size();
	g=periodize_vector(g);
	return g;
}
bool Shape::cross_periodic_boundary(Vector3<double> a){
	bool result=false;
	if(PERIODIC_MODE==1){
		if(a.GetX()>size_system.GetX()/2.0)result=true;
		else if(a.GetX()<-size_system.GetX()/2.0)result=true;
		else if(a.GetY()>size_system.GetY()/2.0)result=true;
		else if(a.GetY()<-size_system.GetY()/2.0)result=true;
	}
	if(PERIODIC_MODE==2){
		if(a.GetY()>size_system.GetY()/2.0)result=true;
		else if(a.GetY()<-size_system.GetY()/2.0)result=true;
	}
	return result;
}
void Shape::CalcCrossBoundaryFlag(){
	for(int i=0;i<p_v.size();i++){
		p_v[i]->SetCrossBoundaryFlag(0);
	}
	for(int i=0;i<p_f.size();i++){
		p_f[i]->SetCrossBoundaryFlag(0);
	}
	for(int i=0;i<p_e.size();i++){//辺が跨いでいるとき
		Vector3<double> a=p_v[p_e[i]->GetVertex(0)]->GetLoc(0)-p_v[p_e[i]->GetVertex(1)]->GetLoc(0);
		if(cross_periodic_boundary(a)){
			p_v[p_e[i]->GetVertex(0)]->SetCrossBoundaryFlag(1);
			p_v[p_e[i]->GetVertex(1)]->SetCrossBoundaryFlag(1);
			for(int j:p_e[i]->GetFace()){
				p_f[j]->SetCrossBoundaryFlag(1);
			}
		}
	}
}
void Shape::CalcFaceNormalVector(int n){
#ifdef _OPENMP
#pragma omp parallel for num_threads(NUM_THREAD)
#endif
	for(int i=0;i<p_f.size();i++){
		p_f[i]->CalcNormalVector(n);
	}
}
void Shape::CalcVertexNormalVector(int n){
#ifdef _OPENMP
#pragma omp parallel for num_threads(NUM_THREAD)
#endif
	for(int i=0;i<p_v.size();i++){
		p_v[i]->CalcNormalVector(n);
	}
}
void Shape::CalcFaceCentroid(int n){
#ifdef _OPENMP
#pragma omp parallel for num_threads(NUM_THREAD)
#endif
	for(int i=0;i<p_f.size();i++){
		p_f[i]->CalcCenter(n);
	}
}
/*
void Shape::CalcVoronoiArea(){
	p_s->CalcVertexNormalVector();
	for(int i=0;i<p_v.size();i++){
		p_v[i]->GetFace()
	}
}
*/

//FLIP
void Shape::FlipOperation(){
#ifdef _OPENMP
#pragma omp parallel for schedule(guided) num_threads(NUM_THREAD)
#endif
	for(int i=0;i<p_e.size();i++){//判定は並列計算で
		if(p_e[i]->GetFace().size()==2){
			if(p_e[i]->GetFlipFlag()>0){
				p_e[i]->SetFlipFlag(p_e[i]->GetFlipFlag()-1);
			}
			else{
				if(p_e[i]->FlipJudgement()){
					p_e[i]->SetFlipFlag(-1);
				}
			}
		}
	}

	for(int i=0;i<p_e.size();i++){
		if(p_e[i]->GetFace().size()==2){
			if(p_e[i]->GetFlipFlag()==-1){
				p_e[i]->FlipOperation_forge();
			}
		}
	}
}
//各座標の最小値を計算(self-collision計算用)
void Shape::CalcMin(int deg){//degの変更
	std::vector<double> vx;
	std::vector<double> vy;
	std::vector<double> vz;
	for(int i=0;i<p_v.size();i++){
		vx.push_back(p_v[i]->GetLoc(deg).GetX());
		vy.push_back(p_v[i]->GetLoc(deg).GetY());
		vz.push_back(p_v[i]->GetLoc(deg).GetZ());
	}
	minx=*min_element(vx.begin(),vx.end());
	miny=*min_element(vy.begin(),vy.end());
	minz=*min_element(vz.begin(),vz.end());
	if(std::isnan(minx)||std::isnan(miny)||std::isnan(minz)){
		std::cout<<"Min_loc has nan."<<std::endl;
		std::cout<<"min_loc="<<minx<<" "<<miny<<" "<<minz<<std::endl;
		exit(1);
	}
	std::cout<<"min_loc="<<minx<<" "<<miny<<" "<<minz<<std::endl;
}
//offsetlocの計算(オフセットした後の座標)-self-collision計算用
void Shape::CalcOffsetLoc(int deg){
	std::vector<int>v1;
	std::vector<int>v2;
	std::vector<int>v3;

	for(int i=0;i<p_v.size();i++){
		p_v[i]->CalcOffsetLoc(deg);
		v1.push_back(p_v[i]->GetOffsetLoc().GetX());
		v2.push_back(p_v[i]->GetOffsetLoc().GetY());
		v3.push_back(p_v[i]->GetOffsetLoc().GetZ());
	}
	maxx=*max_element(v1.begin(),v1.end());//最大位置
	maxy=*max_element(v2.begin(),v2.end());
	maxz=*max_element(v3.begin(),v3.end());
	std::cout<<"max="<<maxx<<" "<<maxy<<" "<<maxz<<std::endl;
}
void Shape::DivideMeshList(int deg){
    CalcMin(deg);
    CalcOffsetLoc(deg);
    if(std::isnan(GetMAXX())||std::isnan(GetMAXY())||std::isnan(GetMAXZ())){
        std::cout<<"Max value has Nan parameter."<<std::endl;
        exit(1);
    }
	//node_list.clear();//多分いらない

	int m=maxx+(maxx+1)*maxy+(maxx+1)*(maxy+1)*maxz+1;
	for(int i=0;i<Layer_NUM;i++){
		node_list[i]=(std::vector<_node> (m));//vectorのコンストラクタ宣言で要素数mのvectorを代入
		std::cout<<"node_list["<<i<<"].size()="<<node_list[i].size()<<std::endl;
	}
	std::cout<<"node_list.size()="<<node_list[0].size()<<std::endl;

	//点をそれぞれのvectorに入れる
    for(int n=0;n<p_v.size();n++){
        node_list[p_v[n]->GetLayerFlag()][p_v[n]->GetOffsetLoc().GetX()+p_v[n]->GetOffsetLoc().GetY()*(maxx+1)+p_v[n]->GetOffsetLoc().GetZ()*(maxx+1)*(maxy+1)].index.push_back(n);//インデックスにpush_backする
    }
}
void Shape::ConnectBoundaryVertex(){//2だけ離れた点のみを繋ぐ
	if(BOUNDARY_FIX_MODE==1)return;
	for(int i=0;i<p_v.size();i++){
		if(p_v[i]->GetBoundaryFlag()==1){
			for(int m:p_v[i]->GetEdge()){
				if(p_e[m]->GetFace().size()==1){
					int j;
					if(p_e[m]->GetVertex(0)!=i) j=p_e[m]->GetVertex(0);
					else if(p_e[m]->GetVertex(1)!=i) j=p_e[m]->GetVertex(1);
					else continue;

					for(int n:p_v[j]->GetEdge()){
						if(p_e[n]->GetFace().size()==1&&n!=m){
							int k;
							if(p_e[n]->GetVertex(0)==j)k=p_e[n]->GetVertex(1);
							else if(p_e[n]->GetVertex(1)==j)k=p_e[n]->GetVertex(0);
							else continue ;
							//iとkが繋げられるかどうかの判定
							bool flag=false;
							int q;
							for(int p:p_v[i]->GetEdge()){
								if(std::max(i,k)==std::max(p_e[p]->GetVertex(0),p_e[p]->GetVertex(1))&&std::min(i,k)==std::min(p_e[p]->GetVertex(0),p_e[p]->GetVertex(1))){//元々すでに繋がれていたとき
									q=p;flag=true;
								}
							}
							if(!flag){
								double dis=(periodize_vector(p_v[k]->GetLoc(0)-p_v[i]->GetLoc(0))).CalcNorm();
								if(dis<BASE_LENGTH){//距離が近い
									int e_index=p_e.size();
									int f_index=p_f.size();
									Edge* _p_e=new Edge(i,k,this);
									PushEdge(_p_e);
									SetBaseConstant(p_e[e_index]);
									p_v[k]->PushNextVertex(i);
									p_v[i]->PushNextVertex(k);
									p_v[i]->PushEdge(e_index);
									p_v[k]->PushEdge(e_index);
									Face* _p_f=new Face(j,p_v[j]->GetNextVertex(p_v[j]->GetNextVertex().size()-1),p_v[j]->GetNextVertex(0),this);
									PushFace(_p_f);
									SetBaseConstant(p_f[f_index]);
									p_v[i]->PushFace(f_index);
									p_v[j]->PushFace(f_index);
									p_v[k]->PushFace(f_index);
									p_v[j]->SetBoundaryFlag(0);
									p_e[e_index]->PushFace(f_index);
									p_f[f_index]->PushEdge(e_index);
									
									p_e[m]->PushFace(f_index);
									p_f[f_index]->PushEdge(m);
									p_e[n]->PushFace(f_index);
									p_f[f_index]->PushEdge(n);
									p_v[i]->sortCounterClockwise();
									p_v[j]->sortCounterClockwise();
									p_v[k]->sortCounterClockwise();
								}
							}
							else{//線だけすでに存在しているとき(面だけ追加)
								if(p_e[q]->GetFace().size()==1){
									int f_index=p_f.size();
									Face* _p_f=new Face(j,p_v[j]->GetNextVertex(p_v[j]->GetNextVertex().size()-1),p_v[j]->GetNextVertex(0),this);
									PushFace(_p_f);
									SetBaseConstant(p_f[f_index]);
									p_v[i]->PushFace(f_index);
									p_v[j]->PushFace(f_index);
									p_v[k]->PushFace(f_index);
									p_v[i]->SetBoundaryFlag(0);
									p_v[j]->SetBoundaryFlag(0);
									p_v[k]->SetBoundaryFlag(0);
									p_e[m]->PushFace(f_index);
									p_e[n]->PushFace(f_index);
									p_e[q]->PushFace(f_index);
									p_f[f_index]->PushEdge(m);
									p_f[f_index]->PushEdge(n);
									p_f[f_index]->PushEdge(q);
									p_v[i]->sortCounterClockwise();
									p_v[j]->sortCounterClockwise();
									p_v[k]->sortCounterClockwise();

								}
							}
						}
					}
				}
			}
		}
	}
}
void Shape::EraseBoundaryFaceAndEdge(){
	if(BOUNDARY_FIX_MODE==1)return;//固定境界では不用
	std::vector<int> elist;//消すedgeのindexのリスト
	for(int i=0;i<p_e.size();i++){
		if(p_e[i]->GetFace().size()==1&&p_v[p_e[i]->GetVertex(0)]->GetFace().size()>2&&p_v[p_e[i]->GetVertex(1)]->GetFace().size()>2)elist.push_back(i);//一つの面しか接していない頂点があると消えるとまずい
	}
	std::sort(elist.begin(),elist.end());//昇順にする(あとで消すときのため(sortされてないとインデックスをずらせない))
	for(int i=0;i<elist.size();i++){
		p_e[elist[i]]->CalcLength(0);
		if(p_e[elist[i]]->GetLength()>p_e[elist[i]]->GetBaseLength()*1.5){//平衡長さの1.5倍になったら切れるようにする
			std::cout<<"---------------------"<<std::endl;
			std::cout<<"Erase a Boundary Edge"<<std::endl;
			int fi=p_e[elist[i]]->GetFace(0);
			int vi1=p_e[elist[i]]->GetVertex(0);//境界の点
			int vi2=p_e[elist[i]]->GetVertex(1);
			std::cout<<vi1<<" and "<<vi2<<" was erased"<<std::endl;
			std::cout<<"---------------------"<<std::endl;
			int vin;//新しい境界上の点(vi1でもvi2でもない点)
			for(int j=0;j<3;j++){
				if(p_f[fi]->GetVertex(j)!=vi1&&p_f[fi]->GetVertex(j)!=vi2){
					vin=p_f[fi]->GetVertex(j);
				}
			}
			p_v[vin]->SetBoundaryFlag(1);//新しい境界点に設定
			//辺と面の削除
			findAndErase_p_e(elist[i]);
			for(int k=i+1;k<elist.size();k++){
				int a=elist[k];//消した分だけ詰める
				elist[k]=a-1;
			}
			p_v[vi1]->sortCounterClockwise();
			p_v[vi2]->sortCounterClockwise();
			p_v[vin]->sortCounterClockwise();
		}
	}
	
}
void Shape::findAndErase_p_e(int index){
	//隣り合うもの同士を消す
	p_v[p_e[index]->GetVertex(0)]->findAndErase_vi_next(p_e[index]->GetVertex(1));//v_nextを消す
	p_v[p_e[index]->GetVertex(1)]->findAndErase_vi_next(p_e[index]->GetVertex(0));
	std::vector<int> fi;//消す面のリスト
	p_e.erase(p_e.begin()+index);
	for(int i=0;i<p_f.size();i++){
		//消した線を含む面は消す必要がある
		for(int j=0;j<p_f[i]->GetEdge().size();j++){
			if(p_f[i]->GetEdge(j)==index){
					fi.push_back(i);//消すべき面を保存しておく
			}
		}//この段階で消してしまうとややこしくなる恐れ

		//p_f[i]->SortEdge();//並べ替える
		for(int j=0;j<p_f[i]->GetEdge().size();j++){
			if(p_f[i]->GetEdge(j)>index){//edgeの中でindexよりも大きいものがあるとき
				p_f[i]->Decrement_ei(j);//後のインデックスは一つずれる
			}
		}
	}
	for(int i=0;i<p_v.size();i++){
		for(int j=0;j<p_v[i]->GetEdge().size();j++){
			if(p_v[i]->GetEdge(j)==index){//等しいものがあるときは消す
				p_v[i]->findAndErase_ei(index);
			}
		}
		//p_v[i]->SortEdge();//並び変える
		for(int j=0;j<p_v[i]->GetEdge().size();j++){
			if(p_v[i]->GetEdge(j)>index){
				p_v[i]->Decrement_ei(j);//後のインデックスは一つずれる
			}
		}
	}
	//面を消す
	std::sort(fi.begin(),fi.end());
	if(fi.size()!=0){
		for(int i=0;i<fi.size();i++){
			findAndErase_p_f(fi[i]);
			for(int j=i+1;j<fi.size();j++){
				fi[j]=fi[j]-1;//インデックスがズレていく
			}
		}
	}
}
    	
void Shape::findAndErase_p_f(int index){
	p_f.erase(p_f.begin()+index);//面を消す
	
	for(int i=0;i<p_e.size();i++){	
		for(int j=0;j<p_e[i]->GetFace().size();j++){
			if(p_e[i]->GetFace(j)==index){
				p_e[i]->findAndErase_fi(index);
			}
		}
		//p_e[i]->SortFace();//並び変えないとおかしくなる
		for(int j=0;j<p_e[i]->GetFace().size();j++){
			if(p_e[i]->GetFace(j)>index){
				p_e[i]->Decrement_fi(j);//後のインデックスは一つずれる
			}
		}
	}
	for(int i=0;i<p_v.size();i++){
		for(int j=0;j<p_v[i]->GetFace().size();j++){
			if(p_v[i]->GetFace(j)==index){
				p_v[i]->findAndErase_fi(index);
			}
		}
		//p_v[i]->SortFace();//並び変え
		for(int j=0;j<p_v[i]->GetFace().size();j++){
			if(p_v[i]->GetFace(j)>index){
				p_v[i]->Decrement_fi(j);//後のインデックスは一つずれる
			}
		}
	}
}

void Shape::updateCellTime(){
#ifdef _OPENMP
#pragma omp parallel for num_threads(NUM_THREAD)
#endif
	for(int i=0;i<p_v.size();i++){
		p_v[i]->SetCellTime(p_v[i]->GetCellTime()+DELTA_TIME);
	}
}
/*-------------------------------------------*/
//--------------------------------
//getter
std::vector<Vertex*> Shape::GetVertex()const{
	return p_v;
}
Vertex* Shape::GetVertex(int i)const{
	return p_v[i];
}
std::vector<Edge*>Shape::GetEdge()const{
	return p_e;
}
Edge* Shape::GetEdge(int i)const{
	return p_e[i];
}
std::vector<Face*>Shape::GetFace()const{
	return p_f;
}
Face* Shape::GetFace(int i)const{
	return p_f[i];
}

double Shape::GetMINX()const{
	return minx;
}
double Shape::GetMINY()const{
	return miny;
}
double Shape::GetMINZ()const{
	return minz;
}

int Shape::GetMAXX()const{
	return maxx;
}
int Shape::GetMAXY()const{
	return maxy;
}
int Shape::GetMAXZ()const{
	return maxz;
}
std::vector<int> Shape::GetNodeIndex(int i,int j)const{
	return node_list[i][j].index;
}

Vector3<double> Shape::Get_SIZE_SYSTEM()const{
	return size_system;
}

//push_back用
void Shape::PushVertex(Vertex* _p_v){
	p_v.push_back(_p_v);
}
void Shape::PushEdge(Edge*_p_e){
	p_e.push_back(_p_e);
}
void Shape::PushFace(Face*_p_f){
	p_f.push_back(_p_f);
}

void Shape::CalcProximityPoint(){
#ifdef _OPENMP
#pragma omp parallel for num_threads(NUM_THREAD)
#endif
	for(int i=0;i<p_v.size();i++){
		p_v[i]->CalcNormalVector(0);//法線ベクトルを前もって計算しておく
	}

	//第4近接点まで出したいので,そのためのリストを作成
	double cos[4];//cosを入れる用(合計)
	int n[4];//足し合わせた個数
	for(int i=0;i<4;i++){//初期化
		cos[i]=0;
		n[i]=0;
	}
#ifdef _OPENMP
#pragma omp parallel for num_threads(NUM_THREAD)
#endif
	for(int i=0;i<p_v.size();i++){
		Vector3<double> n1=p_v[i]->GetNormalVector();//iの法線ベクトル
		std::vector<int> v1;//第2近接点のインデックスのリスト
		std::vector<int> v2;//第3近接点のインデックスのリスト
		std::vector<int> v3;//第4近接点のインデックスのリスト

		//1.第２近接点の計算
		for(int j:p_v[i]->GetNextVertex()){//第0近接点周り
			for(int k:p_v[j]->GetNextVertex()){//第1近接点周り
				v1.push_back(k);//kが第2近接点のインデックス
			}
		}
		std::sort(v1.begin(), v1.end());//並べ替え
		v1.erase(std::unique(v1.begin(), v1.end()), v1.end());//重複要素の削除		
		v1.erase(std::remove(v1.begin(), v1.end(), i), v1.end());//iと同じものを消す
		for(int j:p_v[i]->GetNextVertex()){//第1近接点と同じものを消す
			v1.erase(std::remove(v1.begin(), v1.end(), j), v1.end());
		}//第2近接点のリストv1ができた

		//2.第3近接点の計算
		for(int j:v1){//第1近接点周り
			for(int k:p_v[j]->GetNextVertex()){//第2近接点周り
				v2.push_back(k);//kが第3近接点のインデックス
			}
		}
		std::sort(v2.begin(), v2.end());//並べ替え
		v2.erase(std::unique(v2.begin(), v2.end()), v2.end());//重複要素の削除
		v2.erase(std::remove(v2.begin(), v2.end(), i), v2.end());//iと同じものを消す(第０近接点と同じものを消す)
		for(int j:p_v[i]->GetNextVertex()){//第1近接点と同じものを消す
			v2.erase(std::remove(v2.begin(), v2.end(), j), v2.end());
		}
		for(int j:v1){//第２近接点と同じ物を消す
			v2.erase(std::remove(v2.begin(), v2.end(), j), v2.end());
		}//第3近接点のリストv2ができた

		//3.第4近接点の計算
		for(int j:v2){//第2近接点周り
			for(int k:p_v[j]->GetNextVertex()){//第3近接点周り
				v3.push_back(k);//kは第4近接点のインデックス
			}
		}
		std::sort(v3.begin(), v3.end());//並べ替え
		v3.erase(std::unique(v3.begin(), v3.end()), v3.end());//重複要素の削除
		v3.erase(std::remove(v3.begin(), v3.end(), i), v3.end());//iと同じものを消す(第０近接点と同じものを消す)
		for(int j:p_v[i]->GetNextVertex()){//第1近接点と同じものを消す
			v3.erase(std::remove(v3.begin(), v3.end(), j), v3.end());
		}
		for(int j:v1){//第２近接点と同じ物を消す
			v3.erase(std::remove(v3.begin(), v3.end(), j), v3.end());
		}
		for(int j:v2){//第3近接点と同じ物を消す
			v3.erase(std::remove(v3.begin(), v3.end(), j), v3.end());
		}//第4近接点のリストv3ができた

		//第1近接点
		for(int j:p_v[i]->GetNextVertex()){
			#pragma omp atomic
			cos[0]+=p_v[j]->GetNormalVector()*n1;
			#pragma omp atomic
			n[0]++;
		}
		for(int j:v1){
			#pragma omp atomic
			cos[1]+=p_v[j]->GetNormalVector()*n1;
			#pragma omp atomic
			n[1]++;
		}
		for(int j:v2){
			#pragma omp atomic
			cos[2]+=p_v[j]->GetNormalVector()*n1;
			#pragma omp atomic
			n[2]++;
		}
		for(int j:v3){
			#pragma omp atomic
			cos[3]+=p_v[j]->GetNormalVector()*n1;
			#pragma omp atomic
			n[3]++;
		}
	}
	std::ofstream fout("d_Neighbor.dat",std::ios::app);
	fout<<"0 1"<<"\n";
	for(int i=0;i<4;i++){
		cos[i]/=(double)n[i];
		fout<<i+1<<" "<<cos[i]<<"\n";
	}
	fout.close();
	
}

void Shape::Lloyd_Algorithm(){
	for(int i=0;i<p_f.size();i++){
		p_f[i]->CalcCircumCenter(0);
		p_f[i]->CalcNormalVector(0);
	}
#ifdef _OPENMP
#pragma omp parallel for num_threads(NUM_THREAD)
#endif
	for(int i=0;i<p_v.size();i++){
		p_v[i]->CalcNormalVector2(0);//頂点法線を計算
		if(p_v[i]->GetBoundaryFlag()==0){
			Vector3<double> g=Vector3<double>(0.0,0.0,0.0);
			for(int j=1;j<p_v[i]->GetFace().size();j++){
				Vector3<double> b=periodize_vector(p_f[p_v[i]->GetFace(j)]->GetCircumCenter()-p_f[p_v[i]->GetFace(0)]->GetCircumCenter());
				g+=b/(double)p_v[i]->GetFace().size();
			}
			g+=p_f[p_v[i]->GetFace(0)]->GetCircumCenter();
			g=periodize_vector(g);//voronoi点の重心
			//Vector3<double> x=periodize_vector(((periodize_vector(g-p_v[i]->GetLoc(0))).Projection(p_v[i]->GetNormalVector()))+p_v[i]->GetLoc(0));
			Vector3<double> x=g;
			p_v[i]->SetLoc(x,0);
		}
	}
}
void Shape::ForcedDisplacement(){
	#ifdef _OPENMP
	#pragma omp parallel for num_threads(NUM_THREAD)
	#endif
	for(int i=0;i<p_v.size();i++){
		p_v[i]->SetLoc(Vector3<double>(p_v[i]->GetLoc(0).GetX(),alpha*p_v[i]->GetLoc(0).GetY(),p_v[i]->GetLoc(0).GetZ()),0);
	}
}
void Shape::VoronoiCenterEvaluation(int step){
	double total_dis=0.0;
	int number=0;
	#ifdef _OPENMP
	#pragma omp parallel for num_threads(NUM_THREAD)
	#endif
	for(int i=0;i<p_v.size();i++){
		if(p_v[i]->GetBoundaryFlag()==0){
			p_v[i]->CalcNormalVector(0);
			#ifdef _OPENMP
			#pragma omp atomic
			#endif
			number++;
			Vector3<double> g=Vector3<double>(0.0,0.0,0.0);
			for(int j=1;j<p_v[i]->GetNextVertex().size();j++){
				Vector3<double> b=periodize_vector(p_v[p_v[i]->GetNextVertex(j)]->GetLoc(0)-p_v[p_v[i]->GetNextVertex(0)]->GetLoc(0));
				g+=b/(double)p_v[i]->GetNextVertex().size();
			}
			g+=p_v[p_v[i]->GetNextVertex(0)]->GetLoc(0);
			g=periodize_vector(g);
			#ifdef _OPENMP
			#pragma omp atomic
			#endif
			total_dis+=((periodize_vector(p_v[i]->GetLoc(0)-g)).Projection(p_v[i]->GetNormalVector())).CalcNorm();
		}
	}
	total_dis/=(double)number;
	std::ofstream fout("d_VoronoiCenterEvaluation.dat",std::ios::app);
	fout<<(double)step/(double)CELL_DIVISION_STEP_AVE<<" "<<total_dis<<"\n";
	fout.close();
}
void Shape::CellCenterEvaluation(int step){
double total_dis=0.0;
	int number=0;
	#ifdef _OPENMP
	#pragma omp parallel for num_threads(NUM_THREAD)
	#endif
	for(int i=0;i<p_v.size();i++){
		if(p_v[i]->GetBoundaryFlag()==0){
			p_v[i]->CalcNormalVector(0);
			p_v[i]->CalcLocMapping();
			#ifdef _OPENMP
			#pragma omp atomic
			#endif
			number++;
			Vector3<double> g=p_v[i]->GetLoc(0)-p_v[i]->GetLocMapping();
			g=periodize_vector(g);
			#ifdef _OPENMP
			#pragma omp atomic
			#endif
			total_dis+=g.CalcNorm();
		}
	}
	total_dis/=(double)number;
	std::ofstream fout("d_CellCenterEvaluation.dat",std::ios::app);
	fout<<(double)step/(double)CELL_DIVISION_STEP_AVE<<" "<<total_dis<<"\n";
	fout.close();
}

/*------------------Vertex Class---------------------------*/
//コンストラクタ
Vertex::Vertex(double x,double y,double z,Shape* _p_s){
	Vector3<double> v(x,y,z);
	p_s=_p_s;
	FixFlag=0;
	BoundaryFlag=0;
	apical_flag=0;
	colorFlag=0;
	cross_boundary_flag=0;
	flag_rec=0;
	cell_time=0.0;	
	for(int i=0;i<DEGREE_ACCURACY;i++)	loc[i]=v;
	force[0]=Vector3<double>(0.0,0.0,0.0);
	for(int i=0;i<NUM_THREAD;i++){
		frc_thread[i]=Vector3<double>(0.0,0.0,0.0);
	}
	externalforce=Vector3<double>(0.0,0.0,0.0);
	NormalVector=Vector3<double>(0.0,0.0,0.0);
	gaussian_curvature=0.0;
	mean_curvature=0.0;
	voronoi_area=0.0;
	voronoi_perimeter=0.0;
	layer_flag=0;
	id_x=0;
	id_y=0;
	if(x<-19.3)id_x--;
	else if(x>18.7)id_x++;
	if(y>16)id_y++;
	else if(y<-17)id_y--;
}
//計算
//頂点の法線ベクトルの計算
void Vertex::CalcNormalVector(int n){
	Vector3<double> sumNormalVector(0.0, 0.0, 0.0);
	for(int i = 0; i < fi.size(); i++){
		//p_s->GetFace(fi[i])->CalcArea(n);
		/*
		double alpha;
		for(int j=0;j<3;j++){
			if(p_s->GetVertex(p_s->GetFace(fi[i])->GetVertex(j))==this){
				Vector3<double> a=p_s->GetVertex(p_s->GetFace(fi[i])->GetVertex((j+2)%3))->GetLoc(n)-loc[n];
				a=p_s->periodize_vector(a);
				Vector3<double> b=p_s->GetVertex(p_s->GetFace(fi[i])->GetVertex((j+1)%3))->GetLoc(n)-loc[n];
				b=p_s->periodize_vector(b);
				alpha=acos((a*b)/(a.CalcNorm()*b.CalcNorm()));
				break;
			}
		}
		*/
		/*
		p_s->GetFace(fi[i])->CalcCenter(n);
		//Vector3<double> g=p_s->GetFace(fi[i])->GetCenter(n);
		*/
		//double dis=(p_s->periodize_vector(g-loc[n])).CalcNorm();
		p_s->GetFace(fi[i])->CalcNormalVector(n);
		sumNormalVector += p_s->GetFace(fi[i])->GetNormalVector();
		//sumNormalVector += p_s->GetFace(fi[i])->GetNormalVector()*p_s->GetFace(fi[i])->GetArea();
		//sumNormalVector += p_s->GetFace(fi[i])->GetNormalVector()*p_s->GetFace(fi[i])->GetArea()*alpha;
		//sumNormalVector += alpha*p_s->GetFace(fi[i])->GetNormalVector()/dis;
	}
	NormalVector = sumNormalVector/(sumNormalVector.CalcNorm()+std::pow(10.0,-5));
}
void Vertex::CalcNormalVector2(int n){
	Vector3<double> sumNormalVector(0.0, 0.0, 0.0);
	for(int i = 0; i < fi.size(); i++){
		sumNormalVector += p_s->GetFace(fi[i])->GetNormalVector();
	}
	NormalVector = sumNormalVector/(sumNormalVector.CalcNorm()+std::pow(10.0,-5));
}
void Vertex::CalcVoronoiArea(int n){//頂点での法線方向から見た時
	//CalcNormalVector(n);
	double tmp_area=0.0;
	if(voronoi_area_type==2){//面積を厳密に計算する
		if(GetBoundaryFlag()==1)sortCounterClockwise();//念のため並び替え
		/*
		for(int i=0;i<fi.size();i++){
			p_s->GetFace(fi[i])->CalcCenter(n);
		}
		*/
		for(int i=0;i<fi.size();i++){
			if(GetBoundaryFlag()==1&&i==fi.size())break;
			//1周しないとき
			Vector3<double> gi=p_s->periodize_vector(p_s->GetFace(fi[i])->GetCenter()-loc[n]);
			Vector3<double> gi_1=p_s->periodize_vector(p_s->GetFace(fi[(i+1)%fi.size()])->GetCenter()-loc[n]);
			tmp_area+=0.5*(gi%gi_1)*NormalVector;
		}
	}
	else if(voronoi_area_type==1){
		for(int j:fi){
        	tmp_area+=(p_s->GetFace(j)->GetArea()/3.0)*(p_s->GetFace(j)->GetNormalVector()*GetNormalVector());
        }
	}
	voronoi_area=tmp_area;
}
/*
double Vertex::CalcVoronoiArea(Vector3<double> n1,int n2){
	Vector3<double> n=n1/n1.CalcNorm();
	double tmp_area=0.0;
	if(voronoi_area_type==2){
		if(GetBoundaryFlag()==1){//境界のときは念のため並び替え
			sortCounterClockwise();
			for(int i=0;i<fi.size();i++){
				p_s->GetFace(fi[i])->CalcCenter(n2);
			}
			for(int i=0;i<fi.size()-1;i++){
				//1周しないから
				Vector3<double> gi=p_s->periodize_vector(p_s->GetFace(fi[i])->GetCenter()-loc[n2]);
				Vector3<double> gi_1=p_s->periodize_vector(p_s->GetFace(fi[(i+1)%fi.size()])->GetCenter()-loc[n2]);
				tmp_area+=0.5*(gi%gi_1)*n;
			}
		}
		else{
			for(int i=0;i<fi.size();i++){
				p_s->GetFace(fi[i])->CalcCenter(n2);
			}
			for(int i=0;i<fi.size();i++){
				Vector3<double> gi=p_s->periodize_vector(p_s->GetFace(fi[i])->GetCenter()-loc[n2]);
				Vector3<double> gi_1=p_s->periodize_vector(p_s->GetFace(fi[(i+1)%fi.size()])->GetCenter()-loc[n2]);
				tmp_area+=0.5*(gi%gi_1)*n;
			}//areaを計算
		}
	}
	else if(voronoi_area_type==1){
		for(int j:fi){
        	tmp_area+=(p_s->GetFace(j)->GetArea()/3.0)*(p_s->GetFace(j)->GetNormalVector()*n);
        }
	}
	return tmp_area;
}*/


void Vertex::CalcVoronoiPerimeter(int n){//頂点の法線で射影しない
	//CalcNormalVector(n);
	double tmp_perimeter=0.0;
	/*
	for(int i=0;i<fi.size();i++){
		p_s->GetFace(fi[i])->CalcCenter(n);
	}
	*/
	if(GetBoundaryFlag()==1)sortCounterClockwise();//念のため
	for(int i=0;i<fi.size()-1;i++){
		if(GetBoundaryFlag()==1&&i==fi.size())break;
		Vector3<double> gi=p_s->periodize_vector(p_s->GetFace(fi[i])->GetCenter());
		Vector3<double> gi_1=p_s->periodize_vector(p_s->GetFace(fi[(i+1)%fi.size()])->GetCenter());
		//tmp_perimeter+=((p_s->periodize_vector(gi-gi_1)).Projection(NormalVector)).CalcNorm();
		tmp_perimeter+=((p_s->periodize_vector(gi-gi_1))).CalcNorm();
	}//perimeterの計算
	voronoi_perimeter=tmp_perimeter;
}
void Vertex::CalcMeanCurvature(){//juliche's model
	//if(BoundaryFlag==1)return;
	double total_area=0.0;
	double total_curvature=0.0;
	for(int i:fi){
		total_area+=p_s->GetFace(i)->GetArea()/3.0;
	}
	for(int i:ei){
		total_curvature+=p_s->GetEdge(i)->GetLength()*(p_s->GetEdge(i)->GetDihedralAngle()-M_PI);
	}
	//std::cout<<total_curvature<<std::endl;
	total_curvature/=(4*total_area);
	
	mean_curvature=total_curvature;
}
void Vertex::CalcMeanCurvature(int deg){//using Laplace-Beltrami operator 
	double area_mixed=0.0;//voronoi面積
	Vector3<double> k(0.0,0.0,0.0);
	for(int i:fi){
		for(int j=0;j<3;j++){
			if(p_s->GetVertex(p_s->GetFace(i)->GetVertex(j))==this){	
				Vector3<double> o=p_s->periodize_vector(p_s->GetFace(i)->GetCircumCenter()-loc[deg]);
				Vector3<double> a=p_s->GetVertex(p_s->GetFace(i)->GetVertex((j+1)%3))->GetLoc(deg);
				Vector3<double> b=p_s->GetVertex(p_s->GetFace(i)->GetVertex((j+2)%3))->GetLoc(deg);
				Vector3<double> mid_a=p_s->periodize_vector(a-loc[deg]);//中点
				mid_a=p_s->periodize_vector(loc[deg]+mid_a/2.0);
				Vector3<double> mid_b=p_s->periodize_vector(b-loc[deg]);//中点
				mid_b=p_s->periodize_vector(loc[deg]+mid_b/2.0);
				area_mixed+=0.5*(mid_a%o).CalcNorm()+0.5*(mid_b%o).CalcNorm();

				double cos_alpha=p_s->periodize_vector(loc[deg]-a)*p_s->periodize_vector(b-a)/(p_s->periodize_vector(loc[deg]-a).CalcNorm()*p_s->periodize_vector(b-a).CalcNorm());
				double cos_beta=p_s->periodize_vector(loc[deg]-b)*p_s->periodize_vector(a-b)/(p_s->periodize_vector(loc[deg]-b).CalcNorm()*p_s->periodize_vector(a-b).CalcNorm());
				k+=(p_s->periodize_vector(loc[deg]-b)*cos_alpha/sqrt(1-cos_alpha*cos_alpha))+(p_s->periodize_vector(loc[deg]-a)*cos_beta/sqrt(1-cos_beta*cos_beta));
				break;
			}
		}
	}
	k/=-4.0*area_mixed;
	mean_curvature=k*NormalVector;
}
void Vertex::CalcCurvature(){//using Laplace-Beltrami operator 
	for(int i:fi){
		p_s->GetFace(i)->CalcCircumCenter(0);
		CalcNormalVector(0);
	}
	double area_mixed=0.0;//voronoi面積
	double sum_theta=0.0;//内角の和
	Vector3<double> k(0.0,0.0,0.0);
	for(int i:fi){
		for(int j=0;j<3;j++){
			if(p_s->GetVertex(p_s->GetFace(i)->GetVertex(j))==this){	
				Vector3<double> o=p_s->periodize_vector(p_s->GetFace(i)->GetCircumCenter()-loc[0]);
				Vector3<double> a=p_s->GetVertex(p_s->GetFace(i)->GetVertex((j+1)%3))->GetLoc(0);
				Vector3<double> b=p_s->GetVertex(p_s->GetFace(i)->GetVertex((j+2)%3))->GetLoc(0);
				Vector3<double> mid_a=p_s->periodize_vector(a-loc[0]);//中点
				mid_a=p_s->periodize_vector(loc[0]+mid_a/2.0);
				Vector3<double> mid_b=p_s->periodize_vector(b-loc[0]);//中点
				mid_b=p_s->periodize_vector(loc[0]+mid_b/2.0);
				area_mixed+=0.5*(mid_a%o).CalcNorm()+0.5*(mid_b%o).CalcNorm();

				sum_theta+=acos((p_s->periodize_vector(mid_a-loc[0])*p_s->periodize_vector(mid_b-loc[0]))/(p_s->periodize_vector(mid_a-loc[0]).CalcNorm()*p_s->periodize_vector(mid_b-loc[0]).CalcNorm()));

				double cos_alpha=p_s->periodize_vector(loc[0]-a)*p_s->periodize_vector(b-a)/(p_s->periodize_vector(loc[0]-a).CalcNorm()*p_s->periodize_vector(b-a).CalcNorm());
				double cos_beta=p_s->periodize_vector(loc[0]-b)*p_s->periodize_vector(a-b)/(p_s->periodize_vector(loc[0]-b).CalcNorm()*p_s->periodize_vector(a-b).CalcNorm());
				k+=(p_s->periodize_vector(loc[0]-b)*cos_alpha/sqrt(1-cos_alpha*cos_alpha))+(p_s->periodize_vector(loc[0]-a)*cos_beta/sqrt(1-cos_beta*cos_beta));
				break;
			}
		}
	}
	k/=4.0*area_mixed;
	gaussian_curvature=(2*M_PI-sum_theta)/area_mixed;
	mean_curvature=k*NormalVector;
}
void Vertex::CalcLocMapping(){
	CalcNormalVector(0);
	std::vector<Vector3<double>> v_list(fi.size());
	std::vector<Vector3<double>> alpha_list(fi.size());
	for(int i=0;i<fi.size();i++){
		p_s->GetFace(fi[i])->CalcCenter(0);
		Vector3<double> gi=p_s->GetFace(fi[i])->GetCenter();
		Vector3<double> gi_i=gi-(gi*NormalVector)*NormalVector;
		v_list[i]=gi_i;
	}
	Vector3<double> t=p_s->periodized_center(v_list);
	Vector3<double> a(0.0,0.0,0.0);//=p_s->periodized_center(alpha_list);
	for(int i=0;i<fi.size();i++){
		Vector3<double> alpha=-(p_s->periodize_vector(t-p_s->GetFace(fi[i])->GetCenter())*p_s->GetFace(fi[i])->GetNormalVector())*NormalVector/(NormalVector*p_s->GetFace(fi[i])->GetNormalVector());
		//alpha_list[i]=alpha;
		a+=alpha/(double)fi.size();
	}
	//Vector3<double> a=p_s->periodized_center(alpha_list);
	Vector3<double> x=a+t;
	loc_mapping=p_s->periodize_vector(x);
}

//offsetの計算
void Vertex::CalcOffsetLoc(int deg){
	Vector3<int> v1((int)((loc[deg].GetX()-p_s->GetMINX())/RC),(int)((loc[deg].GetY()-p_s->GetMINY())/RC),(int)((loc[deg].GetZ()-p_s->GetMINZ())/RC));
	offsetloc=v1;
}
void Vertex::sortCounterClockwise(){//反時計回りに並び変え
    if(BoundaryFlag==0){   
        std::vector<int> tmp_ei;
        std::vector<int> tmp_vi;
        std::vector<int> tmp_fi;

        int curr_vidx=GetNextVertex(0);//隣の頂点を一つ指定
        tmp_vi.push_back(curr_vidx);//push_backする
        //p_s->isConsistent();
            
        for(int i=0;i<GetEdge().size();i++){//周囲の辺についてのループ
            bool ans =false;
            for(int fidx:p_s->GetVertex(curr_vidx)->GetFace()){//curr_vidx周りの面について
                for(int k=0;k<3;k++){//三角形のどの頂点か?
                    //三角形の点が反時計回りに登録されているという前提
                    if(p_s->GetFace(fidx)->GetVertex(k)==curr_vidx&&p_s->GetVertex(p_s->GetFace(fidx)->GetVertex((k+2)%3))==this){
                        curr_vidx=p_s->GetFace(fidx)->GetVertex((k+1)%3);
                        tmp_fi.push_back(fidx);//面をpush_back
                        ans=true;
                        if(curr_vidx!=tmp_vi[0]){//一周して
                            tmp_vi.push_back(curr_vidx);//被らないため
                        }
                        break;//更新したのでfor文を一つ抜ける
                    }
                }
                if(ans){
                break;
                }
            }   
        }
        //std::cout<<tmp_fi.size()<<" "<<GetEdge().size()<<std::endl;
        //std::cout<<tmp_vi.size()<<" "<<GetEdge().size()<<std::endl;
        assert(tmp_fi.size()==GetFace().size());
        assert(tmp_vi.size()==GetEdge().size());

        //tmp_eiへの代入(ローカルのインデックスとグローバルのインデックスを区別すること)
        for(int i=0;i<tmp_vi.size();i++){
            for(int j=0;j<GetEdge().size();j++){
                if(p_s->GetEdge(GetEdge(j))->GetVertex(1)==tmp_vi[i]||p_s->GetEdge(GetEdge(j))->GetVertex(0)==tmp_vi[i]){
                    tmp_ei.push_back(GetEdge(j));
                    //std::cout<<i<<" "<<j<<std::endl;
                }
            }
        }  
            //std::cout<<p_s->GetVertex(p_s->GetEdge(GetEdge(3))->GetVertex(1))<<" "<<p_s->GetVertex(p_s->GetEdge(GetEdge(3))->GetVertex(0))<<" "<<this<<std::endl;
            //std::cout<<tmp_ei.size()<<" "<<GetEdge().size()<<std::endl;
        assert(tmp_ei.size()==GetEdge().size());
        //代入(push_back)
        for(int i=0;i<tmp_vi.size();i++){
            SetNextVertex(i,tmp_vi[i]);
            SetEdge(i,tmp_ei[i]);
            SetFace(i,tmp_fi[i]);
        }
    }
    else{
        //std::vector<int> tmp;//境界の頂点のインデックスを保持
        std::vector<int> tmp_ei;//線
        std::vector<int> tmp_vi;//点
        std::vector<int> tmp_fi;//面

		std::vector<int> list1=fi;
		std::sort(list1.begin(),list1.end());
		bool ans=false;
		for(int fidx:fi){
			for(int j=0;j<3;j++){
				if(p_s->GetVertex(p_s->GetFace(fidx)->GetVertex(j))==this){
					std::vector<int>list2=p_s->GetVertex(p_s->GetFace(fidx)->GetVertex((j+1)%3))->GetFace();
					std::sort(list2.begin(),list2.end());
					std::vector<int> list_inter;
                    std::set_intersection(list1.begin(), list1.end(), list2.begin(), list2.end(), std::back_inserter(list_inter));
					if(list_inter.size()==1){
						tmp_vi.push_back(p_s->GetFace(fidx)->GetVertex((j+1)%3));
						tmp_vi.push_back(p_s->GetFace(fidx)->GetVertex((j+2)%3));
						tmp_fi.push_back(fidx);
						ans=true;
					}
					break;
				}
			}
			if(ans){
				break;
			}
		}
		/*
        int n1,n2;
        for(int i=0;i<GetEdge().size();i++){
            if(p_s->GetEdge(GetEdge(i))->GetFace().size()==1){//境界上の線であるとき
                //n1はv0のインデックス,n2は隣のインデックス
                n1=p_s->GetEdge(GetEdge(i))->GetVertex(0);
                n2=p_s->GetEdge(GetEdge(i))->GetVertex(1);
                //std::cout<<n1<<" "<<n2<<std::endl;
                if(p_s->GetVertex(n2)==this){//逆の時は入れ替える
                    int r=n1;
                    //std::cout<<r<<std::endl;
                    n1=n2;
                    //std::cout<<n1<<std::endl;
                    n2=r;    
                }
                tmp.push_back(n2);
                //std::cout<<n1<<" "<<n2<<std::endl;
                    
                //反時計回りに対応しているとき
                for(int k=0;k<3;k++){
                    //std::cout<<p_s->GetFace(p_s->GetEdge(GetEdge(i))->GetFace(0))->GetVertex(k)<<std::endl;
                    if(p_s->GetFace(p_s->GetEdge(GetEdge(i))->GetFace(0))->GetVertex(k)==n1&&p_s->GetFace(p_s->GetEdge(GetEdge(i))->GetFace(0))->GetVertex((k+1)%3)==n2){
                        tmp_fi.push_back(p_s->GetEdge(GetEdge(i))->GetFace(0));
                        tmp_vi.push_back(n2);
                        tmp_vi.push_back(p_s->GetFace(p_s->GetEdge(GetEdge(i))->GetFace(0))->GetVertex((k+2)%3));
                        break;
                    }
                }
            }
        }
        //std::cout<<tmp.size()<<" "<<tmp_vi.size()<<std::endl;
        assert(tmp.size()==2);//境界上の頂点が2つあるはず
        assert(tmp_vi.size()==2);//頂点の1つ目とその隣のみ入ってるはず
		*/
        int curr_vidx=tmp_vi[1];
        while(tmp_vi.size()<GetEdge().size()){//サイズが等しくなるまで
            bool ans=false;
            for(int fidx:p_s->GetVertex(curr_vidx)->GetFace()){
                //curr_vidx周りの面について
                for(int k=0;k<3;k++){
                    //三角形の点が反時計回りに登録されているという前提
                    if(p_s->GetFace(fidx)->GetVertex(k)==curr_vidx&&p_s->GetVertex(p_s->GetFace(fidx)->GetVertex((k+2)%3))==this){
                        curr_vidx=p_s->GetFace(fidx)->GetVertex((k+1)%3);
                        if(tmp_fi.size()<GetFace().size()){
                            tmp_fi.push_back(fidx);//面をpush_back
                        }
                        ans=true;
                        if(curr_vidx!=tmp_vi[0]){//一周して
                            tmp_vi.push_back(curr_vidx);//被らないため
                        }
                        break;//更新したのでfor文を一つ抜ける
                    }
                }
                if(ans){
                   break;//さらに抜ける
                }    
            }
        }
        assert(tmp_vi.size()==GetEdge().size());
        assert(tmp_fi.size()==GetFace().size());
        //tmp_eiへの代入
        for(int i=0;i<tmp_vi.size();i++){
            for(int j=0;j<GetEdge().size();j++){
                if(p_s->GetEdge(GetEdge(j))->GetVertex(1)==tmp_vi[i]||p_s->GetEdge(GetEdge(j))->GetVertex(0)==tmp_vi[i]){
                    tmp_ei.push_back(GetEdge(j));
                }
            }
        }
        assert(tmp_ei.size()==GetEdge().size());
        //代入(push_back)
        for(int i=0;i<tmp_vi.size();i++){
            SetNextVertex(i,tmp_vi[i]);
            SetEdge(i,tmp_ei[i]);
            SetFace(i,tmp_fi[i]);
        }
    }
}

int Vertex::GetFlagRec()const{
	return flag_rec;
}
double Vertex::GetCellTime()const{
	return cell_time;
}
int Vertex::GetFixFlag() const{
	return FixFlag;
}
int Vertex::GetBoundaryFlag() const{
	return BoundaryFlag;
}
int Vertex::GetApicalFlag() const{
	return apical_flag;
}
int Vertex::GetCrossBoundaryFlag() const{
	return cross_boundary_flag;
}
int Vertex::GetColorFlag()const{
	return colorFlag;
}
int Vertex::GetLayerFlag()const{
	return layer_flag;
}

Vector3<double> Vertex::GetNormalVector()const{
	return NormalVector;
}
Vector3<double>Vertex::GetLoc(int n)const{
	return loc[n];
}
Vector3<double>Vertex::GetLocMapping()const{
	return loc_mapping;
}
Vector3<double>Vertex::GetForce(int n)const{
	return force[n];
}
Vector3<double>Vertex::GetFrcThread(int n)const{
	return frc_thread[n];
}
Vector3<int>Vertex::GetOffsetLoc()const{
	return offsetloc;
}
double Vertex::GetGaussianCurvature()const{
	return gaussian_curvature;
}
double Vertex::GetMeanCurvature()const{
	return mean_curvature;
}
double Vertex::GetVoronoiArea()const{
	return voronoi_area;
}
double Vertex::GetVoronoiPerimeter()const{
	return voronoi_perimeter;
}

std::vector<int>Vertex::GetEdge()const{
	return ei;
}
int Vertex::GetEdge(int i)const{
	return ei[i];
}
std::vector<int>Vertex::GetFace()const{
	return fi;
}
int Vertex::GetFace(int i)const{
	return fi[i];
}
std::vector<int>Vertex::GetNextVertex()const{
	return vi_next;
}
int Vertex::GetNextVertex(int i)const{
	return vi_next[i];
}


void Vertex::SetFlagRec(){
	flag_rec=1;
}

void Vertex::SetExternalForce(const Vector3<double>& v){
	externalforce=v;
}

void Vertex::SetLoc(const Vector3<double>& v,int n){
	loc[n]=v;
}
void Vertex::addLoc(const Vector3<double>& v,int n){
	loc[n]+=v;
}
void Vertex::addForce(const Vector3<double>& v,int n){
	force[n]+=v;
}
void Vertex::resetForce(int n){
	force[n]=externalforce;
}
void Vertex::addFrcThread(const Vector3<double>& v,int n){
	frc_thread[n]+=v;
}
void Vertex::resetFrcThread(int n){
	frc_thread[n]=Vector3<double>(0.0,0.0,0.0);
}

void Vertex::SetOffsetLoc(const Vector3<int>& v){
	offsetloc=v;
}

void Vertex::SetFixFlag(int i){
	FixFlag = i;
}
void Vertex::SetBoundaryFlag(int i){
	BoundaryFlag = i;
}
void Vertex::SetApicalFlag(int i){
	apical_flag = i;
}
void Vertex::SetCrossBoundaryFlag(int i){
	cross_boundary_flag=i;
}
void Vertex::SetColorFlag(int i){
	colorFlag=i;
}
void Vertex::SetLayerFlag(int i){
	layer_flag=i;
}
void Vertex::SetCellTime(double a){
	cell_time=a;
}
void Vertex::SetNextVertex(int i,int j){
	vi_next[i]=j;
}
void Vertex::SetEdge(int i,int j){
	ei[i]=j;
}
void Vertex::SetFace(int i,int j){
	fi[i]=j;
}

void Vertex::PushEdge(int e){
	ei.push_back(e);
}
void Vertex::PushFace(int f){
	fi.push_back(f);
}
void Vertex::PushNextVertex(int v){
	vi_next.push_back(v);
}

void Vertex::SortAndErase(){
	std::sort(vi_next.begin(),vi_next.end());
	vi_next.erase(std::unique(vi_next.begin(),vi_next.end()),vi_next.end());
}
void Vertex::SortEdge(){
	std::sort(ei.begin(),ei.end());
}
void Vertex::SortFace(){
	std::sort(fi.begin(),fi.end());
}

void Vertex::findAndErase_ei(int search){
	auto itr = std::find(ei.begin(), ei.end(), search);
	if(itr == ei.end()) {
		std::cerr << "Couldn't find Edge-" << search << std::endl;
      	return;
	}
	ei.erase(itr);	
}
    	
void Vertex::findAndErase_fi(int search){
	auto itr = std::find(fi.begin(), fi.end(), search);
    if(itr == fi.end()) {
      	std::cerr << "Couldn't find Face" << search << std::endl;
		return;
    }
	fi.erase(itr);
}
void Vertex::findAndErase_vi_next(int search){
	auto itr = std::find(vi_next.begin(), vi_next.end(), search);
    if(itr == vi_next.end()) {
      	std::cerr << "Couldn't find Vertex-" << search << std::endl;
      	return;
    }
	vi_next.erase(itr);
}
void Vertex::Replace_ei(int i1,int i2){
	std::replace(ei.begin(),ei.end(),i1,i2);
}
void Vertex::Replace_fi(int i1,int i2){
	std::replace(fi.begin(),fi.end(),i1,i2);
}
void Vertex::Replace_vi_next(int i1,int i2){
	std::replace(vi_next.begin(),vi_next.end(),i1,i2);
}
void Vertex::Decrement_ei(int i){
	ei[i]--;
}
void Vertex::Decrement_fi(int i){
	fi[i]--;
}

/*-----------Edge class------------------*/
Edge::Edge(int v0,int v1,Shape* _p_s){
	vi[0]=v0;
	vi[1]=v1;
	Length=0;//初期化
	apical_flag=0;
	flipFlag=1;
	p_s=_p_s;
	CalcLength(0);
	//CalcDihedralAngle();
}
void Edge::CalcLength(int n){
	Length=(p_s->periodize_vector(p_s->GetVertex(vi[1])->GetLoc(n)-p_s->GetVertex(vi[0])->GetLoc(n))).CalcNorm();
}
void Edge::CalcDihedralAngle(int n){
	assert(fi.size()==2);
	p_s->GetFace(fi[0])->CalcNormalVector(n);
	p_s->GetFace(fi[1])->CalcNormalVector(n);
	double tmp=p_s->GetFace(fi[0])->GetNormalVector()*p_s->GetFace(fi[1])->GetNormalVector();
	if(tmp>1.0){
		tmp=1.0;
	}
	double theta=acos(tmp);
	p_s->GetFace(fi[1])->CalcNormalVector(n);
    p_s->GetFace(fi[0])->CalcNormalVector(n);
    p_s->GetFace(fi[0])->CalcCenter(n);
    p_s->GetFace(fi[1])->CalcCenter(n);

	Vector3<double> r_ij=p_s->periodize_vector(p_s->GetFace(fi[1])->GetCenter()-p_s->GetFace(fi[0])->GetCenter());
    Vector3<double> n_ji=p_s->GetFace(fi[0])->GetNormalVector()-p_s->GetFace(fi[1])->GetNormalVector();//相対位置ベクトルとの関係
    if(n_ji*r_ij<0){
        theta=-theta;
    }
	DihedralAngle=M_PI-theta;
}
/*----------------------------------*/
/*		Flip(shapeとedgeに関数)	     */
/*----------------------------------*/
bool Edge::FlipJudgement(){
	p_s->GetFace(fi[0])->CalcCenter(0);
	p_s->GetFace(fi[1])->CalcCenter(0);
	double leng1=p_s->periodize_vector(p_s->GetFace(fi[0])->GetCenter()-p_s->GetFace(fi[1])->GetCenter()).CalcNorm();
	if(leng1>=FLIP_THRESHOLD_LENGTH){//ダメ
		return false;
	}
	int vi1=vi[0];
	int vi2=vi[1];
	for(int j=0;j<2;j++){
		for(int i=0;i<3;i++){
			if(p_s->GetFace(fi[j])->GetVertex(i)==vi1&&p_s->GetFace(fi[j])->GetVertex((i+1)%3)==vi2){
				fi1=fi[j];
				vi_fi1=p_s->GetFace(fi[j])->GetVertex((i+2)%3);//反時計回り側f1
			}
			if(p_s->GetFace(fi[j])->GetVertex(i)==vi2&&p_s->GetFace(fi[j])->GetVertex((i+1)%3)==vi1){
				fi2=fi[j];
				vi_fi2=p_s->GetFace(fi[j])->GetVertex((i+2)%3);//時計回りが側f2
			}
		}
	}	
	if(vi_fi1==vi_fi2){//同じだったら強制終了
		std::cout<<"triangle"<<vi1<<","<<vi2<<","<<vi_fi1<<std::endl;
		std::cout<<"There are same faces."<<std::endl;
		exit(1);
	}
	for(int i=0;i<p_s->GetVertex(vi_fi1)->GetEdge().size();i++){//フリップ後の辺と同じ辺がすでにあるときはreturn
		//int edge_number1=p_s->GetVertex(vi_fi1)->GetEdge(i);
		for(int j=0;j<p_s->GetVertex(vi_fi2)->GetEdge().size();j++){
			if(p_s->GetVertex(vi_fi1)->GetEdge(i)==p_s->GetVertex(vi_fi2)->GetEdge(j)){
				return false;
			}
		}
	}


	p_s->GetFace(fi1)->CalcNormalVector(0);
	p_s->GetFace(fi2)->CalcNormalVector(0);
	double n1_n2=p_s->GetFace(fi1)->GetNormalVector()*p_s->GetFace(fi2)->GetNormalVector();

	//二つの面の法線ベクトルの内積
	if(n1_n2<=0){//負の時はダメ
		return false;
	}

	//長さの変化
	double dis1=(p_s->periodize_vector(p_s->GetVertex(vi1)->GetLoc(0)-p_s->GetVertex(vi2)->GetLoc(0))).CalcNorm();
	double dis2=(p_s->periodize_vector(p_s->GetVertex(vi_fi1)->GetLoc(0)-p_s->GetVertex(vi_fi2)->GetLoc(0))).CalcNorm();
	if(dis1<=dis2){//ダメ
		return false;
	}
	double cos_v1_f1=p_s->periodize_vector(p_s->GetVertex(vi_fi1)->GetLoc(0)-p_s->GetVertex(vi1)->GetLoc(0))*p_s->periodize_vector(p_s->GetVertex(vi2)->GetLoc(0)-p_s->GetVertex(vi1)->GetLoc(0))/((p_s->periodize_vector(p_s->GetVertex(vi_fi1)->GetLoc(0)-p_s->GetVertex(vi1)->GetLoc(0))).CalcNorm()*(p_s->periodize_vector(p_s->GetVertex(vi2)->GetLoc(0)-p_s->GetVertex(vi1)->GetLoc(0))).CalcNorm());
	double cos_v2_f1=p_s->periodize_vector(p_s->GetVertex(vi_fi1)->GetLoc(0)-p_s->GetVertex(vi2)->GetLoc(0))*p_s->periodize_vector(p_s->GetVertex(vi1)->GetLoc(0)-p_s->GetVertex(vi2)->GetLoc(0))/((p_s->periodize_vector(p_s->GetVertex(vi_fi1)->GetLoc(0)-p_s->GetVertex(vi2)->GetLoc(0))).CalcNorm()*(p_s->periodize_vector(p_s->GetVertex(vi1)->GetLoc(0)-p_s->GetVertex(vi2)->GetLoc(0))).CalcNorm());
	double cos_v1_f2=p_s->periodize_vector(p_s->GetVertex(vi_fi2)->GetLoc(0)-p_s->GetVertex(vi1)->GetLoc(0))*p_s->periodize_vector(p_s->GetVertex(vi2)->GetLoc(0)-p_s->GetVertex(vi1)->GetLoc(0))/((p_s->periodize_vector(p_s->GetVertex(vi_fi2)->GetLoc(0)-p_s->GetVertex(vi1)->GetLoc(0))).CalcNorm()*(p_s->periodize_vector(p_s->GetVertex(vi2)->GetLoc(0)-p_s->GetVertex(vi1)->GetLoc(0))).CalcNorm());
	double cos_v2_f2=p_s->periodize_vector(p_s->GetVertex(vi_fi2)->GetLoc(0)-p_s->GetVertex(vi2)->GetLoc(0))*p_s->periodize_vector(p_s->GetVertex(vi1)->GetLoc(0)-p_s->GetVertex(vi2)->GetLoc(0))/((p_s->periodize_vector(p_s->GetVertex(vi_fi2)->GetLoc(0)-p_s->GetVertex(vi2)->GetLoc(0))).CalcNorm()*(p_s->periodize_vector(p_s->GetVertex(vi1)->GetLoc(0)-p_s->GetVertex(vi2)->GetLoc(0))).CalcNorm());
		
	double theta_v1=acos(cos_v1_f1)+acos(cos_v1_f2);
	double theta_v2=acos(cos_v2_f2)+acos(cos_v2_f1);

	if((theta_v1<2.0*M_PI/3.0)&&(theta_v2<2.0*M_PI/3.0)){//フリップがOKの場合
		return true;
	}
	return false;//ダメの場合
}
void Edge::FlipOperation(){
	p_s->GetFace(fi[0])->CalcCenter(0);
	p_s->GetFace(fi[1])->CalcCenter(0);
	double leng1=p_s->periodize_vector(p_s->GetFace(fi[0])->GetCenter()-p_s->GetFace(fi[1])->GetCenter()).CalcNorm();
	if(leng1<FLIP_THRESHOLD_LENGTH){
		int fi1;
		int fi2;
		int vi1=vi[0];
		int vi2=vi[1];
		int vi_fi1;//辺に関係のない頂点のインデックス(f1側)
		int vi_fi2;//辺に関係のない頂点のインデックス(f2側)
		for(int j=0;j<2;j++){
			for(int i=0;i<3;i++){
				if(p_s->GetFace(fi[j])->GetVertex(i)==vi1&&p_s->GetFace(fi[j])->GetVertex((i+1)%3)==vi2){
					fi1=fi[j];
					vi_fi1=p_s->GetFace(fi[j])->GetVertex((i+2)%3);//反時計回り側f1
				}
				if(p_s->GetFace(fi[j])->GetVertex(i)==vi2&&p_s->GetFace(fi[j])->GetVertex((i+1)%3)==vi1){
					fi2=fi[j];
					vi_fi2=p_s->GetFace(fi[j])->GetVertex((i+2)%3);//時計回りが側f2
				}
			}
		}
		if(vi_fi1==vi_fi2){
			std::cout<<"triangle"<<vi1<<","<<vi2<<","<<vi_fi1<<std::endl;
			std::cout<<"There are same faces."<<std::endl;
			exit(1);
		}
		
		for(int i=0;i<p_s->GetVertex(vi_fi1)->GetEdge().size();i++){//フリップ後の辺と同じ辺がすでにあるときはreturn
			//int edge_number1=p_s->GetVertex(vi_fi1)->GetEdge(i);
			for(int j=0;j<p_s->GetVertex(vi_fi2)->GetEdge().size();j++){
				if(p_s->GetVertex(vi_fi1)->GetEdge(i)==p_s->GetVertex(vi_fi2)->GetEdge(j)){
					return;
				}
			}
		}
		p_s->GetFace(fi1)->CalcNormalVector(0);
		p_s->GetFace(fi2)->CalcNormalVector(0);
		double n1_n2=p_s->GetFace(fi1)->GetNormalVector()*p_s->GetFace(fi2)->GetNormalVector();
		//二つの面の法線ベクトルの内積
		if(n1_n2<=0){//負の時はreturn
			return;
		}

		//長さの変化
		double dis1=(p_s->periodize_vector(p_s->GetVertex(vi1)->GetLoc(0)-p_s->GetVertex(vi2)->GetLoc(0))).CalcNorm();
		double dis2=(p_s->periodize_vector(p_s->GetVertex(vi_fi1)->GetLoc(0)-p_s->GetVertex(vi_fi2)->GetLoc(0))).CalcNorm();
		if(dis1<=dis2){
			return;
		}
		double cos_v1_f1=p_s->periodize_vector(p_s->GetVertex(vi_fi1)->GetLoc(0)-p_s->GetVertex(vi1)->GetLoc(0))*p_s->periodize_vector(p_s->GetVertex(vi2)->GetLoc(0)-p_s->GetVertex(vi1)->GetLoc(0))/((p_s->periodize_vector(p_s->GetVertex(vi_fi1)->GetLoc(0)-p_s->GetVertex(vi1)->GetLoc(0))).CalcNorm()*(p_s->periodize_vector(p_s->GetVertex(vi2)->GetLoc(0)-p_s->GetVertex(vi1)->GetLoc(0))).CalcNorm());
		double cos_v2_f1=p_s->periodize_vector(p_s->GetVertex(vi_fi1)->GetLoc(0)-p_s->GetVertex(vi2)->GetLoc(0))*p_s->periodize_vector(p_s->GetVertex(vi1)->GetLoc(0)-p_s->GetVertex(vi2)->GetLoc(0))/((p_s->periodize_vector(p_s->GetVertex(vi_fi1)->GetLoc(0)-p_s->GetVertex(vi2)->GetLoc(0))).CalcNorm()*(p_s->periodize_vector(p_s->GetVertex(vi1)->GetLoc(0)-p_s->GetVertex(vi2)->GetLoc(0))).CalcNorm());
		double cos_v1_f2=p_s->periodize_vector(p_s->GetVertex(vi_fi2)->GetLoc(0)-p_s->GetVertex(vi1)->GetLoc(0))*p_s->periodize_vector(p_s->GetVertex(vi2)->GetLoc(0)-p_s->GetVertex(vi1)->GetLoc(0))/((p_s->periodize_vector(p_s->GetVertex(vi_fi2)->GetLoc(0)-p_s->GetVertex(vi1)->GetLoc(0))).CalcNorm()*(p_s->periodize_vector(p_s->GetVertex(vi2)->GetLoc(0)-p_s->GetVertex(vi1)->GetLoc(0))).CalcNorm());
		double cos_v2_f2=p_s->periodize_vector(p_s->GetVertex(vi_fi2)->GetLoc(0)-p_s->GetVertex(vi2)->GetLoc(0))*p_s->periodize_vector(p_s->GetVertex(vi1)->GetLoc(0)-p_s->GetVertex(vi2)->GetLoc(0))/((p_s->periodize_vector(p_s->GetVertex(vi_fi2)->GetLoc(0)-p_s->GetVertex(vi2)->GetLoc(0))).CalcNorm()*(p_s->periodize_vector(p_s->GetVertex(vi1)->GetLoc(0)-p_s->GetVertex(vi2)->GetLoc(0))).CalcNorm());
		

		double theta_v1=acos(cos_v1_f1)+acos(cos_v1_f2);
		double theta_v2=acos(cos_v2_f2)+acos(cos_v2_f1);

		if((theta_v1<2.0*M_PI/3.0)&&(theta_v2<2.0*M_PI/3.0)){
			std::cout<<"Flip Operation-"<<vi1<<" and "<<vi2<<std::endl;
			std::cout<<"Changed to Edge-"<<vi_fi1<<","<<vi_fi2<<std::endl;
			p_s->FlipCount();

			
			char fname[100]="d_Flip.dat";
			std::ofstream fout(fname,std::ios::app);
			fout<<"Flip-"<<vi1<<" and "<<vi2<<"\n";
			fout<<"Changed to Edge-"<<vi_fi1<<" , "<<vi_fi2<<"\n";
			fout.close();
			//Flip処理(線を5つ)
			int ei_f1_vi1;
			int ei_f1_vi2;
			int ei_f2_vi1;
			int ei_f2_vi2;
			int ei_vi1_vi2;
			for(int i=0;i<3;i++){//辺の対応をつける
				if(p_s->GetEdge(p_s->GetFace(fi1)->GetEdge(i))->GetVertex(0)!=vi2&&p_s->GetEdge(p_s->GetFace(fi1)->GetEdge(i))->GetVertex(1)!=vi2){
					ei_f1_vi1=p_s->GetFace(fi1)->GetEdge(i);//fi1の辺でどちらもvi2でない
				}
				if(p_s->GetEdge(p_s->GetFace(fi1)->GetEdge(i))->GetVertex(0)!=vi1&&p_s->GetEdge(p_s->GetFace(fi1)->GetEdge(i))->GetVertex(1)!=vi1){
					ei_f1_vi2=p_s->GetFace(fi1)->GetEdge(i);//fi1の辺でどちらもvi1でない
				}
				if(p_s->GetEdge(p_s->GetFace(fi2)->GetEdge(i))->GetVertex(0)!=vi1&&p_s->GetEdge(p_s->GetFace(fi2)->GetEdge(i))->GetVertex(1)!=vi1){
					ei_f2_vi2=p_s->GetFace(fi2)->GetEdge(i);//fi2について、どちらもvi1でない
				}
				if(p_s->GetEdge(p_s->GetFace(fi2)->GetEdge(i))->GetVertex(0)!=vi2&&p_s->GetEdge(p_s->GetFace(fi2)->GetEdge(i))->GetVertex(1)!=vi2){
					ei_f2_vi1=p_s->GetFace(fi2)->GetEdge(i);//fi2について、どちらもvi2でない
				}
				if(p_s->GetEdge(p_s->GetFace(fi1)->GetEdge(i))->GetVertex(0)!=vi_fi1&&p_s->GetEdge(p_s->GetFace(fi1)->GetEdge(i))->GetVertex(1)!=vi_fi1){
					ei_vi1_vi2=p_s->GetFace(fi1)->GetEdge(i);//f1の辺でどちらもvi_fi1でない
				}
			}

			//ei_f1_vi1とei_f2_vi2は変更なし
			//ei_f1_vi2について
			p_s->GetEdge(ei_f1_vi2)->Replace_fi(fi1,fi2);//fi1をfi2に置き換える
			//ei_f2_vi1について
			p_s->GetEdge(ei_f2_vi1)->Replace_fi(fi2,fi1);//fi2をfi1に置き換える
			//ei_vi1_vi2について
			p_s->GetEdge(ei_vi1_vi2)->Replace_vi(vi1,vi_fi1);//vi1をvi_fi1に置き換える
			p_s->GetEdge(ei_vi1_vi2)->Replace_vi(vi2,vi_fi2);//vi2をvi_fi2に置き換える

			//faceについて
			//fi1
			p_s->GetFace(fi1)->Replace_vi(vi2,vi_fi2);//vi2をvi_fi2に置き換える****
			p_s->GetFace(fi1)->Replace_ei(ei_f1_vi2,ei_f2_vi1);//辺を置き換える
			//fi2
			p_s->GetFace(fi2)->Replace_vi(vi1,vi_fi1);//vi1をvi_fi1に置き換える****
			p_s->GetFace(fi2)->Replace_ei(ei_f2_vi1,ei_f1_vi2);//辺を置き換える

			//Vertexについて
			//vi1について
			p_s->GetVertex(vi1)->findAndErase_ei(ei_vi1_vi2);
			p_s->GetVertex(vi1)->findAndErase_fi(fi2);
			p_s->GetVertex(vi1)->findAndErase_vi_next(vi2);
			//vi2について
			p_s->GetVertex(vi2)->findAndErase_ei(ei_vi1_vi2);
			p_s->GetVertex(vi2)->findAndErase_fi(fi1);
			p_s->GetVertex(vi2)->findAndErase_vi_next(vi1);
			//vi_fi1
			p_s->GetVertex(vi_fi1)->PushEdge(ei_vi1_vi2);
			p_s->GetVertex(vi_fi1)->PushFace(fi2);
			p_s->GetVertex(vi_fi1)->PushNextVertex(vi_fi2);
			//vi_fi2			
			p_s->GetVertex(vi_fi2)->PushEdge(ei_vi1_vi2);
			p_s->GetVertex(vi_fi2)->PushFace(fi1);
			p_s->GetVertex(vi_fi2)->PushNextVertex(vi_fi1);

			//並べ替え
			p_s->GetVertex(vi_fi1)->sortCounterClockwise();
			p_s->GetVertex(vi_fi2)->sortCounterClockwise();

			//LineApical(3通り)
			if((p_s->GetEdge(ei_f1_vi1)->GetApicalFlag()==2&&p_s->GetEdge(ei_f1_vi2)->GetApicalFlag()==2)||(p_s->GetEdge(ei_f2_vi1)->GetApicalFlag()==2&&p_s->GetEdge(ei_f2_vi2)->GetApicalFlag()==2)){
				p_s->GetEdge(ei_vi1_vi2)->SetApicalFlag(2);//元々交差していないのときはflagを立てる
			}
			else if((p_s->GetEdge(ei_f1_vi1)->GetApicalFlag()==2&&p_s->GetEdge(ei_f2_vi1)->GetApicalFlag()==2)||(p_s->GetEdge(ei_f1_vi2)->GetApicalFlag()==2&&p_s->GetEdge(ei_f2_vi2)->GetApicalFlag()==2)){
				p_s->GetEdge(ei_vi1_vi2)->SetApicalFlag(0);//元々交差しているときはflagを消す
			}
			//対角に繋がっている場合はそのままでOK
			flipFlag=FLIP_COOLTIME_STEP;
		}
	}
}

void Edge::FlipOperation_forge(){
	int vi1=vi[0];
	int vi2=vi[1];

	std::cout<<"Flip Operation-"<<vi1<<" and "<<vi2<<std::endl;
	std::cout<<"Changed to Edge-"<<vi_fi1<<","<<vi_fi2<<std::endl;
	p_s->FlipCount();
		
	char fname[100]="d_Flip.dat";
	std::ofstream fout(fname,std::ios::app);
	fout<<"Flip-"<<vi1<<" and "<<vi2<<"\n";
	fout<<"Changed to Edge-"<<vi_fi1<<" , "<<vi_fi2<<"\n";
	fout.close();
	//Flip処理(線を5つ)
	int ei_f1_vi1;
	int ei_f1_vi2;
	int ei_f2_vi1;
	int ei_f2_vi2;
	int ei_vi1_vi2;
	for(int i=0;i<3;i++){//辺の対応をつける
		if(p_s->GetEdge(p_s->GetFace(fi1)->GetEdge(i))->GetVertex(0)!=vi2&&p_s->GetEdge(p_s->GetFace(fi1)->GetEdge(i))->GetVertex(1)!=vi2){
			ei_f1_vi1=p_s->GetFace(fi1)->GetEdge(i);//fi1の辺でどちらもvi2でない
		}
		if(p_s->GetEdge(p_s->GetFace(fi1)->GetEdge(i))->GetVertex(0)!=vi1&&p_s->GetEdge(p_s->GetFace(fi1)->GetEdge(i))->GetVertex(1)!=vi1){
			ei_f1_vi2=p_s->GetFace(fi1)->GetEdge(i);//fi1の辺でどちらもvi1でない
		}
		if(p_s->GetEdge(p_s->GetFace(fi2)->GetEdge(i))->GetVertex(0)!=vi1&&p_s->GetEdge(p_s->GetFace(fi2)->GetEdge(i))->GetVertex(1)!=vi1){
			ei_f2_vi2=p_s->GetFace(fi2)->GetEdge(i);//fi2について、どちらもvi1でない
		}
		if(p_s->GetEdge(p_s->GetFace(fi2)->GetEdge(i))->GetVertex(0)!=vi2&&p_s->GetEdge(p_s->GetFace(fi2)->GetEdge(i))->GetVertex(1)!=vi2){
			ei_f2_vi1=p_s->GetFace(fi2)->GetEdge(i);//fi2について、どちらもvi2でない
		}
		if(p_s->GetEdge(p_s->GetFace(fi1)->GetEdge(i))->GetVertex(0)!=vi_fi1&&p_s->GetEdge(p_s->GetFace(fi1)->GetEdge(i))->GetVertex(1)!=vi_fi1){
			ei_vi1_vi2=p_s->GetFace(fi1)->GetEdge(i);//f1の辺でどちらもvi_fi1でない
		}
	}

	//ei_f1_vi1とei_f2_vi2は変更なし
	//ei_f1_vi2について
	p_s->GetEdge(ei_f1_vi2)->Replace_fi(fi1,fi2);//fi1をfi2に置き換える
	//ei_f2_vi1について
	p_s->GetEdge(ei_f2_vi1)->Replace_fi(fi2,fi1);//fi2をfi1に置き換える
	//ei_vi1_vi2について
	p_s->GetEdge(ei_vi1_vi2)->Replace_vi(vi1,vi_fi1);//vi1をvi_fi1に置き換える
	p_s->GetEdge(ei_vi1_vi2)->Replace_vi(vi2,vi_fi2);//vi2をvi_fi2に置き換える
	//faceについて
	//fi1
	p_s->GetFace(fi1)->Replace_vi(vi2,vi_fi2);//vi2をvi_fi2に置き換える****
	p_s->GetFace(fi1)->Replace_ei(ei_f1_vi2,ei_f2_vi1);//辺を置き換える
	//fi2
	p_s->GetFace(fi2)->Replace_vi(vi1,vi_fi1);//vi1をvi_fi1に置き換える****
	p_s->GetFace(fi2)->Replace_ei(ei_f2_vi1,ei_f1_vi2);//辺を置き換える

	//Vertexについて
	//vi1について

	if(p_s->GetEdge(p_s->GetVertex(vi1)->GetEdge(0))==this){//sortの必要をなくすため(faceだけずれてしまう)
		int tmp_num=p_s->GetVertex(vi1)->GetEdge().size()-1;
		int tmp0=p_s->GetVertex(vi1)->GetEdge(0);
		int tmp1=p_s->GetVertex(vi1)->GetEdge(tmp_num);
		
		p_s->GetVertex(vi1)->SetEdge(0,tmp1);
		p_s->GetVertex(vi1)->SetEdge(tmp_num,tmp0);
		
		tmp0=p_s->GetVertex(vi1)->GetNextVertex(0);
		tmp1=p_s->GetVertex(vi1)->GetNextVertex(tmp_num);
		p_s->GetVertex(vi1)->SetNextVertex(0,tmp1);
		p_s->GetVertex(vi1)->SetNextVertex(tmp_num,tmp0);
	}

	p_s->GetVertex(vi1)->findAndErase_ei(ei_vi1_vi2);
	p_s->GetVertex(vi1)->findAndErase_fi(fi2);
	p_s->GetVertex(vi1)->findAndErase_vi_next(vi2);
	//vi2について
	if(p_s->GetEdge(p_s->GetVertex(vi2)->GetEdge(0))==this){//sortの必要をなくすため
		int tmp_num=p_s->GetVertex(vi2)->GetEdge().size()-1;
		int tmp0=p_s->GetVertex(vi2)->GetEdge(0);
		int tmp1=p_s->GetVertex(vi2)->GetEdge(tmp_num);
		
		p_s->GetVertex(vi2)->SetEdge(0,tmp1);
		p_s->GetVertex(vi2)->SetEdge(tmp_num,tmp0);
		
		tmp0=p_s->GetVertex(vi2)->GetNextVertex(0);
		tmp1=p_s->GetVertex(vi2)->GetNextVertex(tmp_num);
		p_s->GetVertex(vi2)->SetNextVertex(0,tmp1);
		p_s->GetVertex(vi2)->SetNextVertex(tmp_num,tmp0);
	}
	p_s->GetVertex(vi2)->findAndErase_ei(ei_vi1_vi2);
	p_s->GetVertex(vi2)->findAndErase_fi(fi1);
	p_s->GetVertex(vi2)->findAndErase_vi_next(vi1);
	//vi_fi1
	p_s->GetVertex(vi_fi1)->PushEdge(ei_vi1_vi2);
	p_s->GetVertex(vi_fi1)->PushFace(fi2);
	p_s->GetVertex(vi_fi1)->PushNextVertex(vi_fi2);
	//vi_fi2			
	p_s->GetVertex(vi_fi2)->PushEdge(ei_vi1_vi2);
	p_s->GetVertex(vi_fi2)->PushFace(fi1);
	p_s->GetVertex(vi_fi2)->PushNextVertex(vi_fi1);
	//並べ替え
	p_s->GetVertex(vi_fi1)->sortCounterClockwise();
	p_s->GetVertex(vi_fi2)->sortCounterClockwise();


	//LineApicalの処理(3通り)
	if((p_s->GetEdge(ei_f1_vi1)->GetApicalFlag()==2&&p_s->GetEdge(ei_f1_vi2)->GetApicalFlag()==2)||(p_s->GetEdge(ei_f2_vi1)->GetApicalFlag()==2&&p_s->GetEdge(ei_f2_vi2)->GetApicalFlag()==2)){
		p_s->GetEdge(ei_vi1_vi2)->SetApicalFlag(2);//元々交差していないのときはflagを立てる
	}
	else if((p_s->GetEdge(ei_f1_vi1)->GetApicalFlag()==2&&p_s->GetEdge(ei_f2_vi1)->GetApicalFlag()==2)||(p_s->GetEdge(ei_f1_vi2)->GetApicalFlag()==2&&p_s->GetEdge(ei_f2_vi2)->GetApicalFlag()==2)){
		p_s->GetEdge(ei_vi1_vi2)->SetApicalFlag(0);//元々交差しているときはflagを消す
	}
	//対角に繋がっている場合はそのままでOK
	flipFlag=FLIP_COOLTIME_STEP;

	if(p_s->GetEdge(ei_f1_vi2)->GetFlipFlag()==-1){
		p_s->GetEdge(ei_f1_vi2)->SetFlipFlag(0);
	}
	if(p_s->GetEdge(ei_f2_vi1)->GetFlipFlag()==-1){
		p_s->GetEdge(ei_f2_vi1)->SetFlipFlag(0);
	}
	if(p_s->GetEdge(ei_f1_vi1)->GetFlipFlag()==-1){
		p_s->GetEdge(ei_f1_vi1)->SetFlipFlag(0);
	}
	if(p_s->GetEdge(ei_f2_vi2)->GetFlipFlag()==-1){
		p_s->GetEdge(ei_f2_vi2)->SetFlipFlag(0);
	}
}

void Edge::SetBaseLength(const double a){
	BaseLength=a;
}
void Edge::SetBaseDihedralAngle(const double a){
	BaseDihedralAngle=a;
}
/*
void Edge::ResetBaseDihedralAngle(){
	BaseDihedralAngle+=p_s->GetApicalConstrictionAngle();
}
void Edge::SetApicalBaseDihedralAngle(){
	BaseDihedralAngle-=p_s->GetApicalConstrictionAngle();
}
*/
void Edge::SetK1_Length(const double a){
	k1_length=a;
}
void Edge::SetK2_Length(const double a){
	k2_length=a;
}
void Edge::SetK_Theta(const double a){
	k_theta=a;
}

int Edge::GetVertex(const int i)const{
	return vi[i];
}
std::vector<int>Edge::GetFace()const{
	return fi;
}
int Edge::GetFace(const int i)const{
	return fi[i];
}
int Edge::GetFlipFlag()const{
	return flipFlag;
}
void Edge::SetFlipFlag(int a){
	flipFlag=a;
}
int Edge::GetApicalFlag()const{
	return apical_flag;
}
void Edge::SetApicalFlag(int i){
	apical_flag=i;
}

double Edge::GetLength()const{
	return Length;
}
double Edge::GetDihedralAngle()const{
	return DihedralAngle;
}
double Edge::GetBaseLength()const{
	return BaseLength;
}
double Edge::GetBaseDihedralAngle()const{
	return BaseDihedralAngle;
}
double Edge::Get_k1_length()const{
	return k1_length;
}
double Edge::Get_k2_length()const{
	return k2_length;
}
double Edge::Get_k_theta()const{
	return k_theta;
}
void Edge::findAndErase_fi(int search){
	auto itr = std::find(fi.begin(), fi.end(), search);
    if(itr == fi.end()) {
      	std::cerr << "Couldn't find Face" << search << std::endl;
		return;
    }
	fi.erase(itr);
}
void Edge::SortFace(){
	std::sort(fi.begin(),fi.end());
}
void Edge::PushFace(int f){
	fi.push_back(f);
}
void Edge::Replace_vi(int i1,int i2){
	if(vi[0]==i1){
		vi[0]=i2;//代入
		return;
	}
	if(vi[1]==i1){
		vi[1]=i2;//代入
		return;
	}
}
void Edge::Replace_fi(int i1,int i2){
	std::replace(fi.begin(),fi.end(),i1,i2);
}
void Edge::Decrement_fi(int i){
	fi[i]--;
}
/*---------------------Face Class--------------------------*/
Face::Face(int v0,int v1,int v2,Shape* _p_s){
	PushVertex(v0);PushVertex(v1);PushVertex(v2);
	//p_v.push_back(v0);p_v.push_back(v1);p_v.push_back(v2);
	p_s=_p_s;
	//CalcCenter();
	//CalcNormalVector();
	//CalcArea();
	cross_boundary_flag=0;
}

void Face::CalcCenter(int n){
	Vector3<double> s(0.0,0.0,0.0);
	for(int i=1;i<3;i++){//i=1,2
		Vector3<double> b=p_s->GetVertex(vi[i])->GetLoc(n)-p_s->GetVertex(vi[0])->GetLoc(n);
		b=p_s->periodize_vector(b);
		s+=b;
	}
	Vector3<double> g=p_s->GetVertex(vi[0])->GetLoc(n)+s/(double)3.0;
	g=p_s->periodize_vector(g);
	center=g;
}
void Face::CalcNormalVector(int n){
	Vector3<double> tmp1=p_s->periodize_vector(p_s->GetVertex(vi[1])->GetLoc(n)-p_s->GetVertex(vi[0])->GetLoc(n));
	Vector3<double> tmp2=p_s->periodize_vector(p_s->GetVertex(vi[2])->GetLoc(n)-p_s->GetVertex(vi[0])->GetLoc(n));
	Vector3<double> tmp_n(tmp1 % tmp2);
	NormalVector=tmp_n/(tmp_n.CalcNorm()+std::pow(10.0,-5));
}
void Face::CalcArea(int n){
	Vector3<double> v0=p_s->periodize_vector(p_s->GetVertex(vi[1])->GetLoc(n)-p_s->GetVertex(vi[0])->GetLoc(n));
	Vector3<double> v1=p_s->periodize_vector(p_s->GetVertex(vi[2])->GetLoc(n)-p_s->GetVertex(vi[0])->GetLoc(n));
	Area=0.5*(v0%v1).CalcNorm();
}
void Face::CalcCircumCenter(int n){
	
	Vector3<double> a=p_s->GetVertex(vi[2])->GetLoc(n);
	Vector3<double> b=p_s->GetVertex(vi[0])->GetLoc(n);
	Vector3<double> c=p_s->GetVertex(vi[1])->GetLoc(n);

	double a_2=(p_s->periodize_vector(b-c)).CalcSqr();
	double b_2=(p_s->periodize_vector(c-a)).CalcSqr();
	double c_2=(p_s->periodize_vector(a-b)).CalcSqr();

	std::vector<std::pair<double,int>> list;
	list.push_back(std::make_pair(a_2,2));
	list.push_back(std::make_pair(b_2,0));
	list.push_back(std::make_pair(c_2,1));
	std::sort(list.begin(),list.end());

	if(list[2].first>=(list[1].first+list[0].first)){//when triangle is obtuse
		int num=list[2].second;
		Vector3<double> g(0.0,0.0,0.0);
		Vector3<double> d=p_s->GetVertex(vi[(num+1)%3])->GetLoc(n);
		Vector3<double> e=p_s->GetVertex(vi[(num+2)%3])->GetLoc(n);
		g=p_s->periodize_vector(d-e);
		g=e+(g/2.0);
		circum_center=p_s->periodize_vector(g);
	}
	else{
		//double a_2=(b-c).CalcSqr();
		//double b_2=(c-a).CalcSqr();
		//double c_2=(a-b).CalcSqr();
		//CalcArea(n);

		Vector3<double> ab=p_s->periodize_vector(b-a);
		Vector3<double> ac=p_s->periodize_vector(c-a);

		Vector3<double> ao=((ab.CalcSqr()*ac.CalcSqr()-ac.CalcSqr()*(ac*ab))/(2*(ab.CalcSqr()*ac.CalcSqr()-(ab*ac)*(ab*ac))))*ab+((ab.CalcSqr()*ac.CalcSqr()-ab.CalcSqr()*(ac*ab))/(2*(ab.CalcSqr()*ac.CalcSqr()-(ab*ac)*(ab*ac))))*ac;

		//Vector3<double> v1=(a_2*(b_2+c_2-a_2)/(16*GetArea()*GetArea()))*a;
		//Vector3<double> v2=(b_2*(c_2+a_2-b_2)/(16*GetArea()*GetArea()))*b;
		//Vector3<double> v3=(c_2*(a_2+b_2-c_2)/(16*GetArea()*GetArea()))*c;

		//Vector3<double> dis1=p_s->periodize_vector(v2-v1);
		//Vector3<double> dis2=p_s->periodize_vector(v3-v1);
		//Vector3<double> dis=3.0*v1+p_s->periodize_vector(dis1+dis2);
		circum_center=p_s->periodize_vector(a+ao);
		//circum_center=v1+v2+v3;//外心
	}
}
void Face::SetCrossBoundaryFlag(int i){
	cross_boundary_flag=i;
}
void Face::SetBaseArea(const double a){
	BaseArea=a;
}
void Face::SetBaseArea(){
	double a=p_s->GetEdge(ei[0])->GetBaseLength();
	double b=p_s->GetEdge(ei[1])->GetBaseLength();
	double c=p_s->GetEdge(ei[2])->GetBaseLength();
	double s=(a+b+c)/2.0;
	if(s*(s-a)*(s-b)*(s-c)>0){
		BaseArea=sqrt(s*(s-a)*(s-b)*(s-c));
	}
}
void Face::SetK_Area(const double a){
	k_area=a;
}
int Face::GetCrossBoundaryFlag()const{
	return cross_boundary_flag;
}
std::vector<int> Face::GetVertex()const{
	return vi;
}
int Face::GetVertex(const int i)const{
	return vi[i];
}
std::vector<int> Face::GetEdge()const{
	return ei;
}
int Face::GetEdge(const int i)const{
	return ei[i];
}
Vector3<double> Face::GetCenter()const{
	return center;
}
Vector3<double> Face::GetCircumCenter()const{
	return circum_center;
}
Vector3<double>Face::GetNormalVector()const{
	return NormalVector;
}
double Face::GetArea()const{
	return Area;
}
double Face::GetBaseArea()const{
	return BaseArea;
}
double Face::Get_k_area()const{
	return k_area;
}

void Face::SortEdge(){
	std::sort(ei.begin(),ei.end());
}
void Face::PushVertex(int v){
	vi.push_back(v);
}
void Face::PushEdge(int e){
	ei.push_back(e);
}
void Face::Replace_vi(int i1,int i2){
	std::replace(vi.begin(),vi.end(),i1,i2);
}
void Face::Replace_ei(int i1,int i2){
	std::replace(ei.begin(),ei.end(),i1,i2);
}
void Face::Decrement_ei(int i){
	ei[i]--;
}
void Face::Decrement_vi(int i){
	vi[i]--;
}