#include"restructure.hpp"
#include"Random.hpp"
#include"class.hpp"
#include"ODE_solver.hpp"
#include<iostream>
#include<fstream>
#include<vector>
#include<cassert>
#include<algorithm>
#include<cmath>
#include<utility>
#include"vector3.hpp"
#include<omp.h>
namespace restructure{
    void cellDivision(Random* p_rnd,Shape* p_s,int step){
        for(int i=0;i<p_s->GetFace().size();i++){
            if(p_s->GetFace(i)->GetVertex().size()!=3){
                std::cout<<"Face-"<<i<<" has "<<p_s->GetFace(i)->GetVertex().size()<<std::endl;
                exit(1);
            }
        }
        if(CELL_TIME_MODE==1){//各ステップでランダムに1つ選ぶ
            p_s->isConsistent();    
            Vertex* p_v1;
            int id;
            do{
                id=p_rnd->RAND_N(p_s->GetVertex().size());
                p_v1=p_s->GetVertex(id);
            }while(p_v1->GetEdge().size()<2||p_v1->GetBoundaryFlag()==1);//とりあえず3は含めない

            std::cout<<"--------------------------------------"<<std::endl;
            std::cout<<"Cell-"<<id<<" Division started"<<std::endl;
            std::cout<<"--------------------------------------"<<std::endl;
            p_s->GetVertex(id)->sortCounterClockwise();//並び変え
            DivideVoronoi(p_rnd,p_s,id,step);
        }
        else if(CELL_TIME_MODE==2){//細胞周期で管理する用
            //p_s->isConsistent();
            std::vector<std::pair<double,int>> list;//分裂する細胞のリスト<細胞周期,インデックス>

#ifdef _OPENMP
#pragma omp parallel num_threads(NUM_THREAD)
{//この括弧は改行しないとダメ
            std::vector<std::pair<double,int>> list_private;
            #pragma omp for nowait// for文をスレッドで分割して,終了段階で待たない
            for(int i=0;i<p_s->GetVertex().size();i++){
                double cell_time=DELTA_TIME*CELL_DIVISION_STEP_AVE;//分裂の基準となる時間
                if(DivisionDistributionMode==1){
                    //分裂に空間分布を与える
                    double r_dis=(p_s->GetVertex(i)->GetLoc(0).GetX()*p_s->GetVertex(i)->GetLoc(0).GetX()+p_s->GetVertex(i)->GetLoc(0).GetY()*p_s->GetVertex(i)->GetLoc(0).GetY());
                    cell_time+=DIV_alpha*DELTA_TIME*CELL_DIVISION_STEP_AVE*std::pow(r_dis,DIV_beta);//分裂の基準となる時間
                }
                if(p_s->GetVertex(i)->GetBoundaryFlag()==0&&p_s->GetVertex(i)->GetFlagRec()==0&&p_s->GetVertex(i)->GetCellTime()>=p_rnd->RAND_Normal(cell_time,CELL_DIVISION_STEP_SD*DELTA_TIME)){    
                    std::pair<double,int>p=std::make_pair(p_s->GetVertex(i)->GetCellTime()-cell_time,i);//境界でない点限定
                    list_private.push_back(p);
                }
            }
            if(list_private.size()!=0){
                #pragma omp critical//クリティカルに設定
                list.insert(list.end(),list_private.begin(),list_private.end());
                //#pragma omp barrier
            }
}
#endif
#ifndef _OPENMP
            for(int i=0;i<p_s->GetVertex().size();i++){
                double cell_time=DELTA_TIME*CELL_DIVISION_STEP_AVE;//分裂の基準となる時間
                if(DivisionDistributionMode==1){
                    //分裂に空間分布を与える
                    double r_dis=(p_s->GetVertex(i)->GetLoc(0).GetX()*p_s->GetVertex(i)->GetLoc(0).GetX()+p_s->GetVertex(i)->GetLoc(0).GetY()*p_s->GetVertex(i)->GetLoc(0).GetY());
                    cell_time+=DIV_alpha*DELTA_TIME*CELL_DIVISION_STEP_AVE*std::pow(r_dis,DIV_beta);//分裂の基準となる時間
                }
                if(p_s->GetVertex(i)->GetBoundaryFlag()==0&&p_s->GetVertex(i)->GetFlagRec()==0&&p_s->GetVertex(i)->GetCellTime()>=p_rnd->RAND_Normal(cell_time,CELL_DIVISION_STEP_SD*DELTA_TIME)){    
                //if(p_s->GetVertex(i)->GetBoundaryFlag()==0&&p_s->GetVertex(i)->GetApicalFlag()==0&&p_s->GetVertex(i)->GetFlagRec()==0&&p_s->GetVertex(i)->GetCellTime()>=cell_time){  
                    std::pair<double,int>p=std::make_pair(p_s->GetVertex(i)->GetCellTime()-cell_time,i);//境界でない点限定
                    list.push_back(p);
                }
            }
#endif
            std::sort(list.begin(),list.end());//ソートする
            std::reverse(list.begin(),list.end());//要素の並びを逆にする(降順)
            for(int i=0;i<list.size();i++){
                std::cout<<"--------------------------------------"<<std::endl;
                std::cout<<"Cell-"<<list[i].second<<" Division started"<<std::endl;
                std::cout<<"--------------------------------------"<<std::endl;
                //p_s->GetVertex(list[i].second)->sortCounterClockwise();
                p_s->DivisionCount(0);
                if(!DivideVoronoi(p_rnd,p_s,list[i].second,step)){
                    continue;//分裂できなかったとき用
                }
                p_s->DivisionCount(1);
                p_s->GetVertex(list[i].second)->SetCellTime(list[i].first);
                if(FlagRecMode==1)p_s->GetVertex(list[i].second)->SetFlagRec();//一度分裂したものは分裂しない設定用
            }
        }
    }
    //取り敢えず境界点以外
    bool DivideVoronoi(Random* p_rnd,Shape* p_s,int num,int step){//numは分裂する頂点のインデックス,stepは現在のステップ
        Vertex* p_v=p_s->GetVertex(num);
        //EvaluateNormalVector(p_s,p_v);
        std::ofstream ofs("d_Division.dat",std::ios::app);
        ofs<<"step-"<<step<<"\n";
        ofs<<"vertex-"<<num<<"\n";
        for(int i=0;i<p_v->GetNextVertex().size();i++){
            //std::cout<<"next_vertex["<<"]="<<p_s->periodize_vector(p_s->GetVertex(p_v->GetNextVertex(i))->GetLoc(0)-p_v->GetLoc(0))<<std::endl;
            ofs<<"next_vertex["<<i<<"]="<<p_s->periodize_vector(p_s->GetVertex(p_v->GetNextVertex(i))->GetLoc(0)-p_v->GetLoc(0))<<"\n";
        }
        p_v->CalcNormalVector(0);
        Vector3<double> n=p_v->GetNormalVector();//p_vでの法線ベクトル
        ofs<<"normal_vector(pre)="<<n<<"\n";
        if(!isInnerVertex(p_v,n)){//内部にないとき
            p_s->DivisionCount(2);
            OutputOFF(p_s,step,num,2);
            return false;
        }
        //std::cout<<"NormalVector="<<n<<std::endl;
        std::vector<Vector3<double>> voronoi_list;//重心のlist(p_vを基準にn方向から投影したもの)
        for(int i=0;i<p_v->GetFace().size();i++){//p_v周りのすべての面に対して
            /*
            for(int j=0;j<3;j++){//中点も入れていたとき用(計算量が増える)
                if(p_s->GetFace(p_v->GetFace(i))->GetVertex(j)==num){
                    voronoi_list.push_back(p_s->periodize_vector(p_s->periodize_vector(p_s->periodize_vector(p_v->GetLoc()+p_s->GetVertex(p_s->GetFace(p_v->GetFace(i))->GetVertex((j+1)%3))->GetLoc())/2.0)-p_v->GetLoc()).Projection(n));
                    break;
                }
            }
            */
            
            p_s->GetFace(p_v->GetFace(i))->CalcCenter(0);
            //voronoi_list.push_back(p_s->GetFace(p_v->GetFace(i))->GetCenter()-p_v->GetLoc()-(n*(p_s->GetFace(p_v->GetFace(i))->GetCenter()-p_v->GetLoc()))*p_s->GetFace(p_v->GetFace(i))->GetNormalVector()/(n*p_s->GetFace(p_v->GetFace(i))->GetNormalVector()));
            voronoi_list.push_back(p_s->periodize_vector(p_s->GetFace(p_v->GetFace(i))->GetCenter()-p_v->GetLoc(0)).Projection(n));
        }  
        //EvaluateCentroid(p_s,p_v,voronoi_list);
        
        for(int i=0;i<voronoi_list.size();i++)ofs<<"voronoi_list="<<voronoi_list[i]<<"\n";
        Vector3<double> axis;//軸のベクトル
        
        if(CELL_DIVISION_MODE==1){//ランダムに分裂軸を選ぶ
            double theta;
            theta=p_rnd->RAND_R(0.0,2*M_PI);//ランダムに選ぶ(0~2PI)
            std::cout<<"axis="<<(int)(theta*180/M_PI)<<"°"<<std::endl;//軸
            axis=(voronoi_list[0]/voronoi_list[0].CalcNorm()).Rotation(n,theta);//単位ベクトル
        }
        else if(CELL_DIVISION_MODE==2){//細胞の最長方向
            //多角形だから、それぞれの頂点を通る軸から交点間距離の最大を選べば十分であろうという前提（ラフに求める）
            std::vector<Vector3<double>> tmp_voronoi;//axisの候補を入れる
            std::vector<double> dis;//交点間の距離
            for(int i=0;i<voronoi_list.size();i++){//軸を決めるためのループ
                Vector3<double> tmp=voronoi_list[i]/voronoi_list[i].CalcNorm();//voronoi頂点方向単位ベクトル
                tmp_voronoi.push_back(tmp);//tmpがaxisの候補
                std::vector<double>tmp_dis;//各頂点を選択したときの交点間距離(1つとは限らないから念のため)
                for(int j=0;j<voronoi_list.size();j++){//軸に対してもう一つの交点を出すためのループ
                    if((tmp%voronoi_list[(j+1)%voronoi_list.size()])*(tmp%voronoi_list[j])<1e-8){
                        double v=(tmp%(voronoi_list[j])).CalcNorm()/((tmp%(voronoi_list[(j+1)%voronoi_list.size()]-voronoi_list[j])).CalcNorm()+1e-5);//交点へのベクトルの比例係数
                        Vector3<double> vv=voronoi_list[j]+v*(voronoi_list[(j+1)%voronoi_list.size()]-voronoi_list[j]);//交点のベクトル
                        tmp_dis.push_back((voronoi_list[i]-vv).CalcNorm());
                    }
                }
                if(tmp_dis.size()==0){//交点が二つないとき(正しく判定できない)
                    std::cout<<"i="<<i<<" tmp_axis="<<tmp<<std::endl;
                    std::cout<<"tmp_axis has only one crosspoint."<<std::endl;
                    p_s->DivisionCount(3);//エラーログ(3)
                    OutputOFF(p_s,step,num,3);//そのときのOFFファイル
                    return false;
                }
                auto itr=std::max_element(tmp_dis.begin(),tmp_dis.end());//それぞれの頂点を通る軸でもう一つの交点との距離が最大のもの(予選)
                dis.push_back(*itr);
                tmp_dis.clear();
            }
            if(dis.size()<tmp_voronoi.size()){
                std::cout<<"restructure:cellDivisionMode 2 error."<<std::endl;
                exit(1);
            }
            auto max_iter=std::max_element(dis.begin(),dis.end());//頂点を通る最大の線の中で最大を選ぶ(本選)
            std::size_t index=std::distance(dis.begin(),max_iter);
            axis=tmp_voronoi[index];
        }
        else if(CELL_DIVISION_MODE==3){//x方向分裂軸
            double nx=n.GetX();
            double nz=n.GetZ();
            if(nx==0&&nz==0){//x軸に垂直のときはy軸方向に与える
                axis=Vector3<double>(0.0,1.0,0.0);
            }
            else{
                axis=Vector3<double>(nz,0,-nx);
                axis/=axis.CalcNorm();
            }
            axis=axis.Projection(n);
        }
        else if(CELL_DIVISION_MODE==4){//円周方向分裂軸
            double nx=n.GetX();
            double ny=n.GetY();
            double nz=n.GetZ();
            double x1;
            double y1;
            double z1;
            if(nz!=0){//x1とy1は座標に依存するように与える
                x1=p_v->GetLoc(0).GetX();
                y1=p_v->GetLoc(0).GetY();
                z1=-(nx*x1+ny*y1)/nz;
                axis=Vector3<double>(x1,y1,z1);
                axis=axis.Projection(n);
                axis/=axis.CalcNorm();
            }
            else{
                double theta=p_rnd->RAND_R(0.0,2*M_PI);//ランダムに選ぶ(0~2PI)
                std::cout<<"axis="<<(int)(theta*180/M_PI)<<"°"<<std::endl;//軸
                axis=(voronoi_list[0]/voronoi_list[0].CalcNorm()).Rotation(n,theta);//単位ベクトル
            }
        }
        std::cout<<"axis="<<axis<<std::endl;
        ofs<<"axis="<<axis<<"\n";

        /*1.ボロノイ領域との交点*/
        //交点のリスト(軸とvoronoi多角形との交点)
        std::cout<<"---------------------------"<<std::endl;
        std::cout<<"voronoi and axis crosspoint"<<std::endl;
        std::cout<<"---------------------------"<<std::endl;
        std::vector<std::pair<int,Vector3<double>>> crosspoint;//2つあるはず(境界かどうかに関わらず)

        for(int i=0;i<voronoi_list.size();i++){
            //int number=i/2;//(voronoi点に中点を入れていたとき用)
            int number=i;//一応そのまま使えるように
            if((axis%voronoi_list[(i+1)%voronoi_list.size()])*(axis%voronoi_list[i])<=0){
                if((voronoi_list[i]-voronoi_list[(i+1)%voronoi_list.size()]).CalcNorm()<1e-4){//voronoi点同士が近すぎるときは飛ばす
                    continue;//同じ点があったとき用
                }
                double s=(axis%voronoi_list[i]).CalcNorm()/((axis%(voronoi_list[i]-voronoi_list[(i+1)%voronoi_list.size()])).CalcNorm());//交点の位置を表す
                Vector3<double> aa=voronoi_list[i]+s*(voronoi_list[(i+1)%voronoi_list.size()]-voronoi_list[i]);//交点のベクトル
                if(crosspoint.size()==0){//0のときはそのまま入れる
                    std::pair<int,Vector3<double>>p1=std::make_pair(number,aa);
                    crosspoint.push_back(p1);
                }
                else{
                    Vector3<double> vv1=crosspoint[0].second;
                    Vector3<double> vv2=crosspoint[crosspoint.size()-1].second;
                    if(((vv1-aa).CalcNorm()>1e-3&&crosspoint[0].first!=number)&&((vv2-aa).CalcNorm()>1e-3&&crosspoint[crosspoint.size()-1].first!=number)){//同じ点があったとき用(voronoi点が近くにあるとき,もしくは,数値誤差でダブっているとき)
                        std::pair<int,Vector3<double>>p1 =std::make_pair(number,aa);
                        crosspoint.push_back(p1);
                    }     
                }
                //p_vに対する相対ベクトルであることに注意
            }
        }//分割点の計算が完了

        //std::cout<<"crosspoint.size()="<<crosspoint.size()<<std::endl;
        ofs<<"crosspoint.size()="<<crosspoint.size()<<"\n";
        for(int i=0;i<crosspoint.size();i++){
            std::cout<<"crosspoint["<<i<<"]=<"<<crosspoint[i].first<<","<<crosspoint[i].second<<">"<<std::endl;
            ofs<<"crosspoint="<<crosspoint[i].first<<","<<crosspoint[i].second<<"\n";
        }        

        if(crosspoint.size()>2){//交点の個数が多いとき
            std::vector<double> s_list;
            for(int i=0;i<crosspoint.size();i++){
                s_list.push_back(crosspoint[i].second*axis);//内積を入れる
            }
            auto itr_max=std::max_element(s_list.begin(),s_list.end());
            auto itr_min=std::min_element(s_list.begin(),s_list.end());

            std::pair<int,Vector3<double>>max_vec=crosspoint[std::distance(s_list.begin(),itr_max)];//軸を表すベクトルとの内積が最大のもの
            std::pair<int,Vector3<double>>min_vec=crosspoint[std::distance(s_list.begin(),itr_min)];//軸のベクトルとの内積が最小のもの

            crosspoint.clear();
            crosspoint.push_back(max_vec);
            crosspoint.push_back(min_vec);
        }
        else if(crosspoint.size()<2){//基本的にどうしようもない(射影したときにvoronoi点が周りを囲んでいないということ)
            p_s->DivisionCount(4);
            OutputOFF(p_s,step,num,4);
            std::cout<<"crosspoint.size()<2"<<std::endl;
            std::cout<<"stop division"<<std::endl;
            std::cout<<"--------------------"<<std::endl;
            return false;
        }
        assert(crosspoint.size()==2);//2でないときはerror

        if((crosspoint[0].second-crosspoint[1].second)*axis<0){//axis方向に1→0と並ぶように
            //axisに対してどちら側であるかの位置関係を決めておく
            using std::swap;
            swap(crosspoint[0],crosspoint[1]);//入れ替え
            /*
            std::pair<int,Vector3<double>> pp=crosspoint[0];
            crosspoint[0]=crosspoint[1];
            crosspoint[1]=crosspoint[0];
            */
        }

        //元々の隣の点(next_vertex)のリスト
        std::vector<Vector3<double>> v_next;//p_vに対する相対ベクトルをnで射影する
        for(int i=0;i<p_v->GetNextVertex().size();i++){
            v_next.push_back(p_s->periodize_vector(p_s->GetVertex(p_v->GetNextVertex(i))->GetLoc(0)-p_v->GetLoc(0)).Projection(n));//push_backする
        }
        for(int i=0;i<v_next.size();i++){
            ofs<<"v_next="<<v_next[i]<<"\n";
        }

        /*-----------------境界でないとき限定------------------------------*/
        
        if(p_v->GetBoundaryFlag()==0){

            //新しい細胞の場所を設定する
            double length1=crosspoint[0].second.CalcNorm();//交点1側の距離
            double length2=crosspoint[1].second.CalcNorm();//交点2側の距離
            double length_c=(length1+length2)/(3.0);//交点間を三分割したときの距離

            //分裂軸を1:1:1に分ける
            Vector3<double> location1=p_s->periodize_vector(crosspoint[0].second+(p_s->periodize_vector(crosspoint[1].second-crosspoint[0].second)/(double)3.0));//0側
            Vector3<double> location2=p_s->periodize_vector(crosspoint[1].second+(p_s->periodize_vector(crosspoint[0].second-crosspoint[1].second)/(double)3.0));//1側

            location1+=p_v->GetLoc(0);//新しい点の位置
            location2+=p_v->GetLoc(0);//(確定)
            location1=p_s->periodize_vector(location1);
            location2=p_s->periodize_vector(location2);
            ofs<<"location1="<<location1<<"\n";
            ofs<<"location2="<<location2<<"\n";
            
            if(p_s->periodize_vector(location1-location2).CalcNorm()<0.1*BASE_LENGTH){
                p_s->DivisionCount(5);
                OutputOFF(p_s,step,num,5);
                std::cout<<"location1 and location2 are too close."<<std::endl;
                std::cout<<"-------------------------------------"<<std::endl;
                return false;
            }
            

            //3.ネットワーク構成用に領域分割を行う(2021/10/21ここから改善)
            std::vector<int> area1;//軸で分割される領域の一方に存在する点のローカルのインデックスを保持
            std::vector<int> area2;//軸の反対側
            std::vector<std::pair<double,int>> area1_network;//条件の判定を行うためのリスト
            std::vector<std::pair<double,int>> area2_network;
            int v_boundary[2];//新しい2つの点の両方に繋がる2つの点を決める
            for(int i=0;i<v_next.size();i++){//v_nextはVector3<double>のlist
                if((axis%v_next[i])*n<=0){//隣の頂点とaxisの外積とnとの内積の正負で判断
                    area1.push_back(i);
                }
                else{
                    area2.push_back(i);
                }
            }
            ofs<<"area1.size="<<area1.size()<<"\n";
            ofs<<"area2.size="<<area2.size()<<"\n";

            //新しくできる三角形で問題があるやつはそもそも候補から除外する(二面角)+角度の差が小さい順に並び変えるようのlistを作成
            for(int i:area1){//iはローカルインデックス
                Vector3<double> b1=p_s->periodize_vector(p_s->GetVertex(p_v->GetNextVertex((i+1)%p_v->GetNextVertex().size()))->GetLoc(0)-location1);
                Vector3<double> b2=p_s->periodize_vector(p_s->GetVertex(p_v->GetNextVertex((i-1+p_v->GetNextVertex().size())%p_v->GetNextVertex().size()))->GetLoc(0)-location2);

                Vector3<double> a1=p_s->periodize_vector(p_s->GetVertex(p_v->GetNextVertex(i))->GetLoc(0)-location1);
                Vector3<double> a2=p_s->periodize_vector(location2-location1);
                Vector3<double> a3=p_s->periodize_vector(p_s->GetVertex(p_v->GetNextVertex(i))->GetLoc(0)-location2);

                Vector3<double> n1=-a1%a2;//新しくできる三角形の法線
                Vector3<double> n2=a1%b1;
                Vector3<double> n3=b2%a3;
                if(!isDihedralAngleCorrect(n1,n2)){
                    p_s->DivisionCount(6);
                    //OutputOFF(p_s,step,num,6);
                    continue;
                }
                if(!isDihedralAngleCorrect(n1,n3)){
                    p_s->DivisionCount(6);
                    //OutputOFF(p_s,step,num,6);
                    continue;
                }
                double cos_1=(a1*a2)/(a1.CalcNorm()*a2.CalcNorm());
                double cos_2=-(a3*a2)/(a3.CalcNorm()*a2.CalcNorm());
                double cos_3=(a1*a3)/(a1.CalcNorm()*a3.CalcNorm());
                /*
                cos150=-0.866
                cos160=-0.94
                cos165=-0.966
                cos170=-0.985
                cos175=-0.996
                */
               
                //if(cos_3<-0.886)continue;
                if(cos_3<-0.866||cos_3>0.986)continue;
                
                std::pair<double,int> pa=std::make_pair(fabs(cos_1-cos_2),i);
                area1_network.push_back(pa);
            }
            for(int i:area2){
                Vector3<double> b1=p_s->periodize_vector(p_s->GetVertex(p_v->GetNextVertex((i+1)%p_v->GetNextVertex().size()))->GetLoc(0)-location2);
                Vector3<double> b2=p_s->periodize_vector(p_s->GetVertex(p_v->GetNextVertex((i-1+p_v->GetNextVertex().size())%p_v->GetNextVertex().size()))->GetLoc(0)-location1);

                Vector3<double> a1=p_s->periodize_vector(p_s->GetVertex(p_v->GetNextVertex(i))->GetLoc(0)-location1);
                Vector3<double> a2=p_s->periodize_vector(location2-location1);
                Vector3<double> a3=p_s->periodize_vector(p_s->GetVertex(p_v->GetNextVertex(i))->GetLoc(0)-location2);

                Vector3<double> n1=a1%a2;//新しくできる三角形の法線
                Vector3<double> n2=a3%b1;
                Vector3<double> n3=b2%a1;
                if(!isDihedralAngleCorrect(n1,n2)){
                    p_s->DivisionCount(6);
                    //OutputOFF(p_s,step,num,6);
                    continue;
                }
                if(!isDihedralAngleCorrect(n1,n3)){
                    p_s->DivisionCount(6);
                    //OutputOFF(p_s,step,num,6);
                    continue;
                }
                
                double cos_1=(a1*a2)/(a1.CalcNorm()*a2.CalcNorm());
                double cos_2=-(a3*a2)/(a3.CalcNorm()*a2.CalcNorm());
                double cos_3=(a1*a3)/(a1.CalcNorm()*a3.CalcNorm());
                //if(cos_3<-0.866)continue;
                if(cos_3<-0.866||cos_3>0.996)continue;
                
                
                std::pair<double,int> pa=std::make_pair(fabs(cos_1-cos_2),i);
                area2_network.push_back(pa);
            }
            //それぞれのfirstを並べ替える(降順)

            //もともとある三角形の妥当性判断------
            
            std::sort(area1_network.begin(),area1_network.end());
            std::sort(area2_network.begin(),area2_network.end());

            /*
            std::vector<int> index_list[area1_network.size()*area2_network.size()];//indexのリスト
            for(int i=0;i<index_list.size();i++){
                index_list[i]=i;
            }
            std::vector<double,std::pair<int,int>> candidate[area1_network.size()*area2_network.size()];
            */
            int list1[p_v->GetEdge().size()]={0};//0で初期化(0,1,2,3の値を取る)
            double theta_thres=-0.866;//(~150°)
            for(int i=0;i<p_v->GetEdge().size();i++){
                Vector3<double> c1=p_s->periodize_vector(p_s->GetVertex(p_v->GetNextVertex((i+1)%p_v->GetEdge().size()))->GetLoc(0)-p_s->GetVertex(p_v->GetNextVertex(i))->GetLoc(0));
                Vector3<double> c2=p_s->periodize_vector(p_s->GetVertex(p_v->GetNextVertex(i))->GetLoc(0)-location1);
                Vector3<double> c3=p_s->periodize_vector(p_s->GetVertex(p_v->GetNextVertex((i+1)%p_v->GetEdge().size()))->GetLoc(0)-location1);
                if((-c2*c1/(c1.CalcNorm()*c2.CalcNorm()))<theta_thres||(c2*c3/(c2.CalcNorm()*c3.CalcNorm()))<theta_thres||(c1*c3/(c1.CalcNorm()*c3.CalcNorm()))<theta_thres){
                    list1[i]+=1;//1がダメ
                }

                Vector3<double> c4=p_s->periodize_vector(p_s->GetVertex(p_v->GetNextVertex(i))->GetLoc(0)-location2);
                Vector3<double> c5=p_s->periodize_vector(p_s->GetVertex(p_v->GetNextVertex((i+1)%p_v->GetEdge().size()))->GetLoc(0)-location2);
                if((-c4*c1/(c1.CalcNorm()*c4.CalcNorm()))<theta_thres||(c4*c5/(c4.CalcNorm()*c5.CalcNorm()))<theta_thres||(c1*c5/(c1.CalcNorm()*c5.CalcNorm()))<theta_thres){
                    list1[i]+=2;//2がダメ
                }
            }
            for(int i=0;i<p_v->GetEdge().size();i++){
                //std::cout<<"list[i]="<<list1[i]<<std::endl;
                if(list1[i]==3){
                    std::cout<<"There is a big angle triangle."<<std::endl;
                    std::cout<<"No choice is left."<<std::endl;
                    p_s->DivisionCount(7);
                    //OutputOFF(p_s,step,num,7);
                    return false;
                }
            }



            bool result1=false;
                //候補
            for(int i=0;i<area1_network.size();i++){//ダメな配置を回避するような選び方ができていればOK
                for(int j=0;j<area2_network.size();j++){
                    bool result3=false;
                    int can1=area1_network[i].second;
                    int can2=area2_network[j].second;
                    //region1
                    int k=can1;
                    while(k%p_v->GetEdge().size()!=can2){
                        //std::cout<<k<<std::endl;
                        if(list1[k%p_v->GetEdge().size()]==1){//region1につなげるのがダメなものがあるとき
                            result3=true;
                            break;
                        }
                        k++;
                    }
                    //std::cout<<result3<<std::endl;
                    if(result3==true){
                        continue;
                    }
                    //region2
                    k=can2;
                    while(k%p_v->GetEdge().size()!=can1){
                        //std::cout<<k<<std::endl;
                        if(list1[k%p_v->GetEdge().size()]==2){//region2につなげるのがダメなものがあるとき
                            result3=true;
                            break;
                        }
                        k++;
                    }
                    if(result3==true){
                        continue;
                    }
                    result1=true;//2つともclearしたときのみresult1はtrueになる
                    v_boundary[0]=can1;
                    v_boundary[1]=can2;
                    break;
                }
                if(result1){
                    break;
                }
            }
            //std::cout<<result1<<std::endl;
            if(!result1){//すべて失敗したときはfalse
                std::cout<<"There is a big angle triangle(2)."<<std::endl;
                p_s->DivisionCount(7);
                //OutputOFF(p_s,step,num,7);
                return false;
            }

            //これより上側を変更(2021/10/21)
            assert(v_boundary[0]!=v_boundary[1]);
            std::cout<<"v_boundary[0]="<<v_boundary[0]<<std::endl;
            std::cout<<"v_boundary[1]="<<v_boundary[1]<<std::endl;
            ofs<<"v_boundary[0]="<<v_boundary[0]<<"\n";
            ofs<<"v_boundary[1]="<<v_boundary[1]<<"\n";
            
            std::cout<<"n1="<<p_v->GetNextVertex(v_boundary[0])<<" n2="<<p_v->GetNextVertex(v_boundary[1])<<std::endl;
            ofs<<"n1="<<p_v->GetNextVertex(v_boundary[0])<<" n2="<<p_v->GetNextVertex(v_boundary[1])<<"\n";

            /*----------v_boundaryを基準に2つのregionに分割する---------*/

            //分割用のvector
            //vi,eiは境界の点を含まないが,fiは分け方が異なる(面の性質上仕方ない)
            std::vector<int> region1vi;
            std::vector<int> region2vi;
            std::vector<int> region1ei;
            std::vector<int> region2ei;
            std::vector<int> region1fi;//region1側の面
            std::vector<int> region2fi;//region2側の面のリスト

            //push_back(v_boundary[0]からv_boundary[1]まで反時計回り)region1
            region1fi.push_back(p_v->GetFace(v_boundary[0]));
            //面だけ余分に二つ必要
            int i=v_boundary[0]+1;
            while(i%p_v->GetEdge().size()!=v_boundary[1]){
                region1vi.push_back(p_v->GetNextVertex(i%p_v->GetEdge().size()));
                region1ei.push_back(p_v->GetEdge(i%p_v->GetEdge().size()));
                region1fi.push_back(p_v->GetFace(i%p_v->GetFace().size()));
                i++;
            }
            //push_back(v_boundary[1]からv_boundary[0]まで反時計回り)region2
            region2fi.push_back(p_v->GetFace(v_boundary[1]));
            i=v_boundary[1]+1;
            while(i%p_v->GetEdge().size()!=v_boundary[0]){
                region2vi.push_back(p_v->GetNextVertex(i%p_v->GetNextVertex().size()));
                region2ei.push_back(p_v->GetEdge(i%p_v->GetEdge().size()));
                region2fi.push_back(p_v->GetFace(i%p_v->GetFace().size()));
                i++;
            }
            assert(region1fi.size()+region2fi.size()==p_v->GetFace().size());
            assert(region1vi.size()+region2vi.size()+2==p_v->GetNextVertex().size());

            Vector3<double> region1Place,region2Place;
            region1Place=p_s->periodize_vector(location1);
            region2Place=p_s->periodize_vector(location2);
            
            int n1=p_v->GetNextVertex(v_boundary[0]);//v_boundaryはローカル,n1,n2はグローバルなインデックス
            int n2=p_v->GetNextVertex(v_boundary[1]);

            std::cout<<"region1="<<region1Place<<std::endl;
            std::cout<<"region2="<<region2Place<<std::endl;
            std::cout<<"p_v("<<num<<")="<<p_v->GetLoc(0)<<std::endl;
            ofs<<"p_v("<<num<<")="<<p_v->GetLoc(0)<<"\n";


            //分裂前の細胞(p_v)からregion2の要素を消す
            for(int vi2:region2vi){
                p_v->findAndErase_vi_next(vi2);
            }
            for(int ei2:region2ei){
                p_v->findAndErase_ei(ei2);
            }
            for(int fi2:region2fi){
                p_v->findAndErase_fi(fi2);
            }

            /*対応関係の変更と追加*/
            //region1は元の細胞を移動させて対応
            p_v->SetLoc(region1Place,0);//座標を設定

            //新しく作る要素のインデックス(点1つ、面2つ、線3つ)
            int vi_itr=p_s->GetVertex().size();
            int li_itr=p_s->GetEdge().size();
            int li_itr2=li_itr+1;
            int li_itr3=li_itr+2;
            int fi_itr=p_s->GetFace().size();//n1側
            int fi_itr2=fi_itr+1;//n2側のインデックス

            //新しい点-コンストラクタを用いて宣言
            Vertex* tmp_v=new Vertex(region2Place.GetX(),region2Place.GetY(),region2Place.GetZ(),p_s);
            p_s->PushVertex(tmp_v);
            p_s->GetVertex(vi_itr)->SetApicalFlag(p_s->GetVertex(num)->GetApicalFlag());
            p_s->GetVertex(vi_itr)->SetColorFlag(p_s->GetVertex(num)->GetColorFlag());
            p_s->GetVertex(vi_itr)->SetLayerFlag(p_s->GetVertex(num)->GetLayerFlag());//layerのフラグを引き継ぐ

            //新しい辺3つ
            Edge* tmp_e1=new Edge(vi_itr,num,p_s);
            Edge* tmp_e2=new Edge(vi_itr,n1,p_s);
            Edge* tmp_e3=new Edge(vi_itr,n2,p_s);
            p_s->PushEdge(tmp_e1);
            p_s->PushEdge(tmp_e2);
            p_s->PushEdge(tmp_e3);
            //定数の設定
            p_s->SetBaseConstant(p_s->GetEdge(li_itr));
            p_s->SetBaseConstant(p_s->GetEdge(li_itr2));
            p_s->SetBaseConstant(p_s->GetEdge(li_itr3));

            //面の頂点を反時計回りに配置する
            int fi1=region1fi[0];//n1を含むregion1側の面
            int fi2=region1fi[region1fi.size()-1];//n2を含むregion1側の面
            //隣の面
            Face* tmp_f1;//n1側
            for(int i=0;i<3;i++){
                if(p_s->GetFace(fi1)->GetVertex(i)==num&&p_s->GetFace(fi1)->GetVertex((i+1)%3)==n1){
                    tmp_f1=new Face(p_s->GetFace(fi1)->GetVertex((i+1)%3),p_s->GetFace(fi1)->GetVertex(i),vi_itr,p_s);//隣り合う面同士は共通する点の順番が逆
                    break;
                }
                else if(p_s->GetFace(fi1)->GetVertex(i)==n1&&p_s->GetFace(fi1)->GetVertex((i+1)%3)==num){
                    tmp_f1=new Face(p_s->GetFace(fi1)->GetVertex((i+1)%3),p_s->GetFace(fi1)->GetVertex(i),vi_itr,p_s);//隣り合う面同士は共通する点の順番が逆
                    //順番が重要
                    break;
                }
            }
            Face* tmp_f2;//n2側
            for(int i=0;i<3;i++){
                if(p_s->GetFace(fi2)->GetVertex(i)==num&&p_s->GetFace(fi2)->GetVertex((i+1)%3)==n2){
                    tmp_f2=new Face(p_s->GetFace(fi2)->GetVertex((i+1)%3),p_s->GetFace(fi2)->GetVertex(i),vi_itr,p_s);//隣り合う面同士は共通する点の順番が逆
                    break;
                }
                else if(p_s->GetFace(fi2)->GetVertex(i)==n2&&p_s->GetFace(fi2)->GetVertex((i+1)%3)==num){
                    tmp_f2=new Face(p_s->GetFace(fi2)->GetVertex((i+1)%3),p_s->GetFace(fi2)->GetVertex(i),vi_itr,p_s);//隣り合う面同士は共通する点の順番が逆
                    break;
                }
            }
            //pushする
            p_s->PushFace(tmp_f1);
            p_s->PushFace(tmp_f2);
            //定数の設定
            p_s->SetBaseConstant(p_s->GetFace(fi_itr));
            p_s->SetBaseConstant(p_s->GetFace(fi_itr2));

            //角度
            AngleRestraint(p_s,p_s->GetFace(fi_itr));
            AngleRestraint(p_s,p_s->GetFace(fi_itr2));
            
            //新しく作った要素の対応関係の設定
            //新しい頂点に辺を対応付ける
            for(int eidx:region2ei){
                p_s->GetVertex(vi_itr)->PushEdge(eidx);
            }
            p_s->GetVertex(vi_itr)->PushEdge(li_itr);
            p_s->GetVertex(vi_itr)->PushEdge(li_itr2);
            p_s->GetVertex(vi_itr)->PushEdge(li_itr3);
            //新しい頂点に面を対応付ける
            for(int fidx:region2fi){
                p_s->GetVertex(vi_itr)->PushFace(fidx);
            }
            p_s->GetVertex(vi_itr)->PushFace(fi_itr);
            p_s->GetVertex(vi_itr)->PushFace(fi_itr2);
            //新しい頂点に隣の頂点を対応付ける
            p_s->GetVertex(vi_itr)->PushNextVertex(num);
            for(int i:region2vi){
                p_s->GetVertex(vi_itr)->PushNextVertex(i);
            }
            p_s->GetVertex(vi_itr)->PushNextVertex(n1);
            p_s->GetVertex(vi_itr)->PushNextVertex(n2);

            //新しい辺に面を対応付ける
            p_s->GetEdge(li_itr)->PushFace(fi_itr);
            p_s->GetEdge(li_itr)->PushFace(fi_itr2);

            p_s->GetEdge(li_itr2)->PushFace(fi_itr);
            p_s->GetEdge(li_itr2)->PushFace(region2fi[region2fi.size()-1]);

            p_s->GetEdge(li_itr3)->PushFace(fi_itr2);
            p_s->GetEdge(li_itr3)->PushFace(region2fi[0]);

            //新しい面に辺を対応付ける
            p_s->GetFace(fi_itr)->PushEdge(li_itr);
            p_s->GetFace(fi_itr)->PushEdge(li_itr2);
            int id1;//edge用のインデックス
            for(int id:p_s->GetFace(fi1)->GetEdge()){
                if((n1==p_s->GetEdge(id)->GetVertex(0)&&num==p_s->GetEdge(id)->GetVertex(1))||(n1==p_s->GetEdge(id)->GetVertex(1)&&num==p_s->GetEdge(id)->GetVertex(0))){
                id1=id;
                }
            }
            p_s->GetFace(fi_itr)->PushEdge(id1);

            p_s->GetFace(fi_itr2)->PushEdge(li_itr);
            p_s->GetFace(fi_itr2)->PushEdge(li_itr3);
            int id2;
            for(int id:p_s->GetFace(fi2)->GetEdge()){
                if((n2==p_s->GetEdge(id)->GetVertex(0)&&num==p_s->GetEdge(id)->GetVertex(1))||(n2==p_s->GetEdge(id)->GetVertex(1)&&num==p_s->GetEdge(id)->GetVertex(0))){
                    id2=id;
                }
            }
            p_s->GetFace(fi_itr2)->PushEdge(id2);

            //既存の要素の対応関係の変更
            //region2の辺
            for(int lidx:region2ei){
                p_s->GetEdge(lidx)->Replace_vi(num,vi_itr);
            }
            //region2の面
            for(int fidx:region2fi){
                p_s->GetFace(fidx)->Replace_vi(num,vi_itr);
                p_s->GetFace(fidx)->Replace_ei(id1,li_itr2);
                p_s->GetFace(fidx)->Replace_ei(id2,li_itr3);
            }
            //region2の点
            for(int vidx:region2vi){
                p_s->GetVertex(vidx)->Replace_vi_next(num,vi_itr);
            }

            //辺id1,id2について
            p_s->GetEdge(id1)->Replace_fi(region2fi[region2fi.size()-1],fi_itr);
            p_s->GetEdge(id2)->Replace_fi(region2fi[0],fi_itr2);

            //vについて
            p_v->PushEdge(li_itr);
            p_v->PushFace(fi_itr2);
            p_v->PushFace(fi_itr);
            p_v->PushNextVertex(vi_itr);
            //n1について
            p_v->p_s->GetVertex(n1)->PushEdge(li_itr2);
            p_v->p_s->GetVertex(n1)->PushFace(fi_itr);
            p_v->p_s->GetVertex(n1)->PushNextVertex(vi_itr);
            
            //n2について
            p_v->p_s->GetVertex(n2)->PushEdge(li_itr3);
            p_v->p_s->GetVertex(n2)->PushFace(fi_itr2);
            p_v->p_s->GetVertex(n2)->PushNextVertex(vi_itr);

            p_s->GetVertex(n1)->sortCounterClockwise();
            p_s->GetVertex(n2)->sortCounterClockwise();
            p_v->sortCounterClockwise();
            p_s->GetVertex(vi_itr)->sortCounterClockwise();

            //二面角の判定を行う(裏返っているものがないか)
            std::vector<int> face_list;
            for(int i:region1fi){
                face_list.push_back(i);
            }
            for(int i:region2fi){
                face_list.push_back(i);
            }
            face_list.push_back(fi_itr);
            face_list.push_back(fi_itr2);

            for(int i:face_list){
                bool result=false;
                for(int j=0;j<3;j++){
                    if(p_s->GetEdge(p_s->GetFace(i)->GetEdge(j))->GetFace().size()!=2){
                        result=true;//境界線であるとき
                        break;
                    }
                    p_s->GetEdge(p_s->GetFace(i)->GetEdge(j))->CalcDihedralAngle(0);
                }
                if(result){
                    continue;
                }
                if(fabs(p_s->GetEdge(p_s->GetFace(i)->GetEdge(0))->GetDihedralAngle()-M_PI)>M_PI/2.0&&fabs(p_s->GetEdge(p_s->GetFace(i)->GetEdge(1))->GetDihedralAngle()-M_PI)>M_PI/2.0&&fabs(p_s->GetEdge(p_s->GetFace(i)->GetEdge(2))->GetDihedralAngle()-M_PI)>M_PI/2.0){
                    //二面角が逆転していると判断
                    /*
                    std::vector<double> leng_list;
                    for(int j=0;j<3;j++){
                        p_s->GetEdge(p_s->GetFace(i)->GetEdge(j))->CalcLength(0);
                        leng_list.push_back(p_s->GetEdge(p_s->GetFace(i)->GetEdge(j))->GetLength());
                    }
                    */
                    //auto max_itera=std::max_element(leng_list.begin(),leng_list.end());
                    //std::size_t inde=std::distance(leng_list.begin(),max_itera);
                    p_s->DivisionCount(8);//エラーログ(5)
                    //OutputOFF(p_s,step,num,8);//そのときのOFFファイル
                    //p_s->GetEdge(p_s->GetFace(i)->GetEdge(inde))->FlipOperation_forge();
                    //std::cout<<"Forge/**********"<<std::endl;
                }
            }

            std::cout<<"CellDivision finished"<<std::endl;
            std::cout<<"Cell-"<<vi_itr<<" has made"<<std::endl;
            std::cout<<"----------------------------"<<std::endl;
            ofs<<"Cell-"<<vi_itr<<" has made"<<std::endl;
            ofs.close();

            //p_s->isConsistent();
        }
        return true;
    }
    bool isInnerVertex(Vertex* p_v,Vector3<double> n){//内部にあるかどうかの判定
        std::vector<Vector3<double>> list;//(nで定まる平面上に投影したも)
        for(int i=0;i<p_v->GetNextVertex().size();i++){//p_v周りのすべての面に対して
            list.push_back(p_v->p_s->periodize_vector(p_v->p_s->GetVertex(p_v->GetNextVertex(i))->GetLoc(0)-p_v->GetLoc(0)).Projection(n));
        }
        double theta=0.0;
        for(int i=0;i<list.size();i++){
            theta+=acos(list[i]*list[(i+1)%list.size()]/(list[i].CalcNorm()*list[(i+1)%list.size()].CalcNorm()));
        }
        if(theta>7*M_PI/4){//角度で判定
            return true;
        }
        else{
            return false;
        }
    }
    bool isDihedralAngleCorrect(Vector3<double> n1,Vector3<double> n2){
        n1/=n1.CalcNorm();
        n2/=n2.CalcNorm();
        if(n1*n2>=0.0){
            return true;
        }
        else{
            return false;
        }
    }

    void OutputOFF(Shape* p_s,int step,int num,int type){//num:頂点のインデックス,type:エラーのタイプ
        if(OUTPUT_DivisionErrorOFFfile==1){
            char fname[100];
            sprintf(fname,"_%u_%u",num,type);
            p_s->OutputOFF(step,fname);
            std::ofstream fout("d_error.dat",std::ios::app);
            fout<<"step-"<<step<<"\n";
            fout<<"vertex-"<<num<<"\n";
            fout<<"type-"<<type<<"\n";
            fout.close();
        }
    }
    
    void EvaluateNormalVector(Shape*p_s,Vertex*p_v){
        std::vector<Vector3<double>> vlist(6,Vector3<double>(0.0,0.0,0.0));
	    for(int i : p_v->GetFace()){
		    p_s->GetFace(i)->CalcArea(0);
		    double alpha;
		    for(int j=0;j<3;j++){
			    if(p_s->GetVertex(p_s->GetFace(i)->GetVertex(j))==p_v){
				    Vector3<double> a= p_s->periodize_vector(p_s->GetVertex(p_s->GetFace(i)->GetVertex((j+2)%3))->GetLoc(0)-p_v->GetLoc(0));
				    Vector3<double> b=p_s->periodize_vector(p_s->GetVertex(p_s->GetFace(i)->GetVertex((j+1)%3))->GetLoc(0)-p_v->GetLoc(0));
				    alpha=acos((a*b)/(a.CalcNorm()*b.CalcNorm()));
				    break;
			    }
		    }
            p_s->GetFace(i)->CalcCenter(0);
		    Vector3<double> g=p_s->GetFace(i)->GetCenter();//重心
	    	double dis=p_s->periodize_vector(g-p_v->GetLoc(0)).CalcNorm();//面の重心との距離
    		p_s->GetFace(i)->CalcNormalVector(0);
            vlist[0] += p_s->GetFace(i)->GetNormalVector();
            vlist[1] += alpha*p_s->GetFace(i)->GetNormalVector();
            vlist[2] += p_s->GetFace(i)->GetArea()*p_s->GetFace(i)->GetNormalVector();
            vlist[3] += p_s->GetFace(i)->GetArea()*p_s->GetFace(i)->GetNormalVector()*alpha;
            vlist[4] += p_s->GetFace(i)->GetNormalVector()/(dis*dis);
            vlist[5] += p_s->GetFace(i)->GetNormalVector()*alpha/(dis);
	    }
        
        for(int k=0;k<vlist.size();k++){
            
            double number=0.0;
            Vector3<double> ren(0.0,0.0,0.0);
            for(int i=0;i<p_v->GetFace().size();i++){//p_v周りのすべての面に対して
                p_s->GetFace(p_v->GetFace(i))->CalcCenter(0);
                ren+=p_s->periodize_vector(p_s->GetFace(p_v->GetFace(i))->GetCenter()-p_v->GetLoc(0)).Projection(vlist[k]);
            }
            ren/=(double)(p_v->GetFace().size());
            number=ren.CalcNorm();
            
            char fname[100];
            sprintf(fname,"d_NormalVector%d.dat",k);
            std::ofstream fout(fname,std::ios::app);
            fout<<number<<"\n";
            fout.close();
        }
    }
    void AngleRestraint(Shape* p_s,Face* p_f){//面の最小角度を分類
        Vector3<double> a=p_s->GetVertex(p_f->GetVertex(0))->GetLoc(0);
        Vector3<double> b=p_s->GetVertex(p_f->GetVertex(1))->GetLoc(0);
        Vector3<double> c=p_s->GetVertex(p_f->GetVertex(2))->GetLoc(0);

        double cos_1=(p_s->periodize_vector(a-c)*p_s->periodize_vector(b-c))/(p_s->periodize_vector(a-c).CalcNorm()*p_s->periodize_vector(b-c).CalcNorm());
        double cos_2=(p_s->periodize_vector(b-a)*p_s->periodize_vector(c-a))/(p_s->periodize_vector(b-a).CalcNorm()*p_s->periodize_vector(c-a).CalcNorm());
        double cos_3=(p_s->periodize_vector(c-b)*p_s->periodize_vector(a-b))/(p_s->periodize_vector(c-b).CalcNorm()*p_s->periodize_vector(a-b).CalcNorm());

        double max_cos=cos_1;
        if(max_cos<cos_2)max_cos=cos_2;
        if(max_cos<cos_3)max_cos=cos_3;

        if(max_cos>0.996)p_s->InternalAngleErrorCount(0);
        else if(max_cos<=0.996&&max_cos>0.985)p_s->InternalAngleErrorCount(1);
        else if(max_cos<=0.985&&max_cos>0.966)p_s->InternalAngleErrorCount(2);
        else if(max_cos<=0.966&&max_cos>0.94)p_s->InternalAngleErrorCount(3);
        else if(max_cos<=0.94&&max_cos>0.866)p_s->InternalAngleErrorCount(4);
    }
    void EvaluateCentroid(Shape*p_s,Vertex*p_v,std::vector<Vector3<double>> voronoi_list){
        //comparison of physical centroid and topological centroid
        //voronoi_listは相対ベクトルのリスト
        Vector3<double> g1(0.0,0.0,0.0);//物理的重心
        Vector3<double> g2(0.0,0.0,0.0);//幾何的重心
        double sum_s=0.0;

        for(int i=1;i<voronoi_list.size();i++){
            g1+=p_s->periodize_vector(voronoi_list[i]-voronoi_list[0])/(double)voronoi_list.size();
        }
        g1=p_s->periodize_vector(g1+voronoi_list[0]);

        Vector3<double> v1=voronoi_list[0];
        for(int i=0;i<voronoi_list.size()-1;i++){
            Vector3<double> gg=(v1+voronoi_list[i]+voronoi_list[(i+1)%voronoi_list.size()])/3.0;
            double s=0.5*((voronoi_list[i]-v1)%(voronoi_list[(i+1)%voronoi_list.size()]-v1))*p_v->GetNormalVector();
            sum_s+=s;
            g2+=s*gg;
        }
        g2/=sum_s;

        double dis_1=g1.CalcNorm();
        double dis_2=g2.CalcNorm();
        double dis_3=(g1*g2)/(dis_1*dis_2);

        std::ofstream fout("d_Centroid1.dat",std::ios::app);
        fout<<dis_1<<"\n";
        fout.close();
        std::ofstream ofs("d_Centroid2.dat",std::ios::app);
        ofs<<dis_2<<"\n";
        ofs.close();
        
        std::ofstream ofst("d_Centroid3.dat",std::ios::app);
        ofst<<dis_3<<"\n";
        ofst.close();
        
    }
}//namespace restructure