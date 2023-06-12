#ifndef CLASS_HPP
#define CLASS_HPP

#include<cassert>
#include<algorithm>
#include<vector>
#include"vector3.hpp"
#include"parameters.hpp"
#include<string>
#include"force.hpp"
#include"Random.hpp"

class Vertex;
class Edge;
class Face;

//構造体の定義(self-collision用)
typedef struct _node{
	std::vector<int> index;
}_node;

//全体のクラス
class Shape{
	private:
		std::vector<Vertex *>p_v;//頂点のクラス
		std::vector<Face*>p_f;//面のクラス
		std::vector<Edge*>p_e;//辺のクラス

		Vector3<double> size_system;//システムの大きさ(x-y平面上)

		std::vector<std::vector<_node> > node_list=std::vector<std::vector<_node> > (Layer_NUM);//indexを保持(各層にて)
		//locの最小値(排除効果を求める際の座標のオフセット分)
		double minx;//x座標の最小値
		double miny;//y座標の最小値
		double minz;//z座標の最小値


		//offsetlocの最大値(オフセット後rcで割った後intでcastしたもの(整数の座標))
		int maxx;
		int maxy;
		int maxz;

		//カウント用
		int flip_count;//flipの回数
		int division_count[9];
		/*cell_division
		0:分裂試行回数
		1:実際の分裂回数
		2:normal_vectorエラー
		3:最長方向分裂軸を与えた時のエラー
		4:crosspointが足りない
		5:2つの点が近すぎるエラー
		6:二面角判定エラー(新しくできる三角形)
		7:接続が不適切判定を受けたエラー回数
		8:二面角不適切判定
		*/
		int dihedral_count;//二面角が不適切なステップの数

		int internal_angle_error_count[5];
		/*4:cos30=0.866,3:cos20=0.94,2:cos15=0.966,1:cos10=0.985,0:cos5=0.996*/

	public:
		Random* p_rnd;

		//エネルギー計算
		double TotalEnergy[15];

		double LineEnergy[NUM_THREAD];
		double AreaEnergy[NUM_THREAD];
		double DihedralAngleEnergy[NUM_THREAD];
		double Z_Energy[NUM_THREAD];
		double Y_Energy[NUM_THREAD];
		double V_Center_Energy[NUM_THREAD];
		double RepulsiveEnergy[NUM_THREAD];
		double ApicalConstrictionEnergy[NUM_THREAD];
		double LateralEnergy[NUM_THREAD];

		double VoronoiAreaEnergy[NUM_THREAD];
		double VoronoiPerimeterEnergy[NUM_THREAD];

		//コンストラクタとデストラクタ
		Shape(Random* p_rnd);
		~Shape();
		void ReadOffFile(std::string fname);//OFFfileからの読み込み
		void EdgeConstruction();//頂点と面の情報から辺を構成
		void AssignNextVertex();//全ての頂点について隣の頂点のインデックスを入れる
		void SetBaseConstant();//定数の設定を行う(平衡、比例定数)
		void SetBaseConstant(Edge*);//辺に関する定数の設定(後から辺を追加した時用)
		void SetBaseConstant(Face*);//面に関する定数の設定(後から面を追加した時用)
		void SetFixAndBoundary();//fixflagとboundaryflagの設定
		void SetLayerFlag();//Layerflagの設定
		void CellTimeInitialize();//細胞周期をランダムに与える
		void CellTimeInitialize(const char*);//細胞周期をランダムに与える
		void sortCounterClockwise();//すべての点を並び変える(並列化用)
		
		//output
		void OutputVTK(int step);//vtkファイルの出力
		void OutputVTK_Voronoi(int step);//voronoi領域をvtkファイルで出力
		void OutputOFF(int);//Offファイルの出力
		void OutputOFF_periodic(int);//周期境界の接続関係を残したOFFファイル出力
		void OutputOFF_Mapping(int step);//loc_mappingを出力
		void OutputOFF(int,char*);
		void Outputarray(const std::vector<int>);//datファイルとして引数の配列の中身を出力
		void OutputFlipCount(int);//カウントデータ
		void OutputDivisionCount(int,int);//分裂エラーデータ
		void OutputParameters();//パラメータ出力
		void OutputCellTime(int step);//細胞時間のtextファイルを出力
		void OutputResult();//結果(カウントデータ)を出力
		void OutputInternalCount(int num,int type);//角度のカウントを出力
		void OutputEnergy(int step);//エネルギーの出力
		void OutputRange(int);//[x,y,z]のrangeを出力
		void OutputDistortion(int);//Rangeの範囲
		
		//カウントを増やす
		void FlipCount();//countを増やす
		void DivisionCount(int);//countを増やす
		void DihedralCount();//二面角のエラー回数
		void InternalAngleErrorCount(int);

		//ちゃんと線、面、点が登録されているかどうか
		bool isConsistent();
		//二面角でおかしいところがないか?
		bool isSuitable(int);

		//計算
		//periodic
		Vector3<double> periodize_vector(Vector3<double>);
		Vector3<double> periodized_center(std::vector<Vector3<double>>);
		bool cross_periodic_boundary(Vector3<double> a);//if vector cross periodic boundary:true
		//cross_periodic_boundary
		void CalcCrossBoundaryFlag();

		//Faceについて
		void CalcFaceNormalVector(int);//すべての面の法線ベクトルを計算
		void CalcVertexNormalVector(int);//すべての頂点の法線ベクトルを計算
		void CalcFaceCentroid(int);//重心を計算(全ての面について)

		void FlipOperation();//flip(すべての辺について)

		//各座標の最小値を計算
		void CalcMin(int);//double
		//offsetlocの計算
		void CalcOffsetLoc(int);//maxも求める(int)
		//collision計算
		void DivideMeshList(int);

		//境界点同士を繋ぐ
		void ConnectBoundaryVertex();
		//境界の面と辺を消す
		void EraseBoundaryFaceAndEdge();
		void findAndErase_p_e(int);
		void findAndErase_p_f(int);
		//細胞周期のアップデート(deltaTimeだけ増やす)
		void updateCellTime();

		//Getter
		std::vector<Vertex*> GetVertex()const;//頂点のlistを得る
		Vertex* GetVertex(int i)const;//頂点のリストのi番目のインデックスを得る
		std::vector<Edge*> GetEdge()const;//辺のリストを得る
		Edge* GetEdge(int i)const;//辺のリストのi番目を得る
		std::vector<Face*> GetFace()const;//面のlistを得る
		Face*GetFace(int i)const;//面のリストのi番目を得る

		//getter
		double GetMINX()const;
		double GetMINY()const;
		double GetMINZ()const;
		
		int GetMAXX()const;
		int GetMAXY()const;
		int GetMAXZ()const;

		std::vector<int> GetNodeIndex(int,int)const;

		Vector3<double> Get_SIZE_SYSTEM()const;//translation vector(when system is periodic)

		//push_back用(std::vectorに追加する)
		void PushVertex(Vertex*);
		void PushEdge(Edge*);
		void PushFace(Face*);

		//近接点のcosを計算し、統計平均を出す
		void CalcProximityPoint();

		//CVT構成アルゴリズムを流用したもの
		void Lloyd_Algorithm();

		//強制変位
		void ForcedDisplacement();
		
		//評価用(誤差関数の計算(分散))
		void VoronoiCenterEvaluation(int step);
		void CellCenterEvaluation(int step);
};
//頂点
class Vertex{
	private:
		
		std::vector<int> ei; //辺のインデックス
		std::vector<int> fi; //細胞のインデックス
		std::vector<int> vi_next; //隣の頂点のインデックス
		
		int FixFlag;//固定する頂点であるかどうか?
		int BoundaryFlag;//境界上の頂点であるかどうか?
		int apical_flag;//SurfaceApicalConstriction用
		int cross_boundary_flag;//周期境界用

		//多層用
		int layer_flag;//どの層であるかのインデックス
	

		int colorFlag;//由来細胞の識別子

		int flag_rec;//分裂したかどうか?
		double cell_time;//細胞周期
		
		Vector3<double> loc[DEGREE_ACCURACY];//座標
		Vector3<double> loc_mapping;//三角形の重心から逆に求めた座標
		Vector3<double> force[DEGREE_ACCURACY];//力
		Vector3<double> frc_thread[NUM_THREAD];//並列計算用

		Vector3<double> externalforce;//外力
		Vector3<double> NormalVector;//面の法線ベクトルの平均で定義する

		double gaussian_curvature;//ガウス曲率
		double mean_curvature;//平均曲率

		double voronoi_area;//voronoiでの面積(法線方向射影)
		double voronoi_perimeter;//voronoiでの周長(法線方向射影)

		Vector3<int> offsetloc;//メッシュに分ける用の座標
	public:

		//識別子
		int id_x;
		int id_y;

		Shape* p_s;//頂点が属するshape
		//コンストラクタ
		Vertex(double x,double y,double z,Shape* _p_s);

		//計算
		void CalcNormalVector(int);//vertexの法線ベクトルを計算
		void CalcNormalVector2(int);//すでに面法線が計算されている前提
		void CalcVoronoiArea(int);//voronoi点の面積を計算
		//double CalcVoronoiArea(Vector3<double> n1,int n);//voronoi点の面積のn1方向射影成分を返す
		void CalcVoronoiPerimeter(int);//voronoi点の周長を計算
		
		//julicher
		void CalcMeanCurvature();//mean-curvature(先に周囲の面積、長さ、二面角を計算しておく)
		//laplace-beltrami-operator
		void CalcMeanCurvature(int);//mean-curvature

		void CalcCurvature();//laplace-beltrami-operator

		void CalcLocMapping();//loc_mappingの計算

		//Offsetの計算
		void CalcOffsetLoc(int);//

		//反時計回りに並び変え
		void sortCounterClockwise();

		//Getter
		int GetFlagRec()const;
		double GetCellTime()const;

		int GetFixFlag()const;
		int GetBoundaryFlag()const;
		int GetApicalFlag()const;
		int GetCrossBoundaryFlag()const;
		int GetColorFlag()const;
		int GetLayerFlag()const;

		Vector3<double> GetNormalVector()const;
		Vector3<double> GetLoc(int)const;
		Vector3<double> GetLocMapping()const;
		Vector3<double> GetForce(int)const;
		Vector3<double> GetFrcThread(int)const;

		Vector3<int> GetOffsetLoc()const;

		double GetGaussianCurvature()const;
		double GetMeanCurvature()const;

		double GetVoronoiArea()const;
		double GetVoronoiPerimeter()const;

		std::vector<int> GetEdge()const;
		int GetEdge(int)const;
		std::vector<int> GetFace()const;
		int GetFace(int)const;
		std::vector<int> GetNextVertex()const;
		int GetNextVertex(int)const;
		
		//setter
		void SetFlagRec();

		void SetLoc(const Vector3<double>& v,int n);
		void addLoc(const Vector3<double>& v,int n);
		
		void addForce(const Vector3<double>& v,int n);
		void resetForce(int);
		void addFrcThread(const Vector3<double>&,int n);
		void resetFrcThread(int);
		

		void SetExternalForce(const Vector3<double>&);

		void SetOffsetLoc(const Vector3<int>& v);//offsetを設定

		void SetFixFlag(int i);//fixflagの設定
		void SetBoundaryFlag(int i);//boundaryflagの設定
		void SetApicalFlag(int i);//apicalflagの設定
		void SetCrossBoundaryFlag(int i);
		void SetColorFlag(int i);
		void SetLayerFlag(int i);

		void SetCellTime(double);
		//値の設定
		void SetNextVertex(int,int);//i番目の値をjにする
		void SetEdge(int,int);
		void SetFace(int,int);
		
		//push_back用
		void PushEdge(int);
		void PushFace(int);
		void PushNextVertex(int);

		void SortAndErase();//vi_nextの重複要素を削除する用

		void SortEdge();//小さい順に並び変える
		void SortFace();

		//要素を消すとき用(同じものを探して消す)
		void findAndErase_ei(int);
    	void findAndErase_fi(int);
		void findAndErase_vi_next(int);
		
		//要素の変更用(i1をi2に置き換える)
		void Replace_fi(int,int);
		void Replace_ei(int,int);
		void Replace_vi_next(int,int);

		void Decrement_ei(int);
		void Decrement_fi(int);

};
//辺
class Edge{
	private:
		//対応関係
		int vi[2];//頂点のインデックス
		std::vector<int> fi;//面のインデックス

		double Length;//長さ
		double DihedralAngle;//二面角

		double BaseLength;//平衡辺長
		double BaseDihedralAngle;//平衡二面角

		double k1_length;//係数
		double k2_length;
		double k_theta;

		int flipFlag;//直前でフリップしたかどうか
		int apical_flag;//LineApicalConstriction用

		//フリップの一時保持
		int fi1;//vi[1],vi[0]の順に回る
		int fi2;//vi[0],vi[1]の順に回る
		int vi_fi1;//辺に関係のない頂点のインデックス(f1側)
		int vi_fi2;//辺に関係のない頂点のインデックス(f2側)
		/*----------------------*/
	public:
		Shape* p_s;//辺の属するshape

		Edge(int v0,int v1,Shape* _p_s);

		void CalcLength(int);//長さの計算
		void CalcDihedralAngle(int);//二面角の計算

		bool FlipJudgement();//判定のみ
		void FlipOperation();//判定つき
		void FlipOperation_forge();//強制的にフリップ
		
		//定数のsetter
		void SetBaseLength(const double a);
		void SetBaseDihedralAngle(const double a);
		//void ResetBaseDihedralAngle();//BaseDihedralAngleをリセット
		//void SetApicalBaseDihedralAngle();//BaseDihedralAngleを変更

		void SetK1_Length(const double a);
		void SetK2_Length(const double a);
		void SetK_Theta(const double a);

		//getter
		int GetVertex(const int i)const;
		std::vector<int>GetFace()const;
		int GetFace(const int i)const;

		int GetFlipFlag()const;
		void SetFlipFlag(int);

		int GetApicalFlag()const;
		void SetApicalFlag(int);

		//定数のgetter
		double GetLength()const;
		double GetDihedralAngle()const;
		double GetBaseLength()const;
		double GetBaseDihedralAngle()const;
		double Get_k1_length()const;
		double Get_k2_length()const;
		double Get_k_theta()const;


		void findAndErase_fi(int);//面を消す
		
		void SortFace();//小さい順
		void PushFace(int);//vectorに入れる用

	   	void Replace_vi(int,int);//int n1をint n2に置き換える
		void Replace_fi(int,int);

		void Decrement_fi(int);

};
//面
class Face{
	private:
		//対応関係
		std::vector<int> vi;//頂点のインデックス(反時計回り固定)
		std::vector<int> ei;//辺のインデックス

		int cross_boundary_flag;

		Vector3<double> center;//面の重心座標
		Vector3<double> NormalVector;//面の法線ベクトル
		double Area;//面積
		Vector3<double> circum_center;//外心

		double BaseArea;//平衡面積
		double k_area;//面積弾性比例係数
	public:
		Shape* p_s;//面の属するShape

		Face(int v0,int v1,int v2,Shape* _p_s);
		//Face(int v0,int v1,int v2);

		void CalcCenter(int);//重心の計算
		void CalcNormalVector(int);//法線ベクトルの計算
		void CalcArea(int);//面積の計算
		void CalcCircumCenter(int);//外心の計算

		//setter
		void SetCrossBoundaryFlag(int i);

		void SetBaseArea(const double a);
		void SetBaseArea();
		void SetK_Area(const double a);

		//getter
		int GetCrossBoundaryFlag()const;

		std::vector<int>GetVertex()const;
		int GetVertex(const int i)const;

		std::vector<int> GetEdge()const;
		int GetEdge(const int i)const;
		

		Vector3<double> GetCenter()const;
		Vector3<double> GetNormalVector()const;
		double GetArea()const;
		Vector3<double> GetCircumCenter()const;
		
		double GetBaseArea()const;
		double Get_k_area()const;

		//sort用
		void SortEdge();//小さい順(要素を消してインデックスがズレるとき用)

		//push_back用
		void PushVertex(int);
		void PushEdge(int);

		//入れ替え用
		void Replace_vi(int,int);
		void Replace_ei(int,int);

		void Decrement_ei(int);
		void Decrement_vi(int);
};

#endif
