#ifndef PARAMETERS_HPP
#define PARAMETERS_HPP

constexpr char OFF_FILE_NAME[] = "InitSquare_periodic40.off";//入力offファイルの名前
constexpr double DELTA_TIME = 2e-4;//時間間隔
constexpr int STEP_END=500000;//繰り返し回数
constexpr int OUTPUT_VTK_STEP=5000;//何回おきにvtkファイルを出力するか
constexpr int CELL_DIVISION_STEP_AVE=500000;//何回おきに分裂するか(平均)
constexpr int CELL_DIVISION_STEP_SD=5000;//細胞周期の標準偏差

constexpr double alpha=1.0;
constexpr int TurnOffStep=0;

//細胞周期
constexpr int CELL_TIME_INITIALIZE_MODE=0;//0:ランダム,1:パラメータファイルから読み込み
constexpr char CELL_TIME_FILE_NAME[]="";
constexpr int Cell_Time_Output_Step=CELL_DIVISION_STEP_AVE/4;


// Trial code for future development of mathematical models for modeling a spatial pattern of the differential growth.
// does not affect monolayer simulations.
/*------------------------------------------------------------------*/
constexpr int DivisionDistributionMode=0;//1:空間分布あり,0:空間分布なし
//細胞周期に空間分布を与えるとき用のパラメータ
constexpr double DIV_alpha=0.1;
constexpr double DIV_beta=0.5;
/*------------------------------------------------------------------*/

constexpr int DEGREE_ACCURACY=1;//計算精度(1次:1,2次:2,3次:3,4次:4)
constexpr int NUM_THREAD=16;//並列計算の使用スレッド数

constexpr double VISCOSITY=0.25;//粘性係数
constexpr double BOUNDARY_VISCOSITY=0.25;//境界での粘性係数

constexpr double K_T_AREA=10;//面積弾性係数
constexpr double K_T_INTERIOR_ANGLE=0;//内角保存(triangle-internal angle)
constexpr double K1_T_LENGTH=0;//1次に比例する辺の弾性係数
constexpr double K2_T_LENGTH=2;//2次に比例する辺の弾性係数

constexpr double K_T_THETA=10;//二面角弾性係数
constexpr double fluctSD=0.02;//平衡二面角の標準偏差

// Trial code for future development of mathematical models for examining bending energies.
// does not affect monolayer simulations.
/*------------------------------------------------------------------*/
constexpr double Bending_Elastic_Constant=0;//曲げ弾性率
constexpr double fluct_curvature=0;//平衡平均曲率への摂動
constexpr int CurvatureType=1;//1:julicher,2:laplace-beltrami
/*------------------------------------------------------------------*/

constexpr int voronoi_area_type=1;//1:三角形の面積の3分の1の射影で重みづけ,2:voronoi_areaで重みづけ
constexpr double K_Z=0.1;//z方向の変形を抑制(-/+)0.0000625,0.000125,0.00025,0.001,0.002,0.004,0.008,0.016
constexpr double KzType=0;//0:両側,1:z>0のみ束縛,2:z<0のみ束縛
constexpr double K_Y=0;//Y方向に圧縮力を加える

constexpr double pressure=0;//圧力項(各点の法線方向)

constexpr double K_APICAL_DIS_S=5;//apical収縮の収縮力(surface)
constexpr double K_APICAL_DIS_L=5;//apical収縮の収縮力(line)
constexpr double M_APICAL_THETA_S=1;//apical収縮の際の二面角エネルギーの比例係数の重みづけ(surface)
constexpr double M_APICAL_THETA_L=1;//apical収縮の際の二面角エネルギーの比例係数の重みづけ(line)

constexpr double ApicalRadius=5.0;//半径(>1)

//こっちは円周指定のときだけ
constexpr double ApicalRadius_Max=5.0;//半径(>1)
constexpr double ApicalRadius_Min=0.0;//半径(>1)

constexpr int SurfaceApicalConstriction_circumferenceMode=0;//1:あり
constexpr int SurfaceApicalConstriction_circleMode=0;//1:あり
constexpr int LineApicalConstrictionMode=0;//1:あり
constexpr int SurfaceBeltApicalConstrictionMode=0;//1:あり
constexpr int LineBeltApicalConstrictionMode=0;//1:あり

constexpr int CellDivisionRegionMode=1;//1:全体で分裂,2:外側で分裂,3:内側で分裂,0:全体でも分裂しない

//排除体積効果
constexpr double K_REP=10;//排除効果
constexpr double RC=0.9;//排除半径

/*----------------------------------------------------------------------*/
/*
// Trial code for future development of mathematical models for multilayering;
// does not affect monolayer simulations.
//2層間の張力
*/
constexpr int Layer_NUM=1;//層の数
constexpr int LayerDivFlag[]={0,1};//1のときにFlagRecを設定(分裂しない設定)
constexpr double K_lateral=10;//比例係数
constexpr double R_LOW=1.5;
constexpr double R_EQ=2.0;
constexpr double R_UP=2.5;
//二層間の点においてR_LOW<=r<=R_UPにあるときにrをR_EQに近づける効果をもつ
/*----------------------------------------------------------------------*/

constexpr double K_V_CENTER=10;//重心に寄せるエネルギー(面内について)
constexpr double K_V_OUT_CENTER=0;//重心に寄せるエネルギー(法線方向について)

constexpr double K_V_AREA=0;//細胞の面積保存(-)
constexpr double K_V_PERIMETER=0;//細胞の周長保存(-)

constexpr double Fluct=0;//(廃止)法線方向の揺らぎ(力)の大きさ(rand(-Fluct,Fluct))

constexpr double BASE_LENGTH=1;//平衡長さ
constexpr double BASE_AREA=0.433;//平衡面積

constexpr double BaseVoronoiArea=2*BASE_AREA;
constexpr double BaseVoronoiPerimeter=2*sqrt(3)*BASE_LENGTH;

constexpr double BASE_DIHEDRAL_ANGLE=3.14159265;//平衡二面角
constexpr double APICAL_CONSTRICTION_ANGLE_S=0.26;//apical収縮の際に変化する二面角量(surface型)
constexpr double APICAL_CONSTRICTION_ANGLE_L=0.26;//apical収縮の際に変化する二面角量(line型)

constexpr double FLIP_THRESHOLD_LENGTH=0.25;//flip判定の閾値(基本的には0.5773)
constexpr int FLIP_COOLTIME_STEP=5;//一度flipした後のcool_timeの長さ

constexpr int PERIODIC_MODE=1;//periodic:1,non-periodic:0,Y-axis periodic:2
//周期領域の計算領域(xy平面)
constexpr double size_x=40;
constexpr double size_y=40*sqrt(3)/2.0;

constexpr int BOUNDARY_FIX_MODE=1;//1:境界点を固定,0:境界点を固定しない(when non-periodic mode)
constexpr int CELL_DIVISION_MODE=3;//0:なし,1:ランダム方向分裂軸,2:最長方向分裂軸,3:x方向分裂軸,4:円周方向分裂軸
constexpr int FLIP_OPERATION_MODE=1;//1でflipする,0でflipしない
constexpr int CELL_TIME_MODE=2;//細胞分裂する細胞の選び方
constexpr int FlagRecMode=0;//0のとき、一度分裂した細胞も分裂する,1のとき一度分裂した細胞は分裂しない

constexpr int OUTPUT_DivisionErrorOFFfile=0;//0:出力しない,1:出力する

#endif
