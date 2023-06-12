#ifndef FORCE_HPP
#define FORCE_HPP
#include<omp.h>
#include"vector3.hpp"
//#include"class.hpp"

class Shape;
class Edge;
class Face;
class Vertex;

namespace force{
    //intは計算精度を示すdegree(ただし使ってない)
    void CalcTotalForce(Shape*,int);

    //三角形領域
    void CalcLineForce(Shape*,int);
    void CalcAreaForce(Shape*,int);

    //曲げ剛性
    void CalcDihedralAngleForce(Shape*,int);
    void CalcCurvatureForce(Shape*,int);//Julicher
    void CalcCurvatureForce2(Shape*,int);//Laplace-Beltrami

    void CalcInteriorAngleForce(Shape*,int);
    void CalcZ_axisRestraintForce(Shape*,int);
    void CalcY_axisRestraintForce(Shape*,int);
    void CalcApicalConstrictionForce(Shape*,int);

    //遠距離相互作用成分(r<rcのときに斥力が働く)
    void CalcRepulsiveForce(Shape* p_s,int deg);//単層間のみ
    //2層間は引力と斥力両方あり、閾値も異なる(lateral方向の話に限る)
    void CalcLateralForce(Shape* p_s,int deg);

    //圧力
    void CalcPressure(Shape* p_s,int deg);
    

    //補正エネルギー
    void CalcVoronoiCenterForce(Shape*,int);//重心に近づける力

    //voronoi領域
    void CalcVoronoiAreaForce(Shape*,int);
    void CalcVoronoiPerimeterForce(Shape*,int);

    //揺らぎ
    void Fluctuation(Shape* p_s,int deg);//揺らぎ

    //リダクション処理
    void OMP_Reduction_Frc(Shape* p_s,int deg);//並列計算で求めたfrc_threadをforceにまとめる


//------------(単独でも使えるもの)--------------------//
    void CalcLineForce(Edge*,int);
    void CalcAreaForce(Face*,int);
    void CalcDihedralAngleForce(Edge*,int);


    //計算用
    Vector3<double> CalcDihedralGlad(Shape*,Vector3<double>,Vector3<double>,Vector3<double>,Vector3<double>,int);//0と1がe,2が左,3が右

    //確認用-E(x+dx)-E(x)
    void CalcCurvatureEnergyDifference(Shape*,int);//julicher
    double CalcCurvatureEnergyDifference(Shape*,int m,int n,double dh);//Julicher
}
#endif