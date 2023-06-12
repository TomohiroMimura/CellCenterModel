#ifndef RESTRUCTION_HPP
#define RESTRUCTION_HPP

#include<vector>
#include"vector3.hpp"
#include"class.hpp"
#include"Random.hpp"


namespace restructure{
    void cellDivision(Random*,Shape*,int step);//細胞分裂
	bool DivideVoronoi(Random*,Shape*,int num,int step);//分裂軸を選ぶ

	//std::vector<std::pair<int,Vector3<double>>>  CalcCrossPoint(Random*,Vertex*,Vector3<double>,std::vector<Vector3<double>>);//使わない
	bool isInnerVertex(Vertex*,Vector3<double>);
    bool isDihedralAngleCorrect(Vector3<double> n1,Vector3<double> n2);//2つの法線ベクトル間
    void OutputOFF(Shape* p_s,int step,int num,int type);//error毎の出力

    void EvaluateNormalVector(Shape*,Vertex*);//法線ベクトルの重みづけ評価
    void AngleRestraint(Shape* p_s,Face* p_f);//新しくできる三角形の最小角度の分類
    void EvaluateCentroid(Shape*,Vertex*,std::vector<Vector3<double>>);//重心に近いかどうかの評価
}

#endif