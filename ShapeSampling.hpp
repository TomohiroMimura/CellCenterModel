#ifndef SHAPESAMPLINGHPP
#define SHAPESAMPLINGHPP

#include<vector>
#include<iostream>
#include<cmath>
#include<string>
#include"parameters.hpp"
#include"class.hpp"
#include"Random.hpp"
#include"restructure.hpp"
class ShapeSampling{
        Shape* p_s;
        Random* p_rnd;
    public:
        ShapeSampling(Random*);
        ~ShapeSampling();
        
        void CalcMDstep();//計算
        //void Fluctuation();//面積が小さい三角形の頂点に摂動を与える
        
        void SurfaceApicalConstriction_circumference();//細胞円状領域(円周)
        void SurfaceApicalConstriction_circle();//細胞円状領域(全体)
        void LineApicalConstriction();//細胞円状領域(円周)
        void SurfaceBeltApicalConstriction();//細胞ベルト領域(細胞)
        void LineBeltApicalConstriction();//細胞ベルト領域(線)

        void Lloyd_Algorithm();
        Shape* GetShape()const;
        Random* GetRandom()const;
 
};
#endif