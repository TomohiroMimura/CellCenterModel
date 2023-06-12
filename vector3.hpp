#ifndef VECTOR3_HPP
#define VECTOR3_HPP
#include<iostream>
#include<cmath>

template<class T>class Vector3{
	private:
		T x;
		T y;
		T z;
	public:
		//コンストラクタ
		Vector3();
		Vector3(const T& a,const T& b,const T&c);
		Vector3(const Vector3<T>& v);
		//セッター(setter)
		void SetX(const T& a);
		void SetY(const T& a);
		void SetZ(const T& a);
		void Set(const T& a,const T& b,const T& c);
		void Set(const Vector3<T>& v);
		//ゲッター(getter)
		T GetX() const;
		T GetY() const;
		T GetZ() const;
		
		T CalcNorm() const;//ノルム計算
		T CalcSqr() const;//二乗和（ノルムの二乗)

		Vector3<T> Projection(const Vector3<T>& n)const;//法線ベクトルnで表される平面上にベクトルを投影
		//ベクトルの回転
		Vector3<T> Rotation(const Vector3<T>& n,const T& theta)const;//nベクトルを回転軸としてthetaだけ回転したベクトルを出力(radian)
		//オペレータのオーバーロード
		Vector3<T>& operator = (const Vector3<T>& v);//代入
		Vector3<T>& operator +=(const Vector3<T>& v);//足して代入
		Vector3<T>& operator -=(const Vector3<T>& v);//引いて代入
		Vector3<T>& operator *=(const T& a);//スカラー倍して代入
		Vector3<T>& operator /=(const T& a);//スカラーで割って代入
		Vector3<T> operator -()const;//-1倍
		bool operator==(const Vector3<T>& v)const;//比較演算子
		bool operator!=(const Vector3<T>& v)const;

		//外部定義オペレーター(エラーになるので削除)
	//	friend Vector3<T> operator*(const T& a,const Vector3<T>& v );//スカラー倍(左から)
	//	friend Vector3<T> operator*(const Vector3<T>& v,const T& a);//スカラー倍(右から)
	//	friend Vector3<T> operator/(const Vector3<T>& v,const T& a);//スカラーで割る
	//	friend Vector3<T> operator+(const Vector3<T>& v1,const Vector3<T>&v2);//ベクトル同士の足し算
	//	friend Vector3<T> operator-(const Vector3<T>& v1,const Vector3<T>&v2);//ベクトル同士の引き算
	//	friend T operator*(const Vector3<T>& v1,const Vector3<T>& v2);//スカラー積(内積)
	//	friend Vector3<T> operator%(const Vector3<T>& v1,const Vector3<T>&v2);//ベクトル積(外積)
};

/*
********************************************************************   
   上の関数の実装

*/
//デフォルトコンストラクタ
template<class T> Vector3<T>::Vector3(){
	x=0.0;y=0.0;z=0.0;
}
template<class T>Vector3<T>::Vector3(const T& a,const T& b,const T& c){
	x=a;y=b;z=c;
}
//同一クラスであれば他のインスタンスからでもアクセスできる
template<class T>Vector3<T>::Vector3(const Vector3<T>& v){
	x=v.x;y=v.y;z=v.z;
}
//セッター
template<class T>void Vector3<T>::SetX(const T& a){
	x=a;
}
template<class T>void Vector3<T>::SetY(const T& a){
	y=a;
}
template<class T>void Vector3<T>::SetZ(const T& a){
	z=a;
}
template<class T>void Vector3<T>::Set(const T& a,const T& b,const T& c){
	x=a;y=b;z=c;
}
template<class T>void Vector3<T>::Set(const Vector3<T>& v){
	x=v.x;
	y=v.y;
	z=v.z;
}
//ゲッター
template<class T> T Vector3<T>::GetX()const{
	return x;
}
template<class T> T Vector3<T>::GetY()const{
	return y;
}
template<class T> T Vector3<T>::GetZ()const{
	return z;
}
//
template<class T> T Vector3<T>::CalcNorm()const{
	return sqrt(x*x+y*y+z*z);
}
template<class T> T Vector3<T>::CalcSqr()const{
	return x*x+y*y+z*z;
}
template<class T> Vector3<T> Vector3<T>::Projection(const Vector3<T>& n)const{
	if(n.CalcNorm()<1e-9){
		std::cout<<"ProjectionError:vector n is too small."<<std::endl;
	}
	Vector3<T>n1=n/n.CalcNorm();
	Vector3<T> v;
	Vector3<T> a(x,y,z);
	v=a-(a*n1)*n1;
	return v;
}
//回転
template<class T> Vector3<T> Vector3<T>::Rotation(const Vector3<T>& n,const T& theta)const{
	if(n.CalcNorm()<1e-9){
		std::cout<<"RotationError:vector n is too small."<<std::endl;
	}
	Vector3<T> v(GetX(),GetY(),GetZ());
	Vector3<T> n1=n/n.CalcNorm();
	Vector3<T> tmp_vec1(v*cos(theta));
	Vector3<T> tmp_vec2((1.0-cos(theta))*(v*n1)*n1);
	Vector3<T> tmp_vec3((n1%v)*sin(theta));
	return tmp_vec1+tmp_vec2+tmp_vec3;
}

//オペレータ
//代入
template<class T> Vector3<T>& Vector3<T>::operator=(const Vector3<T>& v){
	x=v.x;y=v.y;z=v.z;
	return *this;
}
//足して代入
template<class T> Vector3<T>& Vector3<T>::operator+=(const Vector3<T>& v){
	x+=v.x;y+=v.y;z+=v.z;
	return *this;
}
//ひいて代入
template<class T> Vector3<T>& Vector3<T>::operator-=(const Vector3<T>& v){
	x-=v.x;y-=v.y;z-=v.z;
	return *this;
}
//かけて代入
template<class T> Vector3<T>& Vector3<T>::operator*=(const T& a){
	x*=a;y*=a;z*=a;
	return *this;
}
//割って代入
template<class T> Vector3<T>& Vector3<T>::operator/=(const T& a){
	x/=a;y/=a;z/=a;
	return *this;
}
//-1倍(加法の逆元)
template<class T>Vector3<T> Vector3<T>::operator-()const{
	return Vector3<T>(-x,-y,-z);
}
//比較演算
template<class T> bool Vector3<T>::operator==(const Vector3<T>& v)const{
	if((x==v.x)&&(y==v.y)&&(z==v.z)){
		return true;
	}
	return false;
}
template<class T>bool Vector3<T>::operator!=(const Vector3<T>& v)const{
	if((x==v.x)&&(y==v.y)&&(z==v.z)){
		return false;
	}
	return true;
}

//外部定義オペレータ(classのメンバ関数ではない)
//スカラー倍(左から)
template<class T>Vector3<T> operator*(const T& a,const Vector3<T>& v){
	return Vector3<T>(a*v.GetX(),a*v.GetY(),a*v.GetZ());
}
//スカラー倍(右から)
template<class T>Vector3<T> operator*(const Vector3<T>& v,const T& a){
	return Vector3<T>(a*v.GetX(),a*v.GetY(),a*v.GetZ());
}
//ベクトルをスカラーで割る
template<class T>Vector3<T> operator/(const Vector3<T>& v,const T& a){
	return Vector3<T>(v.GetX()/a,v.GetY()/a,v.GetZ()/a);
}
//ベクトル同士の足し算
template<class T>Vector3<T> operator+(const Vector3<T>& v1,const Vector3<T>& v2){
	return Vector3<T>(v1.GetX()+v2.GetX(),v1.GetY()+v2.GetY(),v1.GetZ()+v2.GetZ());
}
//ベクトル同士の引き算
template<class T>Vector3<T> operator-(const Vector3<T>& v1,const Vector3<T>& v2){
	return Vector3<T>(v1.GetX()-v2.GetX(),v1.GetY()-v2.GetY(),v1.GetZ()-v2.GetZ());
}
//スカラー積(内積)
template<class T>T operator*(const Vector3<T>& v1,const Vector3<T>& v2){
	return v1.GetX()*v2.GetX()+v1.GetY()*v2.GetY()+v1.GetZ()*v2.GetZ();
}
//ベクトル積(外積)
template<class T>Vector3<T> operator%(const Vector3<T>& v1,const Vector3<T>& v2){
	return Vector3<T>((v1.GetY()*v2.GetZ()-v1.GetZ()*v2.GetY()),(v1.GetZ()*v2.GetX()-v1.GetX()*v2.GetZ()),(v1.GetX()*v2.GetY()-v1.GetY()*v2.GetX()));
}
//出力用operator
template<class T>std::ostream& operator<<(std::ostream& output,const Vector3<T>& v){
	output<<v.GetX()<<" "<<v.GetY()<<" "<<v.GetZ();
	return output;
}


#endif //VECTOR3_HPP
