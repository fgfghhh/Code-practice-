#pragma once
#include<iostream>
#include"Point.h"
using namespace std;
#include"Line.h"
#include<cmath>
//#include"double_equals.h"
#include<vector>
#include<algorithm>
//构建一个圆弧类，使用 圆弧的圆心，圆弧起始点与x轴正方向的夹角(逆时针为正），
//圆弧起点到终点扫过的角度(逆时针为正），圆弧半径R的方式画圆弧。实现一些圆弧的基本操作
class Arc
{
public:
	Arc();

	//构造函数，初始化圆弧的参数，a赋值给起点与x轴正方向夹角，
	//b赋值给起点到终点扫过的角度，r赋值给圆的半径
	Arc(Point pt1, double a, double b, double r);
	
	//计算圆弧起点并返回起点坐标
	Point Gain_Star();

	//计算圆弧终点并返回终点坐标
	Point Gain_Fin();

	/*计算圆弧某比例处的点坐标，比例值是从起点开始计算
	参数a表示比例，在[0—1]范围内
	返回值为此比例处点的坐标*/
	Point proportionPt(double a);

	//计算并返回圆弧的长度
	double Length(); 

	//输出圆弧的两个端点和圆心及圆弧长度长度测试用
	void cPrintArc(); 

	//克隆圆弧,返回值是新圆弧的指针
	Arc* CopyFun();

	//判断点pt是否在圆弧上,参数表示待判断的点；返回值返回1或-1.
    //点返回值为 - 1时表示点不在圆弧上，返回值为1时，表示点在圆弧上
	int IsInArc(Point pt);

	/*圆弧自我平移，参数传入一个点Dis（x,y)，表示平移的方向和距离。
	数值为正时表示向上（向右）平移，为负时相反*/
	void Translation(Point Dis);

	/*圆弧平移生成新圆弧，参数传入一个点Dis（x,y)，表示平移的方向和距离。
    数值为正时表示向上（向右）平移，为负时相反；返回值为新圆弧的指针*/
	Arc* NewArc(Point Dis);


	/*对圆弧进行镜像，参数l是一条直线，表示圆弧相对直线l镜像
    返回值是镜像后的圆弧指针*/
	Arc *NewMirrorArc(Line l);

	//圆弧自我反向
	void ReverseArc();

	//圆弧自我旋转degree(degree为弧度值）当弧度为负时表示顺时针，为正时表示逆时针；
    //参数pt表示线段绕点pt旋转
	void RotationArc(double degree, Point pt);

	//计算两个圆弧的交点，参数为另一条圆弧；
	//返回值为两个坐标，如果不存在交点，返回空值；
	//如果存在一个交点，返回一个交点坐标；
	//如果存在两个交点，返回两个交点的坐标值
	vector<Point> ArcPoint(Arc c1);

	//将圆弧某端点延伸一段距离
    //参数p为布尔值，如果p为true，延长点为起点；反之为终点；
    //参数a为延长距离
    //返回值为延伸后的圆弧指针
	void Extent(bool p, double a);

	// 判断两个点坐标相同
	int double_equal_c(double a, double b);

	~Arc();

public:
	Point m_center;//圆弧的圆心
	double m_Stangle;//圆弧起始点与x轴正方向的夹角(逆时针为正）
	double m_Swangle;//圆弧起点到终点扫过的角度(逆时针为正）
	double m_R;//圆弧半径R
};