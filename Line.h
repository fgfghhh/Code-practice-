#pragma once
#include<iostream>
#include"Point.h"
using namespace std;
#include<algorithm>
#define _MATH_DEFINES_DEFINED
#include<math.h>
#define PI  3.141592654
//#include"double_equals.h"
//const double EPSINON = 0.00001;

//构造一个直线类，并实现一些对直线的简单操作
class Line
{
public:
	//默认构造函数
	Line();

	//构造函数，用初始化直线的两个端.pt1赋值给起点，pt2赋值给终点
	Line(Point pt1, Point pt2);

	/*计算直线某比例处的点坐标，比例值是从起点开始计算
    参数a表示比例，在[0—1]范围内
	返回值为此比例处点的坐标*/
	Point proportionPt(double a);

	//计算并返回直线的长度
	double Length();

	//输出直线的两个端点和直线长度，做测试用
	void PrintLine();

	//克隆一条新直线，返回值为新的直线
	Line* CopyFun();

	/*判断点是否在线段L上,参数pt表示要确定的点
	返回值为index，当值为-1时，表示点不在线上，当值为1时，表示点在线上*/
	int IsInLine(Point pt)const;

	/*线段自我平移，参数传入一个点Dis（x,y)，表示平移的方向和距离。
	数值为正时表示向上（向右）平移，为负时相反*/
	void Translation(Point Dis);

	/*线段平移后生成一条新线，参数传入一个点Dis（x,y)，表示平移的方向和距离。
    数值为正时表示向上（向右）平移，为负时相反
    返回值是新生成的线*/
	Line* NewLine(Point Dis);

	/*对线段进行镜像，参数l是一条直线，表示线段相对直线l镜像
	返回值是镜像后的直线*/
	Line* NewMirrorLine(Line l);

	//线段自我反向
	void ReverseLine();

	//线段自我旋转degree(degree为弧度值）当弧度为负时表示顺时针，为正时表示逆时针；
	//参数pt表示线段绕点pt旋转
	void RotationLine(double degree, Point pt);

	//判断直线是否相交，参数线段l为另一条直线
	//返回值为bool，如果为false，表示不相，为true表示相交
	bool IsOverlap(Line l);

	//线段求交点，并返回交点坐标。
	//参数是另一条直线，和点pt
	//如果存在交点，pt记录交点的坐标；如果不存在交点，则返回false
	//Point LinePoint(Line l1);
	bool LinePoint(Line l, Point& pt);
	//Point pt;//用来存放两直线的交点,
	
	//线段自我延伸;
	//参数p用来判定延长点为哪个,true表示直线的起点，false表示直线的终点
	//参数a表示延长距离
	void Extent(bool p,double a);

	//线段延伸生成新线;
	//参数p用来判定延长点为哪个,true表示直线的起点，false表示直线的终点
	//参数a表示延长距离
	//返回值为新线
	Line NewExtent(bool p, double a);



	//析构函数
	~Line() {};

public:
	Point m_ptf;//直线起点
	Point m_pts;//直线终点
};
