#include<iostream>
using namespace std;
#include"arc.h"
#include"Point.h"
#include"Line.h"
#include"double_equals.h"
#include<algorithm>

//void PrintOut (vector<Point> p)
//{
//	
//}

int main() {
	//计算点绕某点旋转后的坐标
   //参数pt表示被绕的点；degree表示旋转的弧度，逆时针为正，反正为负
	/*Point p1(0, 0), p2(2, 2);
	Point p = p2.RotationPoint(-PI / 4, p1);
	p.PrintPoint();
*/

	/*--------------------------------------------------------------------------*/
	//返回直线某比例处点的坐标
	/*Point p1(-1,0), p2(2,3);
	Line L(p1, p2 );
	Point p = L.proportionPt(0.333 );
	p.PrintPoint();*/

	//计算并返回直线的长度
	/*Point p2(0, 0), p1(3, 3);
	Line L(p1, p2);
	L.PrintLine();*/

	//判断点pt是否在线段L上
	/*Point p1(0,0), p2(-3.5355339, -3.5355339), pt1(2.5,2.5), pt2(4, 4),pt3(1,-5);
	Line L(p1, p2);
	int ret = L.IsInLine(pt1);
	cout << ret << endl;*/
	//L.IsInLine(pt1);
	//L.IsInLine(pt2);
	//L.IsInLine(pt3);

	/*线段自我平移，参数传入一个点Dis（x,y)，表示平移的方向和距离。
	数值为正时表示向上（向右）平移，为负时相反*/
	/*Point p1(0, 0), p2(3, 3);
	Line L(p1, p2);
	L.Translation(3, 1);
	L.Translation(3, -1);
	L.PrintLine();*/

	//线段镜像后新的线段
	/*Point p1(-3, 3), p2(0, 0),pt1(0,0),pt2(0,2);
	Line L(p1, p2),l(pt1,pt2);
	Line* newLine = L.NewMirrorLine(l);
	newLine->PrintLine();*/

	//线段自我反向
	/*Point p1(1, 1), p2(3, 3);
	Line L(p1, p2);
	L.ReverseLine();
	L.PrintLine();*/

	//判断直线是否相交
	/*Point p1(-5, 0), p2(5, 0), p3(0, 0), p4(5,0),p5(-6,1),p6(7,-1);
	Line L1(p1, p2);
	Line L2(p3, p4);
	Line L3(p5, p6);
    bool ret1 = L1.IsOverlap(L2);
	bool ret2 = L1.IsOverlap(L3);
	cout << ret1 << endl;*/

	//线段求交点，并返回交点.
	//Point p1(4, 1), p2(4, 5), p3(8, 1), p4(3, 0), p5(2, 0), p6(3, 1);
	//Line L1(p1, p2);
	//Line L2(p2, p3);
	////Line L3(p6, p5);
	//bool pt1 = L1.IsOverlap(L2);
	//cout << pt1 << endl;
	
    
	/*Point pt2 = L1.LinePoint(L3);
	pt2.PrintPoint();
	pt2.PrintPoint();*/

	//线段自我延伸,参数a为延长长度,参数p判定延长点为那个，true表示起点，false表示终点
	/*Point p1(0, 0), p2(1.4, 1.4);
	Line L(p1, p2);
	L.Extent(false,1.4);
	L.PrintLine();*/

	/*Point p1(2, 2), p2(0, 0);
	Line L(p1, p2);
	L.Extent(1);
	L.PrintLine();*/

	/*Point p1(0, 0), p2(2, 2);
	Line L(p1, p2);
	L.Extent(-1);
	L.PrintLine();*/

	/*Point p1(2, 2), p2(0, 0);
	Line L(p1, p2);
	L.Extent(-1);
	L.PrintLine();*/

	//线段自我旋转degree(degree为弧度值）当弧度为负时表示顺时针，参数pt表示线段绕点pt旋转
	/*Point p1(0, 0), p2(0, 2),p3(1,1);
	Line L(p2, p3);
	L.RotationLine (PI/4, p1);
	L.PrintLine();*/

	//double e = asin(0.5);
	//cout << e << endl;

/*------------------------------------------------------------------------------------------*/
	//计算圆弧起点并返回起点坐标
	//double f = sqrt(12.5);
	/*Point pf(0, 0), pt1(f, f);
	Arc c(pf, -1.25*PI, PI , 5);*/

	//Point p1(0, 0);
	//Arc c(p1, 0.25*PI, -0.5*PI, 5);
	//Point p = c.Gain_Star();
	//p.PrintPoint();

	//计算圆弧起点并返回终点坐标
    /*Point pf(0, 0);
    Arc c(pf, -0.25*PI, -0.25*PI , 5);
    Point p = c.Gain_Fin();
    p.PrintPoint();*/

	//计算并返回圆弧的长度
	/*Point pf(0, 0), ps(5, 0), cen(0, 0);
	Arc c(pf, 0, PI/2, 5);
	double L = c.Length();
	cout << L << endl;*/

	//判断点pt是否在圆弧上,参数表示待判断的点；返回值返回1或-1.
	//点返回值为 - 1时点不在圆弧上，返回值为1时，点在圆弧上
   /* double f1 = sqrt(12.5),f2 = sqrt(8);
	cout << f2 << endl;*/
	
    /*Point pf(0, 0), ps(4,0), cen(4, 1),ps2(8,1),ps3(4,5);
	Arc c(pf, 0.5*PI, -1 * PI, 4);
	int index1 = c.IsInArc(ps);
	cout << index1 << endl;
	Arc c2(cen, 0.5*PI, -1 * PI, 4);
	int index2 = c2.IsInArc(ps2);
	cout << index2 << endl;
	int index3 = c2.IsInArc(ps3);
	cout << index3 << endl;*/
Point p1(0, 0), p2(3, 3), p3(0, 5), p4(1, 1), p5(0, -5), p6(3.5355339, 3.5355339), p7(3.5355339, 3.5355339);
Arc c2(p1, 0, PI, 5);
int ret = c2.IsInArc(p7);
cout << ret << endl;

	//int index = c.

	//圆弧自我反向
  /*  Point pf(0, 0);
	Arc c(pf, 0.75 * PI, PI / 4, 5);
	c.ReverseArc();
	c.cPrintArc();*/
    

	/*对圆弧进行镜像，参数l是一条直线，表示圆弧相对直线l镜像
	返回值是镜像后的圆弧指针*/
   /* Point cen(0, 0);
	Point p1(-1, 0), p2(0, -1);
	Line l(p1, p2);
	Arc c(cen, 0, PI / 2, 5);
	Arc *Newcircle1 = c.NewMirrorArc(l);
	Newcircle1->cPrintArc();*/


	//圆弧自我旋转degree(degree为弧度值）当弧度为负时表示顺时针，为正时表示逆时针；
	//参数pt表示线段绕点pt旋转
     /*Point cen(0, 0);
	 Arc c(cen, 0, PI/4, 5);
	 c.RotationArc(PI / 4, cen);
	 c.cPrintArc();
*/
	 //将圆弧某端点延伸一段距离
     //参数p为布尔值，如果p为true，延长点为起点；反之为终点；
     //参数a为延长距离
     //返回值为延伸后的圆弧指针
   /*  Point cen(0, 0);
	 Arc c(cen, 0, PI/4, 5);
	 c.Extent(false,3.925);
	 c.cPrintArc();*/

	/*double a = atan(1);
	double b = atan(-1);
	cout << a << b << endl;*/

	//圆弧求交点，并返回交点.
    /*Point cen1(0, 0), cen2(4, 0);
	Arc c1(cen1, 0.5*PI, -1* PI, 4), c2(cen2, 0.5*PI,1 * PI, 4);
	vector<Point> p = c1.ArcPoint(c2);
	for (vector<Point>::iterator it = p.begin(); it != p.end(); it++)
	{
		cout << "(" << (*it).m_x << "," << (*it).m_y << ")" << " ";
	}
	cout << endl;*/


//Point p1(0,0);
//vector<Point> p{ 2,p1 };
//for (vector<Point>::iterator it = p.begin(); it != p.end(); it++)
//{
//	cout << "(" << (*it).m_x << "," << (*it).m_y <<")"<< " ";
//}
//cout << endl; 

    /*计算圆弧某比例处的点坐标，比例值是从起点开始计算
	参数a表示比例，在[0—1]范围内
	返回值为此比例处点的坐标*/
  /* Point p1(0, 0);
   Arc c(p1, 0.25*PI,-0.5*PI, 5);
   Point pt = c.proportionPt(0.5);
   pt.PrintPoint();*/


   //*圆弧自我平移，参数传入一个点Dis（x,y)，表示平移的方向和距离。
   //数值为正时表示向上（向右）平移，为负时相反*/
  /* Point p1(0, 0),p2(-3,1);
   Arc c(p1, 0.25*PI, -0.5*PI, 5);
   cout << "旧圆弧" << endl;
   c.cPrintArc();
   cout << "新圆弧" << endl;
   c.Translation(p2);
   c.cPrintArc();*/
   system("pause");
	return 0;
}