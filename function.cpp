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
	//�������ĳ����ת�������
   //����pt��ʾ���Ƶĵ㣻degree��ʾ��ת�Ļ��ȣ���ʱ��Ϊ��������Ϊ��
	/*Point p1(0, 0), p2(2, 2);
	Point p = p2.RotationPoint(-PI / 4, p1);
	p.PrintPoint();
*/

	/*--------------------------------------------------------------------------*/
	//����ֱ��ĳ�������������
	/*Point p1(-1,0), p2(2,3);
	Line L(p1, p2 );
	Point p = L.proportionPt(0.333 );
	p.PrintPoint();*/

	//���㲢����ֱ�ߵĳ���
	/*Point p2(0, 0), p1(3, 3);
	Line L(p1, p2);
	L.PrintLine();*/

	//�жϵ�pt�Ƿ����߶�L��
	/*Point p1(0,0), p2(-3.5355339, -3.5355339), pt1(2.5,2.5), pt2(4, 4),pt3(1,-5);
	Line L(p1, p2);
	int ret = L.IsInLine(pt1);
	cout << ret << endl;*/
	//L.IsInLine(pt1);
	//L.IsInLine(pt2);
	//L.IsInLine(pt3);

	/*�߶�����ƽ�ƣ���������һ����Dis��x,y)����ʾƽ�Ƶķ���;��롣
	��ֵΪ��ʱ��ʾ���ϣ����ң�ƽ�ƣ�Ϊ��ʱ�෴*/
	/*Point p1(0, 0), p2(3, 3);
	Line L(p1, p2);
	L.Translation(3, 1);
	L.Translation(3, -1);
	L.PrintLine();*/

	//�߶ξ�����µ��߶�
	/*Point p1(-3, 3), p2(0, 0),pt1(0,0),pt2(0,2);
	Line L(p1, p2),l(pt1,pt2);
	Line* newLine = L.NewMirrorLine(l);
	newLine->PrintLine();*/

	//�߶����ҷ���
	/*Point p1(1, 1), p2(3, 3);
	Line L(p1, p2);
	L.ReverseLine();
	L.PrintLine();*/

	//�ж�ֱ���Ƿ��ཻ
	/*Point p1(-5, 0), p2(5, 0), p3(0, 0), p4(5,0),p5(-6,1),p6(7,-1);
	Line L1(p1, p2);
	Line L2(p3, p4);
	Line L3(p5, p6);
    bool ret1 = L1.IsOverlap(L2);
	bool ret2 = L1.IsOverlap(L3);
	cout << ret1 << endl;*/

	//�߶��󽻵㣬�����ؽ���.
	//Point p1(4, 1), p2(4, 5), p3(8, 1), p4(3, 0), p5(2, 0), p6(3, 1);
	//Line L1(p1, p2);
	//Line L2(p2, p3);
	////Line L3(p6, p5);
	//bool pt1 = L1.IsOverlap(L2);
	//cout << pt1 << endl;
	
    
	/*Point pt2 = L1.LinePoint(L3);
	pt2.PrintPoint();
	pt2.PrintPoint();*/

	//�߶���������,����aΪ�ӳ�����,����p�ж��ӳ���Ϊ�Ǹ���true��ʾ��㣬false��ʾ�յ�
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

	//�߶�������תdegree(degreeΪ����ֵ��������Ϊ��ʱ��ʾ˳ʱ�룬����pt��ʾ�߶��Ƶ�pt��ת
	/*Point p1(0, 0), p2(0, 2),p3(1,1);
	Line L(p2, p3);
	L.RotationLine (PI/4, p1);
	L.PrintLine();*/

	//double e = asin(0.5);
	//cout << e << endl;

/*------------------------------------------------------------------------------------------*/
	//����Բ����㲢�����������
	//double f = sqrt(12.5);
	/*Point pf(0, 0), pt1(f, f);
	Arc c(pf, -1.25*PI, PI , 5);*/

	//Point p1(0, 0);
	//Arc c(p1, 0.25*PI, -0.5*PI, 5);
	//Point p = c.Gain_Star();
	//p.PrintPoint();

	//����Բ����㲢�����յ�����
    /*Point pf(0, 0);
    Arc c(pf, -0.25*PI, -0.25*PI , 5);
    Point p = c.Gain_Fin();
    p.PrintPoint();*/

	//���㲢����Բ���ĳ���
	/*Point pf(0, 0), ps(5, 0), cen(0, 0);
	Arc c(pf, 0, PI/2, 5);
	double L = c.Length();
	cout << L << endl;*/

	//�жϵ�pt�Ƿ���Բ����,������ʾ���жϵĵ㣻����ֵ����1��-1.
	//�㷵��ֵΪ - 1ʱ�㲻��Բ���ϣ�����ֵΪ1ʱ������Բ����
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

	//Բ�����ҷ���
  /*  Point pf(0, 0);
	Arc c(pf, 0.75 * PI, PI / 4, 5);
	c.ReverseArc();
	c.cPrintArc();*/
    

	/*��Բ�����о��񣬲���l��һ��ֱ�ߣ���ʾԲ�����ֱ��l����
	����ֵ�Ǿ�����Բ��ָ��*/
   /* Point cen(0, 0);
	Point p1(-1, 0), p2(0, -1);
	Line l(p1, p2);
	Arc c(cen, 0, PI / 2, 5);
	Arc *Newcircle1 = c.NewMirrorArc(l);
	Newcircle1->cPrintArc();*/


	//Բ��������תdegree(degreeΪ����ֵ��������Ϊ��ʱ��ʾ˳ʱ�룬Ϊ��ʱ��ʾ��ʱ�룻
	//����pt��ʾ�߶��Ƶ�pt��ת
     /*Point cen(0, 0);
	 Arc c(cen, 0, PI/4, 5);
	 c.RotationArc(PI / 4, cen);
	 c.cPrintArc();
*/
	 //��Բ��ĳ�˵�����һ�ξ���
     //����pΪ����ֵ�����pΪtrue���ӳ���Ϊ��㣻��֮Ϊ�յ㣻
     //����aΪ�ӳ�����
     //����ֵΪ������Բ��ָ��
   /*  Point cen(0, 0);
	 Arc c(cen, 0, PI/4, 5);
	 c.Extent(false,3.925);
	 c.cPrintArc();*/

	/*double a = atan(1);
	double b = atan(-1);
	cout << a << b << endl;*/

	//Բ���󽻵㣬�����ؽ���.
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

    /*����Բ��ĳ�������ĵ����꣬����ֵ�Ǵ���㿪ʼ����
	����a��ʾ��������[0��1]��Χ��
	����ֵΪ�˱������������*/
  /* Point p1(0, 0);
   Arc c(p1, 0.25*PI,-0.5*PI, 5);
   Point pt = c.proportionPt(0.5);
   pt.PrintPoint();*/


   //*Բ������ƽ�ƣ���������һ����Dis��x,y)����ʾƽ�Ƶķ���;��롣
   //��ֵΪ��ʱ��ʾ���ϣ����ң�ƽ�ƣ�Ϊ��ʱ�෴*/
  /* Point p1(0, 0),p2(-3,1);
   Arc c(p1, 0.25*PI, -0.5*PI, 5);
   cout << "��Բ��" << endl;
   c.cPrintArc();
   cout << "��Բ��" << endl;
   c.Translation(p2);
   c.cPrintArc();*/
   system("pause");
	return 0;
}