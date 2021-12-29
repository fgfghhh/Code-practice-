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

//����һ��ֱ���࣬��ʵ��һЩ��ֱ�ߵļ򵥲���
class Line
{
public:
	//Ĭ�Ϲ��캯��
	Line();

	//���캯�����ó�ʼ��ֱ�ߵ�������.pt1��ֵ����㣬pt2��ֵ���յ�
	Line(Point pt1, Point pt2);

	/*����ֱ��ĳ�������ĵ����꣬����ֵ�Ǵ���㿪ʼ����
    ����a��ʾ��������[0��1]��Χ��
	����ֵΪ�˱������������*/
	Point proportionPt(double a);

	//���㲢����ֱ�ߵĳ���
	double Length();

	//���ֱ�ߵ������˵��ֱ�߳��ȣ���������
	void PrintLine();

	//��¡һ����ֱ�ߣ�����ֵΪ�µ�ֱ��
	Line* CopyFun();

	/*�жϵ��Ƿ����߶�L��,����pt��ʾҪȷ���ĵ�
	����ֵΪindex����ֵΪ-1ʱ����ʾ�㲻�����ϣ���ֵΪ1ʱ����ʾ��������*/
	int IsInLine(Point pt)const;

	/*�߶�����ƽ�ƣ���������һ����Dis��x,y)����ʾƽ�Ƶķ���;��롣
	��ֵΪ��ʱ��ʾ���ϣ����ң�ƽ�ƣ�Ϊ��ʱ�෴*/
	void Translation(Point Dis);

	/*�߶�ƽ�ƺ�����һ�����ߣ���������һ����Dis��x,y)����ʾƽ�Ƶķ���;��롣
    ��ֵΪ��ʱ��ʾ���ϣ����ң�ƽ�ƣ�Ϊ��ʱ�෴
    ����ֵ�������ɵ���*/
	Line* NewLine(Point Dis);

	/*���߶ν��о��񣬲���l��һ��ֱ�ߣ���ʾ�߶����ֱ��l����
	����ֵ�Ǿ�����ֱ��*/
	Line* NewMirrorLine(Line l);

	//�߶����ҷ���
	void ReverseLine();

	//�߶�������תdegree(degreeΪ����ֵ��������Ϊ��ʱ��ʾ˳ʱ�룬Ϊ��ʱ��ʾ��ʱ�룻
	//����pt��ʾ�߶��Ƶ�pt��ת
	void RotationLine(double degree, Point pt);

	//�ж�ֱ���Ƿ��ཻ�������߶�lΪ��һ��ֱ��
	//����ֵΪbool�����Ϊfalse����ʾ���࣬Ϊtrue��ʾ�ཻ
	bool IsOverlap(Line l);

	//�߶��󽻵㣬�����ؽ������ꡣ
	//��������һ��ֱ�ߣ��͵�pt
	//������ڽ��㣬pt��¼��������ꣻ��������ڽ��㣬�򷵻�false
	//Point LinePoint(Line l1);
	bool LinePoint(Line l, Point& pt);
	//Point pt;//���������ֱ�ߵĽ���,
	
	//�߶���������;
	//����p�����ж��ӳ���Ϊ�ĸ�,true��ʾֱ�ߵ���㣬false��ʾֱ�ߵ��յ�
	//����a��ʾ�ӳ�����
	void Extent(bool p,double a);

	//�߶�������������;
	//����p�����ж��ӳ���Ϊ�ĸ�,true��ʾֱ�ߵ���㣬false��ʾֱ�ߵ��յ�
	//����a��ʾ�ӳ�����
	//����ֵΪ����
	Line NewExtent(bool p, double a);



	//��������
	~Line() {};

public:
	Point m_ptf;//ֱ�����
	Point m_pts;//ֱ���յ�
};
