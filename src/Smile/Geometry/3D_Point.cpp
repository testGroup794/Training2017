#include <iostream>
#include <iomanip>
#include <algorithm>
#include <string>
#include <sstream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <ctime>
#include <bitset>
#include <set>
#include <map>
#include <queue>
#include <stack>
#include <deque>
#include <vector>
#include <list>
using namespace std;

const int MAXN=37; 
const int MOD=1000000007; 
const int INF=2147483647;
const double EPS=1e-8;
const double PI=acos(-1.0);

bool zero(double x)
{
	return abs(x)<EPS;
}

double cmp(double a,double b)// 比较实数a和b的大小,返回a-b
{
	if (abs(a-b)<EPS)
	{
		return 0;
	}
	return a-b;
}

bool between(double a,double l,double r)// 判断a是否在[l,r]内
{
	if (cmp(l,r)>0)
	{
		swap(l,r);
	}
	if (cmp(l,a)<=0 && cmp(a,r)<=0)
	{
		return true;
	}
	else
	{
		return false;
	}
}

struct Point // 三维点，三维向量
{
	double x,y,z;
	Point(){}
	Point(double x,double y,double z):x(x),y(y),z(z){}
	Point(const Point &p):x(p.x),y(p.y),z(p.z){}
	bool operator <(const Point &p) const
	{
		return (x<p.x || x==p.x && y<p.y || x==p.x && y==p.y && z<p.z);
	}
	bool operator ==(const Point &p) const
	{
		return (x==p.x && y==p.y && z==p.z);
	}
	const Point& operator =(const Point& p)
	{
		x=p.x;
		y=p.y;
		z=p.z;
		return *this;
	}
	Point operator+(const Point p) 
	{
		return Point(x+p.x,y+p.y,z+p.z);
	}
	Point operator-(const Point p) 
	{
		return Point(x-p.x,y-p.y,z-p.z);
	}
	Point operator*(const double p) 
	{
		return Point(x*p,y*p,z*p);
	}
	Point operator/(const double p) 
	{
		return Point(x/p,y/p,z/p);
	}
	Point std() // 单位向量
	{
		return Point(x/length(),y/length(),z/length());
	}
	double dot(const Point p) // 点积
	{
		return x*p.x+y*p.y+z*p.z;
	}
	Point det(const Point p) // 叉积
	{
		return Point(y*p.z-z*p.y,z*p.x-x*p.z,x*p.y-y*p.x);
	}
	double detLength(Point p) // 叉乘大小
	{
		return det(p).length();
	}
	bool sameDirection(Point p)
	{
		return cmp(dot(p),0)>0;
	}
	double length() // 向量长度
	{
		return sqrt(length2());
	}
	double length2() // 向量长度的平方
	{
		return x*x+y*y+z*z;
	}
	double distance(Point p) // 两点之间的距离
	{
		return sqrt(distance2(p));
	}
	double distance2(Point p) // 两点之间距离的平方
	{
		return (x-p.x)*(x-p.x)+(y-p.y)*(y-p.y)+(z-p.z)*(z-p.z);
	}
	double radian(Point &p) // 向量夹角，弧度
	{
		return acos(dot(p)/length()/p.length());
	}
	double angle(Point &p) // 向量夹角，角度
	{
		return radian(p)*180.0/PI;
	}
	Point rotate(Point b, double alpha) // 向量绕OB轴旋转alpha弧度,沿着旋转轴往原点看旋转的角度为顺时针
	{
		Point n=b.std();
		double s=sin(alpha),c=cos(alpha);
		double tx=x*(n.x*n.x*(1-c)+c)+y*(n.x*n.y*(1-c)-n.z*s)+z*(n.x*n.z*(1-c)+n.y*s);
		double ty=x*(n.x*n.y*(1-c)+n.z*s)+y*(n.y*n.y*(1-c)+c)+z*(n.y*n.z*(1-c)-n.x*s);
		double tz=x*(n.x*n.z*(1-c)-n.y*s)+y*(n.y*n.z*(1-c)+n.x*s)+z*(n.z*n.z*(1-c)+c);
		return Point(tx,ty,tz);
	}
	friend ostream& operator <<(ostream &os, const Point &p)
	{
		if (zero(p.x)) os<<0.0<<" ";
			else os<<p.x<<" ";
		if (zero(p.y)) os<<0.0<<" ";
			else os<<p.y<<" ";
		if (zero(p.z)) os<<0.0;
			else os<<p.z;
		return os;
	}
	friend istream& operator >>(istream &is, Point &p)
	{
		is>>p.x>>p.y>>p.z;
		return is;
	}
};

struct Line
{
	Point s,t;
	Line(){}
	Line(double sx,double sy,double sz,double tx,double ty,double tz):s(sx,sy,sz),t(tx,ty,tz){}
	Line(Point &s,Point &t):s(s),t(t){}
	bool operator <(const Line &p) const
	{
		return (s<p.s || s==p.s && t<p.t);
	}
	bool isSegment() // 是不是线段（两个点重合的时候不是），枚举点来构造线段或者直线的时候一般都需要判（##TODO：如果WA把这里加上EPS试试）
	{
		return !(s==t);
	}
	Line move(Point dist)
	{
		s=s+dist;
		t=t+dist;
		return *this;
	}
	bool inX(double x) // x坐标在线段的x区间内
	{
		return between(x,s.x,t.x);
	}
	bool inY(double y) // y坐标在线段的y区间内
	{
		return between(y,s.y,t.y);
	}
	bool inZ(double z) // y坐标在线段的y区间内
	{
		return between(z,s.z,t.z);
	}
	bool onLine(Point p) // 点在直线上
	{
		if (zero((p-s).det(t-s).length2())) return true;
		else return false;
	}
	bool onSegment(Point p) // 点在线段上
	{
		if (zero((p-s).detLength(t-s)) && between(p.x,s.x,t.x) && between(p.y,s.y,t.y) && between(p.z,s.z,t.z)) return true;
		else return false;
	}
	double distPointToSegment(Point &p) // 点到线段的距离，注意是线段不是直线！
	{
		if (cmp((p-s).dot(t-s),0)<0)
		{
			return (p-s).length();
		}
		if (cmp((p-t).dot(s-t),0)<0)
		{
			return (p-t).length();
		}
		return abs((s-p).detLength(t-p)/s.distance(t));
	}
	double distPointToLine(Point &p) // 点到直线的距离
	{
		return abs((s-p).detLength(t-p)/s.distance(t));
	}
	Point pedal(Point &p) // 点到直线的垂足
	{
		double r=(t-s).dot(p-s)/(t-s).dot(t-s);
		Point res=s+(t-s)*r;
		return res;
	}
	bool perpendicular(Line &p) // 线段是否垂直
	{
		return zero((s-t).dot(p.s-p.t));
	}
	bool parallel(Line &p) // 线段是否平行
	{
		return zero((s-t).detLength(p.s-p.t));
	}
	bool same(Line &p) // 直线重合
	{
		return onLine(p.s) && onLine(p.t);
	}
	bool samePlane(Line &p) // 直线共面
	{
		return zero((s-t).det(t-p.s).dot(p.t-s));
	}
	bool lineAcrossLine(Line &p) // 直线与直线相交
	{
		return samePlane(p) && !parallel(p) && !same(p);
	}
	Point lineIntersectionLine(Line &p) // 直线与直线求交点
	{
		Point p1=(p.t-s).det(p.s-s);
		Point p2=(p.s-t).det(p.t-t);
		double a1,a2;
		a1=p1.length();
		if (p1.sameDirection(p2))
		{
			a2=p2.length();
		}
		else
		{
			a2=-p2.length();
		}
		double x=(s.x*a2+t.x*a1)/(a1+a2);
		double y=(s.y*a2+t.y*a1)/(a1+a2);
		double z=(s.z*a2+t.z*a1)/(a1+a2);
		if (zero(x)) x=0.0; // 防止-0.00的情况出现
		if (zero(y)) y=0.0;
		if (zero(z)) z=0.0;
		return Point(x,y,z);
	}
	double distLineToLine(Line &p) // 直线到直线的距离
	{
		if (samePlane(p))
		{
			if (parallel(p))
			{
				return distPointToLine(p.s);
			}
			else
			{
				return 0.0;
			}
		}
		else
		{
			Point q=(s-t).det(p.s-p.t);
			return abs((s-p.s).dot(q))/q.length();
		}
	}
	Line commonPerpendicularLine(Line &l) // 直线到直线的公垂线
	{
		Point p=t-s;
		Point q=l.t-l.s;
		Point t=p.det(q);
		double a=p.x*p.y*q.y-p.y*p.y*q.x-p.z*p.z*q.x+p.x*p.z*q.z;
		double b=p.x*p.x*q.y-p.x*p.y*q.x-p.y*p.z*q.z+p.z*p.z*q.y;
		double c=p.x*p.z*q.x-p.x*p.x*q.z-p.y*p.y*q.z+p.y*p.z*q.y;
		double d=-s.x*a+s.y*b-s.z*c;
		double k=(-l.s.x*a+l.s.y*b-l.s.z*c-d)/(q.x*a-q.y*b+q.z*c);
		double x=q.x*k+l.s.x;
		double y=q.y*k+l.s.y;
		double z=q.z*k+l.s.z;
		p=Point(x,y,z);
		t=p+t;
		Line m=Line(p,t);
		q=lineIntersectionLine(m);
		return Line(q,p);
	}
	friend ostream& operator <<(ostream &os, const Line &p)
	{
		os<<p.s<<" "<<p.t;
		return os;
	}
	friend istream& operator >>(istream &is, Line &p)
	{
		is>>p.s>>p.t;
		return is;
	}
};

struct Plane
{
	Point a,b,c;
	Plane(){}
	Plane(double ax,double ay,double az,double bx,double by,double bz,double cx,double cy,double cz):a(ax,ay,az),b(bx,by,bz),c(cx,cy,cz){}
	Plane(Point &a,Point &b,Point &c):a(a),b(b),c(c){}
	bool isPlane()
	{
		Line l=Line(a,b);
		return !l.onLine(c);
	}
	Plane move(Point dist)
	{
		a=a+dist;
		b=b+dist;
		c=c+dist;
		return *this;
	}
	Point pvec()
	{
		return (a-b).det(b-c);
	}
	bool pointOnPlane(Point p)
	{
		return zero(pvec().dot(p-a));
	}
	bool pointsAtSameSide(Point p1, Point p2)
	{
		return cmp(pvec().dot(p1-a)*pvec().dot(p2-a),0)>0;
	}
	double distancePointToPlane(Point p)
	{
		return abs(pvec().dot(p-a))/pvec().length();
	}
	bool lineOnPlane(Line l)
	{
		return pointOnPlane(l.s) && pointOnPlane(l.t);
	}
	bool lineParallelPlane(Line l)
	{
		return cmp((l.s-l.t).dot(pvec()),0)<0;
	}
	bool lineAcrossPlane(Line l)
	{
		return !lineOnPlane(l) && !lineParallelPlane(l);
	}
	Point lineIntersectionPlane(Line l)
	{
		double t=pvec().dot(a-l.s)/pvec().dot(l.t-l.s);
		return l.s+(l.t-l.s)*t;
	}
	bool planeParallelPlane(Plane p)
	{
		return zero(pvec().det(p.pvec()).length());
	}
	bool planePerpendicularPlane(Plane p)
	{
		return zero(pvec().dot(p.pvec()));
	}
	bool planeAcrossPlane(Plane p)
	{
		return !planeParallelPlane(p);
	}
	Line planeIntersectionPlane(Plane p)
	{
		Line res;
		if (lineParallelPlane(Line(p.a,p.b)))
		{
			res.s=lineIntersectionPlane(Line(p.b,p.c));
		}
		else
		{
			res.s=lineIntersectionPlane(Line(p.a,p.b));
		}
		if (lineParallelPlane(Line(p.c,p.a)))
		{
			res.t=lineIntersectionPlane(Line(p.b,p.c));
		}
		else
		{
			res.t=lineIntersectionPlane(Line(p.c,p.a));
		}
		return res;
	}
	friend ostream& operator <<(ostream &os, const Plane &p)
	{
		os<<p.a<<" "<<p.b<<" "<<p.c;
		return os;
	}
	friend istream& operator >>(istream &is, Plane &p)
	{
		is>>p.a>>p.b>>p.c;
		return is;
	}
};

Line commonPerpendicularLine(Line l , Line r) // 直线到直线的公垂线
{
	Point a=l.t-l.s;
	Point b=r.t-r.s;
	Point t=a.det(b);
	Point s=l.s+t;
	Plane p=Plane(l.s,l.t,s);
	Point c=p.lineIntersectionPlane(r);
	Point e=c+t;
	Line w=Line(c,e);
	Point d=w.lineIntersectionLine(l);
	return Line(d,c);
}

struct Tetrahedron // 四面体类
{
	Point A,B,C,D;
	double a,b,c,d,e,f;
	Tetrahedron(){}
	Tetrahedron(double a,double b,double c,double d,double e,double f):a(a),b(b),c(c),d(d),e(e),f(f){}
	Tetrahedron(Point a,Point b,Point c,Point d):A(a),B(b),C(c),D(d){}
	double volumeByPoints()
	{
		Point ta=B-A,tb=C-A,tc=D-A;
		return abs(ta.det(tb).dot(tc))/6.0;
	}
	double volumeBySegments() // 6 segments must be sorted as: OA OB OC AB AC BC
	{
		double x,y;
		double aa=a*a,bb=b*b,cc=c*c,dd=d*d,ee=e*e,ff=f*f;
		x=4.0*aa*bb*cc-aa*(bb+cc-ff)*(bb+cc-ff)-bb*(cc+aa-ee)*(cc+aa-ee);
		y=cc*(aa+bb-dd)*(aa+bb-dd)-(aa+bb-dd)*(bb+cc-ff)*(cc+aa-ee);
		return sqrt(x-y)/12.0;
	}
};

int main()
{
	// freopen("s.in","r",stdin);
	// freopen("s.out","w",stdout);
	Point a,b;
	Line l,r;
	Plane p,q;
	double angle;
	while (cin>>a>>b>>angle)
	{
		cout<<a.rotate(b,angle)<<endl;
	}
	return 0;
}