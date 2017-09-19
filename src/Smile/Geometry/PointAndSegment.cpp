#include <bits/stdc++.h>
using namespace std; // 输出注意 -0.00 要不要清掉

const int MAX_N =1007;
const double EPS=1e-7;
const double INF=40007;
const double PI=acos(-1.0);

bool zero(double x) { // 实数判0
	return abs(x)<EPS;
}

double cmp(double a,double b) { // 比较实数a和b的大小,返回a-b
	return abs(a-b)<EPS ? 0 : a-b;
}

bool between(double a,double l,double r) { // 判断a是否在[l,r]内
	if (cmp(l,r)>0) swap(l,r);
	return (cmp(l,a)<=0 && cmp(a,r)<=0);
}

double add(double a,double b) { //用在向量叉积点积上的实数加法，修正0的偏差
	return abs(a+b)<EPS*(abs(a)+abs(b)) ? 0 : a+b;
}

struct Point // 这货是向量,不过也可以当点用
{
	double x,y;
	Point(){}
	Point(double x,double y):x(x),y(y){}
	Point(const Point &p):x(p.x),y(p.y){}
	/*bool operator <(const Point &p) const { // 极角排序用
		if (y * p.y <= 0) {
			if (y > 0 || p.y > 0) return y < p.y;
			if (y == 0 && p.y == 0) return x < p.x;
		}
		return x*p.y-y*p.x > 0;
	}*/
	bool operator <(const Point &p) const {
		return cmp(x,p.x)<0 || (cmp(x,p.x)==0 && cmp(y,p.y)<0);
	}
	bool operator ==(const Point &p) const {
		return cmp(x,p.x)==0 && cmp(y,p.y)==0;
	}
	const Point& operator =(const Point& p) {
		x=p.x; y=p.y;
		return *this;
	}
	Point operator +(Point p) {
		return Point(add(x,p.x),add(y,p.y));
	}
	Point operator -(Point p) {
		return Point(add(x,-p.x),add(y,-p.y));
	}
	Point operator *(double d) {
		return Point(x*d,y*d);
	}
	Point operator /(double d) {
		return Point(x/d,y/d);
	}
	double dot(Point p) { // 内积 点乘
		return add(x*p.x,y*p.y);
	}
	double det(Point p) { // 外积 叉乘, 有向面积
		return add(x*p.y,-y*p.x);
	}
	double length() { // 向量长度
		return sqrt(length2());
	}
	double length2() { // 向量长度的平方
		return x*x+y*y;
	}
	double distance(Point p) { // 两点之间的距离
		return sqrt(add((x-p.x)*(x-p.x),(y-p.y)*(y-p.y)));
	}
	double distance2(Point p) { // 两点之间距离的平方
		return add((x-p.x)*(x-p.x),(y-p.y)*(y-p.y));
	}
	Point rotate(double alpha) { // 向量绕原点逆时针旋转alpha度，sin、cos有精度损失
		return Point(x*cos(alpha)-y*sin(alpha),x*sin(alpha)+y*cos(alpha));
	}
	Point rotate(double cosT, double sinT) { // 向量绕原点逆时针旋转T度，无精度损失
		return Point(x*cosT-y*sinT,x*sinT+y*cosT);
	}
	Point normal() { // 单位法线，即左转90°后长度归一化
		double l=length();
		return Point(-y/l,x/l);
	}
	double arg() { // 极角，值域(-pi,pi]
		return atan2(y,x);
	}
	double radian(Point p) { // 向量夹角，弧度
		if (det(p)>=0) return acos(dot(p)/length()/p.length());
		else return PI*2.0-acos(dot(p)/length()/p.length());
	}
	double angle(Point p) { // 向量夹角，角度
		return radian(p)*180.0/PI;
	}
	friend ostream& operator <<(ostream &os, const Point &p) {
		os<<p.x<<" "<<p.y;
		return os;
	}
	friend istream& operator >>(istream &is, Point &p) {
		is>>p.x>>p.y;
		return is;
	}
};

struct Segment // 线段/直线类， 如果有方向的话从s到t
{
	Point s,t;
	Segment(){}
	Segment(double sx,double sy,double tx,double ty):s(sx,sy),t(tx,ty){}
	Segment(Point s,Point t):s(s),t(t){}
	Segment(const Segment &p):s(p.s),t(p.t){}
	bool operator <(const Segment &p) const {
		return s<p.s || (s==p.s && t<p.t);
	}
	Segment midperpendicular() { // 返回线段的中垂线
		Point vec=(t-s).rotate(0,1);
		Point mid=(s+t)/2;
		return Segment(mid,mid+vec);
	}
	bool isSegment() { // 是不是线段（两个点重合的时候不是），枚举点来构造线段或者直线的时候一般都需要判
		return !(s==t);
	}
	double length() { // 线段的长度
		return (s-t).length();
	}
	bool inX(double x) { // x坐标在线段的x区间内
		return between(x,s.x,t.x);
	}
	bool inY(double y) { // y坐标在线段的y区间内
		return between(y,s.y,t.y);
	}
	double getY(double x) { // 对于线段所在直线，给出x求y, 要求该直线不能竖直
		double k=(s.y-t.y)/(s.x-t.x);
		double b=s.y-k*s.x;
		return k*x+b;
	}
	bool onLine(Point p) { // 点在直线上
		return zero((p-s).det(t-s));
	}
	bool onSegment(Point p) { // 点在线段上
		return zero((p-s).det(t-s)) && inX(p.x) && inY(p.y);
	}
	bool onLeft(Point p) { // 点在s→t向量的左边，在线上不算
		return cmp((t-s).det(p-s),0)>0;
	}
	double distPointToSegment(Point p) { // 点到线段的距离，注意是线段不是直线！
		if (cmp((p-s).dot(t-s),0)<0) return (p-s).length();
		if (cmp((p-t).dot(s-t),0)<0) return (p-t).length();
		return abs((s-p).det(t-p)/s.distance(t));
	}
	double distPointToLine(Point p) { // 点到直线的距离
		return abs((s-p).det(t-p)/s.distance(t));
	}
	Point pedal(Point p) { // 点到直线的垂足
		double r=(t-s).dot(p-s)/(t-s).dot(t-s);
		return s+(t-s)*r;
	}
	bool perpendicular(Segment p) { // 线段是否垂直
		return zero((s-t).dot(p.s-p.t));
	}
	bool parallel(Segment p) { // 线段是否平行
		return zero((s-t).det(p.s-p.t));
	}
	bool same(Segment p) { // 直线重合
		return onLine(p.s) && onLine(p.t);
	}
	bool lineAcrossLine(Segment p) { // 直线与直线相交
		return !parallel(p) && !same(p);
	}
	Point lineIntersectionLine(Segment p) { // 直线与直线求交点
		double a1=(p.t-s).det(p.s-s);
		double a2=(p.s-t).det(p.t-t);
		return (s*a2+t*a1)/(a2+a1);
	}
	bool lineAcrossSegment(Segment p) { // 直线与线段非严格相交, 即线段跨立直线
		return ((p.s-s).det(t-s))*((t-s).det(p.t-s))>=0;
	}
	bool across(Segment p) { // 线段与线段严格相交
		return max(s.x,t.x)>=min(p.s.x,p.t.x) && 
			min(s.x,t.x)<=max(p.s.x,p.t.x) && 
			max(s.y,t.y)>=min(p.s.y,p.t.y) && 
			min(s.y,t.y)<=max(p.s.y,p.t.y) && 
			((s-p.s).det(p.t-p.s))*((p.t-p.s).det(t-p.s))>0 && 
			((p.s-s).det(t-s))*((t-s).det(p.t-s))>0;
	}
	bool touch(Segment p) { // 线段与线段非严格相交
		return max(s.x,t.x)>=min(p.s.x,p.t.x) && 
			min(s.x,t.x)<=max(p.s.x,p.t.x) && 
			max(s.y,t.y)>=min(p.s.y,p.t.y) && 
			min(s.y,t.y)<=max(p.s.y,p.t.y) && 
			((s-p.s).det(p.t-p.s))*((p.t-p.s).det(t-p.s))>=0 && 
			((p.s-s).det(t-s))*((t-s).det(p.t-s))>=0;
	}
	Segment rotate(Point p, double alpha) { // 线段绕点p逆时针旋转alpha度
		Point ps=s-p;
		Point pt=t-p;
		ps=ps.rotate(alpha);
		pt=pt.rotate(alpha);
		return Segment(p+ps,p+pt);
	}
	friend ostream& operator <<(ostream &os, const Segment &p) {
		os<<p.s<<" "<<p.t;
		return os;
	}
	friend istream& operator >>(istream &is, Segment &p) {
		is>>p.s>>p.t;
		return is;
	}
};

struct Triangle // 三角形类
{
	Point a,b,c;
	Triangle(){}
	Triangle(double ax,double ay,double bx,double by,double cx,double cy):
		a(ax,ay),b(bx,by),c(cx,cy){}
	Triangle(Point a,Point b,Point c):a(a),b(b),c(c){}
	Triangle(const Triangle &p):a(p.a),b(p.b),c(p.c){}
	bool operator ==(const Triangle &p) const { // 三角形全等SSS ##TODO
		double l[3],pl[3];
		l[0]=Segment(a,b).length();
		l[1]=Segment(b,c).length();
		l[2]=Segment(c,a).length();
		sort(l,l+3);
		pl[0]=Segment(p.a,p.b).length();
		pl[1]=Segment(p.b,p.c).length();
		pl[2]=Segment(p.c,p.a).length();
		sort(pl,pl+3);
		for (int i=0;i<3;i++) if (cmp(l[i],pl[i])!=0) return false;
		return true;
	}
	bool isTriangle() { // 是三角形，任意两边之和大于第三边
		return cmp(a.distance(b)+a.distance(c),b.distance(c))>0 &&
			cmp(b.distance(a)+b.distance(c),a.distance(c))>0 &&
			cmp(c.distance(b)+c.distance(a),b.distance(a))>0;
	}
	double perimeter() { // 三角形周长
		return (a-b).length()+(a-c).length()+(b-c).length();
	}
	double area() { // 三角形面积
		return abs((b-a).det(c-a))/2;
	}
	Point massCenter() { // 重心，三条中线交点
		return(a+b+c)/3;
	}
	Point circumCenter() { // 外心，外接圆圆心，三条中垂线交点
		Point res;
		double a1=b.x-a.x, b1=b.y-a.y, c1=(a1*a1+b1*b1)/2;
		double a2=c.x-a.x, b2=c.y-a.y, c2=(a2*a2+b2*b2)/2;
		double d=a1*b2-a2*b1;
		res.x=a.x+(c1*b2-c2*b1)/d;
		res.y=a.y+(a1*c2-a2*c1)/d;
		return res;
	}
	Point orthoCenter() { // 垂心， 三条垂线交点, 根据欧拉定理通过重心和外心求
		return massCenter()*3-circumCenter()*2;
	}
	Point innerCenter() { // 内心， 内切圆圆心， 三条角平分线交点
		Point res;
		double la=(b-c).length();
		double lb=(c-a).length();
		double lc=(a-b).length();
		res.x=(la*a.x+lb*b.x+lc*c.x)/(la+lb+lc);
		res.y=(la*a.y+lb*b.y+lc*c.y)/(la+lb+lc);
		return res;
	}
	Point FermatPoint() { // 费马点， 离三个顶点距离和最小的点 ####TODO
		if (cmp((b-a).radian(c-a)/PI,120.0/180.0)<=0) return a;
		if (cmp((c-b).radian(a-b)/PI,120.0/180.0)<=0) return b;
		if (cmp((a-c).radian(b-c)/PI,120.0/180.0)<=0) return c;
		bool flag=cmp((b-a).det(c-a),0)>0;
		Point b1,c1;
		if (flag) {
			c1=a+(b-a).rotate(-PI/3);
			b1=a+(c-a).rotate(PI/3);
		} else {
			c1=a+(b-a).rotate(PI/3);
			b1=a+(c-a).rotate(-PI/3);
		}
		return (Segment(b,b1).lineIntersectionLine(Segment(c,c1)));
	}
	Triangle rotate(Point p, double alpha) { // 三角形绕点p逆时针旋转alpha度
		Point pa=a-p;
		Point pb=b-p;
		Point pc=c-p;
		pa=pa.rotate(alpha);
		pb=pb.rotate(alpha);
		pc=pc.rotate(alpha);
		return Triangle(p+pa,p+pb,p+pc);
	}
	bool pointInTriangle(Point p) { // 点在三角形内（边界不算） ##TODO
		if (cmp(((b-a).det(p-a))*((c-a).det(p-a)),0)>=0) return false;
		if (cmp(((a-b).det(p-b))*((c-b).det(p-b)),0)>=0) return false;
		if (cmp(((b-c).det(p-c))*((a-c).det(p-c)),0)>=0) return false;
		return true;
	}
	friend istream& operator >>(istream &is, Triangle &p) {
		is>>p.a>>p.b>>p.c;
		return is;
	}
	friend ostream& operator <<(ostream &os, const Triangle &p) {
		os<<p.a<<" "<<p.b<<" "<<p.c;
		return os;
	}
};

struct Circle // 圆类
{
	Point p; // 圆心
	double r; // 半径
	Circle(){}
	Circle(double sx,double sy,double r):p(sx,sy),r(r){}
	Circle(Point p,double r):p(p),r(r){}
	Circle(const Circle &c):p(c.p),r(c.r){}
	bool operator <(const Circle &o) const {
		if (cmp(r,o.r)!=0) return cmp(r,o.r)<0;
		if (cmp(p.x,o.p.x)!=0) return cmp(p.x,o.p.x)<0;
		return cmp(p.y,o.p.y)<0;
	}
	bool operator ==(const Circle &o) const {
		return p==o.p && cmp(r,o.r)==0;
	}
	vector<Point> lineCrossCircle(Segment l) { // 直线与圆求交点，res[0]res[1]与传入的直线l.s和l.t同方向
		vector<Point> res;
		res.clear();
		double xa=p.x, ya=p.y;
		double xb=l.s.x, yb=l.s.y;
		double xc=l.t.x, yc=l.t.y;
		double dx=xc-xb, dy=yc-yb;
		double A=dx*dx+dy*dy;
		double B=2*dx*(xb-xa)+2*dy*(yb-ya);
		double C=(xb-xa)*(xb-xa)+(yb-ya)*(yb-ya)-r*r;
		double delta=B*B-A*C*4;
		if (cmp(delta,0)>=0) {
			double ta=(-B-sqrt(delta))/(A*2);
			double tb=(-B+sqrt(delta))/(A*2);
			res.push_back(Point(xb+ta*dx,yb+ta*dy));
			res.push_back(Point(xb+tb*dx,yb+tb*dy));
		}
		return res;
	}
	vector<Point> segmentCrossCircle(Segment l) { // 线段与圆求交点，res[0]res[1]与传入的直线l.s和l.t同方向
		vector<Point> res;
		res.clear();
		double xa=p.x, ya=p.y;
		double xb=l.s.x, yb=l.s.y;
		double xc=l.t.x, yc=l.t.y;
		double dx=xc-xb, dy=yc-yb;
		double A=dx*dx+dy*dy;
		double B=2*dx*(xb-xa)+2*dy*(yb-ya);
		double C=(xb-xa)*(xb-xa)+(yb-ya)*(yb-ya)-r*r;
		double delta=B*B-A*C*4;
		if (cmp(delta,0)>=0) {
			double ta=(-B-sqrt(delta))/(A*2);
			double tb=(-B+sqrt(delta))/(A*2);
			if (cmp(ta,1.0)<=0 && cmp(ta,0.0)>=0)
				res.push_back(Point(xb+ta*dx,yb+ta*dy));
			if (cmp(tb,1.0)<=0 && cmp(tb,0.0)>=0)
				res.push_back(Point(xb+tb*dx,yb+tb*dy));
		}
		return res;
	}
	bool circleCrossCircle(Circle c) { // 圆与圆判断是否相交
		return cmp(p.distance(c.p),r+c.r)<=0;
	}
	pair <Point,Point> circleIntersectionCircle(Circle c) { // 圆与圆求交点
		double d=(p-c.p).length();
		double cosT = (r*r+d*d-c.r*c.r)/(d*r*2);
		double sinT = sqrt(1.0 - cosT*cosT);
		Point v = (c.p-p)/(c.p-p).length()*r;
		return make_pair(p+v.rotate(cosT,-sinT),p+v.rotate(cosT,sinT));
	}
	bool between(Point a, Point b, Point t) { // 圆上的点t在优弧AB上，包含端点
		Point pa=a-p; Point pb=b-p; Point pt=t-p;
		return (cmp(pa.det(pt)*pb.det(pt),0)>0) ? false : Segment(a,b).touch(Segment(p,t));
	}
	bool between(Point a, Point b, Point c, Point t) { // 圆上的点t在弧ABC上，包含端点（AC为弧的两端，B是弧上一点），用线段判交旋转实现
		return between(a,b,t)==between(a,b,c);
	}
	bool betweenABC(Point a, Point b, Point c, Point t) { // 圆上的点t在弧ABC上，包含端点（AC为弧的两端，B是弧上一点），用向量旋转实现
		Point ta=a-p;
		Point tb=b-p;
		Point tc=c-p;
		Point tt=t-p;
		double sinT=ta.y/r;
		double cosT=-ta.x/r;
		Point pa=ta.rotate(cosT,sinT);
		Point pb=tb.rotate(cosT,sinT);
		Point pc=tc.rotate(cosT,sinT);
		Point pt=tt.rotate(cosT,sinT);
		if (cmp(pc.arg(),pb.arg())<0) {
			swap(ta,tc);
			sinT=ta.y/r;
			cosT=-ta.x/r;
			pa=ta.rotate(cosT,sinT);
			pb=tb.rotate(cosT,sinT);
			pc=tc.rotate(cosT,sinT);
			pt=tt.rotate(cosT,sinT);
		}
		return cmp(pt.arg(),pc.arg())<=0;
	}
	friend istream& operator >>(istream &is, Circle &c) {
		is>>c.p>>c.r;
		return is;
	}
	friend ostream& operator <<(ostream &os, const Circle &c) {
		os<<c.p<<" "<<c.r;
		return os;
	}
};

Point getConvexCmpPoint;

bool getConvexCmp(Point pa, Point pb) { // 求凸包用的极角排序，对<重载
	double tmp=(pa-getConvexCmpPoint).det(pb-getConvexCmpPoint); //getConvexCmpPoint为基点
	if(tmp>0) return 1;
	if(tmp==0 && pa.distance(getConvexCmpPoint) < pb.distance(getConvexCmpPoint)) return 1;
	return 0;
}

struct Polygon // 多边形类
{
	int n; // 顶点数
	Point pt[MAX_N]; // 逆时针排列的所有顶点，实际存放n+1个点，pt[n]=pt[0]
	Polygon(){}
	bool isConvex() { // 是否凸多边形（顺时针逆时针输入多边形均可）
		pt[n+1]=pt[1]; // 存放n+2个点，pt[n+1]=pt[1]，方便写循环
		bool posi=false, nega=false;
		for (int i=1;i<=n;i++) {
			double x=(pt[i-1]-pt[i]).det(pt[i+1]-pt[i]);
			if (cmp(x,0)>0) posi=true;
			if (cmp(x,0)<0) nega=true;
		}
		return !(posi && nega);
	}
	double perimeter() { // 周长
		double res=0;
		for (int i=0; i<n; i++) res+=pt[i].distance(pt[i+1]);
		return res;
	}
	double area() { // 面积
		double res=0;
		for (int i=0; i<n; i++) res+=pt[i].det(pt[i+1]);
		return res/2;
	}
	int pointInPolygon(Point p) { // 转角法判断点和多边形的关系，-1在多边形外，1在多边型内，0在多边形上
		pt[n]=pt[0];
		int wn=0;
		for (int i=0;i<n;i++) {
			if (Segment(pt[i],pt[i+1]).onSegment(p)) return 0; // 点在多边形边界上
			double k=cmp((pt[i+1]-pt[i]).det(p-pt[i]),0);
			double d1=cmp(pt[i].y,p.y);
			double d2=cmp(pt[i+1].y,p.y);
			if (k>0 && d1<=0 && d2>0) wn++;
			if (k<0 && d1>0 && d2<=0) wn--;
		}
		return wn?1:-1;
	}
	void graham(int n, Polygon &convex) { // 寻找凸包 graham 扫描法 时间O(n)
		int top=convex.n=0;
		if (n==1) top=0, convex.pt[0]=pt[0];
		if (n>=2) { // 单调栈
			top=1;
			convex.pt[0]=pt[0];
			convex.pt[1]=pt[1];
			for (int i=2; i<n; i++) {
				while (top>=1 && ((convex.pt[top]-convex.pt[top-1])
					.det(pt[i]-convex.pt[top-1]))<=0) top--;
				convex.pt[++top]=pt[i];
			}
		}
		convex.n=++top;
		convex.pt[convex.n]=convex.pt[0];
	}
	int getConvex(Polygon &convex) { // 求凸包
		for(int i=1;i<n;i++) // pt[0] 存储的是 左下方的点
			if ((pt[i].y!=pt[0].y && pt[i].y<pt[0].y) || 
				(pt[i].y==pt[0].y && pt[i].x<pt[0].x)) 
				swap(pt[i],pt[0]);
		getConvexCmpPoint=pt[0];
		sort(pt+1,pt+n,getConvexCmp); // 对 pt[1]-pt[n-1] 进行对 pt[0]的极角排序
		graham(n,convex); // 求凸包
		return convex.n;
	}
	Point massCenter() { // 多边形的重心(凹凸都可以) 
		double sum=0.0, sumx=0, sumy=0;
		Point p1=pt[0], p2=pt[1], p3;
		for (int i=2;i<=n-1;i++) {
			p3=pt[i];
			double area=(p2-p1).det(p3-p2)/2.0;
			sum+=area;
			sumx+=(p1.x+p2.x+p3.x)*area;
			sumy+=(p1.y+p2.y+p3.y)*area;
			p2=p3;
		}
		return Point (sumx/(3.0*sum), sumy/(3.0*sum));
	}
	friend istream& operator >>(istream &is, Polygon &p) {
		is>>p.n;
		for (int i=0;i<p.n;i++) is>>p.pt[i];
		if (p.n==0) is.clear(ios::failbit);
		p.pt[p.n]=p.pt[0];
		return is;
	}
	friend ostream& operator <<(ostream &os, const Polygon &p) {
		os<<p.n<<endl;
		for (int i=0;i<p.n;i++) os<<p.pt[i]<<" ";
		os<<endl;
		return os;
	}
};

struct HalfPlane // 半平面向量，左边区域即半平面
{
	Point p;
	Point v;
	double angle;
	HalfPlane(){}
	HalfPlane(Point p, Point v):p(p),v(v){angle=atan2(v.y,v.x);}
	bool operator <(const HalfPlane& l) const {
		return angle<l.angle;
	}
	Segment toSegment() { // 半平面的线段表示
		return Segment(p,p+v);
	}
	bool onLeft(Point t) { // 点在半平面上（注意半平面的边界算不算，如果不算把||后面部分去掉）
		return toSegment().onLeft(t) || toSegment().onLine(t);
	}
	Point getIntersection(HalfPlane l) { // 两个半平面向量的交点
		return toSegment().lineIntersectionLine(l.toSegment());
	}
};

struct HalfPlaneIntersection // 半平面交
{
	HalfPlane l[MAX_N];
	Polygon poly;
	int n;
	HalfPlaneIntersection(){}
	HalfPlaneIntersection(const Polygon &pl) { // 由逆时针排列的多边形构造，用于计算多边形的核
		n=pl.n;
		for (int i=0;i<n;i++) {
			l[i]=HalfPlane(pl.pt[i],(Point)pl.pt[i+1]-pl.pt[i]);
		}
	}
	int init() {
		sort (l,l+n);
	}
	int getResult() { // 返回半平面交上点的个数，结果保存在poly里，边界也算的时候结果可能会有重点
		int fst=0,lst=0; // 双端队列首尾下标
		Point p[n]; // p[i]为q[i]和q[i+1]的交点
		HalfPlane q[n]; // 双端队列，保存当前在半平面交边界上的半平面向量
		q[fst]=l[0]; // 初始状态双端队列只有一个半平面
		for (int i=1;i<n;i++) {
			while (fst<lst && !l[i].onLeft(p[lst-1])) lst--;
			while (fst<lst && !l[i].onLeft(p[fst])) fst++;
			q[++lst]=l[i];
			if (zero(q[lst].v.det(q[lst-1].v))) { // 平行且同向的两个半平面，取里面那个
				lst--;
				if (q[lst].onLeft(l[i].p)) q[lst]=l[i];
			}
			if (fst<lst) p[lst-1]=q[lst].getIntersection(q[lst-1]);
		}
		while (fst<lst && !q[fst].onLeft(p[lst-1])) lst--; // 删除无用的半平面
		if (lst-fst<=1) return 0; // 结果是空集
		p[lst]=q[lst].getIntersection(q[fst]); // 首尾两个半平面的交点
		poly.n=0; // 返回结果
		for (int i=fst;i<=lst;i++) poly.pt[poly.n++]=p[i];
		poly.pt[poly.n]=poly.pt[0];
		return poly.n;
	}
};

struct CircleUnionArea // 圆面积并
{
	// const int MAX_N=307;
	struct Node
	{
		Point p;
		double a;
		int d;
		Node(){}
		Node(double sx,double sy,double a,int d):p(sx,sy),a(a),d(d){}
		Node(Point p,double a,int d):p(p),a(a),d(d){}
		Node(const Node &c):p(c.p),a(c.a),d(c.d){}
		bool operator <(const Node &o) const {
			return a<o.a;
		}
	};
	Circle c[MAX_N],tc[MAX_N]; // tc 存放待求的圆， m 为待求圆的个数
	int n,m;
	double solve() { // 圆面积并
		n=0;
		sort(tc,tc+m);
		m=unique(tc,tc+m)-tc;
		for (int i=m-1;i>=0;i--) {
			bool ok=true;
			for (int j=i+1;j<m;j++) {
				double d=(tc[i].p-tc[j].p).length();
				if (cmp(d,abs(tc[i].r-tc[j].r))<=0) {
					ok=false;
					break;
				}
			}
			if (ok) c[n++]=tc[i];
		}
		double ans=0;
		for (int i=0;i<n;i++) {
			vector<Node> ev;
			ev.clear();
			Point boundary=c[i].p+Point(-c[i].r,0);
			ev.push_back(Node(boundary,-PI,0));
			ev.push_back(Node(boundary,PI,0));
			for (int j=0;j<n;j++) {
				if (i==j) continue;
				double d=(c[i].p-c[j].p).length();
				if (cmp(d,(c[i].r+c[j].r))<0) {
					pair<Point,Point> ret=c[i].circleIntersectionCircle(c[j]);
					double x=(ret.first-c[i].p).arg();
					double y=(ret.second-c[i].p).arg();
					if (cmp(x,y)>0) {
						ev.push_back(Node(ret.first,x,1));
						ev.push_back(Node(boundary,PI,-1));
						ev.push_back(Node(boundary,-PI,1));
						ev.push_back(Node(ret.second,y,-1));
					} else {
						ev.push_back(Node(ret.first,x,1));
						ev.push_back(Node(ret.second,y,-1));
					}
				}
			}
			sort(ev.begin(),ev.end());
			int sum=ev[0].d;
			for (int j=1;j<ev.size();j++) {
				if (sum==0) {
					ans+=ev[j-1].p.det(ev[j].p)/2;
					double x=ev[j-1].a;
					double y=ev[j].a;
					double area=c[i].r*c[i].r*(y-x)/2;
					Point v1=ev[j-1].p-c[i].p;
					Point v2=ev[j].p-c[i].p;
					area-=v1.det(v2)/2;
					ans+=area;
				}
				sum+=ev[j].d;
			}
		}
		return ans;
	}
};


int i,j,k,n,m,x,y,z,ans,cnt,tcase,xcase;
double alpha;
Polygon pg,cov;

int main()
{
	// freopen("s.in","r",stdin);
	// freopen("s.out","w",stdout);
	ios::sync_with_stdio(false);
	while (cin>>pg) {
		cout<<fixed<<setprecision(2)<<pg.massCenter()<<endl;
	}
	return 0;
}
