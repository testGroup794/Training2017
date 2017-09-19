#include <bits/stdc++.h>
using namespace std;

const int MAXN =1007;
const double EPS=1e-10, PI=acos(-1.0);

double add(double a,double b) { // 实数加法，修正0的偏差
	return abs(a+b)<EPS*(abs(a)+abs(b)) ? 0 : a+b;
}

struct Point {
	double x,y;
	Point(){}
	Point(double x,double y):x(x),y(y){}
	Point(const Point &p):x(p.x),y(p.y){}
	Point operator-(Point p) { // 重载减号
		return Point(add(x,-p.x),add(y,-p.y));
	}
	double det(Point p) { // 外积 叉乘, 有向面积
		return add(x*p.y,-y*p.x);
	}
	double distance(Point p) { // 两点之间的距离
		return sqrt(add((x-p.x)*(x-p.x),(y-p.y)*(y-p.y)));
	}
	friend istream& operator >>(istream &is, Point &p) {
		is>>p.x>>p.y;
		return is;
	}
	friend ostream& operator <<(ostream &os, const Point &p) {
		os<<p.x<<" "<<p.y;
		return os;
	}
};

Point ori[MAXN];	// 输入的原始点集
Point convex[MAXN];	// 按逆时针方向排列在凸包中的各个顶点
int top;			// 栈顶

bool cmp(Point pa, Point pb) { // 极角排序对<重载
	double tmp=(pa-ori[0]).det(pb-ori[0]); //ori[0]为基点
	if(tmp>0) return 1;
	if(tmp==0 && pa.distance(ori[0]) < pb.distance(ori[0])) return 1;
	return 0;
}

void graham(int n) { // 寻找凸包 graham 扫描法 时间O(n)
	if (n==1) top=0, convex[0]=ori[0];
	if (n>=2) { // 单调栈
		top=1;
		convex[0]=ori[0];
		convex[1]=ori[1];
		for (int i=2; i<n; i++) {
			while (top>=1 && ((convex[top]-convex[top-1]).det(ori[i]-convex[top-1]))<=0) top--;
			convex[++top]=ori[i];
		}
	}
}

int main()
{
	int n;
	double r;
	ios::sync_with_stdio(false);
	while (cin>>n>>r) {
		for(int i=0;i<n;i++) cin>>ori[i];
		for(int i=1;i<n;i++) // ori[0] 存储的是 左下方的点
			if ((ori[i].y!=ori[0].y && ori[i].y<ori[0].y) || (ori[i].y==ori[0].y && ori[i].x<ori[0].x)) 
				swap(ori[i],ori[0]);
		sort(ori+1,ori+n,cmp); // 对 ori[1]-ori[n-1] 进行对 ori[0]的极角排序
		graham(n); // 求凸包
		cout<<top<<endl;
		double ans=convex[0].distance(convex[top]);
		for (int i=0; i<top; i++) ans+=convex[i].distance(convex[i+1]);
		cout<<fixed<<setprecision(0)<<ans+2*PI*r<<endl;
	}
	return 0;
}
