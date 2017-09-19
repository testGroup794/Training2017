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
#define PR 1e-9
#define N 1007
struct TPoint
{
    double x,y,z;
    TPoint(){}
    TPoint(double _x,double _y,double _z):x(_x),y(_y),z(_z){}
    TPoint operator *(double d) {return TPoint(x*d,y*d,z*d);}
    TPoint operator / (double d) {return TPoint(x/d,y/d,z/d);}
    TPoint operator+(const TPoint p) {return TPoint(x+p.x,y+p.y,z+p.z);}
    TPoint operator-(const TPoint p) {return TPoint(x-p.x,y-p.y,z-p.z);}
    TPoint operator*(const TPoint p) {return TPoint(y*p.z-z*p.y,z*p.x-x*p.z,x*p.y-y*p.x);}//叉积
    double operator^(const TPoint p) {return x*p.x+y*p.y+z*p.z;}//点积
};
struct face//
{
    int a,b,c;//凸包一个面上的三个点的编号
    bool ok;//该面是否是最终凸包中的面
};
struct T3dhull
{
    int n;//初始点数
    TPoint ply[N];//初始点
    int trianglecnt;//凸包上三角形数
    face tri[N];//凸包三角形
    int vis[N][N];//点i到点j是属于哪个面
    double dist(TPoint a){return sqrt(a.x*a.x+a.y*a.y+a.z*a.z);}//两点长度
    double area(TPoint a,TPoint b,TPoint c){return dist((b-a)*(c-a));}//三角形面积*2
    double volume(TPoint a,TPoint b,TPoint c,TPoint d){return (b-a)*(c-a)^(d-a);}//四面体有向体积*6
    double ptoplane(TPoint &p,face &f)//正：点在面同向
    {
        TPoint m=ply[f.b]-ply[f.a],n=ply[f.c]-ply[f.a],t=p-ply[f.a];
        return (m*n)^t;
    }
    void deal(int p,int a,int b)
    {
        int f=vis[a][b];
        face add;
        if(tri[f].ok)
        {
            if((ptoplane(ply[p],tri[f]))>PR) dfs(p,f);
            else
            {
                add.a=b,add.b=a,add.c=p,add.ok=1;
                vis[p][b]=vis[a][p]=vis[b][a]=trianglecnt;
                tri[trianglecnt++]=add;
            }
        }
    }
    void dfs(int p,int cnt)//维护凸包，如果点p在凸包外更新凸包
    {
        tri[cnt].ok=0;
        deal(p,tri[cnt].b,tri[cnt].a);
        deal(p,tri[cnt].c,tri[cnt].b);
        deal(p,tri[cnt].a,tri[cnt].c);
    }
    bool same(int s,int e)//判断两个面是否为同一面
    {
        TPoint a=ply[tri[s].a],b=ply[tri[s].b],c=ply[tri[s].c];
        return fabs(volume(a,b,c,ply[tri[e].a]))<PR
            &&fabs(volume(a,b,c,ply[tri[e].b]))<PR
            &&fabs(volume(a,b,c,ply[tri[e].c]))<PR;
    }
    void construct()//构建凸包
    {
        int i,j;
        trianglecnt=0;
        if(n<4) return ;
        bool tmp=true;
        for(i=1;i<n;i++)//前两点不共点
        {
            if((dist(ply[0]-ply[i]))>PR)
            {
                swap(ply[1],ply[i]); tmp=false; break;
            }
        }
        if(tmp) return;
        tmp=true;
        for(i=2;i<n;i++)//前三点不共线
        {
            if((dist((ply[0]-ply[1])*(ply[1]-ply[i])))>PR)
            {
                swap(ply[2],ply[i]); tmp=false; break;
            }
        }
        if(tmp) return ;
        tmp=true;
        for(i=3;i<n;i++)//前四点不共面
        {
            if(fabs((ply[0]-ply[1])*(ply[1]-ply[2])^(ply[0]-ply[i]))>PR)
            {
                swap(ply[3],ply[i]); tmp=false; break;
            }
        }
        if(tmp) return ;
        face add;
        for(i=0;i<4;i++)//构建初始四面体
        {
            add.a=(i+1)%4,add.b=(i+2)%4,add.c=(i+3)%4,add.ok=1;
            if((ptoplane(ply[i],add))>0) swap(add.b,add.c);
            vis[add.a][add.b]=vis[add.b][add.c]=vis[add.c][add.a]=trianglecnt;
            tri[trianglecnt++]=add;
        }
        for(i=4;i<n;i++)//构建更新凸包
        {
            for(j=0;j<trianglecnt;j++)
            {
                if(tri[j].ok&&(ptoplane(ply[i],tri[j]))>PR)
                {
                    dfs(i,j); break;
                }
            }
        }
        int cnt=trianglecnt;
        trianglecnt=0;
        for(i=0;i<cnt;i++)
        {
            if(tri[i].ok)
                tri[trianglecnt++]=tri[i];
        }
    }
    double area()//表面积
    {
        double ret=0;
        for(int i=0;i<trianglecnt;i++)
            ret+=area(ply[tri[i].a],ply[tri[i].b],ply[tri[i].c]);
        return ret/2;
    }
    double volume()//体积
    {
        double res=0;
        TPoint tmp(0,0,0);
        for(int i=0;i<trianglecnt;i++)
            res+=volume(tmp,ply[tri[i].a],ply[tri[i].b],ply[tri[i].c]);
        return fabs(res/6);
    }
    int facetri() {return trianglecnt;}//表面三角形数
    int facepolygon()//表面多边形数
    {
        int ans=0,i,j,k;
        for(i=0;i<trianglecnt;i++)
        {
            for(j=0,k=1;j<i;j++)
            {
                if(same(i,j)) {k=0;break;}
            }
            ans+=k;
        }
        return ans;
    }
    double getMinDis(TPoint Q)//the minimum distance from the point to the hull
    {
        double ans = 0x7fffffff;
        for (int i = 0; i < trianglecnt; i++)
        {
            ans = min(ans,
                    volume(Q, ply[tri[i].a], ply[tri[i].b], ply[tri[i].c])
                            / area(ply[tri[i].a], ply[tri[i].b], ply[tri[i].c]));
        }
        return ans;
    }
    TPoint barycenter()//三维凸包重心
    {
        TPoint ans(0,0,0),o(0,0,0);
        double all=0;
        for(int i=0; i<trianglecnt; i++)
        {
            double vol=volume(o,ply[tri[i].a],ply[tri[i].b],ply[tri[i].c]);
            ans=ans+(o+ply[tri[i].a]+ply[tri[i].b]+ply[tri[i].c])/4.0*vol;
            all+=vol;
        }
        ans=ans/all;
        return ans;
    }
    double ptoface(TPoint p,int i)//点到面的距离
    {
        return fabs(volume(ply[tri[i].a],ply[tri[i].b],ply[tri[i].c],p)/dist((ply[tri[i].b]-ply[tri[i].a])*(ply[tri[i].c]-ply[tri[i].a])));
    }
}hull;
int main()
{
    // freopen("c.in", "r", stdin);
    // freopen("c-antb.out", "w", stdout);
    while(~scanf("%d",&hull.n))
    {
        if (hull.n==0)
        {
            break;
        }
        int i;
        for(i=0;i<hull.n;i++)
            scanf("%lf%lf%lf",&hull.ply[i].x,&hull.ply[i].y,&hull.ply[i].z);
        hull.construct();
        //printf("%d\n",hull.facepolygon());
        //printf("%.3lf\n", hull.area());
        /*int q;
        scanf("%d",&q);
        for (int i = 0; i < q; ++i)
        {
            TPoint Q;
            scanf("%lf%lf%lf", &Q.x, &Q.y, &Q.z);
            printf("%.4lf\n", hull.getMinDis(Q));
        }*/
        TPoint p=hull.barycenter();
        double ans=1e8;
        for (i=0;i<hull.trianglecnt;i++)
        {
            ans=min(ans,hull.ptoface(p,i));
        }
        printf("%.3f\n",ans);
    }
    return 0;
}