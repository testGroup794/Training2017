#include <cstdio>
#include <iostream>
#include <algorithm>

#define REP(I,A,B) for (int I=(A),I##_END_=(B);I<=I##_END_;I++)
#define RI(X) scanf("%d",&X)

using namespace std;

void open()
{
	freopen("CF61E.in","r",stdin);
	freopen("CF61E.out","w",stdout);
}
void close()
{
	fclose(stdin);
	fclose(stdout);
}

const int MAXN = 1010101;

int t[MAXN<<1];
long long s[MAXN];

int n;

struct num{
	int val;
	int id;
	int rk;
}a[MAXN];

bool cmpvl(num a,num b){
	return a.val<b.val;
}
bool cmpid(num a,num b){
	return a.id<b.id;
}

void init(){
	RI(n);
	REP(i,1,n)
		RI(a[i].val);
	REP(i,1,n)
		a[i].id=i;
	sort(a+1,a+1+n,cmpvl);
	REP(i,1,n)
		a[i].rk=i;
	sort(a+1,a+1+n,cmpid);
}

void add(int l,int r,int id){
	t[ l+r | (l!=r) ]++;
	if (l==r) return ;
	int mid=(l+r)/2;
	if (id<=mid) add(l,mid,id);
	else add(mid+1,r,id);
}
void sub(int l,int r,int id){
	t[ l+r | (l!=r) ]--;
	if (l==r) return ;
	int mid=(l+r)/2;
	if (id<=mid) sub(l,mid,id);
	else sub(mid+1,r,id);
}

int get(int l,int r,int ll,int rr){
	if (l==ll && r==rr)
		return t[ l+r | (l!=r)];
	int mid=(l+r)/2;
	if (rr<=mid)
		return get(l,mid,ll,rr);
	else if (mid<ll)
		return get(mid+1,r,ll,rr);
	else
		return get(l,mid,ll,mid)+get(mid+1,r,mid+1,rr);
}

void work(){
	REP(i,1,n)
	{
		if (a[i].rk<n)
			s[i]=get(1,n,a[i].rk+1,n);
		else
			s[i]=0;
		add(1,n,a[i].rk);
	}
	REP(i,1,n)
	{
		if (a[i].rk>1)
			s[i]*=get(1,n,1,a[i].rk-1);
		else
			s[i]=0;
		sub(1,n,a[i].rk);
	}
	long long ans=0;
	REP(i,1,n)
		ans+=s[i];
	cout << ans << endl;
}

int main()
{
	open();
	init();
	work();
	close();
	return 0;
}
