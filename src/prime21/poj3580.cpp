#include <stdio.h>
#include <algorithm>

using namespace std;

const int N= 101010;
const int M= 101010;
const int inf= 0x3f3f3f3f;

int tot;

char opt[20];

struct splay_node{
	int fa;
	int c[2];
	bool rev;
	int mn;
	int val;
	int add;
	int sz;
}t[N+M];
int root;

inline void setc(const int& p,const int &x,const int &kd){
	if (p) t[p].c[kd]=x;
	if (x) t[x].fa=p;
}
inline void maintain(const int &p){
	t[p].mn=t[p].val;
	if (t[p].c[0]) t[p].mn=min(t[p].mn,t[t[p].c[0]].mn);
	if (t[p].c[1]) t[p].mn=min(t[p].mn,t[t[p].c[1]].mn);
	if (p) t[p].sz=1+t[t[p].c[0]].sz+t[t[p].c[1]].sz;
}
inline void pushdown(const int &p){
	if (t[p].rev)
	{
		swap(t[p].c[0],t[p].c[1]);
		if (t[p].c[0]) t[t[p].c[0]].rev^=1;
		if (t[p].c[1]) t[t[p].c[1]].rev^=1;
		t[p].rev=false;
	}
	if (t[p].add)
	{
		if (t[p].c[0]) t[t[p].c[0]].mn+=t[p].add,t[t[p].c[0]].add+=t[p].add,t[t[p].c[0]].val+=t[p].add;
		if (t[p].c[1]) t[t[p].c[1]].mn+=t[p].add,t[t[p].c[1]].add+=t[p].add,t[t[p].c[1]].val+=t[p].add;
		t[p].add=0;
	}
}

void Insert(int p,int kd){
	if (t[p].c[kd])
		Insert(t[p].c[kd],kd);
	else
	{
		tot++;
		t[tot].mn=0x3f3f3f3f;
		t[tot].add=t[tot].c[0]=t[tot].c[1]=t[tot].fa=t[tot].rev=0;
		t[tot].sz=1;
		setc(p,tot,kd);
	}
	maintain(p);
}

int n;
int a[N];
int m;

void dig(int p){
	if (p==0) return;
	pushdown(p);
	dig(t[p].c[0]);
	fprintf(stderr,"%d\n",t[p].val);
	dig(t[p].c[1]);
}

int prepare(int l,int r){
	if (l>r) return 0;
	if (l==r){
		tot++;
		t[tot].val=t[tot].mn=a[l];
		t[tot].add=t[tot].c[0]=t[tot].c[1]=t[tot].fa=t[tot].rev=0;
		t[tot].sz=1;
		return tot;
	}
	int mid =(l+r)/2;

	int ret=++tot;
	t[ret].val=t[ret].mn=a[mid];
	t[ret].add=t[ret].c[0]=t[ret].c[1]=t[ret].fa=t[ret].rev=0;
	setc(ret,prepare(l,mid-1),0);
	setc(ret,prepare(mid+1,r),1);
	maintain(ret);

	return ret;
}

void rotate(int x){
	int p=t[x].fa;
	if (p==0) return ;
	int kd=t[p].c[1]==x;

	pushdown(p);
	pushdown(x);

	setc(t[p].fa,x,t[t[p].fa].c[1]==p);
	setc(p,t[x].c[kd^1],kd);
	setc(x,p,kd^1);

	maintain(p);
	maintain(x);
}
void splay(int x,int p=0){
	while (t[x].fa!=p)
	{
		if (t[t[x].fa].fa==p)
			rotate(x);
		else
		{
			if ((t[t[t[x].fa].fa].c[1]==t[x].fa)==(t[t[x].fa].c[1]==x))
				rotate(t[x].fa);
			else
				rotate(x);
			rotate(x);
		}
	}

	if (p==0) root=x;
}
int get_rk(int p,int k){
	pushdown(p);
	if (t[t[p].c[0]].sz>=k) return get_rk(t[p].c[0],k);
	else if (t[t[p].c[0]].sz+1==k) return p;
	else return get_rk(t[p].c[1],k-t[t[p].c[0]].sz-1);
}
int get_seg(int l,int r){
	int t1=get_rk(root,l-1);
	splay(t1);
	int t2=get_rk(root,r+1);
	splay(t2,t1);
	return t[t2].c[0];
}
void up(int p){
	while (t[p].fa)
	{
		p=t[p].fa;
		maintain(p);
	}
}
void cut(int x){
	int p=t[x].fa;
	int kd=t[p].c[1]==x;
	t[p].c[kd]=0;
	t[x].fa=0;
	maintain(p);
	up(p);
}
void add(int l,int r,int d){
	int T=get_seg(l,r);
	t[T].add+=d;
	t[T].mn+=d;
	t[T].val+=d;
	up(T);
}
void ins(int x,int d){
	int T=get_seg(x,x);

	++tot;
	t[tot].mn=t[tot].val=d;
	t[tot].add=t[tot].c[0]=t[tot].c[1]=t[tot].fa=t[tot].rev=0;
	t[tot].sz=1;
	pushdown(T);
	setc(T,tot,1);
	maintain(T);
	up(T);
}
void del(int x){
	int T=get_seg(x,x);
	cut(T);
}
int get(int x,int y){
	int T=get_seg(x,y);
	return t[T].mn;
}
void rev(int x,int y){
	int T=get_seg(x,y);
	t[T].rev^=1;
}
void fly(int x,int y,int d){
	if (t<=0) return ;
	int T=get_seg(y-d+1,y);
	cut(T);
	int T1=get_seg(x-1,x-1);
	pushdown(T1);
	setc(T1,T,1);
	maintain(T1);
	up(T1);
}

void work(){
	scanf("%d",&m);
	int x,y,d;
	for (int i=1;i<=m;i++)
	{
		scanf("%s",opt+1);
		if (opt[1]=='A') //add
		{
			scanf("%d%d%d",&x,&y,&d);
			add(x+2,y+2,d);
		}
		else if (opt[1]=='I') //insert
		{
			scanf("%d%d",&x,&d);
			ins(x+2,d);
		}
		else if (opt[1]=='D') //delete
		{
			scanf("%d",&x);
			del(x+2);
		}
		else if (opt[1]=='M') //min
		{
			scanf("%d%d",&x,&y);
			printf("%d\n",get(x+2,y+2));
		}
		else if (opt[1]=='R' && opt[5]=='R') //reverse
		{
			scanf("%d%d",&x,&y);
			rev(x+2,y+2);
		}
		else // revolve
		{
			scanf("%d%d%d",&x,&y,&d);
			fly(x+2,y+2,d%(y-x+1));
		}
		//dig(root);
	}
}

int main(){
	freopen("3580.in","r",stdin);
	freopen("3580.out","w",stdout);

	scanf("%d",&n);
	for (int i=1;i<=n;i++) scanf("%d",a+i);
	root=prepare(1,n);
	
	// avoid the discussion -general way
	Insert(root,0);
	Insert(root,0);
	Insert(root,1);
	Insert(root,1);

	work();

	fclose(stdin);
	fclose(stdout);
	return 0;
}
