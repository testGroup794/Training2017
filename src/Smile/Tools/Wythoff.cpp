#include <bits/stdc++.h>
using namespace std;

int n,m,t,x,y,z;

bool WythoffGame(int x, int y) // 威佐夫博弈，返回是否先手必胜
{
	double a = (1.0 + sqrt(5.0)) / 2.0;  
	double b = (3.0 + sqrt(5.0)) / 2.0;  
	if (x>y) swap(x,y);
	int k=ceil(y/b);
	return x!=int(a*k) || y!=int(b*k);
}

int main()
{
	while (cin>>x>>y) cout<<(WythoffGame(x,y)?1:0)<<endl;
	return 0;
}

