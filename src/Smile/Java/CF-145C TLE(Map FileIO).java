import java.util.*;
import java.io.*;
import java.math.*;
public class Main {
	public static void main(String[] args) throws Exception {
		
		FileInputStream fis=new FileInputStream("h.in");  
        System.setIn(fis);
		PrintStream ps=new PrintStream(new FileOutputStream("h.out"));
		System.setOut(ps);

        Scanner cin = new Scanner(System.in);

        int MAXN=100007;
        int MAXM=1024;
        int MOD=1000000007;
        int [] luckyNumber=new int[MAXM];
        int i,j,k,n,m,t,x,y,tn,tk;
        long ans;
        int [] a=new int[MAXN];
        long [] cnm=new long[MAXN];
        long [] luckyCnt=new long[MAXM];
        long [][] dp=new long[MAXM][MAXM];
        Queue <Integer> qa=new LinkedList <Integer>();
        Queue <Integer> qb=new LinkedList <Integer>();
        Map <Integer,Integer> luckyMap= new HashMap<Integer,Integer>();
        BigInteger tmp=BigInteger.valueOf(1);

        int cnt=0;
        qa.add(4);
        luckyNumber[cnt++]=4;
        luckyMap.put(4,cnt-1);
        qa.add(7);
        luckyNumber[cnt++]=7;
        luckyMap.put(7,cnt-1);
        for (i=1;i<9;i++) {
        	while (qa.peek()!=null) {
        		x=qa.poll();
        		qb.add(x*10+4);
        		luckyNumber[cnt++]=x*10+4;
        		luckyMap.put(x*10+4,cnt-1);
        		qb.add(x*10+7);
        		luckyNumber[cnt++]=x*10+7;
        		luckyMap.put(x*10+7,cnt-1);
        	}
        	while (qb.peek()!=null) {
        		qa.add(qb.poll());
        	}
        }
        
        while (cin.hasNext()) {
        	for (i=0;i<MAXM;i++)
        	{
        		luckyCnt[i]=0;
        	}
        	for (i=0;i<MAXM;i++)
        	{
        		for (j=0;j<MAXM;j++)
        		{
        			dp[i][j]=0;
        		}
        	}
        	n=cin.nextInt(); 
			k=cin.nextInt(); 
			tn=n;
			for (i=0; i<n; i++) {
				a[i]=cin.nextInt();
				if (luckyMap.containsKey(a[i]))
				{
					luckyCnt[luckyMap.get(a[i])]++;
					tn--;
				}
			}

			cnm[0]=1;
			cnm[1]=tn;
			tmp=BigInteger.valueOf(tn);
			for (i=2; i<=(tn+1)/2; i++) {
				tmp=tmp.multiply(BigInteger.valueOf(tn-i+1));
				tmp=tmp.divide(BigInteger.valueOf(i));
				cnm[i]=tmp.mod(BigInteger.valueOf(MOD)).intValue();
			}
			for (i=(tn+1)/2+1; i<=tn; i++)
			{
				cnm[i]=cnm[tn-i];
			}
			for (i=tn+1; i<MAXN; i++)
			{
				cnm[i]=0;
			}

			dp[0][0]=1;
			if (luckyCnt[0]>0)
			{
				dp[0][1]=luckyCnt[0];
			}
			for (i=1;i<cnt;i++)
			{
				if (luckyCnt[i]>0)
				{
					dp[i][0]=1;
					for (j=1;j<=i+1;j++)
					{
						dp[i][j]=dp[i-1][j-1]*luckyCnt[i]+dp[i-1][j];
						dp[i][j]%=MOD;
					}
				}
				else
				{
					for (j=0;j<=i+1;j++)
					{
						dp[i][j]=dp[i-1][j];
					}
				}
			}

			// for (j=0;j<10;j++)
			// {
			// 	System.out.print(""+a[j]+" ");
			// }
			// System.out.println();

			// for (i=0;i<10;i++)
			// {
			// 	for (j=0;j<10;j++)
			// 	{
			// 		System.out.print(""+dp[i][j]+" ");
			// 	}
			// 	System.out.println();
			// }

			ans=0;
			for (i=0;i<=cnt;i++)
			{
				if (k>=i)
				{
					ans+=dp[cnt-1][i]*cnm[k-i];
					ans%=MOD;
				}
			}
			System.out.println(ans);
        }
    }
}