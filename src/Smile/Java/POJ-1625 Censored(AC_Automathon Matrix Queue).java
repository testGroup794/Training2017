import java.util.*;
import java.io.*;
import java.math.*;

public class Main {
	static int		N	= 57;
	static int		LI	= 12;
	static int[]	mp	= new int[300];

	public static class Trie {
		public static int		CN, N, LI;
		public static int[][]	nxt;
		public static int[]		fail;
		public static boolean[]	ed;
		public static int		rt, l;

		Trie() {}
		Trie(int cn, int n, int li) {
			CN = cn;
			N = n;
			LI = li;
			nxt = new int[N * LI][CN];
			fail = new int[N * LI];
			ed = new boolean[N * LI];
		}

		int getBuf(byte c) {
			int t = c;
			return mp[t + 128];
		}
		int newNode() {
			for (int i = 0; i < CN; i++) nxt[l][i] = -1;
			ed[l++] = false;
			return l - 1;
		}
		public void init() {
			l = 0;
			rt = newNode();
		}

		public void insert(byte buf[]) {
			int len = buf.length;
			int now = rt;
			for (int i = 0; i < len; i++) {
				if (nxt[now][getBuf(buf[i])] == -1) {
					nxt[now][getBuf(buf[i])] = newNode();
				}
				now = nxt[now][getBuf(buf[i])];
			}
			ed[now] = true;
		}

		public void build() {
			Queue<Integer> q = new LinkedList<Integer>();
			fail[rt] = rt;
			for (int i = 0; i < CN; i++) {
				if (nxt[rt][i] == -1) {
					nxt[rt][i] = rt;
				} else {
					fail[nxt[rt][i]] = rt;
					q.add(nxt[rt][i]);
				}
			}
			while (!q.isEmpty()) {
				int now = q.remove();
				if (ed[fail[now]]) {
					ed[now] = true;
				}
				for (int i = 0; i < CN; i++) {
					if (nxt[now][i] == -1) {
						nxt[now][i] = nxt[fail[now]][i];
					} else {
						fail[nxt[now][i]] = nxt[fail[now]][i];
						q.add(nxt[now][i]);
					}
				}
			}
		}

		int[][] getMatrix() {
			int[][] dp = new int[l][l];
			for (int i = 0; i < l; i++) {
				for (int j = 0; j < CN; j++) {
					if (!ed[nxt[i][j]]) {
						dp[i][nxt[i][j]]++;
					}
				}
			}
			return dp;
		}
	}

	public static void main(String[] args) throws Exception {
		int n, m, p;
		int[][] gm;
		BigInteger ans;
		BigInteger[][] dp;
		Scanner cin = new Scanner(System.in);
		while (cin.hasNext()) {
			n = cin.nextInt();
			m = cin.nextInt();
			p = cin.nextInt();
			Trie ac = new Trie(n, N, LI);
			Arrays.fill(mp, 0);
			byte[] s = cin.next().getBytes();
			for (int i = 0; i < n; i++) {
				int t = s[i];
				mp[t + 128] = i;
			}
			ac.init();
			for (int i = 0; i < p; i++) {
				byte[] tmp = cin.next().getBytes();
				ac.insert(tmp);
			}
			ac.build();
			gm = ac.getMatrix();
			dp = new BigInteger[N][ac.l];
			for (int i = 0; i < dp.length; i++) {
				for (int j = 0; j < dp[i].length; j++) {
					dp[i][j] = BigInteger.ZERO;
				}
			}
			// System.out.println(ac.l);
			dp[0][0] = BigInteger.ONE;
			for (int i = 1; i <= m; i++) {
				for (int j = 0; j < ac.l; j++) {
					for (int k = 0; k < ac.l; k++) {
						if (gm[j][k] > 0) {
							dp[i][k] = dp[i][k].add(dp[i - 1][j].multiply(BigInteger.valueOf(gm[j][k])));
						}
					}
				}
			}
			ans = BigInteger.ZERO;
			for (int i = 0; i < dp[m].length; i++) {
				ans = ans.add(dp[m][i]);
			}
			System.out.println(ans.toString());
		}
	}
}