import java.util.*;
import java.math.*;
import java.io.*;

public class Main {
    static PrintWriter cout = new PrintWriter(System.out, true);
    static Scanner cin = new Scanner(System.in);

    public static class Point {
        BigDecimal x, y;
        BigDecimal distence2(Point p) {
            return x.subtract(p.x).multiply(x.subtract(p.x)).add(y.subtract(p.y).multiply(y.subtract(p.y)));
        }
    }

    public static Point easyCircumCenter(Point a, Point b, Point c) { // 外心，外接圆圆心，三条中垂线交点
        Point res = new Point();
        BigDecimal a1 = b.x.subtract(a.x), b1 = b.y.subtract(a.y), c1 = (a1.multiply(a1).add(b1.multiply(b1))).divide(BigDecimal.valueOf(2.0));
        BigDecimal a2 = c.x.subtract(a.x), b2 = c.y.subtract(a.y), c2 = (a2.multiply(a2).add(b2.multiply(b2))).divide(BigDecimal.valueOf(2.0));
        BigDecimal d = a1.multiply(b2).subtract(a2.multiply(b1));
        res.x = a.x.add((c1.multiply(b2).subtract(c2.multiply(b1))).divide(d, 30, BigDecimal.ROUND_HALF_UP));
        res.y = a.y.add((a1.multiply(c2).subtract(a2.multiply(c1))).divide(d, 30, BigDecimal.ROUND_HALF_UP));
        return res;
    }

    public static void main(String[] args) {

        int T;
        T = cin.nextInt();
        for (int ca = 1; ca <= T; ca++) {
            Point p1 = new Point();
            Point p2 = new Point();
            Point p3 = new Point();
            p1.x = cin.nextBigDecimal();
            p1.y = cin.nextBigDecimal();
            p2.x = cin.nextBigDecimal();
            p2.y = cin.nextBigDecimal();
            p3.x = cin.nextBigDecimal();
            p3.y = cin.nextBigDecimal();
            Point cir = easyCircumCenter(p1, p2, p3);
            BigDecimal r = cir.distence2(p1);
            Point p0 = new Point();
            p0.x = cin.nextBigDecimal();
            p0.y = cin.nextBigDecimal();
            BigDecimal dis = cir.distence2(p0);
            if (dis.compareTo(r) > 0) cout.println("Accepted");
            else cout.println("Rejected");
        }
    }
}
