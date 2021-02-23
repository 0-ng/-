
### [0.	基础](#0)
### [1.	几何](#1)
### [2.	直线与凸多边形交点](#2)
### [3.	凸包Andrew](#3)
### [4.	凸包Graham](#4)
### [5.	判断线段相交](#5)
### [6.	半平面交](#6)
### [7.	四面体内切球](#7)
### [8.	平面最近点对](#8)
### [9.	旋转卡壳求最大三角形面积](#9)
### [10. 最小圆覆盖](#10)
### [11. 模拟退火求矩阵内所有点最小距离最大](#11)
### [12. 矩阵面积凸包最小矩形覆盖](#12)
### [13. 稳定凸包](#13)
### [14. simpson积分](#14)
### [15. 三维几何凸包](#15)
### [16. 平面最近点对](#16)
### [17. 旋转卡壳求最大四边形面积](#17)
### [18. 可加点的动态凸包](#18)
### [19. KD-Tree](#19)


<span id="0"><h4>0. 基础</h4></span>
1. 判断是否输出负0
2. 判断生成凸包时有无等于
acos： [0,pi] radians.
asin： [-pi/2,+pi/2] radians.
atan： [-pi/2,+pi/2] radians.
平面问题可以先算左下角的子问题，将平面或点旋转后再算一次以此类推


<span id="1"><h4>1. 几何</h4></span>
```cpp
const double EPS=1e-8;
const double PI=acos(-1);
inline int sgn(double a){ return a < -EPS ? -1 : a > EPS; }
inline int cmp(double a, double b){ return sgn(a-b); }
struct Point;
struct Line;
typedef Point Vector;

struct Point{
    double x,y;
    Point(){}
    Point(double a, double b):x(a),y(b){}
    double len(){return sqrt(x*x+y*y);}
    double len2(){return x*x+y*y;}
    double disToPoint(Point p){return sqrt((x-p.x)*(x-p.x)+(y-p.y)*(y-p.y));}
    Point rotLeft(){
        //`逆时针旋转90度`
        return Point(-y,x);
    }
    Point rotRight(){
        //`顺时针旋转90度`
        return Point(y,-x);
    }
    Point trunc(double r){
        //`化为长度为r的向量`
        double l = len();
        if(!sgn(l))return *this;
        r /= l;
        return Point(x*r,y*r);
    }
    double rad(Point a,Point b){//返回弧度角
        Point p = *this;
        return fabs(atan2( fabs((a-p)^(b-p)),(a-p)*(b-p) ));
        //acos(((a-p)*(b-p))/(a-p).len()/(b-p).len()));
        //`计算pa  和  pb 的夹角`
        //`就是求这个点看a,b 所成的夹角`
        //`测试 LightOJ1203`
    }
    double rad(Vector b){//返回弧度角
        return fabs(atan2( fabs((*this)^(b)),(*this)*(b) ));
        //`计算*this和  b 的夹角`
    }
    Point rotate(Point p,double angle){
        //`绕着p点逆时针旋转angle`
        Point v = (*this) - p;
        double c = cos(angle), s = sin(angle);
        return Point(p.x + v.x*c - v.y*s,p.y + v.x*s + v.y*c);
    }
    Vector rotate(double angle){
        //~向量逆时针旋转angle`
        double c = cos(angle), s = sin(angle);
        return Vector(x*c-y*s, x*s+y*c);
    }
    void read(){scanf("%lf%lf",&x,&y);}
    Point operator+(const Vector v)const{return {x+v.x,y+v.y};}
    Vector operator-(const Point p)const{return {x-p.x,y-p.y};}
    double operator^(const Vector v)const{return x*v.y-y*v.x;}//叉乘
    double operator*(const Vector v)const{return x*v.x+y*v.y;}//点乘
    Vector operator*(const double d)const{return {x*d,y*d};}
    Vector operator/(const double d)const{return {x/d,y/d};}
    bool operator==(const Point p)const{return cmp(x,p.x)==0&&cmp(y,p.y)==0;}
    bool operator<(const Point p)const{if(cmp(x,p.x)==0) return y<p.y;return x<p.x;}
};
struct Line{
    Line(){}
    Line(Point a,Point b){s=a;e=b;}
    Line(Point p,double angle){
        //`根据一个点和倾斜角angle确定直线,0<=angle<pi`
        s = p;
        if(sgn(angle-PI/2) == 0)
            e = (s + Point(0,1));
        else
            e = (s + Point(1,tan(angle)));
    }
    Line(double a,double b,double c) {
        //ax+by+c=0
        if (sgn(a) == 0) {
            s = Point(0, -c / b);
            e = Point(1, -c / b);
        } else if (sgn(b) == 0) {
            s = Point(-c / a, 0);
            e = Point(-c / a, 1);
        } else {
            s = Point(0, -c / b);
            e = Point(1, (-c - a) / b);
        }
    }
    double getAngle(){
        //`返回直线倾斜角 0<=angle<pi`
        double k = atan2(e.y-s.y,e.x-s.x);
        if(sgn(k) < 0)k += PI;
        if(sgn(k-PI) == 0)k -= PI;
        return k;
    }
    double len(){return s.disToPoint(e);}
    int relation(Point p){
        //`点和直线关系`
        int c = sgn((p-s)^(e-s));
        if(c < 0)return 1;
        else if(c > 0)return 2;
        else return 3;
        //`1 在左侧`//`2 在右侧`//`3 在直线上`
    }
    bool pointOnSeg(Point p){
        // 点在线段上的判断
        return sgn((p-s)^(e-s)) == 0 && sgn((p-s)*(p-e)) <= 0;
    }
    bool parallel(Line v){//checked
        //`两向量平行(对应直线平行或重合)`
        return sgn((e-s)^(v.e-v.s)) == 0;
    }
    Point point(double t){return s + (e - s)*t;}//返回点P = v + (p - v)*t
    int segCrossSeg(Line v){
        //`两线段相交判断`
        int d1 = sgn((e-s)^(v.s-s));
        int d2 = sgn((e-s)^(v.e-s));
        int d3 = sgn((v.e-v.s)^(s-v.s));
        int d4 = sgn((v.e-v.s)^(e-v.s));
        if( (d1^d2)==-2 && (d3^d4)==-2 )return 2;
        return (d1==0 && sgn((v.s-s)*(v.s-e))<=0) ||
               (d2==0 && sgn((v.e-s)*(v.e-e))<=0) ||
               (d3==0 && sgn((s-v.s)*(s-v.e))<=0) ||
               (d4==0 && sgn((e-v.s)*(e-v.e))<=0);
        //`2 规范相交`//`1 非规范相交`//`0 不相交`
    }
    int lineCrossSeg(Line v){
        //`直线和线段相交判断`
        //`-*this line   -v seg`
        int d1 = sgn((e-s)^(v.s-s));
        int d2 = sgn((e-s)^(v.e-s));
        if((d1^d2)==-2) return 2;
        return (d1==0||d2==0);
        //`2 规范相交`//`1 非规范相交`//`0 不相交`
    }
    int lineCrossLine(Line v){
        //`两直线关系`
        if((*this).parallel(v))
            return v.relation(s)==3;
        return 2;
        //`0 平行`//`1 重合`//`2 相交`
    }
    Point crossPoint(Line v){//checked
        //`求两直线的交点`
        //`要保证两直线不平行或重合`
        double a1 = (v.e-v.s)^(s-v.s);
        double a2 = (v.e-v.s)^(e-v.s);
        return Point((s.x*a2-e.x*a1)/(a2-a1),(s.y*a2-e.y*a1)/(a2-a1));
    }
    double disPointToLine(Point p){//checked
        //点到直线的距离
        return fabs((p-s)^(e-s))/len();
    }
    double disPointToSeg(Point p){
        //点到线段的距离
        if(sgn((p-s)*(e-s))<0 || sgn((p-e)*(s-e))<0)
            return min(p.disToPoint(s),p.disToPoint(e));
        return disPointToLine(p);
    }
    double disSegToSeg(Line v){
        //`返回线段到线段的距离` //`前提是两线段不相交，相交距离就是0了`
        return min(min(disPointToSeg(v.s),disPointToSeg(v.e)),
                    min(v.disPointToSeg(s),v.disPointToSeg(e)));
    }
    Point lineProg(Point p){//checked
        //`返回点p在直线上的投影`
        return s + ( ((e-s)*((e-s)*(p-s)))/((e-s).len2()) );
    }
    Point symmetryPoint(Point p){
        //`返回点p关于直线的对称点`
        Point q = lineProg(p);
        return Point(2*q.x-p.x,2*q.y-p.y);
    }
    Point s,e;
    double ang;
};
struct Circle{
    Circle(){}
    Circle(Point _p, double _r):p(_p), r(_r) {}
    Circle(double x,double y,double _r){
        p=Point(x,y);
        r=_r;
    }
    Circle(Point a,Point b,Point c){
        //`三角形的外接圆`
        Line u = Line((a+b)/2,((a+b)/2)+((b-a).rotLeft()));
        Line v = Line((b+c)/2,((b+c)/2)+((c-b).rotLeft()));
        p = u.crossPoint(v);
        r = p.disToPoint(a);
        //`需要Point的+ /  rotate()  以及Line的crosspoint()`
        //`利用两条边的中垂线得到圆心`
        //`测试：UVA12304`
    }
    Circle(Point a,Point b,Point c,bool t){
//        Line u,v;
//        double m = atan2(b.y-a.y,b.x-a.x), n = atan2(c.y-a.y,c.x-a.x);
//        u.s = a;
//        u.e = u.s + Point(cos((n+m)/2),sin((n+m)/2));
//        v.s = b;
//        m = atan2(a.y-b.y,a.x-b.x) , n = atan2(c.y-b.y,c.x-b.x);
//        v.e = v.s + Point(cos((n+m)/2),sin((n+m)/2));
//        p = u&v;
//        r = p.disToSegment(Line(a,b));
        //`三角形的内切圆`
        //`参数bool t没有作用，只是为了和上面外接圆函数区别`
        //`测试：UVA12304`
        if(((b-a)^(c-a))<0)
            swap(b,c);
        r=((a^b)+(b^c)+(c^a))/((b-a).len()+(c-b).len()+(a-c).len());
        Vector v1=((b-a).rotLeft()).trunc(r);
        Vector v2=((c-b).rotLeft()).trunc(r);
        Line l1={a+v1,b+v1};Line l2={b+v2,c+v2};
        p=l1.crossPoint(l2);
    }
    double area(){return PI*r*r;}
    double circumference(){return 2*PI*r;}
    int relation(Point b){
        //`点和圆的关系`
        double dis = b.disToPoint(p);
        if(sgn(dis-r) < 0)return 2;
        else if(sgn(dis-r)==0)return 1;
        return 0;
        //`0 圆外`//`1 圆上`//`2 圆内`
    }
    int relationSeg(Line v){
        //`线段和圆的关系`
        double dis = v.disPointToSeg(p);
        if(sgn(dis-r) < 0)return 2;
        else if(sgn(dis-r) == 0)return 1;
        return 0;
        //`比较的是圆心到线段的距离和半径的关系`
        //`0 圆外`//`1 圆上`//`2 圆内`
    }
    int relationLine(Line v){//checked
        //`直线和圆的关系`
        double dis = v.disPointToLine(p);
        if(sgn(dis-r) < 0)return 2;
        else if(sgn(dis-r) == 0)return 1;
        return 0;
        //`比较的是圆心到直线的距离和半径的关系`
        //`0 圆外`//`1 圆上`//`2 圆内`
    }
    int relationCircle(Circle v){//`两圆的关系`//checked
        double d = p.disToPoint(v.p);
        if(sgn(d-r-v.r) > 0)return 5;
        if(sgn(d-r-v.r) == 0)return 4;
        double l = fabs(r-v.r);
        if(sgn(d-r-v.r)<0 && sgn(d-l)>0)return 3;
        if(sgn(d-l)==0)return 2;
        if(sgn(d-l)<0)return 1;
        //`5 相离`//`4 外切`//`3 相交`//`2 内切`//`1 内含`
        //`需要Point的distance`
        //`测试：UVA12304`
    }
    int pointCrossCircle(Circle v,Point &p1,Point &p2){//checked
        //`求两个圆的交点，返回0表示没有交点，返回1是一个交点，2是两个交点`
        int rel = relationCircle(v);
        if(rel == 1 || rel == 5)return 0;
        double d = p.disToPoint(v.p);
        double l = (d*d+r*r-v.r*v.r)/(2*d);
        double h = sqrt(r*r-l*l);
        Point tmp = p + (v.p-p).trunc(l);
        p1 = tmp + ((v.p-p).rotLeft().trunc(h));
        p2 = tmp + ((v.p-p).rotRight().trunc(h));
        if(rel == 2 || rel == 4)
            return 1;
        return 2;
        //`需要relationcircle`
        //`测试：UVA12304`
    }
    int pointCrossLine(Line v,Point &p1,Point &p2){//checked
        //`求直线和圆的交点，返回交点个数`
        if(!(*this).relationLine(v))return 0;
        Point a = v.lineProg(p);
        double d = v.disPointToLine(p);
        d = sqrt(r*r-d*d);
        if(sgn(d) == 0){
            p1 = a;
            p2 = a;
            return 1;
        }
        p1 = a + (v.e-v.s).trunc(d);
        p2 = a - (v.e-v.s).trunc(d);
        return 2;
    }
    int getCircle(Point a,Point b,double r1,Circle &c1,Circle &c2){
        //`得到过a,b两点，半径为r1的两个圆`
        Circle x(a,r1),y(b,r1);
        int t = x.pointCrossCircle(y,c1.p,c2.p);
        if(!t)return 0;
        c1.r = c2.r = r;
        return t;
    }
    int getCircle(Line u,Point q,double r1,Circle &c1,Circle &c2){//checked
        //`得到与直线u相切，过点q,半径为r1的圆`
        double dis = u.disPointToLine(q);
        if(sgn(dis-r1*2)>0)return 0;
        if(sgn(dis) == 0){
            c1.p = q + ((u.e-u.s).rotLeft().trunc(r1));
            c2.p = q + ((u.e-u.s).rotRight().trunc(r1));
            c1.r = c2.r = r1;
            return 2;
        }
        Line u1 = Line((u.s + (u.e-u.s).rotLeft().trunc(r1)),(u.e + (u.e-u.s).rotLeft().trunc(r1)));
        Line u2 = Line((u.s + (u.e-u.s).rotRight().trunc(r1)),(u.e + (u.e-u.s).rotRight().trunc(r1)));
        Circle cc = Circle(q,r1);
        Point p1,p2;
        if(!cc.pointCrossLine(u1,p1,p2))cc.pointCrossLine(u2,p1,p2);
        c1 = Circle(p1,r1);
        if(p1 == p2){
            c2 = c1;
            return 1;
        }
        c2 = Circle(p2,r1);
        return 2;
        //`测试：UVA12304`
    }
    int getCircle(Line u,Line v,double r1,Circle &c1,Circle &c2,Circle &c3,Circle &c4){//checked
        //`同时与直线u,v相切，半径为r1的圆`
        if(u.parallel(v))return 0;//两直线平行
        Line u1 = Line(u.s + (u.e-u.s).rotLeft().trunc(r1),u.e + (u.e-u.s).rotLeft().trunc(r1));
        Line u2 = Line(u.s + (u.e-u.s).rotRight().trunc(r1),u.e + (u.e-u.s).rotRight().trunc(r1));
        Line v1 = Line(v.s + (v.e-v.s).rotLeft().trunc(r1),v.e + (v.e-v.s).rotLeft().trunc(r1));
        Line v2 = Line(v.s + (v.e-v.s).rotRight().trunc(r1),v.e + (v.e-v.s).rotRight().trunc(r1));
        c1.r = c2.r = c3.r = c4.r = r1;
        c1.p = u1.crossPoint(v1);
        c2.p = u1.crossPoint(v2);
        c3.p = u2.crossPoint(v1);
        c4.p = u2.crossPoint(v2);
        return 4;
        //`测试：UVA12304`
    }
    int getCircle(Circle cx,Circle cy,double r1,Circle &c1,Circle &c2){//checked
        //`同时与不相交圆cx,cy相切，半径为r1的圆`
        Circle x(cx.p,r1+cx.r),y(cy.p,r1+cy.r);
        int t = x.pointCrossCircle(y,c1.p,c2.p);
        if(!t)return 0;
        c1.r = c2.r = r1;
        return t;
        //`测试：UVA12304`
    }
    int tangentLine(Point q,Line &u,Line &v){
        //`过一点作圆的切线(先判断点和圆的关系)`
        int x = relation(q);
        if(x == 2)return 0;
        if(x == 1){
            u = Line(q,q + (q-p).rotLeft());
            v = u;
            return 1;
        }
        double d = p.disToPoint(q);
        double l = r*r/d;
        double h = sqrt(r*r-l*l);
        u = Line(q,p + ((q-p).trunc(l) + (q-p).rotLeft().trunc(h)));
        v = Line(q,p + ((q-p).trunc(l) + (q-p).rotRight().trunc(h)));
        return 2;
        //`测试：UVA12304`
    }
    double areaCircle(Circle v){
        //`求两圆相交的面积`
        int rel = relationCircle(v);
        if(rel >= 4)return 0.0;
        if(rel <= 2)return min(area(),v.area());
        double d = p.disToPoint(v.p);
        double hf = (r+v.r+d)/2.0;
        double ss = 2*sqrt(hf*(hf-r)*(hf-v.r)*(hf-d));
        double a1 = acos((r*r+d*d-v.r*v.r)/(2.0*r*d));
        a1 = a1*r*r;
        double a2 = acos((v.r*v.r+d*d-r*r)/(2.0*v.r*d));
        a2 = a2*v.r*v.r;
        return a1+a2-ss;
    }
    double areaTriangle(Point a,Point b) {
        //`求圆和三角形pab的相交面积`
        if (sgn((p - a) ^ (p - b)) == 0)return 0.0;
        Point q[5];
        int len = 0;
        q[len++] = a;
        Line l(a, b);
        Point p1, p2;
        if (pointCrossLine(l, q[1], q[2]) == 2) {
            if (sgn((a - q[1]) * (b - q[1])) < 0)q[len++] = q[1];
            if (sgn((a - q[2]) * (b - q[2])) < 0)q[len++] = q[2];
        }
        q[len++] = b;
        if (len == 4 && sgn((q[0] - q[1]) * (q[2] - q[1])) > 0)swap(q[1], q[2]);
        double res = 0;
        for (int i = 0; i < len - 1; i++) {
            if (relation(q[i]) == 0 || relation(q[i + 1]) == 0) {
                double arg = p.rad(q[i], q[i + 1]);
                res += r * r * arg / 2.0;
            } else {
                res += fabs((q[i] - p) ^ (q[i + 1] - p)) / 2.0;
            }
        }
        return res;
        //`测试：POJ3675 HDU3982 HDU2892`
    }
    double areaCircle(Point* poly,int n){
        //`多边形和圆交的面积`
        double ans = 0;
        for(int i = 1;i <= n;i++){
            int j = (i%n+1);
            if(sgn( (poly[j]-p)^(poly[i]-p) ) >= 0)
                ans += areaTriangle(poly[i],poly[j]);
            else ans -= areaTriangle(poly[i],poly[j]);
        }
        return fabs(ans);
        //`测试：POJ3675 HDU3982 HDU2892`
    }

    Point point(double d){ return Point(p.x + cos(d)*r, p.y + sin(d)*r);}//通过圆心角求坐标
    bool contain(Point p){return cmp(r,p.disToPoint(p))>=0;}

    void read(){scanf("%lf%lf%lf",&p.x,&p.y,&r);}
    Point p;
    double r;
};
struct Polygon {
    int n;
    Point p[MAXN];
    Line l[MAXN];
    bool isConvex(){
        //`判断是不是凸的`
        bool s[3];
        memset(s,false,sizeof(s));
        for(int i = 0;i < n;i++){
            int j = (i+1)%n;
            int k = (j+1)%n;
            s[sgn((p[j]-p[i])^(p[k]-p[i]))+1] = true;
            if(s[0] && s[2])return false;
        }
        return true;
    }
    int relationPoint(Point point){
        //判断点是否在多边形内，若点在多边形内返回1，在多边形外部返回0，在多边形上返回-1
        int wn = 0;
        for(int i = 1; i <= n; ++i){
            if(Line(p[i],p[i%+1]).pointOnSeg(point)) return -1;
            int k = sgn((p[i%n+1] - p[i])^(point - p[i]));
            int d1 = sgn(p[i].y - point.y);
            int d2 = sgn(p[i%n+1].y - point.y);
            if(k > 0 && d1 <= 0 && d2 > 0) wn++;
            if(k < 0 && d2 <= 0 && d1 > 0) wn--;
        }
        return (wn != 0);
    }
    Point getBaryCentre(){
        //`得到重心`
        Point ret(0,0);
        double area = 0;
        for(int i = 1;i < n-1;i++){
            double tmp = (p[i]-p[0])^(p[i+1]-p[0]);
            if(sgn(tmp) == 0)continue;
            area += tmp;
            ret.x += (p[0].x+p[i].x+p[i+1].x)/3*tmp;
            ret.y += (p[0].y+p[i].y+p[i+1].y)/3*tmp;
        }
        if(sgn(area)) ret = ret/area;
        return ret;
    }
    int relationcircle(Circle c){
        //`多边形和圆关系`
        int x = 2;
        if(relationPoint(c.p) != 1)return 0;//圆心不在内部
        for(int i = 0;i < n;i++){
            if(c.relationSeg(l[i])==2)return 0;
            if(c.relationSeg(l[i])==1)x = 1;
        }
        return x;
        //` 2 圆完全在多边形内`
        //` 1 圆在多边形里面，碰到了多边形边界`
        //` 0 其它`
    }
};
```
<span id="2"><h4>2.	直线与凸多边形交点</h4></span>
```cpp
vector<Point> convexCut(const vector<Point> &ps,Point q1,Point q2){
    //`直线切凸多边形`
    //`多边形是逆时针的，在q1q2的左侧`
    //`测试:HDU3982`
    vector<Point>qs;
    int n = ps.size();
    for(int i = 0;i < n;i++){
        Point p1 = ps[i], p2 = ps[(i+1)%n];
        int d1 = sgn((q2-q1)^(p1-q1)), d2 = sgn((q2-q1)^(p2-q1));
        if(d1 >= 0)
            qs.push_back(p1);
        if(d1 * d2 < 0)
            qs.push_back(Line(p1,p2).crossPoint(Line(q1,q2)));
    }
    return qs;
}
```
<span id="3"><h4>3.	凸包Andrew</h4></span>
```cpp
double Andrew(){
    sort(p+1,p+1+n);
    int len=0;
    for (int i=1;i<=n;i++){
        while (len>1&&sgn((stk[len]-stk[len-1])^(p[i]-stk[len-1]))<0) len--;
        stk[++len]=p[i];
    }
    int k=len;
    for (int i=n-1;i>=1;i--){
        while (len>k&&sgn((stk[len]-stk[len-1])^(p[i]-stk[len-1]))<0) len--;
        stk[++len]=p[i];
    }
    len--;
//    return len;
    double sum=0;
    for(int i=1;i<=len;i++){
        sum+=stk[i].disToPoint(stk[i%len+1]);
    }
    return sum;//凸包的边长,点数为2时特判
}
```
<span id="4"><h4>4.	凸包Graham</h4></span>
```cpp
Point base;
bool grahamCmp(Point& p1,Point& p2){
    if(sgn((p1-base)^(p2-base))==0)
        return (p1-base).len()<(p2-base).len();
    return ((p1-base)^(p2-base))>0;
}
Point stk[MAXN];
int graham(){
    for(int i=2;i<=n;i++){
        if(p[i]<p[1])
            swap(p[i],p[1]);
    }
    base=p[1];
    sort(p+2,p+1+n,grahamCmp);
    int k=0;
    stk[++k]=p[1];stk[++k]=p[2];
    for(int i=3;i<=n;i++){
        while(k>1&&sgn((stk[k]-stk[k-1])^(p[i]-stk[k-1]))<0)
            k--;
        stk[++k]=p[i];
    }
    return k;
}
```
<span id="5"><h4>5.	判断线段相交</h4></span>
```cpp
struct Line{
    Point s;
    Point e;
};

//线段的两端接触也算相交
bool inter(Line l1,Line l2){
    return
    max(l1.s.x,l1.e.x)>=min(l2.s.x,l2.e.x)&&
    max(l2.s.x,l2.e.x)>=min(l1.s.x,l1.e.x)&&
    max(l1.s.y,l1.e.y)>=min(l2.s.y,l2.e.y)&&
    max(l2.s.y,l2.e.y)>=min(l1.s.y,l1.e.y)&&
    sgn(Cross(l2.s-l1.s,l1.e-l1.s))*sgn(Cross(l1.e-l1.s,l2.e-l1.s))>=0&&
    sgn(Cross(l1.s-l2.s,l2.e-l2.s))*sgn(Cross(l2.e-l2.s,l1.e-l2.s))>=0;
}


//线段两端接触不算相交
bool inter(Line l1,Line l2){
    return
    max(l1.s.x,l1.e.x)>=min(l2.s.x,l2.e.x)&&
    max(l2.s.x,l2.e.x)>=min(l1.s.x,l1.e.x)&&
    max(l1.s.y,l1.e.y)>=min(l2.s.y,l2.e.y)&&
    max(l2.s.y,l2.e.y)>=min(l1.s.y,l1.e.y)&&
    sgn(Cross(l2.s-l1.s,l1.e-l1.s))*sgn(Cross(l1.e-l1.s,l2.e-l1.s))>=0&&
    sgn(Cross(l1.s-l2.s,l2.e-l2.s))*sgn(Cross(l2.e-l2.s,l1.e-l2.s))>=0;
}
```
<span id="6"><h4>6.	半平面交</h4></span>
```cpp
bool HPIcmp(Line a,Line b) {
    if (fabs(a.ang - b.ang) > EPS)return a.ang < b.ang;
    else return ((a.s - b.s) ^ (b.e - b.s)) < 0;
}
struct halfplanes{
    int n;
    Line line[MAXN];
    Point p[MAXN];
    Line Q[MAXN];
    int head,tail;
    bool HPI() {
        sort(line+1, line+n+1, HPIcmp);
        int tot = 1;
        for (int i = 2; i <= n; i++)
            if (fabs(line[i].ang - line[i - 1].ang) > EPS) //去掉斜率重复的
                line[++tot] = line[i];
        head = 1, tail = 2;
        Q[1] = line[1];Q[2] = line[2];
        for (int i = 3; i <= tot; i++) {
            if (fabs((Q[tail].e - Q[tail].s) ^ (Q[tail - 1].e - Q[tail - 1].s)) < EPS ||
                fabs((Q[head].e - Q[head].s) ^ (Q[head + 1].e - Q[head + 1].s)) < EPS)
                return false;
            while (head < tail && (((Q[tail] & Q[tail - 1]) -
                                    line[i].s) ^ (line[i].e - line[i].s)) > EPS)
                tail--;
            while (head < tail && (((Q[head] & Q[head + 1]) -
                                    line[i].s) ^ (line[i].e - line[i].s)) > EPS)
                head++;
            Q[++tail] = line[i];
        }
        while (head < tail && (((Q[tail] & Q[tail - 1]) -
                                Q[head].s) ^ (Q[head].e - Q[head].s)) > EPS)
            tail--;
        while (head < tail && (((Q[head] & Q[head + 1]) -
                                Q[tail].s) ^ (Q[tail].e - Q[tail].s)) > EPS)
            head++;
        if (tail <= head + 1) return false;
        return true;
    }
    void getConvex(Point res[], int &resn){
        resn=0;
        for(int i = head; i < tail; i++)
            res[++resn] = Q[i]&Q[i+1];
        if(head < tail - 1)
            res[++resn] = Q[head]&Q[tail];
    }
}hpi;

//2021-02-22
bool HPIcmp(Line a,Line b) {
    if (fabs(a.ang - b.ang) > EPS)return a.ang < b.ang;
    else return ((a.s - b.s) ^ (b.e - b.s)) < 0;
}
struct halfplanes{
    int n;
    Line line[100050];
    Point p[100050];
    Line Q[100050];
    int head,tail;
    double HPI() {
        rep(i,1,n){
            line[i].ang=atan2(line[i].e.y-line[i].s.y,line[i].e.x-line[i].s.x);
        }
        sort(line+1, line+n+1, HPIcmp);
        int tot = 1;
        rep(i,2,n){
            if (fabs(line[i].ang - line[i - 1].ang) > EPS) //去掉斜率重复的
                line[++tot] = line[i];
        }
        head = 1, tail = 2;
        Q[1] = line[1];Q[2] = line[2];
        rep(i,3,tot){
            if (fabs((Q[tail].e - Q[tail].s) ^ (Q[tail - 1].e - Q[tail - 1].s)) < EPS ||
                fabs((Q[head].e - Q[head].s) ^ (Q[head + 1].e - Q[head + 1].s)) < EPS)
                return false;
            while (head < tail && sgn((Q[tail].crossPoint(Q[tail - 1]) - line[i].s) ^ (line[i].e - line[i].s)) > 0)
                tail--;
            while (head < tail && sgn((Q[head].crossPoint(Q[head + 1]) - line[i].s) ^ (line[i].e - line[i].s)) > 0)
                head++;
            Q[++tail] = line[i];
        }
        while (head < tail && sgn((Q[tail].crossPoint(Q[tail - 1]) -
                                Q[head].s) ^ (Q[head].e - Q[head].s)) > 0)
            tail--;
        while (head < tail && sgn((Q[head].crossPoint(Q[head + 1]) -
                                Q[tail].s) ^ (Q[tail].e - Q[tail].s)) > 0)
            head++;
        if (tail <= head + 1) return false;
        int num=0;
        getConvex(p,num);
        double ans=0;
        rep(i,1,num){
            ans+=(p[i]-p[i%num+1]).len();
        }
        return ans;
    }
    void getConvex(Point res[], int &resn){
        resn=0;
        for(int i = head; i < tail; i++)
            res[++resn] = Q[i].crossPoint(Q[i+1]);
        if(head < tail - 1)
            res[++resn] = Q[head].crossPoint(Q[tail]);
    }
}hpi;
-------------------------------------------------多项式
            double a=0,b=0,c=0;
            a=2*(X[i]-X[j]);
            b=2*(Y[i]-Y[j]);
            c=XX[j]+YY[j]-XX[i]-YY[i];

            Point pp;
            if(sgn(a)==0){
                pp={0,-c/b};
            }else{
                pp={-c/a,0};
            }
            Vector v={b,-a};
            hpi.line[++hpi.n]={pp,pp+v};

//            Vector v=Vector(b,-a);
//            Point p;
//            if(fabs(a)>fabs(b)) p=Point(-c/a,0);
//            else p=Point(0,-c/b);
//            hpi.line[++hpi.n]=Line(p,p+v);
```
<span id="7"><h4>7.	四面体内切球</h4></span>
https://blog.csdn.net/helloiamclh/article/details/51971951
↑注意求xyz那里不是除6而是除(s1+s2+s3+s4)
I - tetrahedron HDU - 5733
```cpp
    Vector3 AB=p[2]-p[1];
    Vector3 AC=p[3]-p[1];
    Vector3 AD=p[4]-p[1];
    double V=fabs((AB^AC)*AD/6);
    if(sgn(V)==0){
        printf("0 0 0 0\n");
        return ;
    }
    double s1=((p[2]-p[4])^(p[3]-p[4])).len()/2;
    double s2=((p[1]-p[4])^(p[3]-p[4])).len()/2;
    double s3=((p[1]-p[4])^(p[2]-p[4])).len()/2;
    double s4=((p[1]-p[3])^(p[2]-p[3])).len()/2;
    double x=p[1].x*s1+p[2].x*s2+p[3].x*s3+p[4].x*s4;
    double y=p[1].y*s1+p[2].y*s2+p[3].y*s3+p[4].y*s4;
    double z=p[1].z*s1+p[2].z*s2+p[3].z*s3+p[4].z*s4;
    x/=(s1+s2+s3+s4);y/=(s1+s2+s3+s4);z/=(s1+s2+s3+s4);
    double r=Point3(x,y,z).point_to_plane(p[1],p[2],p[3]);//r=3*V/(s1+s2+s3+s4)
    printf("%.4f %.4f %.4f %.4f\n",x,y,z,r);
```
<span id="7"><h4>8.	平面最近点对</h4></span>
https://blog.csdn.net/GGN_2015/article/details/80785621
https://www.cnblogs.com/zyxStar/p/4591897.html
C - Quoit Design HDU - 1007
```cpp
bool cmpY(int& a,int& b){
    return p[a].y<p[b].y;
}
int mpt[MAXN];
double Closest_Pair(int left, int right){
    double d = 0x3f3f3f3f;
    if (left == right)
        return d;
    if (left + 1 == right)
        return p[left].disToPoint(p[right]);
    int mid = (left + right) >> 1;
    double d1 = Closest_Pair(left, mid);
    double d2 = Closest_Pair(mid + 1, right);
    d = min(d1, d2);
    int i, j, k = 0;
    //分离出宽度为d的区间
    for (i = left; i <= right; i++){
        if (fabs(p[mid].x - p[i].x) <= d)
            mpt[k++] = i;
    }
    sort(mpt, mpt + k, cmpY);
    //线性扫描
    for (i = 0; i < k; i++){
        for (j = i + 1; j < k && p[mpt[j]].y - p[mpt[i]].y<d; j++){
            double d3 = p[mpt[i]].disToPoint(p[mpt[j]]);
            if (d > d3)    d = d3;
        }
    }
    return d;
}
void solve(){
    sort(p+1,p+n+1);
    printf("%.2f\n",Closest_Pair(1,n)/2);
}
```
<span id="9"><h4>9.	旋转卡壳求最大三角形面积</h4></span>
```cpp
double rotating_calipers_S(Point* stk,int mn){//最大三角形面积
    stk[mn+1]=stk[1];
    int cur=2;
    double ret=0;
    for(int i=1;i<=mn;i++){
        for(int j=i+1;j<=mn;j++){
            Vector v = stk[i]-stk[j];
            while((v^(stk[cur+1]-stk[cur])) < 0)
                cur = (cur%mn+1);
            double tmp=max((stk[cur]-stk[i])^v,(stk[cur+1]-stk[j])^v);
            ret = max(ret,tmp);
        }
    }
    return fabs(ret)/2;
}

-------------------------------↑rubbish
double rota(int n){
    double ret=0;
    int j=2,k=3;
    Vector v;
    for(int i=1;i<=n;i++){
        v=stk[i]-stk[j];
        while((v^(stk[k+1]-stk[k]))<0)
            k=k%n+1;
        ret=max(ret,(stk[k]-stk[j])^v);
        v=stk[i]-stk[k];
        while((v^(stk[j+1]-stk[j]))>0)
            j=j%n+1;
        ret=max(ret,v^(stk[j]-stk[k]));
    }
    return fabs(ret)/2;
}
void solve(){
    if(n<=2){
        printf("0.00\n");
        return;
    }
    int len=Andrew();
    printf("%.2f\n",rota(len));
}
```
<span id="10"><h4>10. 最小圆覆盖</h4></span>
```cpp
Point circumCenter(Point a, Point b, Point c){ //返回三角形的外心
    Point ret;
    double a1 = b.x-a.x,b1 = b.y-a.y,c1 = (a1*a1+b1*b1)/2;
    double a2 = c.x-a.x,b2 = c.y-a.y,c2 = (a2*a2+b2*b2)/2;
    double d = a1*b2-a2*b1;
    ret.x=a.x+(c1*b2-c2*b1)/d;
    ret.y=a.y+(a1*c2-a2*c1)/d;
    return ret;
}
void min_cover_circle(Point p[],int n,Circle &c) {
    c.c = p[1], c.r = 0;
    for (int i = 2; i <= n; i++) {
        if (cmp(p[i].disToPoint(c.c),c.r) > 0) {
            c.c = p[i];
            c.r = 0;
            for (int j = 1; j <= i; j++) {
                if (cmp(p[j].disToPoint(c.c),c.r) > 0) {
                    c.c = Point((p[i].x + p[j].x) / 2, (p[i].y + p[j].y) / 2);
                    c.r = p[j].disToPoint(c.c);
                    for (int g = 1; g <= j; g++) {
                        if (cmp(p[g].disToPoint(c.c), c.r) > 0) {
                            c.c = circumCenter(p[i], p[j], p[g]);
                            c.r = p[i].disToPoint(c.c);
                        }
                    }
                }
            }
        }
    }
}
```
<span id="11"><h4>11. 模拟退火求矩阵内与所有点最小距离最大</h4></span>
```cpp
Point p[2005];
double getDis(Point a){
    double dis=1e9;
    for(int i=1;i<=n;i++)
        dis=min(dis,a.disToPoint(p[i]));
    return dis;
}
Point ans[30];
double dis[30];
Point sa(){
    for(int i=1;i<=20;i++){
        ans[i].x=(rand()%1000+1)/1000.0*x;
        ans[i].y=(rand()%1000+1)/1000.0*y;
        dis[i]=getDis(ans[i]);
    }
    srand(time(0));
    double t=300;
    while(t>1e-7){
        for(int i=1;i<=20;i++){
            for(int j=1;j<=20;j++){
                Point tmp=ans[i];
                double angle=rand()%1000/1000.0*2*PI;
                tmp.x+=t*cos(angle)*(rand()%1000/1000.0);
                tmp.y+=t*sin(angle)*(rand()%1000/1000.0);
                if(tmp.x<0||tmp.x>x||tmp.y<0||tmp.y>y)continue;
                double tmpdis=getDis(tmp);
                if(tmpdis>dis[i]){
                    dis[i]=tmpdis;
                    ans[i]=tmp;
                }
            }
        }
        t*=0.96;
    }
    double dd=0;
    int pp=0;
    for(int i=1;i<=20;i++){
        if(getDis(ans[i])>dd){
            dd=getDis(ans[i]);
            pp=i;
        }
    }
    return ans[pp];
}
```
<span id="12"><h4>12. 矩阵面积凸包最小矩形覆盖</h4></span>
```cpp
double rota(int n){
    Vector v;
    int up=2;
    int left=2;
    int right=2;
    double ret=1e15;
    for(int down=1;down<=n;down++) {

        v = stk[down + 1] - stk[down];
        while (sgn(v ^ (stk[up + 1] - stk[up])) >= 0)
            up = up % n + 1;

        while(cmp(v*(stk[right]-stk[down]),
                  v*(stk[right+1]-stk[down]))<=0)
            right=right%n+1;

        left=up%n+1;
        while(cmp((stk[left]-stk[down+1])*(stk[down]-stk[down+1]),
                  (stk[left+1]-stk[down+1])*(stk[down]-stk[down+1]))<=0)
            left=left%n+1;

        double h=fabs(v^(stk[up]-stk[down]))/v.len();
        double w1=fabs((stk[left]-stk[down])*(stk[down]-stk[down+1]))/v.len();
        double w2=fabs(v*(stk[right]-stk[down]))/v.len();
        double w=(w1+w2);
        double tmp=h*w;
        if(tmp<ret){
            ret=tmp;
            ans[0]=stk[down]-(v*w1/v.len());//矩形四个点
            ans[1]=stk[down]+(v*w2/v.len());
            Vector vvv=v.rotLeft()/v.len();
            ans[2]=ans[1]+vvv*h;
            ans[3]=ans[0]+vvv*h;
        }
    }
    return ret;
}
void solve(){
    int len=Andrew();
    printf("%.0f\n",rota(len));
}
```
<span id="13"><h4>13. 稳定凸包</h4></span>
```cpp
Point p[1010],stk[1010];
int Andrew(){
    sort(p+1,p+n+1);
    int len=0;
    for (int i=1;i<=n;i++){
        while (len>1&&sgn((stk[len]-stk[len-1])^(p[i]-stk[len-1]))<0) len--;
        stk[++len]=p[i];
    }
    int k=len;
    for (int i=n-1;i>=1;i--){
        while (len>k&&sgn((stk[len]-stk[len-1])^(p[i]-stk[len-1]))<0) len--;
        stk[++len]=p[i];
    }
    len--;
    return len;
}
void solve(){
    if(n<6){
        printf("NO\n");
        return ;
    }
    int len=Andrew();
    stk[0]=stk[len];
    stk[len+1]=stk[1];
    stk[len+2]=stk[2];
    for(int i=1;i<=len;i++){
        if(sgn((stk[i+2]-stk[i])^(stk[i+1]-stk[i]))==0||
            sgn((stk[i+1]-stk[i-1])^(stk[i]-stk[i-1]))==0)
            continue;
        printf("NO\n");
        return;
    }
    printf("YES\n");
}
```
<span id="14"><h4>14. simpson积分</h4></span>
```cpp
int w,t,v;
struct unb{
    double x,l,v;
    void read(){
        scanf("%lf%lf%lf",&x,&l,&v);
    }
}ub[15];
double ans;
double d[15];
Point line[15];
double f(double now){
    for(int i=1;i<=n;i++){
        if(cmp(ub[i].l,w)>=0){
            line[i].x=ub[i].x;
            line[i].y=ub[i].x+ub[i].l;
            continue;
        }
        double dis=fmod(now*ub[i].v,2*(w-ub[i].l));
        double tmp=ub[i].x+dis;
        if(tmp<=0)
            tmp+=2*(w-ub[i].l);
        tmp=fmod(tmp+ub[i].x,2*(w-ub[i].l));
        if(tmp>=w-ub[i].l)
            line[i].x=2*(w-ub[i].l)-tmp;
        else
            line[i].x=tmp;
        line[i].y=line[i].x+ub[i].l;
    }
    sort(line+1,line+1+n);
    double len=0;
    double pre=-1;
    double tmplen=0;
    for(int i=1;i<=n;i++){
        if(line[i].x>pre+tmplen){
            len+=tmplen;
            pre=line[i].x;
            tmplen=line[i].y-line[i].x;
        }else{
            tmplen=max(tmplen,line[i].y-pre);
        }
    }
    return w-tmplen-len;
}
double simpson(double a,double b){
    return (b-a)*(f(a)+f(b)+4*f((a+b)/2))/6;
}
double simpson(double l,double r,double now){
    double mid=(l+r)/2;
    double ls=simpson(l,mid);
    double rs=simpson(mid,r);
    if(cmp(rs+ls,now)==0) return ls+rs;
    return simpson(l,mid,ls)+simpson(mid,r,rs);
}
void solve(){

    printf("%.2f\n",simpson(0,t,simpson(0,t))*v);
}
```
<span id="15"><h4>15. 三维几何凸包</h4></span>
```cpp
struct Point3{
    double x,y,z;
    Point3(){}
    Point3(double a, double b, double c):x(a),y(b),z(c){}

    double len(){return sqrt(x*x+y*y+z*z);}
    double norm(){return x*x+y*y+z*z;}
    double disToPoint(Point3 p){return sqrt((x-p.x)*(x-p.x)+(y-p.y)*(y-p.y)+(z-p.z)*(z-p.z));}
//    double disToLine(Line l);
//    double disToSegment(Line l);
//    bool onLine(Line l);
//    bool onSegment(Line l);
//    Point project(Line l);//投影
//    Point symmetryPoint(Line l);
//    int inPolygon(Point poly[]);//点是否在多边形内部
    bool dots_inline(Point3 p2,Point3 p3){
        return sgn(((*this-p2)^(p2-p3)).len())==0;
    }//三点共线
    double point_to_plane(Point3 p1,Point3 p2,Point3 p3){
        return fabs(((p1-p2)^(p2-p3))*(*this-p1))/((p1-p2)^(p2-p3)).len();
    }
    void read(){scanf("%lf%lf%lf",&x,&y,&z);}
    Point3 operator+(Vector3 v){return {x+v.x,y+v.y,z+v.z};}
    Vector3 operator-(Point3 p){return {x-p.x,y-p.y,z-p.z};}
    Vector3 operator^(Vector3 v){return {y*v.z-z*v.y,z*v.x-x*v.z,x*v.y-y*v.x};}//叉乘
    double operator*(Vector3 v){return x*v.x+y*v.y+z*v.z;}//点乘
    Vector3 operator*(double d){return {x*d,y*d,z*d};}
    Vector3 operator/(double d){return {x/d,y/d,z/d};}
    bool operator==(Point3 p){return cmp(x,p.x)==0&&cmp(y,p.y)==0&&cmp(z,p.z)==0;}
    bool operator<(Point3 p){if(cmp(x,p.x)==0&&cmp(y,p.y))return z<p.z;
                            else if(cmp(x,p.x))return y<p.y;return x<p.x;}
};
//Point3 pt[MAXN];
inline double area(Point3 o, Point3 a, Point3 b) {
    return ((o-a)^(o-b)).len();
}//三角形面积*2
inline double volume(Point3 o, Point3 a, Point3 b, Point3 c) {
    return ((o-a)^(o-b))*(c-o) / 6.0;
}

struct CH3D {
    struct face {
        int a, b, c;//表示凸包一个面上的三个点的编号
        bool ok;//表示该面是否属于最终凸包上的面
    };
    int n;//初始顶点数
    Point3 P[MAXN];//初始顶点
    int num;//凸包表面的三角形数
    face F[8 * MAXN];//凸包表面的三角形
    int g[MAXN][MAXN];//凸包表面的三角形
//    double vlen(Point3 a) { return sqrt(a.x * a.x + a.y * a.y + a.z * a.z); }//向量长度
    Point3 cross(const Point3 &a, const Point3 &b, const Point3 &c) {
        return {(b.y - a.y) * (c.z - a.z) - (b.z - a.z) * (c.y - a.y),
                (b.z - a.z) * (c.x - a.x) - (b.x - a.x) * (c.z - a.z),
                (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x)};
    } //叉乘
//    double area(Point3 a, Point3 b, Point3 c) { return ((b - a) ^ (c - a)).len(); }
//    double volume(Point3 a, Point3 b, Point3 c, Point3 d) { 
        //return ((b - a) ^ (c - a)) * (d - a); }//四面体有向体积*6
    double dblcmp(Point3 &p, face &f) {
        Point3 m = P[f.b] - P[f.a];
        Point3 n = P[f.c] - P[f.a];
        Point3 t = p - P[f.a];
        return (m ^ n) * t;
    }//正：点在面同向

    void deal(int p, int a, int b) {
        int f = g[a][b];//搜索与该边相邻的另一个平面
        face add;
        if (F[f].ok) {
            if (dblcmp(P[p], F[f]) > EPS)
                dfs(p, f);
            else {
                add.a = b;
                add.b = a;
                add.c = p;//这里注意顺序，要成右手系
                add.ok = true;
                num++;
                g[p][b] = g[a][p] = g[b][a] = num;
                F[num] = add;
            }
        }
    }

    void dfs(int p, int now) {
        F[now].ok = 0;
        deal(p, F[now].b, F[now].a);
        deal(p, F[now].c, F[now].b);
        deal(p, F[now].a, F[now].c);
    }//递归搜索所有应该从凸包内删除的面
    bool same(int s, int t) {
        Point3 &a = P[F[s].a];
        Point3 &b = P[F[s].b];
        Point3 &c = P[F[s].c];
        return fabs(volume(a, b, c, P[F[t].a])) < EPS &&
               fabs(volume(a, b, c, P[F[t].b])) < EPS &&
               fabs(volume(a, b, c, P[F[t].c])) < EPS;
    }

    void create() {
        num = 0;
        if (n < 4)return;
        bool flag = true;
        for (int i = 2; i <= n; i++) {
            if ((P[1] - P[i]).len() > EPS) {
                swap(P[2], P[i]);
                flag = false;
                break;
            }
        }
        if (flag)return;
        flag = true;
        for (int i = 3; i <= n; i++) {
            if (((P[1] - P[2]) ^ (P[2] - P[i])).len() > EPS) {
                swap(P[3], P[i]);
                flag = false;
                break;
            }
        }
        if (flag)return;
        flag = true;
        for (int i = 4; i <= n; i++) {
            if (fabs(((P[1] - P[2]) ^ (P[2] - P[3])) * (P[1] - P[i])) > EPS) {
                swap(P[4], P[i]);
                flag = false;
                break;
            }
        }
        if (flag)return;

        for (int i = 1; i <= 4; i++) {
            face add;
            add.a = (i) % 4 + 1;
            add.b = (i + 1) % 4 + 1;
            add.c = (i + 2) % 4 + 1;
            add.ok = true;
            if (dblcmp(P[i], add) > 0)swap(add.b, add.c);
            num++;
            g[add.a][add.b] = g[add.b][add.c] = g[add.c][add.a] = num;
            F[num] = add;
        }
        for (int i = 5; i <= n; i++) {
            for (int j = 1; j <= num; j++) {
                if (F[j].ok && dblcmp(P[i], F[j]) > EPS) {
                    dfs(i, j);
                    break;
                }
            }
        }
        int tmp = num;
        num = 0;
        for (int i = 1; i <= tmp; i++)
            if (F[i].ok)
                F[++num] = F[i];
    }

    double CH_area() {//表面积
        double res = 0;
        if (n == 3) {
            Point3 p = cross(P[1], P[2], P[3]);
            res = p.len() / 2.0;
            return res;
        }
        for (int i = 1; i <= num; i++)
            res += area(P[F[i].a], P[F[i].b], P[F[i].c]);
        return res / 2.0;
    }

    double CH_volume() {
        double res = 0;
        Point3 tmp(0, 0, 0);
        for (int i = 1; i <= num; i++)
            res += volume(tmp, P[F[i].a], P[F[i].b], P[F[i].c]);
        return fabs(res / 6.0);
    }

    int triangle() {
        return num;
    }//表面三角形个数
    int polygon() {
        int res = 0;
        for (int i = 1, flag; i <= num; i++) {
            flag = 1;
            for (int j = 1; j < i; j++)
                if (same(i, j)) {
                    flag = 0;
                    break;
                }
            res += flag;
        }
        return res;
    }

    Point3 barycenter() {
        Point3 ans(0, 0, 0), o(0, 0, 0);
        double all = 0;
        for (int i = 1; i <= num; i++) {
            double vol = volume(o, P[F[i].a], P[F[i].b], P[F[i].c]);
            ans = ans + (o + P[F[i].a] + P[F[i].b] + P[F[i].c]) / 4.0 * vol;
            all += vol;
        }
        ans = ans / all;
        return ans;
    }//三维凸包重心
    double ptoface(Point3 p, int i) {
        return fabs(
                volume(P[F[i].a], P[F[i].b], P[F[i].c], p) / 
                ((P[F[i].b] - P[F[i].a]) ^ (P[F[i].c] - P[F[i].a])).len()
                );
    }//点到面的距离
}ch;
void solve(){
    ch.create();
    printf("%d\n",ch.polygon());

}
```
<span id="16"><h4>16. 平面最近点对</h4></span>
```cpp
https://blog.csdn.net/GGN_2015/article/details/80785621
https://www.cnblogs.com/zyxStar/p/4591897.html
C - Quoit Design HDU - 1007

bool cmpY(int& a,int& b){
    return p[a].y<p[b].y;
}
int mpt[MAXN];
double Closest_Pair(int left, int right){
    double d = 0x3f3f3f3f;
    if (left == right)
        return d;
    if (left + 1 == right)
        return p[left].disToPoint(p[right]);
    int mid = (left + right) >> 1;
    double d1 = Closest_Pair(left, mid);
    double d2 = Closest_Pair(mid + 1, right);
    d = min(d1, d2);
    int i, j, k = 0;
    //分离出宽度为d的区间
    for (i = left; i <= right; i++){
        if (fabs(p[mid].x - p[i].x) <= d)
            mpt[k++] = i;
    }
    sort(mpt, mpt + k, cmpY);
    //线性扫描
    for (i = 0; i < k; i++){
        for (j = i + 1; j < k && p[mpt[j]].y - p[mpt[i]].y<d; j++){
            double d3 = p[mpt[i]].disToPoint(p[mpt[j]]);
            if (d > d3)    d = d3;
        }
    }
    return d;
}
void solve(){
    sort(p+1,p+n+1);
    printf("%.2f\n",Closest_Pair(1,n)/2);
}
```


<span id="17"><h4>17. 旋转卡壳求最大四边形面积</h4></span>
```cpp
double rota(int n){
    double ret=0;
    int j=2,k=3,l=2;
    Vector v;
    for(int i=1;i<=n;i++){
        v=stk[i+1]-stk[i];
        while((v^(stk[k+1]-stk[k]))>0)
            k=k%n+1;
        v=stk[i]-stk[k];
//        j=i%n+1;
        while((v^(stk[j+1]-stk[j]))>0)
            j=j%n+1;
        double sum=v^(stk[j]-stk[k]);
        v=stk[k]-stk[i];
//        l=k%n+1;
        while((v^(stk[l+1]-stk[l]))>0)
            l=l%n+1;
        sum+=v^(stk[l]-stk[i]);
        ret=max(ret,sum);
    }
    return fabs(ret)/2;
}


long long rota(int n){
    long long ret=0;
    int j=2,k=3,l=2;
    Vector v;
    for(int i=1;i<=n;i++){//n^2应该一定能过
        l=i+3;
        j=i+1;
        for(int k=i+2;k<=n;k++){
            v=stk[i]-stk[k];
            while((v^(stk[j+1]-stk[j]))>0)
                j=j%n+1;

            v=stk[k]-stk[i];
            while((v^(stk[l+1]-stk[l]))>0)
                l=l%n+1;

            long long sum=(stk[i]^stk[j])+(stk[j]^stk[k])+(stk[k]^stk[l])+(stk[l]^stk[i]);
            ret=max(ret,abs(sum));
        }
    }
//    bool first=true;
//    for(int i=1;i<=n;i++){//nlogn不会证明
//        v=stk[i+1]-stk[i];
//        while((v^(stk[k+1]-stk[k]))>0)
//            k=k%n+1;
//        v=stk[i]-stk[k];
//        while((v^(stk[j+1]-stk[j]))>0)
//            j=j%n+1;
//        long long sum=abs(v^(stk[j]-stk[k]));
//        v=stk[k]-stk[i];
//        if(first){
//            l=k%n+1;
//            first=false;
//        }
//        while((v^(stk[l+1]-stk[l]))>0)
//            l=l%n+1;
//        sum+=abs(v^(stk[l]-stk[i]));
//        ret=max(ret,sum);
//    }
    return ret;
}
```


<span id="18"><h4>18. 可加点的动态凸包</h4></span>
```cpp
//bzoj2300【HAOI2011】防线修建  每次删除某个点或求边长

#include <bits/stdc++.h>
#define pb push_back
#define memarray(array,val) memset(array,val,sizeof(array))
using namespace std;
const int mod=1e9+7;
const double PI=acos(-1.0);
const double EPS=1e-6;
inline int sgn(double a){
    if(a<-EPS)return -1;
    return a>EPS;
}
inline int cmp(double a,double b){
    return sgn(a-b);
}
const int MAXN=2e5+10;
const int INF=0x3f3f3f3f;
int n,m,q;
int ty[MAXN],cancel[MAXN];
bool no[MAXN];
double ans;
struct Point{
    double x,y;
    Point(double _x=0.0,double _y=0.0){
        x=_x;y=_y;
    }
    double len(){
        return sqrt(len2());
    }
    double len2(){
        return x*x+y*y;
    }
    void read(){
        scanf("%lf%lf",&x,&y);
    }
    bool operator<(const Point p)const{
        if(cmp(x,p.x)==0)return y<p.y;
        return x<p.x;
    }
    Point operator+(const Point p)const{
        return (Point){x+p.x,y+p.y};
    }
    Point operator-(const Point p)const{
        return (Point){x-p.x,y-p.y};
    }
    double operator*(const Point p)const{
        return x*p.x+y*p.y;
    }
    double operator^(const Point p)const{
        return x*p.y-y*p.x;
    }
}p[MAXN],ptmp[MAXN];
Point p1,p2,p3;
set<Point>stk;
int Andrew(){
    int num=0;
    for(int i=1;i<=n;i++){
        if(no[i])continue;
        ptmp[++num]=p[i];
    }
    ptmp[++num]=p1;
    ptmp[++num]=p2;
    ptmp[++num]=p3;
    sort(ptmp+1,ptmp+1+num);
    stk.insert(ptmp[num]);
    int len=1;
    for(int i=num-1;i>=1;i--){
        while(len>1){
            set<Point>::iterator sss=stk.begin();
            Point t1=*sss;
            sss++;
            Point t2=*sss;
            if(sgn(((t1-t2)^(ptmp[i]-t2)))<0){
                stk.erase(stk.begin());
                len--;
            }else{
                break;
            }
        }
        stk.insert(ptmp[i]);
        len++;
    }
    return len;
}
void insert(int idx){
    Point point=p[idx];
    set<Point>::iterator it,it2;
    it=lower_bound(stk.begin(),stk.end(),point);
    it2=it;
    it2--;
    Point point1=*it2,point2=*it;
    if(sgn(((point1-point2)^(point-point2)))>=0)return;
    ans-=(point1-point2).len();
    while(it2!=stk.begin()){
        it2--;
        point2=*it2;
        if(sgn(((point2-point)^(point1-point)))>0){
            stk.erase(point1);
            ans-=(point1-point2).len();
            point1=point2;
        }else{
            break;
        }
    }
    Point leftPoint=point1;
    it2=it;
    point1=*it2;
    it2++;
    while(it2!=stk.end()){
        point2=*it2;
        if(sgn(((point1-point)^(point2-point)))>0){
            stk.erase(point1);
            ans-=(point1-point2).len();
            point1=point2;
            it2++;
        }else{
            break;
        }
    }
    if(sgn(((point-point1)^(leftPoint-point1)))>=0){
        stk.insert(point);
        ans+=(point-leftPoint).len();
        ans+=(point-point1).len();
    }
}
void solve() {
    Andrew();
    Point pre={0,0};
    for(auto i:stk){
        ans+=(i-pre).len();
        pre=i;
    }
    stack<double>two;
    for(int i=q;i>=1;i--){
        if(ty[i]==2){
            two.push(ans);
        }else{
            insert(cancel[i]);
        }
    }
    while(!two.empty()){
        double tmp=two.top();
        two.pop();
        printf("%.2f\n",tmp);
    }
}
void init(){
    double x,y;
    scanf("%lf",&x);
    p1=(Point){0,0};
    p2=(Point){x,0};
    scanf("%lf%lf",&x,&y);
    p3=(Point){x,y};
    scanf("%d",&n);
    for(int i=1;i<=n;i++)
        p[i].read();
    scanf("%d",&q);
    for(int i=1;i<=q;i++){
        scanf("%d",ty+i);
        if(ty[i]==1){
            scanf("%d",cancel+i);
            no[cancel[i]]=true;
        }
    }
}
int main()
{
    int T=1;
//    scanf("%d",&T);
    while(T--)
    {
        init();
        solve();
    }
    return 0;
}

```

<span id="19"><h4>19. KD-Tree</h4></span>
```cpp
const int DG=2;
struct Point{
    double x[DG];
    void read(){
        rep(i,0,DG-1){
            scanf("%lf",&x[i]);
        }
    }
    double len(){
        return sqrt(len2());
    }
    double disToPoint(const Point p){
        return (*this-p).len2();
    }
    double len2(){
        double ret=0;
        rep(i,0,DG-1){
            ret+=x[i]*x[i];
        }
        return ret;
    }
    Point operator-(const Point p)const{
        Point ret;
        rep(i,0,DG-1){
            ret.x[i]=x[i]-p.x[i];
        }
        return ret;
    }
}p[MAXN];
struct KDT{
private:
    double min_ans,max_ans;
    Point p[MAXN],mi[MAXN],ma[MAXN],_p;
    int tot,rt,ls[MAXN],rs[MAXN];
    int nw(const Point &_p){
        p[++tot]=_p;
        rep(i,0,DG-1){
            mi[tot].x[i]=ma[tot].x[i]=_p.x[i];
        }
        ls[tot]=rs[tot]=0;
        return tot;
    }
    void upd(int u,const Point &_p){
        rep(i,0,DG-1){
            mi[u].x[i]=min(mi[u].x[i],_p.x[i]);
            ma[u].x[i]=max(ma[u].x[i],_p.x[i]);
        }
    }
    void ins(int &u,bool d){
        if(!u){
            u=nw(_p);
            return;
        }
        ins(_p.x[d]<=p[u].x[d]?ls[u]:rs[u],d^1);
        upd(u,_p);
    }
    double min_dis(int u,const Point &_p){
        if(!u)return 1e15;
        double dis=0;
        rep(i,0,DG-1){
            double tmp=max(_p.x[i]-ma[u].x[i],0.0)+max(mi[u].x[i]-_p.x[i],0.0);
            dis+=tmp*tmp;
        }
        return dis;
    }
    double max_dis(int u,const Point &_p){
        if(!u)return 0;
        double dis=0;
        rep(i,0,DG-1){
            double tmp=max(fabs(_p.x[i]-ma[u].x[i]),fabs(mi[u].x[i]-_p.x[i]));
            dis+=tmp*tmp;
        }
        return dis;
    }
    void query_min(int u){
        if(!u)return;
        min_ans=min(min_ans,p[u].disToPoint(_p));
        double lv=min_dis(ls[u],_p),rv=min_dis(rs[u],_p);
        if(lv>rv){
            if(rv<min_ans)query_min(rs[u]);
            if(lv<min_ans)query_min(ls[u]);
        }else{
            if(lv<min_ans)query_min(ls[u]);
            if(rv<min_ans)query_min(rs[u]);
        }
    }
    void query_max(int u){
        if(!u)return;
        max_ans=max(max_ans,p[u].disToPoint(_p));
        double lv=max_dis(ls[u],_p),rv=max_dis(rs[u],_p);
        if(lv<rv){
            if(rv>max_ans)query_max(rs[u]);
            if(lv>max_ans)query_max(ls[u]);
        }else{
            if(lv>max_ans)query_max(ls[u]);
            if(rv>max_ans)query_max(rs[u]);
        }
    }
public:
    void ins(const Point &P){
        _p=P;
        ins(rt,0);
    }
    double ask_min(const Point &P){
        _p=P;
        min_ans=1e15;
        query_min(rt);
        return min_ans;
    }
    double ask_max(const Point &P){
        _p=P;
        max_ans=0;
        query_max(rt);
        return max_ans;
    }
}tr;
void solve(){
    double mn=1e15;
    double mx=0;
    rep(i,1,n){
        mn=min(mn,tr.ask_min(p[i]));
        mx=max(mx,tr.ask_max(p[i]));
        tr.ins(p[i]);
    }
    printf("%.4f %.4f\n",sqrt(mn),sqrt(mx));
}
void init() {
    scanf("%d",&n);
    rep(i,1,n){
        p[i].read();
    }
    random_shuffle(p+1,p+1+n);
}
```
