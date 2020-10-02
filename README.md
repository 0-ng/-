####1.	基础模板
####2.	快读快写
####3.	线性筛素数
####4.	exgcd求逆元
####5.	求解ax+by=c
####6.	矩阵
####7.	矩阵快速幂
####8.	矩阵快速幂求斐波那契
####9.	漂浮法
####10.	归并排序求逆序对
####11.	RMQ求区间最值
####12.	中国剩余问题
####13.	二维ST表
####14.	线段树
####15.	线段树区间最值单点修改
####16.	Dijstra最短路
####17.	Floyd最短路
####18.	Jhonson全源最短路
####19.	LCA
####20.	Tarjan找割点
####21.	Tarjan找强连通分量个数
####22.	Tarjan缩点求路径最大点权和
####23.	二分图最大匹配匈牙利算法
####24.	拓扑排序
####25.	最大流
####26.	最小费用最大流
####27.	欧拉路径
####28.	KMP
####29.	几何
####30.	直线与凸多边形交点
####31.	凸包Andrew
####32.	凸包Graham
####33.	判断线段相交
####34.	半平面交
####35.	四面体内切球
####36.	平面最近点对
####37.	旋转卡壳求最大三角形面积
####38.	最小圆覆盖
####39.	模拟退火求矩阵内所有点最小距离最大
####40.	矩阵面积凸包最小矩形覆盖
####41.	稳定凸包大整数模板

---------------------
####1.	基础模板
```cpp
#include<map>
#include<set>
#include<cmath>
#include<queue>
#include<stack>
#include<ctime>
#include<vector>
#include<cstdio>
#include<vector>
#include<string>
#include<bitset>
#include<cstdlib>
#include<iomanip>
#include<cstring>
#include<iostream>
#include<algorithm>
using namespace std;
#define pb push_back
#define memarray(array, value) memset(array, value, sizeof(array))
const double EPS=1e-8;
const double PI=acos(-1);
const long long mod=1e9+7;
const int INF=0x3f3f3f3f; 
const int MAXN=1e5+10;
const int N=1e5+10;
```
####2.	快读快写
```cpp
inline __int128 read(){//输入模板
    __int128 x=0,f=1;
    char ch=getchar();
    while(ch<'0'||ch>'9'){
        if(ch=='-') f=-1;
        ch=getchar();
    }
    while(ch>='0'&&ch<='9'){
        x=x*10+ch-'0';
        ch=getchar();
    }
    return x*f;
}
inline void print(__int128 x){//输出模板
    if(x<0){
        putchar('-');
        x=-x;
    }
    if(x>9) print(x/10);
    putchar(x%10+'0');
}
```
####3.	线性筛素数
```cpp
int prime[MAXN],minprime[MAXN];//minprime相当于原来的vis，只是存放的是i的最小因子
int euler(int n){
    int c=0;
    for(int i=2;i<=n;i++){
        if(!minprime[i]) prime[++c]=i,minprime[i]=i;
        for(int j=1;j<=c&&i*prime[j]<=n;j++){
            minprime[i*prime[j]]=prime[j];
            if(i%prime[j]==0) break;
        }
    }
    return c;
}
```
####4.	exgcd求逆元
```cpp
long long exgcd(long long a,long long b,long long &x,long long &y)//扩展欧几里得算法
{
    if(b==0)
    {
        x=1,y=0;
        return a;
    }
    long long ret=exgcd(b,a%b,y,x);
    y-=a/b*x;
    return ret;
}
long long getInv(long long a,long long mod)//求a在mod下的逆元，不存在逆元返回-1
{
    long long x,y;
    long long d=exgcd(a,mod,x,y);
    return d==1?(x%mod+mod)%mod:-1;
}
```
####5.	求解ax+by=c
```cpp
template<class T> void exgcd(T a,T b,T &d,T &x,T &y){
    if(!b) {d=a;x=1;y=0;}
    else {exgcd(b,a%b,d,y,x);y-=x*(a/b);}
}
//求解二元一次方程 a*x+b*y=c,一组解为x,y,无解则返回false
template<class T> bool Solve_equation(T a,T b,T c,T &x,T& y){
    T gcd;
    exgcd(a,b,gcd,x,y);
    if(c%gcd) return false;   //无解
    T k=c/gcd;
    x*=k;y*=k;
    T xplus=b/gcd,yplus=a/gcd;
    if(xplus<0) xplus*=-1;if(yplus<0) yplus*=-1;
    //此时求出的x,y即为一组解，该方程的通解形式为X=x+t*(b/gcd),Y=y-t*(a/gcd) t为任意正整数
    //根据题目要求我们需要构造特殊解
    //x=(x%xplus+xplus)%xplus;y=c-a*x; //x的最小正整数解
    //y=(y%yplus+yplus)%yplus;x=c-b*y; //y的最小正整数解
    return true;
}
```
####6.	矩阵
```cpp
struct Matrix {
    int m, n;
    long long A[210][210];
    long long A_inv[210][210];
    Matrix(){}
    Matrix(int _m, int _n) {
        m = _m;
        n = _n;
        memset(A, 0, sizeof(A));
    }
    void init(int _m,int _n) {
        m = _m;
        n = _n;
    }
    void read(){
        for(int i=1;i<=m;i++)
            for(int j=1;j<=n;j++)
                scanf("%lld",A[i]+j);
    }
    Matrix operator+(const Matrix mat) const {
        Matrix ans{m, n};
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++) {
                ans.A[i][j] = A[i][j] + mat.A[i][j];
                if (ans.A[i][j] >= mod)
                    ans.A[i][j] -= mod;
            }
        return ans;
    }

    Matrix operator-(const Matrix mat) const {
        Matrix ans{m, n};
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++) {
                ans.A[i][j] = A[i][j] - mat.A[i][j];
                if (ans.A[i][j] < 0)
                    ans.A[i][j] += mod;
            }
        return ans;
    }

    Matrix operator*(const Matrix mat) const {
        Matrix ans{m, mat.n};
        for (int i = 0; i < m; i++)
            for (int j = 0; j < mat.n; j++)
                for (int k = 0; k < n; k++)
                    (ans.A[i][j] += A[i][k] * mat.A[k][j]) %= mod;
        return ans;
    }

    long long q_pow(long long a, long long x) {
        long long cnt = 1;
        while (x) {
            if (x & 1)
                cnt = cnt * a % mod;
            a = a * a % mod;
            x >>= 1;
        }
        return cnt;
    }

    long long inv(long long a, long long x) {
        return q_pow(a,x-2);
    }

    bool Gauss(int n) {
        memarray(A_inv, 0);
        for(int i=1;i<=n;i++)
            A_inv[i][i]=1;
        for (int i = 1; i <= n; i++) {
            for (int j = i; j <= n; j++) {
                if (A[j][i]) {
                    for (int k = 1; k <= n&& !A[i][i]; k++){
                        swap(A[i][k], A[j][k]);
                        swap(A_inv[i][k], A_inv[j][k]);
                    }
                }
            }
            if (!A[i][i]) {
                puts("No Solution");
                return false;
            }
            long long inv_ = inv(A[i][i], mod);

            for (int j = i; j <= n; j++)
                A[i][j] = A[i][j] * inv_ % mod;
            for (int j = 1; j <= n; j++)
                A_inv[i][j] = A_inv[i][j] * inv_ % mod;

            for (int j = 1; j <= n; j++) {
                if (j != i) {
                    long long m = A[j][i];
                    for (int k = i; k <= n; k++){
                        A[j][k] = (A[j][k] - m * A[i][k] % mod + mod) % mod;
                    }
                    for (int k = 1; k <= n; k++){
                        A_inv[j][k] = (A_inv[j][k] - m * A_inv[i][k] % mod + mod) % mod;
                    }
                }
            }
        }
        return true;
    }

    Matrix T() {
        Matrix ans{n, m};
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                ans.A[j][i] = A[i][j];
        return ans;
    }

    long long det() {
        long long ans = 1;
        Matrix o = (*this);
        for (int i = 0; i < n; i++) {
            if (!o.A[i][i]) {
                int pos = 0;
                for (int j = i + 1; j < n; j++)
                    if (o.A[j][i]) {
                        pos = j;
                        break;
                    }
                if (!pos)return 0;
                for (int j = i; j < n; j++)
                    swap(o.A[i][j], o.A[pos][j]);
                ans = -ans;
            }
            for (int j = 0; j < n; j++)
                if (i != j) {
                    long long tmp = o.A[j][i] * q_pow(o.A[i][i], mod - 2) % mod;
                    for (int k = i; k < n; k++) {
                        o.A[j][k] = (o.A[j][k] - o.A[i][k] * tmp) % mod + mod;
                        if (o.A[j][k] > mod)o.A[j][k] -= mod;
                    }
                }
            (ans *= o.A[i][i]) %= mod;
        }
        if (ans < 0)ans += mod;
        return ans;
    }
}matrix;
```
####7.	矩阵快速幂
```cpp
struct Matrix{
    int n;
    long long a[N][N];
    Matrix(int _n,bool identity=false){
        n=_n;
        memarray(a,0);
        if(identity){
            for(int i=1;i<=n;i++)
                a[i][i]=1;
        }
    }
    Matrix operator*(Matrix matrix){
        Matrix ret(matrix.n);
        for(int i=1;i<=n;i++)
            for(int k=1;k<=n;k++)
                for(int j=1;j<=n;j++)
                 ret.a[i][j]=(ret.a[i][j]+a[i][k]*matrix.a[k][j])%mod;
        return ret;
    }
};
long long ksm(Matrix a,long long x){
    Matrix ret(a.n,true);
    while(x){
        if(x&1)ret=ret*a;
        a=a*a;
        x>>=1;
    }
    return ?;
}
Matrix mat(N);
```
####8.	矩阵快速幂求斐波那契
//斐波那契前n项和为F(n+2)-1
//斐波那契前n项平方和和为(F(n)+F(n-1))*F(n)//几何画图证明
```
struct Matrix{
    long long a[2][2];
    Matrix operator*(Matrix matrix){
        Matrix ret={(a[0][0]*matrix.a[0][0]%mod+a[0][1]*matrix.a[1][0]%mod)%mod,
                    (a[0][0]*matrix.a[0][1]%mod+a[0][1]*matrix.a[1][1]%mod)%mod,
                    (a[1][0]*matrix.a[0][0]%mod+a[1][1]*matrix.a[1][0]%mod)%mod,
                    (a[1][0]*matrix.a[0][1]%mod+a[1][1]*matrix.a[1][1]%mod)%mod};
        return ret;
    }
};
long long ksm(Matrix a,long long x){
    Matrix ret={1,0,0,1};
    while(x){
        if(x&1)ret=ret*a;
        a=a*a;
        x>>=1;
    }
    return ret.a[0][1];
}
void solve(){
    Matrix mat={1,1,1,0};
    printf("%lld\n",ksm(mat,n));
}
```
####9.	漂浮法
```cpp
void cov(int l,int r,int idx,int now){
    if(l>r) return;
    while(idx<=n&&(l>line[idx].r||r<line[idx].l))
        idx++;
    if(idx>n){
        ans[now]+=(r-l+1);
        return;
    }
    if(l>=line[idx].l&&r<=line[idx].r)
        return;
    if(l<=line[idx].r){
        cov(line[idx].r+1,r,idx+1,now);
    }
    if(r>=line[idx].l){
        cov(l,line[idx].l-1,idx+1,now);
    }
}
```
####10.	归并排序求逆序对
```cpp
long long a[MAXN],b[MAXN];
long long merge_sort(int l,int r){
    if(l>=r) return 0;
    int mid=(l+r)>>1;
    long long ret=merge_sort(l,mid);
    ret+=merge_sort(mid+1,r);
    int i=l,j=mid+1;
    int num=l;
    while(i<=mid&&j<=r){
        if(a[i]<=a[j])
            b[num++]=a[i++];
        else{
            b[num++]=a[j++];
            ret+=mid-i+1;
        }
    }
    while(i<=mid)
        b[num++]=a[i++];
    while(j<=r)
        b[num++]=a[j++];
    for(i=l;i<=r;i++)
        a[i]=b[i];
    return ret;
}
```
####11.	RMQ求区间最值
```cpp
int dpmx[MAXN][21];
int dpmn[MAXN][21];
int lg[MAXN];
void RMQ(){
    for(int j=1;j<=lg[n];j++){
        for(int i=1;i+(1<<j)-1<=n;i++){
            dpmx[i][j]=max(dpmx[i][j-1],dpmx[i+(1<<(j-1))][j-1]);
            dpmn[i][j]=min(dpmn[i][j-1],dpmn[i+(1<<(j-1))][j-1]);
        }
    }
}
void solve(){
    RMQ();
    int l,r;
    while(m--){
        scanf("%d%d",&l,&r);
        int i=lg[r-l+1];
        int mx=max(dpmx[l][i],dpmx[r-(1<<i)+1][i]);
        int mn=min(dpmn[l][i],dpmn[r-(1<<i)+1][i]);
        printf("%d\n",mx-mn);
    }
}
void init(){
    scanf("%d%d",&n,&m);
    lg[0]=-1;
    for(int i=1;i<=n;i++){
        scanf("%d",&dpmx[i][0]);
        dpmn[i][0]=dpmx[i][0];
        lg[i]=lg[i/2]+1;
    }
}
```
####12.	中国剩余问题
```cpp
int a[10],b[10];
long long qmod(long long a,long long b,long long mod){
    long long ret=0;
    while(b){
        if(b&1)
            ret=(ret+a)%mod;
        a=(a+a)%mod;
        b>>=1;
    }
    return ret;
}
long long exgcd(long long a,long long b,long long &x,long long &y){
    if(b==0){
        x=1,y=0;
        return a;
    }
    long long gcd=exgcd(b,a%b,x,y);
    long long t=y;
    y=x-(a/b)*y;
    x=t;
    return gcd;
}
long long china(){
    long long m=a[1],ans=b[1],x,y;
    for(int i=2;i<=n;i++){
        long long c=(b[i]-ans%a[i]+a[i])%a[i];
        long long gcd=exgcd(m,a[i],x,y);
        long long ag=a[i]/gcd;
        if(c%gcd!=0) return -1;
        x=qmod(x,c/gcd,ag);
        ans+=x*m;
        m*=ag;
        ans=(ans%m+m)%m;
    }
    return ans>0?ans:ans+m;
}
void solve(){
    printf("%lld\n",china());
}
void init(){
    scanf("%d",&n);
    for(int i=1;i<=n;i++)
        scanf("%d",&a[i]);
    for(int i=1;i<=n;i++)
        scanf("%d",&b[i]);
}
------------------------------------------
int m[10],a[10];
int gcd(int a,int b){
    return b?gcd(b,a%b):a;
}
void exgcd(int a,int b,int &d,int &x,int &y){
    if(!b) d=a,x=1,y=0;
    else exgcd(b,a%b,d,y,x),y-=x*(a/b);
}
int china(){
    int x,y,d,A=a[1],M=m[1];
    for(int i=2;i<=n;i++){
        exgcd(M,m[i],d,x,y);
        if((a[i]-A)%d) return -1;
        x=(a[i]-A)/d*x%(m[i]/d);
        A+=x*M;
        M=M/d*m[i];
        A%=M;
    }
    return A>0?A:A+M;
}
void solve(){
    printf("%d\n",china());
}
void init(){
    scanf("%d",&n);
    for(int i=1;i<=n;i++)
        scanf("%d",&m[i]);
    for(int i=1;i<=n;i++)
        scanf("%d",&a[i]);
}
```
####13.	二维ST表
```cpp
const int N = 1010;
int a[N][N],n,m,len,Log[N];
int st[3][N][N][15];//0最小，1最大值 
int max(int a,int b,int c,int d){
   int mx = a;if(mx < b) mx = b;if(mx < c) mx = c;if(mx < d) mx = d;
   return mx;
}
int min(int a,int b,int c,int d){
   int mi = a;if(mi > b) mi = b;if(mi > c) mi = c;if(mi > d) mi = d;
   return mi;
}
void init(){
   for(int i = 2;i < N;i++) Log[i] = Log[i/2]+1;
   for(int i = 1;i <= n;i++)
      for(int j = 1;j <= m;j++) 
      st[0][i][j][0] = st[1][i][j][0] = a[i][j];
   for(int k = 1;k <= 12;k++){
      for(int i = 1;i + (1<<k)-1 <= n;i++){
         for(int j = 1;j + (1<<k)-1 <= m;j++){
            int t1 = st[0][i][j][k-1];
            int t2 = st[0][i+(1<<(k-1))][j][k-1];
            int t3 = st[0][i][j+(1<<(k-1))][k-1];
            int t4 = st[0][i+(1<<k-1)][j+(1<<k-1)][k-1];
            st[0][i][j][k] = min(t1,t2,t3,t4);
            t1 = st[1][i][j][k-1];
            t2 = st[1][i+(1<<(k-1))][j][k-1];
            t3 = st[1][i][j+(1<<(k-1))][k-1];
            t4 = st[1][i+(1<<k-1)][j+(1<<k-1)][k-1];
            st[1][i][j][k] = max(t1,t2,t3,t4);
         }
      }
   }
}
int ask(int r,int c,int len){
   int k = Log[len];
   int t1 = st[0][r][c][k];
   int t2 = st[0][r+len-(1<<k)][c][k];
   int t3 = st[0][r][c+len-(1<<k)][k];
   int t4 = st[0][r+len-(1<<k)][c+len-(1<<k)][k];
   int mi = min(t1,t2,t3,t4);
   t1 = st[1][r][c][k];
   t2 = st[1][r+len-(1<<k)][c][k];
   t3 = st[1][r][c+len-(1<<k)][k];
   t4 = st[1][r+len-(1<<k)][c+len-(1<<k)][k];
   int mx = max(t1,t2,t3,t4);
   //printf("%d %d\n",mx,mi);
   return mx - mi;
}
int main(){
   scanf("%d%d%d",&n,&m,&len);
   for(int i = 1;i <= n;i++)
      for(int j = 1;j <= m;j++) scanf("%d",&a[i][j]);
   init();
   int ans = 0x3f3f3f3f;
   for(int i = 1;i <= n-len+1;i++){
      for(int j = 1;j <= m-len+1;j++){
         int tmp = ask(i,j,len);
         ans = ans < tmp?ans:tmp;
      }
   }
   printf("%d\n",ans);
   return 0;
}
```
####14.	线段树
```cpp
int tree[MAXN<<2];
int tag[MAXN<<2];
void push_up(int rt){
    tree[rt]=tree[rt<<1]+tree[rt<<1|1];
}
void push_down(int l,int r,int rt){
    if(tag[rt]){
        int mid=(l+r)>>1;
        tree[rt<<1]=tag[rt]*(mid-l+1);
        tree[rt<<1|1]=tag[rt]*(r-mid);
        tag[rt<<1]=tag[rt];
        tag[rt<<1|1]=tag[rt];
        tag[rt]=0;
    }
}
void build(int l,int r,int rt){
    if(l==r){
        scanf("%d",&tree[rt]);
        return;
    }
    int mid=(l+r)>>1;
    build(lson);
    build(rson);
    push_up(rt);
}
void update(int L,int R,int val,int l,int r,int rt){
    if(L<=l&&R>=r){
        tree[rt]=val*(r-l+1);
        tag[rt]=val;
        return;
    }
    push_down(l,r,rt);
    int mid=(l+r)>>1;
    if(L<=mid)update(L,R,val,lson);
    if(R>mid)update(L,R,val,rson);
    push_up(rt);
}
int query(int L,int R,int l,int r,int rt){
    if(L<=l&&R>=r){
        return tree[rt];
    }
    push_down(l,r,rt);
    int ret=0;
    int mid=(l+r)>>1;
    if(L<=mid)ret+=query(L,R,lson);
    if(R>mid)ret+=query(L,R,rson);
    return ret;
}
```
####15.	线段树区间最值单点修改
```cpp
long long chushi[MAXN], Min[MAXN * 4], Max[MAXN * 4];//记得开4倍空间
void pushup(int rt) {
    Min[rt] = min(Min[rt<<1], Min[rt<<1|1]);
    Max[rt] = max(Max[rt<<1], Max[rt<<1|1]);
}
void build(int l, int r, int rt) {
    if (l == r){
        Min[rt] = chushi[l];
        Max[rt] = chushi[l];
        return;
    }
    int mid = (l + r)>>1;
    build(l, mid, rt<<1);
    build(mid + 1, r, rt<<1|1);
    pushup(rt);
}
int qurrymax(int x, int y, int l, int r, int rt) {
    if (x <= l && y >= r) {
        return Max[rt];
    }
    int mid = (l + r)>>1;
    int ret = -1e9;
    if (x <= mid) ret = max(ret, qurrymax(x, y, l, mid, rt<<1));//如果这个区间的左儿子和目标区间有交集那么搜索左儿子
    if (y > mid) ret = max(ret, qurrymax(x, y, mid + 1, r, rt<<1|1));//如果这个区间的右儿子和目标区间有交集那么搜索右儿子
    return ret;
}
int qurrymin(int x, int y, int l, int r, int rt) {
    if (x <= l && y >= r) {
        return Min[rt];
    }
    int mid = (l + r)>>1;
    int ret=1e9;
    if (x <= mid) ret = min(ret, qurrymin(x, y, l, mid, rt<<1));//如果这个区间的左儿子和目标区间有交集那么搜索左儿子
    if (y > mid) ret = min(ret, qurrymin(x, y, mid + 1, r, rt<<1|1));//如果这个区间的右儿子和目标区间有交集那么搜索右儿子
    return ret;
}
void update(int x, int c, int l, int r, int rt) {
    if (l == r) {
        Min[rt] = c;
        Max[rt] = c;
        return;
    }
    int mid = (l + r)>>1;
    if (x <= mid)update(x, c, l, mid, rt<<1);
    else update(x, c, mid + 1, r, rt<<1|1);
    pushup(rt);
}
void solve(){

    long long op,l,r;
    while(m--){
        scanf("%lld%lld %lld", &op, &l, &r);
        if(op==2){
            if(qurrymax(l,r,1,n,1)-qurrymin(l,r,1,n,1)==r-l){
                printf("YES\n");
            }else{
                printf("NO\n");
            }
        }
        else
            update(l, r, 1, n, 1);
    }
}
void init(){
    scanf("%d%d",&n,&m);
    for(int i=1;i<=n;i++)
        scanf("%d",chushi+i);
    build(1,n,1);
}
```
####16.	Dijkstra最短路
```cpp
struct Edge{
    int to;
    long long dis;
    bool operator<(Edge e)const{
        return dis>e.dis;
    }
};
long long dis[MAXN];
bool vis[MAXN];
void dijkstra(int s){
    for(int i=1;i<=n;i++)
        dis[i]=1e9;
    memarray(vis,false);
    priority_queue<Edge>pq;
    pq.push({s,0});
    dis[s]=0;
    while(!pq.empty()){
        Edge now=pq.top();pq.pop();
        int u=now.to;
        if(vis[u])continue;
        vis[u]=true;
        for(int v=1;v<=n;v++){
            if(vis[v])continue;
            if(dis[v]>dis[u]+Map[u][v]){
                dis[v]=dis[u]+Map[u][v];
                pq.push({v,dis[v]});
            }
        }
    }
}
```
####17.	Floyd最短路
```cpp
for(int k=1;k<=n;k++)
    for(int i=1;i<=n;i++)
        for(int j=1;j<=n;j++)
            Map[i][j]=min(Map[i][j],Map[i][k]+Map[k][j]);
18.	Jhonson全源最短路
//先用SPFA求出虚点0到n个点的距离，若有负环则退出，否则更新每条边的权值使其为正，然后对n个点进行dijkstra求dis

long long Map[3005][3005];
long long dis[MAXN];
long long h[MAXN];
int change[MAXN];
bool SPFA(int s){
    for(int i=0;i<=n;i++){
        h[i]=1e9;
        change[i]=0;
    }
    queue<int>q;
    q.push(s);
    h[s]=0;
    while(!q.empty()){
        int now=q.front();q.pop();
        for(int i=1;i<=n;i++){
            if(Map[now][i]==1e9)continue;
            if(h[i]>h[now]+Map[now][i]){
                h[i]=h[now]+Map[now][i];
                q.push(i);
                change[i]++;
                if(change[i]>n)return true;
            }
        }
    }
    return false;
}
struct Edge{
    int to;
    long long dis;
    bool operator<(Edge e)const{
        return dis>e.dis;
    }
};
vector<Edge>vec[MAXN];
bool vis[MAXN];
void dij(int s){
    for(int i=1;i<=n;i++){
        dis[i]=1e9;
        vis[i]=false;
    }
    priority_queue<Edge>pq;
    pq.push({s,0});
    dis[s]=0;
    while(!pq.empty()){
        Edge now=pq.top();pq.pop();
        int u=now.to;
        if(vis[u])continue;
        vis[u]=true;
        for(auto edge:vec[u]){
            if(vis[edge.to])continue;
            if(dis[edge.to]>dis[u]+edge.dis){
                dis[edge.to]=dis[u]+edge.dis;
                pq.push({edge.to,dis[edge.to]});
            }
        }
    }
}
void solve(){
    if(SPFA(0)){
        printf("-1\n");
        return;
    }
    for(int i=1;i<=n;i++){
        for(int j=1;j<=n;j++){
            if(i==j)continue;
            if(Map[i][j]!=1e9)
                vec[i].push_back({j,Map[i][j]+h[i]-h[j]});
        }
    }
    for(int i=1;i<=n;i++){
        dij(i);
        long long sum=0;
        for(int j=1;j<=n;j++){
            if(dis[j]==1e9)sum+=j*dis[j];
            else sum+=j*(dis[j]+h[j]-h[i]);
        }
        printf("%lld\n",sum);
    }
}
void init(){
    scanf("%d%d",&n,&m);
    for(int i=1;i<=n;i++)
        for(int j=1;j<=n;j++)
            Map[i][j]=1e9;
    for(int i=1;i<=n;i++)
        Map[i][i]=0;
    for(int i=1,u,v,w;i<=m;i++){
        scanf("%d%d%d",&u,&v,&w);
        Map[u][v]=min(Map[u][v],1LL*w);
    }
}
```
####19.	LCA
```cpp
vector<int>vec[MAXN];
int depth[MAXN];
int fa[MAXN][22];
struct LCA{
    void init(){
        depth[1]=1;
        dfs(1,0);
        for(int i=1;i<=20;i++)
            for(int j=1;j<=n;j++)
                fa[j][i]=fa[fa[j][i-1]][i-1];
    }
    void dfs(int now,int pre){
        for(auto to:vec[now]){
            if(to==pre)continue;
            depth[to]=depth[now]+1;
            fa[to][0]=now;
            dfs(to,now);
        }
    }
    int lca(int x,int y){
        if(depth[x]>depth[y])swap(x,y);
        for(int i=20;i>=0;i--)
            if(depth[fa[y][i]]>=depth[x])
                y=fa[y][i];
        if(x==y)return depth[x];
        for(int i=20;i>=0;i--)
            if(fa[x][i]!=fa[y][i]){
                x=fa[x][i];
                y=fa[y][i];
            }
        return depth[fa[x][0]];
    }
}lca;
```
####20.	Tarjan找割点
```cpp
vector<int>vec[MAXN];
int dfn[MAXN],low[MAXN];
int cuo;
bool vis[MAXN];
void dfs(int now,int fa){
    dfn[now]=low[now]=++cuo;
    int son=0;
    for(auto i:vec[now]){
        if(!dfn[i]){
            dfs(i,fa);
            low[now]=min(low[now],low[i]);
            if(low[i]>=dfn[now]&&now!=fa)
                vis[now]=true;
            if(now==fa)
                son++;
        }
        low[now]=min(low[now],dfn[i]);
    }
    if(son>=2&&now==fa)
        vis[now]=true;
}
void solve(){
    for(int i=1;i<=n;i++){
        if(dfn[i])continue;
        dfs(i,i);
    }
    int ans=0;
    for(int i=1;i<=n;i++)
        if(vis[i])ans++;
    printf("%d\n",ans);
    for(int i=1;i<=n;i++)
        if(vis[i])
            printf("%d ",i);
}
void init(){
    scanf("%d%d",&n,&m);
    for(int i=1,u,v;i<=m;i++){
        scanf("%d%d",&u,&v);
        vec[u].push_back(v);
        vec[v].push_back(u);
    }
}
```
####21.	Tarjan找强连通分量个数
```cpp
vector<int>vec[MAXN];
int dfn[MAXN],low[MAXN];
int cuo;
stack<int>s;
bool vis[MAXN];
int ans;
int num[MAXN];
void dfs(int now){
    dfn[now]=low[now]=++cuo;
    s.push(now);
    vis[now]=true;
    for(auto i:vec[now]){
        if(!dfn[i]){
            dfs(i);
            low[now]=min(low[now],low[i]);
        }else if(vis[i]){
            low[now]=min(low[now],dfn[i]);
        }
    }
    if(dfn[now]==low[now]){
        ans++;
        while(s.top()!=now){
            int t=s.top();s.pop();
            num[ans]++;
            vis[t]=false;
        }
        num[ans]++;
        vis[now]=false;
        s.pop();
    }
}
void solve(){
    for(int i=1;i<=n;i++){
        if(dfn[i])continue;
        dfs(i);
    }
    int a=0;
    for(int i=1;i<=ans;i++){
        if(num[i]>1)
            a++;
    }
    printf("%d\n",a);
}
void init(){
    scanf("%d%d",&n,&m);
    for(int i=1,u,v;i<=m;i++){
        scanf("%d%d",&u,&v);
        vec[u].push_back(v);
    }
}
```
####22.	Tarjan缩点求路径最大点权和
```cpp
vector<int>vec[MAXN];
int stac[MAXN];
int dfn[MAXN];
int low[MAXN];
bool vis[MAXN];
int w[MAXN],sd[MAXN];
int top,tim;
void Tarjan(int now){
    dfn[now]=low[now]=++tim;
    stac[++top]=now;
    vis[now]=true;
    for(auto to:vec[now]){
        if(!dfn[to]){
            Tarjan(to);
            low[now]=min(low[now],low[to]);
        }else if(vis[to]){
            low[now]=min(low[now],low[to]);
        }
    }
    if(dfn[now]==low[now]){
        int y=stac[top--];
        while(y){
            sd[y]=now;
            vis[y]=false;
            if(now==y)break;
            w[now]+=w[y];
            y=stac[top--];
        }
    }
}
vector<int>tp[MAXN];
int in[MAXN];
int dis[MAXN];
int topo(){
    queue<int>q;
    for(int i=1;i<=n;i++){
        if(sd[i]==i&&in[i]==0)
            q.push(i);
        dis[i]=w[i];
    }
    int ret=0;
    while(!q.empty()){
        int now=q.front();q.pop();
        ret=max(ret,dis[now]);
        for(auto to:tp[now]){
            in[to]--;
            dis[to]=max(dis[to],dis[now]+w[to]);
            if(in[to]==0)q.push(to);
        }
    }
    return ret;
}
void solve(){
    for(int i=1;i<=n;i++)
        if(!dfn[i])Tarjan(i);
    for(int i=1;i<=n;i++){
        for(auto to:vec[i]){
            if(sd[i]==sd[to])continue;
            tp[sd[i]].push_back(sd[to]);
            in[sd[to]]++;
        }
    }
    printf("%d\n",topo());
}
void init(){
    scanf("%d%d",&n,&m);
    for(int i=1;i<=n;i++)
        scanf("%d",w+i);
    for(int i=1,u,v;i<=m;i++){
        scanf("%d%d",&u,&v);
        vec[u].push_back(v);
    }
}
```
####23.	二分图最大匹配匈牙利算法
```cpp
int pei[505];
vector<int>vec[505];
int vis[505];
bool dfs(int nowboy,int mainboy){
    if(vis[nowboy]==mainboy)return false;
    vis[nowboy]=mainboy;
    for(auto i:vec[nowboy]){
        if(!pei[i]||dfs(pei[i],mainboy)){
            pei[i]=nowboy;
            return true;
        }
    }
    return false;
}
void solve(){
    int ans=0;
    for(int i=1;i<=n;i++){
        if(dfs(i,i))
            ans++;
    }
    printf("%d\n",ans);
}
```
####24.	拓扑排序
```cpp
int tp(){
    queue<int>q;
    for(int i=1;i<=n;i++)
        if(in[i]==0){
            q.push(i);
            ans[i]=1;
        }
    while(!q.empty()){
        int x=q.front();q.pop();
        for(auto i:to[x]){
            in[i]--;
            ans[i]=(ans[i]+ans[x])%mod;
            if(in[i]==0){
                q.push(i);
            }
        }
    }
    long long ret=0;
    for(int i=1;i<=n;i++)
        if(to[i].size()==0)
            ret=(ret+ans[i])%mod;
    return ret;
}
```
####25.	最大流
```cpp
char a[MAXN],b[MAXN];
struct node{
    int e;
    int v;
    int nxt;
}edge[MAXN];
int cnt;
int head[MAXN],dis[MAXN];
int cur[MAXN];
void add(int from,int to,int val) {
    edge[cnt].e = to;
    edge[cnt].v = val;
    edge[cnt].nxt = head[from];
    head[from] = cnt;
    cnt++;
}
bool bfs(){
    memset(dis,-1, sizeof(dis));
    dis[1]=0;
    queue<int>q;
    q.push(1);
    while(!q.empty()){
        int r=q.front();q.pop();
        for(int i=head[r];i!=-1;i=edge[i].nxt){
            int j=edge[i].e;
            if(dis[j]==-1&&edge[i].v){
                dis[j]=dis[r]+1;
                q.push(j);
            }
        }
    }
    return dis[n]!=-1;
}
int dfs(int u,int flo){
    if(u==n) return flo;
    int delta=flo;
    for(int i=cur[u];i!=-1;i=edge[i].nxt){
        cur[u]=edge[i].nxt;
        int v=edge[i].e;
        if(dis[v]==(dis[u]+1)&&edge[i].v>0){
            int d=dfs(v,min(delta,edge[i].v));
            edge[i].v-=d;
            edge[i+1].v+=d;
            delta-=d;
            if(delta==0)break;
        }
    }
    return flo-delta;
}
int dini() {//求最大流量
    int ans=0;
    while(bfs()){
        for(int i=1;i<=n;i++)
            cur[i]=head[i];
        ans+=dfs(1,INF);
    }
    return ans;
}
void solve(){
    int a=dini();
    if(a){
        printf("%d %d\n",a,(1LL*x+a-1)/a);
    }else{
        printf("Orz Ni Jinan Saint Cow!");
    }
}
void init(){
    scanf("%d%d%d",&n,&m,&x);
    memset(head,-1, sizeof(head));
    for(int i=1,u,v,w;i<=m;i++){
        scanf("%d%d%d",&u,&v,&w);
        add(u,v,w);
        add(v,u,0);
    }
}
```
####26.	最小费用最大流
```cpp
const int MAXN=2e2+10;
int INF=0x3f3f3f3f;
struct Edge
{
    int from,to;
    long long cap,flow,cost;
    Edge(int u,int v,int ca,int f,int co):from(u),to(v),cap(ca),flow(f),cost(co){};
};
struct MCMF {
    int n,m,s,t;
    vector<Edge> edges;
    vector<int> G[MAXN];
    int inq[MAXN];//是否在队列中
    long long d[MAXN];//距离
    long long p[MAXN];//上一条弧
    long long a[MAXN];//可改进量
//
    struct rec{
        long long flow,cost,percost;
    };
    vector<rec> r;//增广路
    void init(int n) {
        this->n = n;
        for(int i = 0; i <= n; i++) G[i].clear();
        edges.clear();
        r.clear();
        r.push_back({0,0,0});
    }
    void add_edge(int from,int to, int cap,int cost) {
        edges.push_back((Edge){from,to,cap,0,cost});
        edges.push_back((Edge){to,from,0,0,-cost});
        m = edges.size();
        G[from].push_back(m-2);
        G[to].push_back(m-1);
    }

    bool spfa(int s,int t,int &flow,int &cost) {
        for(int i = 0; i <= n; i++) d[i] = INF;
        for(int i = 0; i <= n; i++) inq[i] = 0;
        d[s] = 0; inq[s] = 1; p[s] = 0; a[s] = INF;
        queue<int> Q;
        Q.push(s);
        while(!Q.empty()) {
            int u = Q.front(); Q.pop();
            inq[u] = 0;
            for(int i = 0; i < G[u].size(); i++) {
                Edge & e = edges[G[u][i]];
                if(e.cap > e.flow && d[e.to] > d[u] + e.cost) {
                    d[e.to] = d[u] + e.cost;
                    p[e.to] = G[u][i];
                    a[e.to] = min(a[u], e.cap - e.flow);
                    if(!inq[e.to]) {Q.push(e.to); inq[e.to] = 1;}
                }
            }
        }
        if(d[t] == INF) return false;
        flow += a[t];
        cost += d[t] * a[t];
        r.push_back({flow,cost,d[t]});
        int u = t;
        while(u != s) {
            edges[p[u]].flow += a[t];
            edges[p[u]^1].flow -= a[t];
            u = edges[p[u]].from;
        }
        return true;
    }

    int Mincost(int s,int t) {
        int cost = 0;
        int flow=0;
        while(1){
            if(!spfa(s,t,flow,cost))
                break;
        }
        return flow;
    }
} mcmf;

void init(){
    mcmf.init(n);
    for(int i=1,u,v,w,f;i<=m;i++){
        scanf("%d%d%d",&u,&v,&w);
        mcmf.add_edge(u,v,1,w);
    }
}
```
####27.	欧拉路径
```cpp
void fleury(int start) {
    int u = start;
    top = 0; path.clear();
    S[top++] = u;
    while (top) {
        u = S[--top];
        if (!G[u].empty())
            DFS(u);
        else path.push_back(u);
    }
}
```
####28.	KMP
```cpp
// s[]是长文本，p[]是模式串，n是s的长度，m是p的长度
/*求模式串的Next数组：*/
for (int i = 2, j = 0; i <= m; i ++ )
{
    while (j && p[i] != p[j + 1]) j = ne[j];
    if (p[i] == p[j + 1]) j ++ ;
    ne[i] = j;
}

// 匹配
for (int i = 1, j = 0; i <= n; i ++ )
{
    while (j && s[i] != p[j + 1]) j = ne[j];
    if (s[i] == p[j + 1]) j ++ ;
    if (j == m)
    {
        j = ne[j];
        // 匹配成功后的逻辑
    }
}
```
####29.	几何
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
        return min(min(disPointToSeg(v.s),disPointToSeg(v.e)),min(v.disPointToSeg(s),v.disPointToSeg(e)));
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
####30.	直线与凸多边形交点
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
####31.	凸包Andrew
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
####32.	凸包Graham
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
####33.	判断线段相交
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
####34.	半平面交
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
####35.	四面体内切球
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
####36.	平面最近点对
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
####37.	旋转卡壳求最大三角形面积
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
####38.	最小圆覆盖
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
####39.	模拟退火求矩阵内与所有点最小距离最大
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
####40.	矩阵面积凸包最小矩形覆盖
```cpp
double rota(int n){
    Vector v;
    int up=2;
    int left=2;
    int right=2;
    double ret=1e9;
    for(int down=1;down<=n;down++) {

        v = stk[down] - stk[down + 1];
        while (sgn(v ^ (stk[up + 1] - stk[up])) <= 0)
            up = up % n + 1;

//        right=down;
        while(cmp((stk[down+1]-stk[down])*(stk[right]-stk[down]),
                (stk[down+1]-stk[down])*(stk[right+1]-stk[down]))<=0)
            right=right%n+1;

        left=up%n+1;
        while(cmp((stk[left]-stk[down+1])*(stk[down]-stk[down+1]),
                (stk[left+1]-stk[down+1])*(stk[down]-stk[down+1]))<=0)
            left=left%n+1;


        double h=((stk[down+1]-stk[down])^(stk[up]-stk[down]));
        double w=(((stk[left]-stk[down])*(stk[down]-stk[down+1]))+((stk[down+1]-stk[down])*(stk[right]-stk[down])));
        ret=min(ret,h*w/(stk[down+1]-stk[down]).norm());
    }
    return ret;
}
void solve(){
    int len=Andrew();
    printf("%.0f\n",rota(len));
}
```
####41.	稳定凸包
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
####42.	大整数模板
```cpp
struct BigInteget{ // 非负整数范围运行
    int digit[N];
    int length;
    BigInteget();
    BigInteget(int x);
    BigInteget(string str);
    BigInteget(const BigInteget& b);
    BigInteget operator =(int x);
    BigInteget operator =(string str);
    BigInteget operator =(const BigInteget& b);
    bool operator <=(const BigInteget& b);
    bool operator ==(const BigInteget& b);
    bool operator >(const BigInteget& b);
    BigInteget operator +(const BigInteget& b);
    BigInteget operator -(const BigInteget& b);
    BigInteget operator *(const BigInteget& b);
    BigInteget operator /(const BigInteget& b);
    BigInteget operator %(const BigInteget& b);
    friend istream& operator>>(istream& in, BigInteget& x);
    friend ostream& operator<<(ostream& out, const BigInteget& x);
    void show();
};
istream& operator >> (istream& in, BigInteget& x){
    string str;
    in >> str;
    x = str;
    return in;
}
ostream& operator << (ostream& out, const BigInteget& x){
    for(int i = x.length - 1; i >= 0; i--)
        out << x.digit[i];
    return out;
}
BigInteget::BigInteget(){
    memset(digit, 0, sizeof(digit));
    length = 0;
}
BigInteget::BigInteget(int x){
    memset(digit, 0, sizeof(digit));
    length = 0;
    if(x == 0)
        digit[length++] = x;
    while(x) {
        digit[length++] = x % 10;
        x /= 10;
    }
}
BigInteget::BigInteget(string str){
    memset(digit, 0, sizeof(digit));
    length = str.size();
    for(int i = 0; i < length; i++)
        digit[i] = str[length - i - 1] - '0';
}
BigInteget::BigInteget(const BigInteget& b){
    memset(digit, 0, sizeof(digit));
    length = b.length;
    for(int i = 0; i < length; i++)
        digit[i] = b.digit[i];
}
BigInteget BigInteget::operator = (int x){
    memset(digit, 0, sizeof(digit));
    length = 0;
    if(x == 0)
        digit[length++] = x;
    while(x){
        digit[length++] = x % 10;
        x /= 10;
    }
    return *this;
}
BigInteget BigInteget::operator = (string str){
    memset(digit, 0, sizeof(digit));
    length = str.size();
    for(int i = 0; i < length; i++)
        digit[i] = str[length - i - 1] - '0';
    return *this;
}
BigInteget BigInteget::operator = (const BigInteget& b){
    memset(digit, 0, sizeof(digit));
    length = b.length;
    for(int i = 0; i < length; i++)
        digit[i] = b.digit[i];
    return *this;
}
bool BigInteget::operator > (const BigInteget &b) {
    if(length > b.length) return true;
    else if(b.length > length) return false;
    else {
        for(int i = length - 1; i >= 0; i--) {
            if(digit[i] == b.digit[i]) continue;
            else return digit[i] > b.digit[i];
        }
    }
    return false;
}
bool BigInteget::operator <= (const BigInteget& b) {
    if(length < b.length) return true;
    else if (b.length < length) return false;
    else {
        for(int i = length - 1; i >= 0; i--) {
            if(digit[i] == b.digit[i]) continue;
            else return digit[i] < b.digit[i];
        }
    }
    return true;
}

bool BigInteget::operator == (const BigInteget& b){
    if(length != b.length) return false;
    else{
        for(int i = length -1; i >= 0; i--){
            if(digit[i] != b.digit[i])
                return false;
        }
    }
    return true;
}
BigInteget BigInteget::operator + (const BigInteget& b){
    BigInteget answer;
    int carry = 0;
    for(int i = 0; i < length || i < b.length; i++){
        int current = carry + digit[i] + b.digit[i];
        carry = current /10;
        answer.digit[answer.length++] = current % 10;
    }
    if(carry){
        answer.digit[answer.length++] = carry;
    }
    return answer;
}
BigInteget BigInteget::operator - (const BigInteget& b){
    BigInteget answer;
    int carry = 0;
    for(int i = 0; i < length; i++){
        int current = digit[i] - b.digit[i] - carry;
        if(current < 0) {
            current += 10;
            carry = 1;
        } else carry  = 0;
        answer.digit[answer.length++] = current;
    }
    while(answer.digit[answer.length - 1] == 0 && answer.length > 1){//书上在这里写得是answer.digit[answer.length]
        answer.length--;
    }
    return answer;
}
BigInteget BigInteget::operator * (const BigInteget& b){
    BigInteget answer;
    answer.length = length + b.length;
    for(int i = 0; i < length; i++){
        for(int j = 0; j < b.length; j++)
            answer.digit[i+j] += digit[i] * b.digit[j];
    }
    for(int i = 0; i < answer.length; i++){
        answer.digit[i+1] += answer.digit[i] / 10;
        answer.digit[i] %= 10;
    }
    while(answer.digit[answer.length - 1] == 0 && answer.length > 1){ //书上在这里写得是answer.digit[answer.length]
        answer.length--;
    }
    return answer;
}

BigInteget BigInteget::operator / (const BigInteget& b){
    BigInteget answer;
    answer.length = length;
    BigInteget remainder = 0;
    BigInteget temp = b;
    for(int i = length - 1; i >= 0; i--){
        if(!(remainder.length == 1 && remainder.digit[0] == 0)){
            for(int j = remainder.length -1; j >= 0; j--)
                remainder.digit[j + 1] = remainder.digit[j];
            remainder.length++;
        }
        remainder.digit[0] = digit[i];
        while(temp <= remainder){
            remainder = remainder - temp;
            answer.digit[i]++;
        }
    }
    while(answer.digit[answer.length - 1] == 0 && answer.length > 1){//书上在这里写得是answer.digit[answer.length]
        answer.length--;
    }
    return answer;
}
BigInteget BigInteget::operator % (const BigInteget &b){
    BigInteget remainder = 0;
    BigInteget temp = b;
    for(int i = length - 1; i >= 0; i--) {
        if(!(remainder.length == 1 && remainder.digit[0] == 0)){
            for(int j = remainder.length - 1; j >= 0; j--)
                remainder.digit[j + 1] = remainder.digit[j];
            remainder.length++;
        }
        remainder.digit[0] = digit[i];
        while(temp <= remainder){
            remainder = remainder - temp;
        }
    }
    return remainder;
}
void BigInteget::show(){
    for(int i = length - 1; i >= 0; i--)
        printf("%d", digit[i]);
}
```
