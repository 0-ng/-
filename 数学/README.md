
### [1.大素数模板](#1)
### [2.1-n因子和或个数](#2)
### [3.	线性筛素数](#3)
### [4.	exgcd求逆元](#4)
### [5.	求解ax+by=c](#5)
### [6.	矩阵](#6)
### [7.	矩阵快速幂](#7)
### [8.	矩阵快速幂求斐波那契](#8)
### [9.	中国剩余问题](#9)
### [10. min25求前n素数和](#10)
### [11. FFT](#11)
### [12. 判断大素数](#12)


<span id="1"><h4>1. 大素数模板</h4></span>
```cpp
#include <bits/stdc++.h>
#define ll long long
using namespace std;
ll f[340000],g[340000],n;
void init()
{
    ll i,j,m;
    for(m=1; m*m<=n; ++m)f[m]=n/m-1;
    for(i=1; i<=m; ++i)g[i]=i-1;
    for(i=2; i<=m; ++i)
    {
        if(g[i]==g[i-1])continue;
        for(j=1; j<=min(m-1,n/i/i); ++j)
        {
            if(i*j<m)f[j]-=f[i*j]-g[i-1];
            else f[j]-=g[n/i/j]-g[i-1];
        }
        for(j=m; j>=i*i; --j)g[j]-=g[j/i]-g[i-1];
    }
}
int main()
{
    while(scanf("%I64d",&n)!=EOF)
    {
        init();
        cout<<f[1]<<endl;
    }
    return 0;
}
----------------------------------------------------
#include<cstdio>
#include<cmath>
using namespace std;
#define LL long long
const int N = 5e6 + 2;
bool np[N];
int prime[N], pi[N];
int getprime()
{
    int cnt = 0;
    np[0] = np[1] = true;
    pi[0] = pi[1] = 0;
    for(int i = 2; i < N; ++i)
    {
        if(!np[i]) prime[++cnt] = i;
        pi[i] = cnt;
        for(int j = 1; j <= cnt && i * prime[j] < N; ++j)
        {
            np[i * prime[j]] = true;
            if(i % prime[j] == 0)   break;
        }
    }
    return cnt;
}
const int M = 7;
const int PM = 2 * 3 * 5 * 7 * 11 * 13 * 17;
int phi[PM + 1][M + 1], sz[M + 1];
void init()
{
    getprime();
    sz[0] = 1;
    for(int i = 0; i <= PM; ++i)  phi[i][0] = i;
    for(int i = 1; i <= M; ++i)
    {
        sz[i] = prime[i] * sz[i - 1];
        for(int j = 1; j <= PM; ++j) phi[j][i] = phi[j][i - 1] - phi[j / prime[i]][i - 1];
    }
}
int sqrt2(LL x)
{
    LL r = (LL)sqrt(x - 0.1);
    while(r * r <= x)   ++r;
    return int(r - 1);
}
int sqrt3(LL x)
{
    LL r = (LL)cbrt(x - 0.1);
    while(r * r * r <= x)   ++r;
    return int(r - 1);
}
LL getphi(LL x, int s)
{
    if(s == 0)  return x;
    if(s <= M)  return phi[x % sz[s]][s] + (x / sz[s]) * phi[sz[s]][s];
    if(x <= prime[s]*prime[s])   return pi[x] - s + 1;
    if(x <= prime[s]*prime[s]*prime[s] && x < N)
    {
        int s2x = pi[sqrt2(x)];
        LL ans = pi[x] - (s2x + s - 2) * (s2x - s + 1) / 2;
        for(int i = s + 1; i <= s2x; ++i) ans += pi[x / prime[i]];
        return ans;
    }
    return getphi(x, s - 1) - getphi(x / prime[s], s - 1);
}
LL getpi(LL x)
{
    if(x < N)   return pi[x];
    LL ans = getphi(x, pi[sqrt3(x)]) + pi[sqrt3(x)] - 1;
    for(int i = pi[sqrt3(x)] + 1, ed = pi[sqrt2(x)]; i <= ed; ++i) ans -= getpi(x / prime[i]) - i + 1;
    return ans;
}
LL lehmer_pi(LL x)
{
    if(x < N)   return pi[x];
    int a = (int)lehmer_pi(sqrt2(sqrt2(x)));
    int b = (int)lehmer_pi(sqrt2(x));
    int c = (int)lehmer_pi(sqrt3(x));
    LL sum = getphi(x, a) +(LL)(b + a - 2) * (b - a + 1) / 2;
    for (int i = a + 1; i <= b; i++)
    {
        LL w = x / prime[i];
        sum -= lehmer_pi(w);
        if (i > c) continue;
        LL lim = lehmer_pi(sqrt2(w));
        for (int j = i; j <= lim; j++) sum -= lehmer_pi(w / prime[j]) - (j - 1);
    }
    return sum;
}
int main()
{
    init();
    LL n;
    while(~scanf("%lld",&n))
    {
        printf("%lld\n",lehmer_pi(n));
    }
    return 0;
}
```
<span id="2"><h4>2. 1-n因子和或个数</h4></span>
```cpp
long long getSum(long long n){//求1-n的因子和（不包括1和本身）
    long long ans=0;
    for(long long i=2;i<=sqrt(n);i++)
    {
        ans+=(n/i-1)*i;
        if(n/i>sqrt(n))
        {
            long long  q=n/i;long long p=sqrt(n)+1;
            ans+=(q-p+1)*(q+p)/2;
        }
    }
    return ans;
}

////////////////////////
long long getSum(long long n){//求1-n的因子个数和
    if(n==0)return 0;
    long long ans=0;
    for(long long i=2;i<=sqrt(n);i++)
    {
        ans+=n/i-1;
        if(n/i>sqrt(n))
        {
            long long  q=n/i;
            long long p=sqrt(n)+1;
            ans+=(q-p+1);
        }
    }
    return ans+n+n-1;//（1和自身）
}
```
<span id="3"><h4>3.	线性筛素数</h4></span>
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
<span id="4"><h4>4.	exgcd求逆元</h4></span>
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
<span id="5"><h4>5.	求解ax+by=c</h4></span>
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
<span id="6"><h4>6.	矩阵</h4></span>
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
<span id="7"><h4>7.	矩阵快速幂</h4></span>
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
<span id="8"><h4>8.	矩阵快速幂求斐波那契</h4></span>
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
<span id="9"><h4>9.	中国剩余问题</h4></span>
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
<span id="10"><h4>10. min25求前n素数和</h4></span>
```cpp

################################素数和↓
const int N = 1000010;
typedef long long LL;
namespace Min25 {
    int prime[N], id1[N], id2[N], flag[N], ncnt, m;
    LL g[N], sum[N], a[N], T, n;
    inline int ID(LL x) {
        return x <= T ? id1[x] : id2[n / x];
    }
    inline LL calc(LL x) {
        return x * (x + 1) / 2 - 1;
    }
    inline LL f(LL x) {
        return x;
    }
    inline void init() {
        T = sqrt(n + 0.5);
        for (int i = 2; i <= T; i++) {
            if (!flag[i]) prime[++ncnt] = i, sum[ncnt] = sum[ncnt - 1] + i;
            for (int j = 1; j <= ncnt && i * prime[j] <= T; j++) {
                flag[i * prime[j]] = 1;
                if (i % prime[j] == 0) break;
            }
        }
        for (LL l = 1; l <= n; l = n / (n / l) + 1) {
            a[++m] = n / l;
            if (a[m] <= T) id1[a[m]] = m; else id2[n / a[m]] = m;
            g[m] = calc(a[m]);
        }
        for (int i = 1; i <= ncnt; i++)
            for (int j = 1; j <= m && (LL)prime[i] * prime[i] <= a[j]; j++)
                g[j] = g[j] - (LL)prime[i] * (g[ID(a[j] / prime[i])] - sum[i - 1]);
    }
    inline LL solve(LL x) {
        if (x <= 1) return x;
        return n = x, init(), g[ID(n)];
    }
}
int main() {
    LL n; scanf("%lld", &n);
    printf("%lld\n", Min25::solve(n));
}



################################素数个数↓
long long n;
int len;
int lim;
int lenp;
long long a[7000000];
long long g[7000000];
inline int get_index(long long x)
{
    if(x<=lim)
        return x;
    else
        return len-n/x+1;
}
int main()
{
    cin>>n;
    lim=sqrt(n);
    for(long long i=1;i<=n;i=a[len]+1)
    {
        a[++len]=n/(n/i);
        g[len]=a[len]-1;
    }
    for(int i=2;i<=lim;++i)
    {
        if(g[i]!=g[i-1])
        {
            ++lenp;
            for(int j=len;a[j]>=1ll*i*i;--j)
                g[j]=g[j]-g[get_index(a[j]/i)]+lenp-1;
        }

    }
    cout<<g[len];
    return 0;
}
```

<span id="11"><h4>11. FFT</h4></span>
```cpp
const double PI=acos(-1.0);
const int MAXN=3e6+10;
struct Complex{
    double x,y;
    Complex(double _x=0.0,double _y=0.0){
        x=_x;y=_y;
    }
    Complex operator-(const Complex &b)const{
        return Complex(x-b.x,y-b.y);
    }
    Complex operator+(const Complex &b)const{
        return Complex(x+b.x,y+b.y);
    }
    Complex operator*(const Complex &b)const{
        return Complex(x*b.x-y*b.y,x*b.y+y*b.x);
    }
};
void change(Complex y[],int len){
    for(int i=1,j=len/2,k;i<len-1;i++){
        if(i<j)swap(y[i],y[j]);
        k=len/2;
        while(j>=k){
            j-=k;
            k/=2;
        }
        if(j<k)j+=k;
    }
}
void fft(Complex y[],int len,int on){
    change(y,len);
    for(int h=2;h<=len;h<<=1){
        Complex wn(cos(-on*2*PI/h),sin(-on*2*PI/h));
        for(int j=0;j<len;j+=h){
            Complex w(1,0);
            for(int k=j;k<j+h/2;k++){
                Complex u=y[k];
                Complex t=w*y[k+h/2];
                y[k]=u+t;
                y[k+h/2]=u-t;
                w=w*wn;
            }
        }
    }
    if(on==-1){
        for(int i=0;i<len;i++)
            y[i].x/=len;
    }
}
char a[MAXN],b[MAXN];
Complex x1[MAXN],x2[MAXN];
int sum[MAXN];
int main(){
    scanf("%s",a);
    scanf("%s",b);
    int len1=strlen(a);
    int len2=strlen(b);
    int len=1;
    while(len<len1*2||len<len2*2)len<<=1;
    for(int i=0;i<len1;i++)x1[i]=Complex(a[len1-i-1]-'0',0);
    for(int i=len1;i<len;i++)x1[i]=Complex(0,0);

    for(int i=0;i<len2;i++)x2[i]=Complex(b[len2-i-1]-'0',0);
    for(int i=len2;i<len;i++)x2[i]=Complex(0,0);

    fft(x1,len,1);fft(x2,len,1);
    for(int i=0;i<len;i++)x1[i]=x1[i]*x2[i];
    fft(x1,len,-1);
    for(int i=0;i<len;i++)
        sum[i]=(int)(x1[i].x+0.5);
    for(int i=0;i<len;i++){
        sum[i+1]+=sum[i]/10;
        sum[i]%=10;
    }
    len=len1+len2+1;
    while(sum[len]<=0&&len>0)len--;
    for(int i=len;i>=0;i--)printf("%d",sum[i]);
    printf("\n");
    return 0 ;
}
-----------------------------------------

-------------------------------------
/*
1 2
1 2
1 2 1
*/
/*
1 4 5 2
*/
const double PI=acos(-1.0);
struct Complex{
    double x,y;
    Complex(double _x=0.0,double _y=0.0){
        x=_x;y=_y;
    }
    Complex operator-(const Complex &b)const{
        return Complex(x-b.x,y-b.y);
    }
    Complex operator+(const Complex &b)const{
        return Complex(x+b.x,y+b.y);
    }
    Complex operator*(const Complex &b)const{
        return Complex(x*b.x-y*b.y,x*b.y+y*b.x);
    }
}a[MAXN],b[MAXN];
int N,M;
int l,r[MAXN];
void fast_fast_tle(Complex *A,int limit,int type)
{
    for(int i=0;i<limit;i++)
        if(i<r[i]) swap(A[i],A[r[i]]);//求出要迭代的序列
    for(int mid=1;mid<limit;mid<<=1)//待合并区间的中点
    {
        Complex Wn( cos(PI/mid) , type*sin(PI/mid) ); //单位根
        for(int R=mid<<1,j=0;j<limit;j+=R)//R是区间的右端点，j表示前已经到哪个位置了
        {
            Complex w(1,0);//幂
            for(int k=0;k<mid;k++,w=w*Wn)//枚举左半部分
            {
                 Complex x=A[j+k],y=w*A[j+mid+k];//蝴蝶效应
                A[j+k]=x+y;
                A[j+mid+k]=x-y;
            }
        }
    }
}
int main()
{
    scanf("%d%d",&N,&M);
    for(int i=0;i<=N;i++) scanf("%lf",&a[i].x);
    for(int i=0;i<=M;i++) scanf("%lf",&b[i].x);
    int limit=1;
    while(limit<=N+M) limit<<=1,l++;
    for(int i=0;i<limit;i++)
        r[i]= ( r[i>>1]>>1 )| ( (i&1)<<(l-1) ) ;
    // 在原序列中 i 与 i/2 的关系是 ： i可以看做是i/2的二进制上的每一位左移一位得来
    // 那么在反转后的数组中就需要右移一位，同时特殊处理一下复数
    fast_fast_tle(a,limit,1);
    fast_fast_tle(b,limit,1);
    for(int i=0;i<=limit;i++) a[i]=a[i]*b[i];
    fast_fast_tle(a,limit,-1);
    for(int i=0;i<=N+M;i++)
        printf("%d ",(int)(a[i].x/limit+0.5));
    return 0;
}
```

<span id="12"><h4>12. 判断大素数</h4></span>
```cpp
/*枚举两个奇数加数，直接套用米勒罗宾素数测试方法模板判断是否是素数。模板地址在另一个博客里。*/
typedef unsigned long long ull;
typedef unsigned long long ULL;

ULL prime[6] = {2, 3, 5, 233, 331};
ULL qmul(ULL x, ULL y, ULL mod) { // 乘法防止溢出， 如果p * p不爆LL的话可以直接乘； O(1)乘法或者转化成二进制加法
    return (x * y - (long long)(x / (long double)mod * y + 1e-3) *mod + mod) % mod;
    /*
    LL ret = 0;
    while(y) {
        if(y & 1)
            ret = (ret + x) % mod;
        x = x * 2 % mod;
        y >>= 1;
    }
    return ret;
    */
}
ULL ksm(ULL a, ULL n, ULL mod) {
    ULL ret = 1;
    while(n) {
        if(n & 1) ret = qmul(ret, a, mod);
        a = qmul(a, a, mod);
        n >>= 1;
    }
    return ret;
}
bool Miller_Rabin(ULL p) {
    if(p < 2) return 0;
    if(p != 2 && p % 2 == 0) return 0;
    ULL s = p - 1;
    while(! (s & 1)) s >>= 1;
    for(int i = 0; i < 5; ++i) {
        if(p == prime[i]) return 1;
        ULL t = s, m = ksm(prime[i], s, p);
        while(t != p - 1 && m != 1 && m != p - 1) {
            m = qmul(m, m, p);
            t <<= 1;
        }
        if(m != p - 1 && !(t & 1)) return 0;
    }
    return 1;
}
```
