### [0.	啥也行](#0)
### [1.	基础模板](#1)
### [2.	快读快写](#2)
### [3.	漂浮法](#3)
### [4.	归并排序求逆序对](#4)
### [5.	大整数模板](#5)
### [6. 快速约瑟夫环求最后一个](#6)
### [7. 高速1\sqrt](#7)
### [8. 麻将\sqrt](#8)
### [9. 0/1背包bitset](#9)
### [10. codeblock复制粘贴](#10)
### [11. bitset黑科技](#11)
### [12. __黑科技](#12)
### [13. 输出黑科技](#13)
### [14. 枚举子集](#14)
### [15. 宏](#15)

---------------------
<span id="0"><h4>0.	啥也行</h4></span>
最小化s到t之间路径长度最大值与最小值的比值可以考虑枚举下界  推广同理
一堆点建图可能只需要把各个方向相邻的点连起来即可
一个交点由两条直线相交得到，一条直线由两个点组成



一、int scanf(const char * restrict format,…)
（一）函数介绍
返回值位成功读入的数据项数
format是一个或者多个 {%[*] [width] [{h | l | I64 | L}]type | ’ ’ | ‘\t’ | ‘\n’ | 非%符号}
用户在键盘上输入的数据首先进入输入缓冲区，scanf标准函数从输入缓冲区里获得数据并记录到存储区里，先进入的数据必须被先处理。如果format位%d，而缓冲区中是3.14，则只会将3存入，并且，下一个%d也不能正确的读取，因为.14无法被%d识别所以一直卡在缓冲区

（二）常用简单的正则表达
scanf("%ns", str)；
表示读取长度为n的字符串
输入：123456
输出str：123 ( 以scanf("%3s", str);为例 )

scanf("%[a-z]", str);
表示读取a-z的小写字母，出现非a-z的小写字母，立即停止读取。
输入：abcd123
输出str：abcd

scanf("%*[a-z]%s", str); 注：该语句一定要加%s，因为前边有*
%*[ ]表示过滤掉满足括号内条件的字符串 %[ ]表示读取
输入：abcd123
输出str：123

scanf("%[^0-9]", str);
^表示非，^0-9表示非0-9的一切字符，因此是遇到0-9就立即停止读取。
输入：abcd123
输出str：abcd

PS：%[^\n] 表示 读取回车符以前的所有字符，常用于读取含空格的字符串。
　　%[^ ] 表示 读取空格符以前的所有字符。

scanf("%*[^\n]%d", &num);
表示过滤掉回车以前所有的字符，并将回车符的下一个数字赋给num。
输入：abcd\n12
输出num：12

PS：%*[^  ]表示过滤空格以前的所有字符

scanf("%*[^\n]");
scanf("%*c");
可以使用该组合语句清空缓冲区中\n及其以前的内容，用在其他的scanf();语句之后，切记不可把两句合二为一，否则当输入缓冲区只有一个\n时，合成一个之后，第一个%匹配失败（没能匹配到一个非\n字符用来跳过）就直接结束匹配了，然后就执行下一条语句了，\n还是留在缓冲区里



<span id="1"><h4>1.	基础模板</h4></span>
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
<span id="2"><h4>2.	快读快写</h4></span>
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
<span id="3"><h4>3.	漂浮法</h4></span>
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
<span id="4"><h4>4.	归并排序求逆序对</h4></span>
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
<span id="5"><h4>5.	大整数模板</h4></span>
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
    while(answer.digit[answer.length - 1] == 0 && answer.length > 1){
    //书上在这里写得是answer.digit[answer.length]
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
    while(answer.digit[answer.length - 1] == 0 && answer.length > 1){ 
    //书上在这里写得是answer.digit[answer.length]
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
    while(answer.digit[answer.length - 1] == 0 && answer.length > 1){
    //书上在这里写得是answer.digit[answer.length]
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
<span id="6"><h4>6. 快速约瑟夫环求最后一个</h4></span>
```cpp
long long joser(long long n,long long k,long long s=1)
{
    if(k==1) return (n-1+s)%n;
    long long ans=0;
    //ans=(ans+k)%i
    for(long long i=2;i<=n;)
    {
        if(ans+k<i) //快速跳跃
        {
            long long leap;
            if((i-ans-1)%(k-1)==0) leap=(i-ans-1)/(k-1)-1;
            else leap=(i-ans-1)/(k-1);
            if(i+leap>n) return ((ans+(n+1-i)*k)+s)%n;
            i+=leap;
            ans+=leap*k;
        }
        else
        {
            ans=(ans+k)%i;
            i++;
        }
    }
    return (ans+s)%n;
}
void solve(){
    long long ans=joser(n,m);
    if(ans==0)printf("%lld\n",n);
    else printf("%lld\n",ans);
}
```
<span id="7"><h4>7. 高速1\sqrt</h4></span>
```cpp
float Q_rsqrt( float number ){
   int i;
   float x2, y;
   const float threehalfs = 1.5F;

   x2 = number * 0.5F;
   y   = number;
   i   = * ( int * ) &y;   // evil doubleing point bit level hacking
   i   = 0x5f3759df - ( i >> 1 ); // what the fuck?
   y   = * ( float * ) &i;
   y   = y * ( threehalfs - ( x2 * y * y ) ); 
   // 1st iteration
   y   = y * ( threehalfs - ( x2 * y * y ) ); 
   // 2nd iteration, this can be removed

   return y;
}
```


<span id="8"><h4>8. 麻将</h4></span>
input
5
1w2w3w4b5b6b7s8s9s1b1b1z2z6z
1w2w3w4b5b6b7s8s9s1b1b2z2z6z
1w2w3w4b5b6b7s8s9s1b1b2z2z2z
1b2b3b4b5b6b2s4s5s5s5s6s7s8s
1b1b1b2b3b4b5b6b7b8b9b9b9b1s

output
0
1
6z 1b2z
Tsumo!
4
2s 3s4s6s9s
4s 2s
5s 3s
8s 3s
4
2b 1s
5b 1s
8b 1s
1s 1b2b3b4b5b6b7b8b9b
```cpp
#include <cstdio>
#include <cmath>
#include <vector>
#include <algorithm>
#include <cstring>
#include <iostream>
#include <stack>
#include <ctime>
#include <set>
#define lson l,mid,rt<<1
#define rson mid+1,r,rt<<1|1
#define pb push_back
#define memarray(array,val) memset(array,val,sizeof(array))
#define rep(i,a,n) for(int i=a;i<=n;i++)
#define f(i,a,n) for(int i=a;i<=n;i++)
#define per(i,n,a) for(int i=n;i>=a;i--)
#define ff(i,n,a) for(int i=n;i>=a;i--)
using namespace std;
const int INF=0x3f3f3f3f;
const double EPS=1e-6;
const int mod=1000000007;
const int MAXN=1e6+10;
const double PI=acos(-1);
int n,m;
char str[50];
int wan[10];
int suo[10];
int bing[10];
int zi[10];
int tywan[10];
int tysuo[10];
int tybing[10];
int tyzi[10];
bool check2(){
    memcpy(tywan,wan,sizeof wan);
    memcpy(tysuo,suo,sizeof wan);
    memcpy(tybing,bing,sizeof wan);
    memcpy(tyzi,zi,sizeof wan);
    rep(i,1,9){
        if(!tywan[i])continue;
        if(tywan[i]<0)return false;
        while(tywan[i]>=3)tywan[i]-=3;
        while(tywan[i]){
            if(i+2>9)return false;
            tywan[i]--;
            tywan[i+1]--;
            tywan[i+2]--;
        }
    }
    rep(i,1,9){
        if(!tybing[i])continue;
        if(tybing[i]<0)return false;
        while(tybing[i]>=3)tybing[i]-=3;
        while(tybing[i]){
            if(i+2>9)return false;
            tybing[i]--;
            tybing[i+1]--;
            tybing[i+2]--;
        }
    }
    rep(i,1,9){
        if(!tysuo[i])continue;
        if(tysuo[i]<0)return false;
        while(tysuo[i]>=3)tysuo[i]-=3;
        while(tysuo[i]){
            if(i+2>9)return false;
            tysuo[i]--;
            tysuo[i+1]--;
            tysuo[i+2]--;
        }
    }
    rep(i,1,7){
        if(tyzi[i]%3)return false;
    }
    return true;
}
bool check(){
    rep(i,1,9){
        if(wan[i]>=2){
            wan[i]-=2;
            if(check2()){
                wan[i]+=2;
                return true;
            }
            wan[i]+=2;
        }
    }
    rep(i,1,9){
        if(bing[i]>=2){
            bing[i]-=2;
            if(check2()){
                bing[i]+=2;
                return true;
            }
            bing[i]+=2;
        }
    }
    rep(i,1,9){
        if(suo[i]>=2){
            suo[i]-=2;
            if(check2()){
                suo[i]+=2;
                return true;
            }
            suo[i]+=2;
        }
    }
    rep(i,1,7){
        if(zi[i]>=2){
            zi[i]-=2;
            if(check2()){
                zi[i]+=2;
                return true;
            }
            zi[i]+=2;
        }
    }
    return false;
}
vector<pair<int,int>>vec[5][10];
void go(int t,int id){
    rep(i,1,9){
        wan[i]++;
        if(check()){
            vec[t][id].pb({1,i});
        }
        wan[i]--;
    }
    rep(i,1,9){
        bing[i]++;
        if(check()){
            vec[t][id].pb({2,i});
        }
        bing[i]--;
    }
    rep(i,1,9){
        suo[i]++;
        if(check()){
            vec[t][id].pb({3,i});
        }
        suo[i]--;
    }
    rep(i,1,7){
        zi[i]++;
        if(check()){
            vec[t][id].pb({4,i});
        }
        zi[i]--;
    }
}
void print(int i,int j){
    char ch;
    if(i==1)ch='w';
    if(i==2)ch='b';
    if(i==3)ch='s';
    if(i==4)ch='z';
    printf("%d%c",j,ch);
}
void solve(){
    if(check()){
        printf("Tsumo!\n");
        return;
    }
    rep(i,1,9){
        if(!wan[i])continue;
        wan[i]--;
        go(1,i);
        wan[i]++;
    }
    rep(i,1,9){
        if(!bing[i])continue;
        bing[i]--;
        go(2,i);
        bing[i]++;
    }
    rep(i,1,9){
        if(!suo[i])continue;
        suo[i]--;
        go(3,i);
        suo[i]++;
    }
    rep(i,1,7){
        if(!zi[i])continue;
        zi[i]--;
        go(4,i);
        zi[i]++;
    }
    int ans=0;
    rep(i,1,4){
        int j=i==4?7:9;
        rep(k,1,j){
            if(vec[i][k].size()){
                ans++;
            }
        }
    }
    printf("%d\n",ans);
    rep(i,1,4){
        int j=i==4?7:9;
        rep(k,1,j){
            if(vec[i][k].size()){
                print(i,k);
                printf(" ");
                for(auto p:vec[i][k]){
                    print(p.first,p.second);
                }
                printf("\n");
            }
        }
    }
}
void init(){
    memarray(bing,0);
    memarray(suo,0);
    memarray(wan,0);
    memarray(zi,0);
    rep(i,1,4){
        int j=i==4?7:9;
        rep(k,1,j){
            vec[i][k].clear();
        }
    }
    scanf("%s",str+1);
    n=strlen(str+1);
    rep(i,1,n){
        int num=str[i]-'0';
        if(str[i+1]=='w'){
            wan[num]++;
        }else if(str[i+1]=='s'){
            suo[num]++;
        }else if(str[i+1]=='b'){
            bing[num]++;
        }else{
            zi[num]++;
        }
        i++;
    }
//    rep(i,1,9){
//        printf("%d ",wan[i]);
//    }
//    printf("\n");
//    rep(i,1,9){
//        printf("%d ",suo[i]);
//    }
//    printf("\n");
//    rep(i,1,9){
//        printf("%d ",bing[i]);
//    }
//    printf("\n");
//    rep(i,1,6){
//        printf("%d ",zi[i]);
//    }
//    printf("\n");
}
int main()
{
    freopen("data.in","r",stdin);
    int T=1;
//    int T;
    scanf("%d",&T);
    while(T--)
    {
        init();
        solve();
    }
    return 0;
}

```


<span id="9"><h4>9. 0/1背包bitset</h4></span>
```cpp
bitset<2000001>d;
d[0]=1;
cin>>n;
for(int&x:A){
    cin>>x;
}
for(int x:A){
    d|=d<<x;
}
```



<span id="10"><h4>10. codeblock复制粘贴</h4></span>
```cpp
ubuntu下codeblocks解决运行窗口的复制粘贴数据问题
1.）在Code::blocks中，点击Settings -> Environment。 
2.）将Terminal to launch console programs选项改成gnome-terminal -t $TITLE -x。（原来是xterm -T $TITLE -e）。 
    gnome-terminal和xterm的参数表示方法不一样。 
注意：在Terminal中复制是Ctrl+Shift+C（注意不要按Ctrl+C，Ctrl+C是强制退出），粘贴是Ctrl+Shift+V。
```

<span id="11"><h4>11. bitset黑科技</h4></span>
```cpp
int Q;
int n;
bitset<MAXN> dp;
int main() {
//    freopen("data.in", "r", stdin);
    cin >> n;
    dp[0] = 1;
    f(i, 1, n){
        int x;
        cin >> x;
        dp |= dp << x;
    }
    cout << dp.count() - 1 << endl;
    dp[0] = 0;
    for(int i = dp._Find_first(); i != dp.size(); i = dp._Find_next(i)){
        cout << i << ' ';
    }
 
    return 0;
}
```

<span id="12"><h4>12. __黑科技</h4></span>
```cpp
•int __builtin_ffs (unsigned int x)
返回x的最后一位1的是从后向前第几位，比如7368（1110011001000）返回4。
•int __builtin_clz (unsigned int x)
返回前导的0的个数。
•int __builtin_ctz (unsigned int x)
返回后面的0个个数，和__builtin_clz相对。
•int __builtin_popcount (unsigned int x)
返回二进制表示中1的个数。
•int __builtin_parity (unsigned int x)
返回x的奇偶校验位，也就是x的1的个数模2的结果。

此外，这些函数都有相应的usigned long和usigned long long版本，只需要在函数名后面加上l或ll就可以了，比如int __builtin_clzll。


状压dp拿出每一个1
 int tmp = (1 << n) - 1 - S;
    while(tmp){
        int j = __builtin_ffs(tmp) - 1;
        tmp ^= (1 << j);
    }
    
```


<span id="13"><h4>13. 输出黑科技</h4></span>
```cpp
 
template <typename... T>
void write(T &&...args) {
    ((cout << args), ...);
}

```


<span id="14"><h4>14. 枚举子集</h4></span>
```cpp
    int S = (1 << m) - 1 - sta;
    int sub = S;
    do{
        if(!st[sta | sub]) (v += dfs(cur + 1, sub)) %= MOD;
        sub = (sub - 1) & S;
    }while(sub != S);
```

<span id="15"><h4>15. 宏</h4></span>
```cpp
#define deb(x) cout << #x << " is " << x << "\n"

// bitwise ops
// also see https://gcc.gnu.org/onlinedocs/gcc/Other-Builtins.html
constexpr int pct(int x) { return __builtin_popcount(x); } // # of bits set
constexpr int bits(int x) { // assert(x >= 0); // make C++11 compatible until USACO updates ...
	return x == 0 ? 0 : 31-__builtin_clz(x); } // floor(log2(x))
constexpr int p2(int x) { return 1<<x; }
constexpr int msk2(int x) { return p2(x)-1; }
```
