### [0.	啥也行](#0)
### [1.	基础模板](#1)
### [2.	快读快写](#2)
### [3.	漂浮法](#3)
### [4.	归并排序求逆序对](#4)
### [5.	大整数模板](#5)
### [6. 快速约瑟夫环求最后一个](#6)
### [7. 高速1\sqrt](#7)

---------------------
<span id="0"><h4>0.	啥也行</h4></span>
最小化s到t之间路径长度最大值与最小值的比值可以考虑枚举下界  推广同理
一堆点建图可能只需要把各个方向相邻的点连起来即可
一个交点由两条直线相交得到，一条直线由两个点组成


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
