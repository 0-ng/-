dfs(数的最后若干位,各种限制条件,当前第几位)
	if 最后一位
    	return 各种限制条件下的返回值
    局部变量 ct=当前位的数字
    局部变量 sum=0;
    for i=0 to ct-1
    	sum+=当前位取i时一定无无限制的合法状态数
        sum+=当前位取i时满足当前限制的合法状态数
    根据ct更新限制条件 不再满足则return sum
    return sum+dfs(当前位后的若干位,更新后的限制条件,下一位)

slv(当前数)
	if(只有一位) return 对应的贡献
    局部变量 ct;
    for ct=可能最高位 to 1
    	if 当前位有数字 break
    局部变量 nw=当前位数字
    局部变量 sum=0
    for i=1 to nw-1
    	sum+=当前位取i后合法情况任意取的贡献
    for i=1 to ct-1
    	for j=1 to 9
        	sum+=第i位取j后合法情况任意取的贡献
    sum+=dfs(去掉第一位后的若干位,限制条件,第二位)
    return sum

main
	预处理当前位取i的各种条件各种限制的贡献
    读入 L R
    --L
    输出 slv(R)-slv(L)
    return 0
    
    
   
   

题目描述
人们选择手机号码时都希望号码好记、吉利。比如号码中含有几位相邻的相同数字、不含谐音不吉利的数字等。手机运营商在发行新号码时也会考虑这些因素，从号段中选取含有某些特征的号码单独出售。为了便于前期规划，运营商希望开发一个工具来自动统计号段中满足特征的号码数量。

工具需要检测的号码特征有两个：号码中要出现至少 3 个相邻的相同数字；号码中不能同时出现 8 和 4。号码必须同时包含两个特征才满足条件。满足条件的号码例如：13000988721、23333333333、14444101000。而不满足条件的号码例如：1015400080、10010012022。

手机号码一定是 11 位数，前不含前导的 0。工具接收两个数 L 和 R，自动统计出 [L,R] 区间内所有满足条件的号码数量。L 和 R 也是 11 位的手机号码。

输入格式
输入文件内容只有一行，为空格分隔的 22 个正整数 L,RL,R。

输出格式
输出文件内容只有一行，为 11 个整数，表示满足条件的手机号数量。
```cpp
long long l,r;
long long dp[15][15][15][2][2][2][2];
int num[15];
long long f(int last,int pre1,int pre2,int c,int d,bool _4,bool _8){
    if(_4&&_8)return 0;
    if(last<=0) return c;
    if(~dp[last][pre1][pre2][c][d][_4][_8])return dp[last][pre1][pre2][c][d][_4][_8];
    long long ret=0;
    int lim=!d?num[last]:9;
    rep(i,0,lim){
        ret+=f(last-1,i,pre1,c||(i==pre2&&i==pre1),d||(i<lim),_4||(i==4),_8||(i==8));
    }
    return dp[last][pre1][pre2][c][d][_4][_8]=ret;
}
long long cal(long long x){
    if(x<1e10)return 0;
    memarray(dp,-1);
    int len=11;
    rep(i,1,len){
        num[i]=x%10;
        x/=10;
    }
    long long ret=0;
    rep(i,1,num[len]){
        ret+=f(10,i,0,0,i<num[len],i==4,i==8);
    }
    return ret;
}
void solve(){
    printf("%lld\n",cal(r)-cal(l));
}
void init(){
    scanf("%lld%lld",&l,&r);
    l--;
}
```
