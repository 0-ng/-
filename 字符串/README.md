
### [1.	KMP](#1)
### [2.	马拉车](#2)
### [3.	AC自动机](#3)

<span id="1"><h4>1.	KMP</h4></span>
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

------------------2021/4/10--------------------
struct KMP{
    char s[MAXN],t[MAXN];
    int Next[MAXN];
    void s_read(){
        scanf("%s",s+1);
    }
    void t_read(){
        scanf("%s",t+1);
    }
    void get_next(){
        int len=strlen(t+1);
        int j=0;
        Next[1]=j;
        rep(i,2,len){
            while(j&&t[i]!=t[j+1])
                j=Next[j];
            if(t[i]==t[j+1])j++;
            Next[i]=j;
        }
        rep(i,1,len){//[1-i]最小循环节
            printf("%d ",i-Next[i]);
        }
        printf("\n");
    }
    void run(){
        int len1=strlen(s+1);
        int len2=strlen(t+1);
        int match_num=0;
        int j=0;
        rep(i,1,len1){
            while(j&&s[i]!=t[j+1])
                j=Next[j];
            if(s[i]==t[j+1])j++;
            if(j==len2){
                match_num++;
                j=Next[j];
            }
        }
        printf("%d\n",match_num);
    }
}kmp;
void solve(){
    kmp.get_next();
    kmp.run();
}
void init(){
    kmp.s_read();
    kmp.t_read();
}


---------
使用第一行字符串的字母组成n个字符的字符串，字符串内不能包含第二行的字符串
用kmp把第i个字母转到j状态记录，快速幂不计到n即可
Sample Input
3
3
ab
ab
4
acd
ca
5
ab
aaa
Sample Output
Case 1: 4
Case 2: 55
Case 3: 24

struct KMP{
    char s[MAXN];
    char t[MAXN];
    int Next[MAXN];
    int len;
    void s_read(){
        scanf("%s",s+1);
    }
    void t_read(){
        scanf("%s",t+1);
    }
    void get_next(char *str){
        int j=0;
        int len= strlen(str+1);
        rep(i,2,len){
            while(j&&str[i]!=str[j+1])
                j=Next[j];
            if(str[i]==str[j+1])
                j++;
            Next[i]=j;
        }
    }
    struct Matrix{
        int n;
        unsigned mat[55][55];
        Matrix(){
            memarray(mat,0);
        }
        Matrix operator*(const Matrix matrix)const{
            Matrix ret;
            ret.n=n;
            rep(i,0,n-1){
                rep(k,0,n-1){
                    rep(j,0,n-1){
                        ret.mat[i][j]+=mat[i][k]*matrix.mat[k][j];
                    }
                }
            }
            return ret;
        }
        void print(){
            printf("\nmatrix:\n");
            rep(i,0,n){
                rep(j,0,n){
                    printf("%u ",mat[i][j]);
                }
                printf("\n");
            }
        }
    };
    void ksm(Matrix a,int x){
        Matrix ret;
        len= strlen(t+1);
        rep(i,0,len){
            ret.mat[i][i]=1;
        }
        ret.n=len;
        a.n=len;
        while(x){
            if(x&1)
                ret=ret*a;
            a=a*a;
            x>>=1;
        }
        unsigned int sum=0;
        rep(i,0,len-1){
            sum+=ret.mat[0][i];
        }
        printf("%u\n",sum);
    }
    void run(){
        get_next(t);
        int len1= strlen(s+1);
        int len2= strlen(t+1);
        Matrix ret;
        ret.mat[0][1]=1;
        ret.mat[0][0]=len1-1;
        rep(i,1,len2){
            rep(j,1,len1){
                int tmp=i;
                while(tmp&&s[j]!=t[tmp+1])
                    tmp=Next[tmp];
                if(s[j]==t[tmp+1])
                    tmp++;
                ret.mat[i][tmp]++;
            }
        }
        ret.n=len2;
        ksm(ret,n);
    }
}kmp;

void solve(){
    kmp.run();
}
void init(){
    scanf("%d",&n);
    kmp.s_read();
    kmp.t_read();
}
```

<span id="2"><h4>2.	马拉车</h4></span>
```cpp
char beginn[MAXN],endd[MAXN];
int len[MAXN];
int  tralation()
{
	int gg=strlen(beginn);
	int noww=0;
	endd[0]='&';
	for(int i=0;i<gg;i++)
	{
		endd[++noww]='#';
		endd[++noww]=beginn[i];
	}
	endd[++noww]='#';
	return noww;
}
int Manacher(int total)
{
	int maxx = 0;
	int mx = 0,id = 0; 
	//id表示的最大回文的中心点是哪一个 而mx表示的是最大回文的中心点的最远的边界是哪一个
	for(int i =1;i<=total;i++)
	{
		if(i<mx) //如果此刻 i的点比mx还要小的话说明 在mx-i处这边是回文 然后在比较一下
			len[i] = min(mx - i,len[2*id-i]); 
		//因为2*id - i 和 i 他们是相对于id对称的 所以说吧 就是要比较mx-i和len[2*id-i] ;
		
		else len[i] = 1; //如果此刻i的点比边界还要大的话 那就需要从一开始加了
		while(endd[i+len[i]]==endd[i-len[i]])
			len[i]++;
		if(i+len[i]>mx)
		{
			mx = i+len[i];
			id = i;
		}
		/*能到结尾最长的回文串 murderforajarof=6*/
		if(i+len[i]-1==total){
            maxx = max(maxx,len[i]);
		}
		/**/

		/*最长的回文串*/
		maxx = max(maxx,len[i]);
		/**/
	}
	return maxx-1;
}
void solve(){

    int total = tralation();
    int ans = Manacher(total);
    printf("%d\n",n-ans);

}
```


<span id="3"><h4>3.	AC自动机</h4></span>
```cpp
struct Aho_Corasick_Automaton{
    int c[MAXN][26],val[MAXN],fail[MAXN],cnt;
    void ins(char *s){
        int len=strlen(s);int now=0;
        for(int i=0;i<len;i++){
            int v=s[i]-'a';
            if(!c[now][v])
                c[now][v]=++cnt;
            now=c[now][v];
        }
        val[now]++;
    }
    void build(){
        queue<int>q;
        for(int i=0;i<26;i++){
            if(c[0][i]){
                fail[c[0][i]]=0;
                q.push(c[0][i]);
            }
        }
        while(!q.empty()){
            int u=q.front();q.pop();
            for(int i=0;i<26;i++){
                if(c[u][i]){
                    fail[c[u][i]]=c[fail[u]][i];
                    q.push(c[u][i]);
                }else{
                    c[u][i]=c[fail[u]][i];
                }
            }
        }
    }
    int query(char *s){
        int len=strlen(s);
        int now=0,ans=0;
        for(int i=0;i<len;i++){
            now=c[now][s[i]-'a'];
            for(int t=now;t&&~val[t];t=fail[t]){
                ans+=val[t];
                val[t]=-1;
            }
        }
        return ans;
    }
}AC;
char p[MAXN];
void solve() {
    int ans=AC.query(p);
    printf("%d\n",ans);
}
void init() {
    scanf("%d",&n);
    for(int i=1;i<=n;i++)
        scanf("%s",p),AC.ins(p);
    AC.build();
    scanf("%s",p);
}
```
