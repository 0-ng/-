
### [1.	KMP](#1)
### [2.	马拉车](#2)
### [3.	AC自动机](#3)
### [4. 字符串转移n\sqrt(n)](#4)
### [5. trie树](#5)

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

struct AutoACMachine{
    int trie[MAXN][30];
    int fail[MAXN];
    int num[MAXN];
    int cnt;
    void init(){
        memarray(trie,0);
        memarray(fail,0);
        memarray(num,0);
        cnt=0;
    }
    void insert(char *str){
        int len=strlen(str);
        int now=0;
        rep(i,0,len-1){
            if(!trie[now][str[i]-'a']){
                trie[now][str[i]-'a']=++cnt;
            }
            now=trie[now][str[i]-'a'];
        }
        num[now]++;
    }
    void getFail(){
        queue<int>q;
        rep(i,0,25){
            if(trie[0][i]){
                fail[trie[0][i]]=0;
                q.push(trie[0][i]);
            }
        }
        while(!q.empty()){
            int now=q.front();q.pop();
            rep(i,0,25){
                if(trie[now][i]){
                    fail[trie[now][i]]=trie[fail[now]][i];
                    q.push(trie[now][i]);
                }else{
                    trie[now][i]=trie[fail[now]][i];
                }
            }
        }
    }
    int query(char *str){
        int len=strlen(str);
        int now=0;
        int ret=0;
        rep(i,0,len-1){
            now=trie[now][str[i]-'a'];
            for(int j=now;j&&~num[j];j=fail[j]){
                ret+=num[j];
                num[j]=-1;
            }
        }
        return ret;
    }
}ac;
void solve(){
    ac.getFail();
    scanf("%s",str);
    printf("%d\n",ac.query(str));
}
void init(){
    scanf("%d",&n);
    ac.init();
    rep(i,1,n){
        scanf("%s",str);
        ac.insert(str);
    }
}
```

<span id="4"><h4>4.	字符串转移n\sqrt(n)</h4></span>
```cpp
#include <bits/stdc++.h>
#include <ext/rope>
using namespace std;
using namespace __gnu_cxx;
int main()
{
    int n,m,p,s;
    cin>>n>>m;
    rope<int>r;
    for(int i=0;i<n;i++)
           r.push_back(i+1);
    while(m--)
    {
        scanf("%d%d",&p,&s);
        r=r.substr(p-1,s)+r.substr(0,p-1)+r.substr(p+s-1,n-(p+s-1));
        //r=r.substr(1,2);
    }
    for(int i=0;i<r.size();i++)
        printf("%d ",r[i]);
}
```

<span id="5"><h4>5.	trie树</h4></span>
```cpp
struct Trie{
    int tree[MAXN][30];
    int num[MAXN*30];
    int cnt;
    void insert(char *str){
        int len=strlen(str);
        int now=0;
        rep(i,0,len-1){
            if(!tree[now][str[i]-'a']){
                tree[now][str[i]-'a']=++cnt;
            }
            now=tree[now][str[i]-'a'];
            num[now]++;
        }
    }
    int check(char *str){
        int len=strlen(str);
        int now=0;
        rep(i,0,len-1){
            if(!tree[now][str[i]-'a']){
                return 0;
            }
            now=tree[now][str[i]-'a'];
        }
        return num[now];
    }
}trie;
```
