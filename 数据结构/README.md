
### [1. RMQ](#1)
### [2. 二维ST表max min](#2)
### [3. 带权并查集](#3)
### [4. 线段树](#4)
### [5. 线段树区间最大最小单点修改](#5)
### [6. 树状数组](#6)
### [7. 二维树状数组求矩阵异或](#7)
### [8. 线段树区间赋值单点查询](#8)
### [9. CDQ分治](#9)
### [10. 主席树区间第k小](#10)
### [11. ST表](#11)
### [12. 莫队](#12)
### [13. 红黑树黑科技](#13)

<span id="1"><h4>1.	RMQ</h4></span>
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
<span id="2"><h4>2. 二维ST表max min</h4></span>
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
<span id="3"><h4>3. 带权并查集</h4></span>
```cpp
int fa[MAXN];
int dis[MAXN];
int Find(int x){
    if(fa[x]==x)return x;
    int t=Find(fa[x]);
    dis[x]+=dis[fa[x]];
    return fa[x]=t;
}
void Union(int x,int y,int w){
    int xx=Find(x);
    int yy=Find(y);
    if(xx!=yy){
        dis[yy]=dis[x]-dis[y]-w;
        fa[yy]=xx;
    }
}
void solve(){
    int flag=0;
    for(int i=1,u,v,w;i<=m;i++){
        scanf("%d%d%d",&u,&v,&w);
        if(Find(u)!=Find(v+1)){
            Union(u,v+1,w);
        }else if(dis[u]-dis[v+1]!=w){
            flag=1;
        }
    }
    if(flag)printf("false\n");
    else printf("true\n");
}
void init(){
    scanf("%d%d",&n,&m);
    for(int i=1;i<=n+1;i++){
        fa[i]=i;
        dis[i]=0;
    }
}
```
<span id="4"><h4>4. 线段树</h4></span>
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
//////////////////////////////////////////////////////////
//树链剖分线段树区间加区间乘区间赋值求区间立方和
long long W[MAXN];
vector<int>vec[MAXN];
int son[MAXN],bigson[MAXN];
int fa[MAXN],depth[MAXN];
void dfs1(int now,int pre){
    son[now]=1;
    int mx=0;
    for(auto to:vec[now]){
        if(to==pre)continue;
        fa[to]=now;
        depth[to]=depth[now]+1;
        dfs1(to,now);
        son[now]+=son[to];
        if(son[to]>mx){
            bigson[now]=to;
            mx=son[to];
        }
    }
}
int cnt;
int dfn[MAXN],top[MAXN],pos[MAXN];
void dfs2(int now,int tp){
    dfn[now]=++cnt;
    top[now]=tp;
    pos[cnt]=now;
    if(son[now]==1)return;
    dfs2(bigson[now],tp);
    for(auto to:vec[now]){
        if(to==fa[now]||to==bigson[now])continue;
        dfs2(to,to);
    }
}
long long sum1[MAXN<<2],sum2[MAXN<<2],sum3[MAXN<<2];
long long lazyMul[MAXN<<2],lazyAdd[MAXN<<2];
void solve(int l,int r,int rt,long long a,long long b){
    if(a!=1){
        long long w1=a;
        long long w2=(a*a)%mod;
        long long w3=(a*a)%mod*a%mod;

        sum1[rt]=(sum1[rt]*w1)%mod;
        sum2[rt]=(sum2[rt]*w2)%mod;
        sum3[rt]=(sum3[rt]*w3)%mod;

        lazyMul[rt]=(lazyMul[rt]*a)%mod;
        lazyAdd[rt]=(lazyAdd[rt]*a)%mod;
    }
    if(b!=0){
        long long w1=b;
        long long w2=(b*b)%mod;
        long long w3=(b*b)%mod*b%mod;

        sum3[rt]=(sum3[rt]+3*w1*sum2[rt]%mod)%mod;
        sum3[rt]=(sum3[rt]+3*w2*sum1[rt]%mod)%mod;
        sum3[rt]=(sum3[rt]+(r-l+1)*w3%mod)%mod;

        sum2[rt]=(sum2[rt]+2*w1*sum1[rt]%mod)%mod;
        sum2[rt]=(sum2[rt]+(r-l+1)*w2)%mod;

        sum1[rt]=(sum1[rt]+(r-l+1)*w1)%mod;

        lazyAdd[rt]=(lazyAdd[rt]+b)%mod;

    }
}
void push_up(int rt){
    sum1[rt]=(sum1[rt<<1]+sum1[rt<<1|1])%mod;
    sum2[rt]=(sum2[rt<<1]+sum2[rt<<1|1])%mod;
    sum3[rt]=(sum3[rt<<1]+sum3[rt<<1|1])%mod;
}
void push_down(int l,int r,int rt){
    int x=lazyMul[rt];
    int y=lazyAdd[rt];

    int mid=(l+r)>>1;
    solve(l,mid,rt<<1,x,y);
    solve(mid+1,r,rt<<1|1,x,y);

    lazyMul[rt]=1;
    lazyAdd[rt]=0;
}
void build(int l,int r,int rt){
    lazyMul[rt]=1;
    lazyAdd[rt]=0;
    if(l==r){
        long long w=W[pos[l]];
        sum1[rt]=w;
        sum2[rt]=(w*w)%mod;
        sum3[rt]=(w*w)%mod*w%mod;
        return;
    }
    int mid=(l+r)>>1;
    build(l,mid,rt<<1);
    build(mid+1,r,rt<<1|1);
    push_up(rt);
}
long long getsum3(int l,int r,int L,int R,int rt){
    if(L<=l&&r<=R){
        return sum3[rt];
    }
    push_down(l,r,rt);
    int mid=(l+r)>>1;
    long long ret=0;
    if(L<=mid) ret+=getsum3(l,mid,L,R,rt<<1);
    if(R>mid) ret+=getsum3(mid+1,r,L,R,rt<<1|1);
    push_up(rt);
    return ret%mod;
}
void update(int l,int r,int L,int R,long long a,long long b,int rt){
    if(L<=l&&r<=R){
        solve(l,r,rt,a,b);
        return;
    }
    int mid=(l+r)>>1;
    push_down(l,r,rt);
    if(L<=mid)update(l,mid,L,R,a,b,rt<<1);
    if(R>mid)update(mid+1,r,L,R,a,b,rt<<1|1);
    push_up(rt);
}
void update(int x,int y,long long a,long long b){//ax+b
    while(top[x]!=top[y]){
        if(depth[top[x]]<depth[top[y]])swap(x,y);
        update(1,n,dfn[top[x]],dfn[x],a,b,1);
        x=fa[top[x]];
    }
    if(depth[x]<depth[y])swap(x,y);
    update(1,n,dfn[y],dfn[x],a,b,1);
}
long long query(int x,int y){
    long long ret=0;
    while(top[x]!=top[y]){
        if(depth[top[x]]<depth[top[y]])swap(x,y);
        ret=(ret+getsum3(1,n,dfn[top[x]],dfn[x],1))%mod;
        x=fa[top[x]];
    }
    if(depth[x]<depth[y])swap(x,y);
    ret=(ret+getsum3(1,n,dfn[y],dfn[x],1))%mod;
    return ret;
}
void solve(){
    dfs1(1,0);
    dfs2(1,1);
    build(1,n,1);
    int q;scanf("%d",&q);
    int op,u,v;
    long long w;
    while(q--){
        scanf("%d",&op);
        if(op!=4){
            scanf("%d%d%lld",&u,&v,&w);
            if(op==1) update(u,v,0,w);
            else if(op==2) update(u,v,1,w);
            else update(u,v,w,0);
        }else{
            scanf("%d%d",&u,&v);
            printf("%lld\n",query(u,v));
        }
    }
}
void init(){
    scanf("%d",&n);
    cnt=0;
    for(int i=1;i<=n;i++)
        vec[i].clear();
    for(int i=1,u,v;i<n;i++){
        scanf("%d%d",&u,&v);
        vec[u].pb(v);
        vec[v].pb(u);
    }
    for(int i=1;i<=n;i++)
        scanf("%lld",W+i);
}

```
<span id="5"><h4>5. 线段树区间最大最小单点修改</h4></span>
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
    if (x <= mid) ret = max(ret, qurrymax(x, y, l, mid, rt<<1));
    //如果这个区间的左儿子和目标区间有交集那么搜索左儿子
    if (y > mid) ret = max(ret, qurrymax(x, y, mid + 1, r, rt<<1|1));
    //如果这个区间的右儿子和目标区间有交集那么搜索右儿子
    return ret;
}
int qurrymin(int x, int y, int l, int r, int rt) {
    if (x <= l && y >= r) {
        return Min[rt];
    }
    int mid = (l + r)>>1;
    int ret=1e9;
    if (x <= mid) ret = min(ret, qurrymin(x, y, l, mid, rt<<1));
    //如果这个区间的左儿子和目标区间有交集那么搜索左儿子
    if (y > mid) ret = min(ret, qurrymin(x, y, mid + 1, r, rt<<1|1));
    //如果这个区间的右儿子和目标区间有交集那么搜索右儿子
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

<span id="6"><h4>6. 树状数组</h4></span>
```cpp
struct BIT{//单点修改区间查询
    int n;
    long long c[MAXN];
    void update(int x,long long val){//单点修改
        for(x;x<=n;x+=lowbit(x))
            c[x]+=val;
    }
    long long query(int x){//区间查询
        int ret=0;
        for(x;x;x-=lowbit(x))
            ret+=c[x];
        return ret;
    }
    /*
    int sum=query(idx);//前idx个数的和
     
     差分[1-idx]区间加1
    update(1,1);????
    update(idx+1,-1);????
     */
};
struct BIT{//区间修改区间查询
    int n;
    long long sum1[MAXN];    //(D[1] + D[2] + ... + D[n])
    long long sum2[MAXN];    //(1*D[1] + 2*D[2] + ... + n*D[n])
    void updata(int i,long long k){//修改前i个数
        int x = i;    //因为x不变，所以得先保存i值
        while(i <= n){
            sum1[i] += k;
            sum2[i] += k * (x-1);
            i += lowbit(i);
        }
    }
    long long getsum(int i){        //求前缀和 前i个数的和
        long long res = 0, x = i;
        while(i > 0){
            res += x * sum1[i] - sum2[i];
            i -= lowbit(i);
        }
        return res;
    }
    /*
    cin>>n;
    for(int i = 1; i <= n; i++){
        cin>>a[i];
        updata(i,a[i] - a[i-1]);   //输入初值的时候，也相当于更新了值
    }
    //[x,y]区间内加上k
    updata(x,k);    //A[x] - A[x-1]增加k
    updata(y+1,-k);        //A[y+1] - A[y]减少k

    //求[x,y]区间和
    int sum = getsum(y) - getsum(x-1);
     */
};
```

<span id="7"><h4>7.	二维树状数组求矩阵异或</h4></span>
//初始全0，区间修改区间查询
```cpp
long long tree[2][2][MAXN][MAXN];
int lowbit(int x){return x&-x;}
void add(int x,int y,long long v){
    for(int i=x;i<=n;i+=lowbit(i))
        for(int j=y;j<=n;j+=lowbit(j))
            tree[x&1][y&1][i][j]^=v;
}
long long sum(int x,int y){
    long long ret=0;
    for(int i=x;i>0;i-=lowbit(i))
        for(int j=y;j>0;j-=lowbit(j))
            ret^=tree[x&1][y&1][i][j];
    return ret;
}
void solve(){
    int ty,x_2,y_2,x_1,y_1;
    long long v;
    while(m--){
        scanf("%d%d%d%d%d",&ty,&x_1,&y_1,&x_2,&y_2);
        if(ty==1){
            long long ans=0;
            ans^=sum(x_2,y_2);
            ans^=sum(x_1-1,y_2);
            ans^=sum(x_2,y_1-1);
            ans^=sum(x_1-1,y_1-1);
            printf("%lld\n",ans);
        }else{
            scanf("%lld",&v);
            add(x_2+1,y_2+1,v);
            add(x_1,y_2+1,v);
            add(x_2+1,y_1,v);
            add(x_1,y_1,v);
        }
    }
}
void init(){
    scanf("%d%d",&n,&m);
}
```



<span id="8"><h4>8.    线段树区间赋值单点查询</h4></span>
```cpp
struct SGT
{
    int sum[MAXN<<2],lazy[MAXN<<2];
    void pushup(int rt){
        sum[rt]=sum[rt<<1]+sum[rt<<1|1];
    }
    void pushdown(int l,int r,int rt){
        if(lazy[rt]!=-1){
            int mid=(l+r)>>1;
            lazy[rt<<1]=lazy[rt];
            lazy[rt<<1|1]=lazy[rt];
            sum[rt<<1]=lazy[rt]*(mid-l+1);
            sum[rt<<1|1]=lazy[rt]*(r-mid);
            lazy[rt]=-1;
        }
    }
    void build(int l,int r,int rt){
        lazy[rt]=-1;
        sum[rt]=0;
        if(l==r)return;
        int mid=(l+r)>>1;
        build(l,mid,rt<<1);
        build(mid+1,r,rt<<1|1);
        pushup(rt);
    }
    void update(int l,int r,int rt,int L,int R,int c)
    {
        if(L<=l&&r<=R){
            sum[rt]=c*(r-l+1);
            lazy[rt]=c;
            return;
        }
        pushdown(l,r,rt);
        int mid=(l+r)>>1;
        if(L<=mid)update(l,mid,rt<<1,L,R,c);
        if(R>mid) update(mid+1,r,rt<<1|1,L,R,c);
        pushup(rt);
    }
    int query(int l,int r,int rt,int p){
        if(l==r)return sum[rt];
        pushdown(l,r,rt);
        int mid=(l+r)>>1;
        if(p<=mid)return query(l,mid,rt<<1,p);
        else return query(mid+1,r,rt<<1|1,p);
    }
}tree;

tree.update(1,n,1,l,r,val);
tree.query(1,n,1,id);
```




<span id="9"><h4>9.    CDQ分治</h4></span>
题目描述
有 n 个元素，第 i 个元素有 a_i,b_i,c_i 三个属性，设 f(i) 表示满足 a_j <= a_i 且 b_j <= b_i 且 c_j <= c_i 且 j != i 的 j 的数量。
对于 d \in [0, n)d∈[0,n)，求 f(i) = d 的数量。
输入格式
第一行两个整数 n,k，表示元素数量和最大属性值。
接下来 n 行，每行三个整数 a_i ,b_i,c_i，分别表示三个属性值。
输出格式
n 行，第 d+1 行表示 f(i) = d 的 i 的数量。

```cpp
struct Point{
    int x,y,z;
    int id,val,ans;
    void read(){
        scanf("%d%d%d",&x,&y,&z);
    }
    bool operator<(const Point p)const{
        if(x==p.x){
            if(y==p.y)return z<p.z;
            return y<p.y;
        }
        return x<p.x;
    }
}p[MAXN],B[MAXN];
bool cmp1(Point p1,Point p2){
    if(p1.y==p2.y)return p1.z<p2.z;
    return p1.y<p2.y;
}
int lowbit(int x){return x&-x;}
int sum[MAXN<<2];
void add(int id,int val){
    for(id;id<=m;id+=lowbit(id))
        sum[id]+=val;
}
int query(int id){
    int ret=0;
    for(id;id;id-=lowbit(id))
        ret+=sum[id];
    return ret;
}
void CDQ(int l,int r){
    if(l==r)return;
    int mid=(l+r)>>1;
    CDQ(l,mid);CDQ(mid+1,r);
    int i,j;
    for(j=mid+1,i=l;j<=r;j++){
        while(i<=mid&&p[i].y<=p[j].y)
            add(p[i].z,p[i].val),i++;
        p[j].ans+=query(p[j].z);
    }
    for(j=l;j<i;j++){
        add(p[j].z,-p[j].val);
    }
    for(int i=l,l1=l,l2=mid+1;i<=r;i++){
        if(l2>r||l1<=mid&&p[l1].y<=p[l2].y)
            B[i]=p[l1++];
        else
            B[i]=p[l2++];
    }
    for(int i=l;i<=r;i++)
        p[i]=B[i];
}
int num[MAXN];
void solve() {
    sort(p+1,p+1+n);
    int _n=1;
    for(int i=2;i<=n;i++){
        if(p[_n].x==p[i].x&&p[_n].y==p[i].y&&p[_n].z==p[i].z)
            p[_n].val++;
        else
            p[++_n]=p[i];
    }
    CDQ(1,_n);
    for(int i=1;i<=_n;i++){
        num[p[i].ans+p[i].val-1]+=p[i].val;
    }
    for(int i=0;i<n;i++)
        printf("%d\n",num[i]);
}
void init(){
    scanf("%d%d",&n,&m);
    for(int i=1;i<=n;i++){
        p[i].read();
        p[i].id=i;
        p[i].val=1;
        p[i].ans=0;
    }
}
```


<span id="10"><h4>10.    主席树区间第k小</h4></span>
```cpp
int n, m;
struct HJT{
    int n;
    int node_cnt;
    int sum[MAXN<<5],rt[MAXN],lc[MAXN<<5],rc[MAXN<<5];
    void init(int _n){
        n=_n;
        build(rt[0],1,n);
    }
    void build(int &rt,int l,int r){
        rt=++node_cnt;
        if(l==r)return;
        int mid=(l+r)>>1;
        build(lc[rt],l,mid);
        build(rc[rt],mid+1,r);
    }
    void modify(int id,int now){
        rt[now]=modify(id,rt[now-1],1,n);
    }
    int modify(int id,int o,int l,int r){
        int oo=++node_cnt;
        lc[oo]=lc[o];
        rc[oo]=rc[o];
        sum[oo]=sum[o]+1;
        if(l==r)return oo;
        int mid=(l+r)>>1;
        if(id<=mid)lc[oo]=modify(id,lc[oo],l,mid);
        else rc[oo]=modify(id,rc[oo],mid+1,r);
        return oo;
    }
    int query(int l,int r,int k){
        return query(rt[l-1],rt[r],1,n,k);
    }
    int query(int u,int v,int l,int r,int k){
        if(l==r)return l;
        int mid=(l+r)>>1;
        int x=sum[lc[v]]-sum[lc[u]];
        if(x>=k)return query(lc[u],lc[v],l,mid,k);
        else return query(rc[u],rc[v],mid+1,r,k-x);
    }
    
    
    //建图优化
    int modify(int l,int r,int pre,int pos,int point){
        int oo=++node_cnt;
        if(l==r&&pre>0){
            sap.addedge(oo,pre,INF);
        }
        if(l==r){
            sap.addedge(oo,point,INF);
            return oo;
        }
        lc[oo]=lc[pre];
        rc[oo]=rc[pre];
        int mid=(l+r)>>1;
        if(pos<=mid)lc[oo]=modify(l,mid,lc[pre],pos,point);
        else rc[oo]=modify(mid+1,r,rc[pre],pos,point);
        if(lc[oo]){
            sap.addedge(oo,lc[oo],INF);
        }
        if(rc[oo]){
            sap.addedge(oo,rc[oo],INF);
        }
        return oo;
    }
    void query(int l,int r,int x,int L,int R,int point){
        if(!x||L>r||R<l)return;
        if(L<=l&&R>=r){
            sap.addedge(point,x,INF);
            return;
        }
        int mid=(l+r)>>1;
        if(L<=mid)query(l,mid,lc[x],L,R,point);
        if(R>mid)query(mid+1,r,rc[x],L,R,point);
    }
    hjt.init(sum);
    rep(i,1,n){
        if(i>1)hjt.query(1,sum,hjt.rt[i-1],p[i].l,p[i].r,i+n);
        hjt.rt[i]=hjt.modify(1,sum,hjt.rt[i-1],p[i].a,i);
    }
    
}hjt;
int a[MAXN],b[MAXN];
void solve(){
    sort(b+1,b+1+n);
    int num=unique(b+1,b+n+1)-(b+1);
    hjt.init(num);
    rep(i,1,n){
        int id=lower_bound(b+1,b+num+1,a[i])-b;
        hjt.modify(id,i);
    }
    int l,r,k;
    while(m--){
        scanf("%d%d%d",&l,&r,&k);
        printf("%d\n",b[hjt.query(l,r,k)]);
    }
}
void init(){
    scanf("%d%d",&n,&m);
    rep(i,1,n){
        scanf("%d",a+i);
        b[i]=a[i];
    }
}


-----------------------区间求和
#include<bits/stdc++.h>

using namespace std;

const int N = 1001000, inf = 1e9;
typedef long long ll;
ll sum[N * 60];
int lc[N * 60], rc[N * 60];
int cnt;

int ins(int o, int l, int r, int x) {
    int oo = ++cnt;
    sum[oo] = sum[o] + x;
    lc[oo] = lc[o];
    rc[oo] = rc[o];
    if (l == r)return oo;
    int mid = (l + r) / 2;
    if (x <= mid)lc[oo] = ins(lc[o], l, mid, x);
    else rc[oo] = ins(rc[o], mid + 1, r, x);
    return oo;
}

ll query(int o1, int o2, int l, int r, ll ql, ll qr) {
    if (ql <= l && r <= qr)return sum[o1] - sum[o2];
    int mid = (l + r) / 2;
    ll ret = 0;
    if (ql <= mid)ret += query(lc[o1], lc[o2], l, mid, ql, qr);
    if (qr > mid)ret += query(rc[o1], rc[o2], mid + 1, r, ql, qr);
    return ret;
}

int rt[N], a[N], n, Q;
ll lans;

int main() {
    freopen("data.in","r",stdin);
    ios::sync_with_stdio(0);
    cin.tie(0);
    cin >> n >> Q;
    for (int i = 1; i <= n; ++i)
        cin >> a[i], rt[i] = ins(rt[i - 1], 1, inf, a[i]);
    while (Q-- > 0) {
        int l, r;
        cin >> l >> r;
        l = (l + lans) % n + 1;
        r = (r + lans) % n + 1;
        if (l > r)swap(l, r);
        ll ans = 0, now_max = 0;
        for (; ans < inf;) {
            ll tmp = query(rt[r], rt[l - 1], 1, inf, 1, ans + 1);
            if (tmp <= ans)break;
            ans=tmp;
        }
        lans = ans + 1;
        cout << (ans + 1) << '\n';
    }
    return 0;
}

```



<span id="11"><h4>11.    ST表</h4></span>
```cpp
    rep(i,1,n){
        scanf("%d",a+i);
        lg[i]=lg[i/2]+1;
    }
    rep(i,1,lg[n]){
        rep(j,1,n-((1<<i)-1)){
            ST[j][i]= ST[ST[j][i - 1]][i - 1];
        }
    }
---------------------------------
int f[MAXN][21]; // 第二维的大小根据数据范围决定，不小于log(MAXN)
for (int i = 1; i <= n; ++i)
    f[i][0] = read(); // 读入数据
for (int i = 1; i <= 20; ++i)
    for (int j = 1; j + (1 << i) - 1 <= n; ++j)
        f[j][i] = max(f[j][i - 1], f[j + (1 << (i - 1))][i - 1]);
```


<span id="12"><h4>12.    莫队</h4></span>
```cpp
#include<iostream>
#include<cstdio>
#include<cmath>
#include<algorithm>
#include<cctype>
using namespace std;
const int MAXN=50005,MAXM=200005;
struct A{
    int l,r,id;
}q[MAXM];
int n,m,a[MAXN],num[MAXN],ans[MAXM],tim,answer=0;
bool cmp(A a,A b)
{
    return a.l/tim==b.l/tim?a.r<b.r:a.l<b.l;
}
inline int read()
{
    char ch=getchar();
    int x=0;
    while(!isdigit(ch)) ch=getchar();
    while(isdigit(ch)){
        x=x*10+ch-'0';
        ch=getchar();
    }
    return x;
}
void add(int pos)
{
    if(++num[a[pos]]==1) answer++;
}
void remove(int pos)
{
    if(!--num[a[pos]])  answer--;
}
int main()
{
    int i,curl=1,curr=0;
    n=read();
    for(i=1;i<=n;i++){
        a[i]=read();
    }
    m=read();
    tim=sqrt(m);
    for(i=1;i<=m;i++){
        q[i].l=read();
        q[i].r=read();
        q[i].id=i;
    }
    sort(q+1,q+1+m,cmp);
    for(i=1;i<=m;i++){
        while(curl<q[i].l){
            remove(curl++);
        }
        while(curl>q[i].l){
            add(--curl);
        }
        while(curr<q[i].r){
            add(++curr);
        }
        while(curr>q[i].r){
            remove(curr--);
        }
        ans[q[i].id]=answer;
    }
    for(i=1;i<=m;i++){
        printf("%d\n",ans[i]);
    }
    return 0;
}
```		      


<span id="13"><h4>13.    红黑树黑科技</h4></span>
/*
template<class T> using Tree = tree<T,null_type,less<T>,rb_tree_tag,tree_order_statistics_node_update>;
定义一颗红黑树
int 关键字类型
null_type无映射(低版本g++为null_mapped_type)
less<int>从小到大排序
rb_tree_tag 红黑树（splay_tree_tag）
tree_order_statistics_node_update结点更新
插入t.insert();
删除t.erase();
Rank:t.order_of_key();
第K值:t.find_by_order(); t.find_by_order(0)第1个数
前驱:t.lower_bound();
后继t.upper_bound();
a.split(v,b)key小于等于v的元素属于a，其余的属于b //b之前的会被清空
T.lower_bound(x)   >=x的min的迭代器
T.upper_bound((x)  >x的min的迭代器
T.find_by_order(k) 有k个数比它小的数


Tree<int> t;
t.insert(1);
t.insert(2);
t.insert(2);
t.insert(3);
debug(t.size());
debug(t.order_of_key(0));
debug(t.order_of_key(1));
debug(t.order_of_key(2));
debug(t.order_of_key(3));
debug(t.order_of_key(4));
debug(*t.find_by_order(0));
debug(*t.find_by_order(1));
debug(*t.find_by_order(2));
debug(*t.find_by_order(3));
debug(*t.find_by_order(4));
debug(*t.lower_bound(0));
debug(*t.lower_bound(1));
debug(*t.lower_bound(2));
debug(*t.lower_bound(3));
debug(*t.lower_bound(4));
debug(t.lower_bound(4)==t.end());

t.size() 3
t.order_of_key(0) 0
t.order_of_key(1) 0
t.order_of_key(2) 1
t.order_of_key(3) 2
t.order_of_key(4) 3
*t.find_by_order(0) 1
*t.find_by_order(1) 2
*t.find_by_order(2) 3
*t.find_by_order(3) 0
*t.find_by_order(4) 0
*t.lower_bound(0) 1
*t.lower_bound(1) 1
*t.lower_bound(2) 2
*t.lower_bound(3) 3
t.lower_bound(4)==t.end() 1

```cpp
#include <bits/stdc++.h>
#include <ext/pb_ds/assoc_container.hpp>

using namespace std;
using namespace __gnu_pbds;

template<class T> using Tree = tree<T,null_type,less<T>,rb_tree_tag,tree_order_statistics_node_update>;

const int MX = 1e5+5;
#define sz(x) (int)(x).size()

int N, a[MX], ind[MX], ans[MX], ret;
vector<int> child[MX];
Tree<int> d[MX];

void comb(int a, int b) {
	if (sz(d[a]) < sz(d[b])) d[a].swap(d[b]);
	for (int i: d[b]) d[a].insert(i);
}

void dfs(int x) {
	ind[x] = x;
	for (int i: child[x]) {
		dfs(i);
		comb(x,i);
	}
	ans[x] = sz(d[x])-d[x].order_of_key(a[x]);
	d[x].insert(a[x]);
}

int main() {
	freopen("promote.in","r",stdin);
	freopen("promote.out","w",stdout);
	cin >> N; for (int i = 1; i <= N; ++i) cin >> a[i];
	for (int i = 2; i <= N; ++i) {
		int p; cin >> p;
		child[p].push_back(i);
	}
	dfs(1);
	for (int i = 1; i <= N; ++i) cout << ans[i] << "\n";
}
```
