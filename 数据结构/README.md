
### [1. RMQ](#1)
### [2. 二维ST表max min](#2)
### [3. 带权并查集](#3)
### [4. 线段树](#4)
### [5. 线段树区间最大最小单点修改](#5)
### [6. 树状数组](#6)
### [7. 二维树状数组求矩阵异或](#7)

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
    update(1,1);
    update(idx+1,-1);
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



<span id="8"><h4>7.    线段树区间赋值单点查询</h4></span>
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
