
# 图论
### [0.性质](#0)
### [1.	Johnson 全源最短路](#1)
### [2.	Dijstra最短路](#2)
### [3.	Floyd最短路](#3)
### [4.	LCA](#4)
### [5.	Tarjan找割点](#5)
### [6.	Tarjan找强连通分量个数](#6)
### [7.	Tarjan缩点求路径最大点权和](#7)
### [8.	二分图最大匹配匈牙利算法](#8)
### [9.	拓扑排序](#9)
### [10.	最大流](#10)
### [11.	最小费用最大流](#11)
### [12.	欧拉路径](#12)
### [13.	树链剖分](#13)
### [14.  2-sat输出可行解](#14)
### [15. 最小割](#15)
### [16. 最小路径覆盖](#16)
### [17. 最长不降序列数量](#17)
### [18. 二分图的最佳完美匹配KM算法](#18)
### [19. 最小费用最大流（逐步炒菜）](#19)

<span id="0"><h4>0. 性质</h4></span>
树的重心性质：
1. 树中所有点到某个点的距离和中，到重心的距离和是最小的，如果有两个距离和，他们的距离和一样。
2. 把两棵树通过一条边相连，新的树的重心在原来两棵树重心的连线上。
3. 一棵树添加或者删除一个节点，树的重心最多只移动一条边的位置。
4. 一棵树最多有两个重心，且相邻。

网络流:
这里首先说一个技巧，无向图的网络流在建边时反向弧直接建成权值为v的边即可，因为这样的边一开始就是可以增广的。


<span id="1"><h4>1. Johnson 全源最短路</h4></span>
```cpp
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
<span id="2"><h4>2. Dijkstra最短路</h4></span>
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
<span id="3"><h4>3. Floyd最短路</h4></span>
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
<span id="4"><h4>4. LCA</h4></span>
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
<span id="5"><h4>5. Tarjan找割点</h4></span>
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
<span id="6"><h4>6. Tarjan找强连通分量个数</h4></span>
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
<span id="7"><h4>7. Tarjan缩点求路径最大点权和</h4></span>
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
<span id="8"><h4>8. 二分图最大匹配匈牙利算法</h4></span>
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
<span id="9"><h4>9. 拓扑排序</h4></span>
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
<span id="10"><h4>10. 最大流</h4></span>
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
<span id="11"><h4>11. 最小费用最大流</h4></span>
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
<span id="12"><h4>12. 欧拉路径</h4></span>
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
<span id="13"><h4>13. 树链剖分</h4></span>
```cpp
//已知一棵包含 NN 个结点的树（连通且无环），每个节点上包含一个数值，需要支持以下操作：
//操作 1： 格式： 1 x y z 表示将树从 x 到 y 结点最短路径上所有节点的值都加上 z。
//操作 2： 格式： 2 x y 表示求树从 x 到 y 结点最短路径上所有节点的值之和。
//操作 3： 格式： 3 x z 表示将以 x 为根节点的子树内所有节点值都加上 z。
//操作 4： 格式： 4 x 表示求以 x 为根节点的子树内所有节点值之和
int n,m;
long long val[MAXN];
int r;
long long p;
vector<int>vec[MAXN];
int fa[MAXN],depth[MAXN];
int siz[MAXN],bigson[MAXN];
void dfs1(int now,int pre){
    int mxson=-1;
    siz[now]=1;
    for(auto to:vec[now]){
        if(to==pre)continue;
        fa[to]=now;
        depth[to]=depth[now]+1;
        dfs1(to,now);
        siz[now]+=siz[to];
        if(mxson<siz[to]){
            bigson[now]=to;
            mxson=siz[to];
        }
    }
}
int top[MAXN],dfn[MAXN],pos[MAXN],cnt[MAXN];
int cuo;
void dfs2(int now,int tp){
    top[now]=tp;
    dfn[now]=++cuo;
    pos[cuo]=now;
    if(bigson[now])
        dfs2(bigson[now],tp);
    for(auto to:vec[now]){
        if(to!=fa[now]&&to!=bigson[now])
            dfs2(to,to);
    }
    cnt[now]=cuo;
}
struct Node{
    long long c;
    long long f;
}t[MAXN<<2];
void push_up(int rt){
    t[rt].c=(t[rt<<1].c+t[rt<<1|1].c)%p;
}
void build(int l,int r,int rt){
    if(l==r){
        t[rt].c=val[pos[l]];
        return;
    }
    int mid=(l+r)>>1;
    build(l,mid,rt<<1);
    build(mid+1,r,rt<<1|1);
    push_up(rt);
}
void down(int l,int r,int rt){
    t[rt<<1].f+=t[rt].f;
    t[rt<<1|1].f+=t[rt].f;
    int mid=(l+r)>>1;
    t[rt<<1].c+=t[rt].f*(mid-l+1)%p;
    t[rt<<1|1].c+=t[rt].f*(r-mid)%p;
    t[rt].f=0;
}
void change(int l,int r,int rt,int ls,int rs,long long z){
    if(ls<=l&&r<=rs){
        t[rt].c+=(r-l+1)*z;
        t[rt].f+=z;
        return;
    }
    if(t[rt].f)down(l,r,rt);
    int mid=(l+r)>>1;
    if(ls<=mid)change(l,mid,rt<<1,ls,rs,z);
    if(rs>mid)change(mid+1,r,rt<<1|1,ls,rs,z);
    push_up(rt);
}
long long getsum(int l,int r,int rt,int ls,int rs){
    if(ls<=l&&r<=rs){
        return t[rt].c;
    }
    if(t[rt].f)down(l,r,rt);
    int mid=(l+r)>>1;
    long long ans=0;
    if(ls<=mid)ans+=getsum(l,mid,rt<<1,ls,rs);
    if(rs>mid)ans+=getsum(mid+1,r,rt<<1|1,ls,rs);
    return ans%p;
}
void change_xtoy(int x,int y,long long z){
    while(top[x]!=top[y]){
        if(depth[top[x]]>depth[top[y]])swap(x,y);
        change(1,cuo,1,dfn[top[y]],dfn[y],z);
        y=fa[top[y]];
    }
    if(depth[x]>depth[y])swap(x,y);
    change(1,cuo,1,dfn[x],dfn[y],z);
}
void getson_xtoy(int x,int y){
    long long ans=0;
    while(top[x]!=top[y]){
        if(depth[top[x]]>depth[top[y]])swap(x,y);
        ans=(ans+getsum(1,cuo,1,dfn[top[y]],dfn[y]))%p;
        y=fa[top[y]];
    }
    if(depth[x]>depth[y])swap(x,y);
    ans=(ans+getsum(1,cuo,1,dfn[x],dfn[y]))%p;
    printf("%lld\n",ans);
}
void change_sontree(int x,long long z){
    change(1,cuo,1,dfn[x],cnt[x],z);
}
void getsum_sontree(int x){
    printf("%lld\n",getsum(1,cuo,1,dfn[x],cnt[x])%p);
}
void solve(){
    dfs1(r,0);
    dfs2(r,r);
    build(1,cuo,1);
    int op,x,y,z;
    while(m--){
        scanf("%d",&op);
        if(op==1){
            scanf("%d%d%d",&x,&y,&z);
            change_xtoy(x,y,z);
        }else if(op==2){
            scanf("%d%d",&x,&y);
            getson_xtoy(x,y);
        }else if(op==3){
            scanf("%d%d",&x,&z);
            change_sontree(x,z);
        }else{
            scanf("%d",&x);
            getsum_sontree(x);
        }
    }
}
void init(){
    scanf("%d%d%d%lld",&n,&m,&r,&p);
    for(int i=1;i<=n;i++)
        scanf("%lld",val+i);
    for(int i=1,u,v;i<n;i++){
        scanf("%d%d",&u,&v);
        vec[u].pb(v);
        vec[v].pb(u);
    }
}
```
<span id="14"><h4>14. 2-sat输出可行解</h4></span>
```cpp
/*
2
08:00 09:00 30
08:15 09:00 20

YES
08:00 08:30
08:40 09:00
*/
struct Edge{
    int l,r;
}edge[MAXN];
int dfn[MAXN],low[MAXN],fa[MAXN];
int cnt,fanum;
stack<int>s;
vector<int>vec[MAXN];
int ans[MAXN];
void tarjan(int now){
    dfn[now]=low[now]=++cnt;
    s.push(now);
    for(int i=0,len=vec[now].size();i<len;i++){
        int to=vec[now][i];
        if(!dfn[to]){
            tarjan(to);
            low[now]=min(low[now],low[to]);
        }else if(!fa[to]){
            low[now]=min(low[now],low[to]);
        }
    }
    if(dfn[now]==low[now]){
        fanum++;
        while(1){
            int t=s.top();s.pop();
            fa[t]=fanum;
            if(t==now)break;
        }
    }
}
void print(){
    for(int i=2;i<=(n<<1|1);i+=2){
        if(fa[i]<fa[i^1])ans[i>>1]=i;
        else ans[i>>1]=i^1;
    }
    printf("YES\n");
    for(int i=1;i<=n;i++){
        int l=edge[ans[i]].l;
        int r=edge[ans[i]].r;
        printf("%02d:%02d %02d:%02d\n",l/60,l%60,r/60,r%60);
    }
}
void solve(){
    for(int i=2;i<=(n<<1|1);i++){
        if(!dfn[i])
            tarjan(i);
    }
    for(int i=2;i<=(n<<1|1);i++){
        if(fa[i]==fa[i^1]){
            printf("NO\n");
            return;
        }
    }
    print();
}
bool judge(int l1,int r1,int l2,int r2){
    return (l1 < r2 && l2 < r1);
}
void init(){
    scanf("%d",&n);
    for(int i=1,h,m,l,r;i<=n;i++){
        scanf("%d:%d",&h,&m);
        l=h*60+m;
        scanf("%d:%d",&h,&m);
        r=h*60+m;
        scanf("%d",&m);
        edge[i<<1]={l,l+m};
        edge[i<<1|1]={r-m,r};
    }
    for(int i=2;i<=(n<<1);i++){
        for(int j=i+1;j<=(n<<1|1);j++){
            if((i>>1)==(j>>1))continue;
            if(judge(edge[i].l,edge[i].r,edge[j].l,edge[j].r)){
                if(judge(edge[i].l,edge[i].r,edge[j^1].l,edge[j^1].r)){
                    vec[i].pb(i^1);
                }else{
                    vec[i].pb(j^1);
                    vec[j].pb(i^1);
                }
            }
        }
    }
}
```


<span id="15"><h4>15. 最小割</h4></span>
```cpp
struct Edge{
    int nx;
    int to;
    int w;
}edge[MAXN];
int cnt=1;
int s,t;
int head[MAXN];
void add(int u,int v,int w){
    cnt++;
    edge[cnt].to=v;
    edge[cnt].w=w;
    edge[cnt].nx=head[u];
    head[u]=cnt;
}
int dis[MAXN];
int dfs(int now,int flow){
    if(now==t)return flow;
    int f=0;
    for(int &i=cur[now];i;i=edge[i].nx){
        int to=edge[i].to;
        if(dis[to]!=dis[now]+1||edge[i].w==0)continue;
        int tmp=dfs(to,min(edge[i].w,flow));
        flow-=tmp;
        edge[i].w-=tmp;edge[i^1].w+=tmp;
        f+=tmp;
        if(flow<=0)break;
    }
    return f;
}
bool bfs(){
    memcpy(cur,head,sizeof(head));
    memarray(dis,0);
    dis[s]=1;
    queue<int>q;
    q.push(s);
    while(!q.empty()){
        int now=q.front();q.pop();
        for(int i=head[now];i;i=edge[i].nx){
            int to=edge[i].to;
            if(dis[to]||edge[i].w==0)continue;
            dis[to]=dis[now]+1;
            q.push(to);
        }
    }
    return dis[t];
}
int dinic(){
    int ans=0;
    while(bfs()){
        ans+=dfs(s,INF);
    }
    return ans;
}
void solve(){
    int ans=dinic();
    printf("%d\n",ans);
}
void init(){
    scanf("%d%d%d%d",&n,&m,&s,&t);
    for(int i=1;i<=n;i++){
        add(i,i+n,1);add(i+n,i,0);
    }
    s+=n;
    for(int i=1,u,v;i<=m;i++){
        scanf("%d%d",&u,&v);
        add(u+n,v,INF);add(v,u+n,0);
        add(v+n,u,INF);add(u,v+n,0);
    }
}
```


<span id="16"><h4>16. 最小路径覆盖</h4></span>
```cpp
struct Edge{
    int nx,to,w;
}edge[MAXN];
int cnt=1;
int head[MAXN];
void add(int u,int v,int w){
    cnt++;
    edge[cnt].to=v;
    edge[cnt].w=w;
    edge[cnt].nx=head[u];
    head[u]=cnt;

    cnt++;
    edge[cnt].to=u;
    edge[cnt].w=0;
    edge[cnt].nx=head[v];
    head[v]=cnt;
}
int dis[MAXN];
int cur[MAXN];
int s,t;
bool bfs(){
    memcpy(cur,head,sizeof head);
    memarray(dis,0);
    queue<int>q;
    q.push(s);
    dis[s]=1;
    while(!q.empty()){
        int now=q.front();q.pop();
        for(int i=head[now];i;i=edge[i].nx){
            int to=edge[i].to;
            if(dis[to]||edge[i].w==0)continue;
            dis[to]=dis[now]+1;
            q.push(to);
        }
    }
    return dis[t];
}
int ansto[MAXN];
bool tag[MAXN];
int dfs(int now,int flow){
    if(now==t)return flow;
    int ret=0;
    for(int i=head[now];~i;i=edge[i].nx){
        int to=edge[i].to;
        if(dis[to]!=dis[now]+1||edge[i].w==0)continue;
        int tmp=dfs(to,min(flow,edge[i].w));
        if(tmp<=0)continue;
        ansto[now]=to;
        if(now!=s)tag[to-n]=1;

        flow-=tmp;
        edge[i].w-=tmp;edge[i^1].w+=tmp;
        ret+=tmp;
//        return ret;
        if(flow<=0)break;
    }
    return ret;
}

int dinic(){
    int ret=0;
    while(bfs()){
        ret+=dfs(s,INF);
    }
    return ret;
}
void solve(){
    int ans=dinic();
    for(int i=1;i<=n;i++){
        if(!tag[i]){
            int now=i;
            printf("%d",now);
            while(ansto[now]&&ansto[now]!=t){
                printf(" %d",ansto[now]-n);
                now=ansto[now]-n;
            }
            printf("\n");
        }
    }
    printf("%d\n",n-ans);
}
void init(){
    memarray(head,-1);
    scanf("%d%d",&n,&m);
    s=2*n+5,t=s+1;
    for(int i=1,u,v;i<=m;i++){
        scanf("%d%d",&u,&v);
        add(u,v+n,1);
    }
    for(int i=1;i<=n;i++)
        add(s,i,1);
    for(int i=1;i<=n;i++)
        add(i+n,t,1);
}
```


<span id="17"><h4>17. 最长不降序列数量</h4></span>
给定正整数序列
计算其最长不下降子序列的长度 s.
如果每个元素只允许使用一次，计算从给定的序列中最多可取出多少个长度为 s 的不下降子序列。
如果允许在取出的序列中多次使用	
 （其他元素仍然只允许使用一次），则从给定序列中最多可取出多少个不同的长度为 s 的不下降子序列。
```cpp
int a[505],dp[505];
void init(){
    scanf("%d",&n);
    for(int i=1;i<=n;i++){
        scanf("%d",a+i);
        dp[i]=1;
    }
    if(n==1){
        printf("1\n");
        printf("1\n");
        printf("1\n");
        return;
    }
    int len=1;
    for(int i=2;i<=n;i++){
        for(int j=1;j<i;j++){
            if(a[i]>=a[j])
                dp[i]=max(dp[i],dp[j]+1);
        }
        len=max(len,dp[i]);
    }
    printf("%d\n",len);
    s=0,t=2*n+1;

    for(int i=1;i<=n;i++){
        if(dp[i]==1)
            add(s,i,1);
    }
    for(int i=1;i<=n;i++){
        if(dp[i]==len)
            add(i+n,t,1);
    }
    for(int i=1;i<=n;i++)
        add(i,i+n,1);
    for(int i=1;i<=n;i++){
        for(int j=i+1;j<=n;j++){
            if(a[i]>a[j]||dp[j]!=dp[i]+1)continue;
            add(i+n,j,1);
        }
    }
    int ans=dinic();
    printf("%d\n",ans);

    add(s,1,INF);
    add(1,1+n,INF);
    if(dp[n]==len){
        add(n+n,t,INF);
        add(n,n+n,INF);
    }
    ans+=dinic();
    printf("%d\n",ans);
}
```

<span id="18"><h4>18. 二分图的最佳完美匹配KM算法</h4></span>
```cpp
bool visx[105],visy[105];
int linker[105];
int slack[105];
int ly[105],lx[105];
int g[105][105];
int nx,ny;
bool dfs(int x)
{
    visx[x]=true;
    for(int y=1; y<=ny; y++)
    {
        if(visy[y])continue;
        int tmp=lx[x]+ly[y]-g[x][y];
        if(tmp==0)
        {
            visy[y]=true;
            if(linker[y]==-1||dfs(linker[y]))
            {
                linker[y]=x;
                return true;
            }
        }
        else if(slack[y]>tmp)
        {
            slack[y]=tmp;
        }
    }
    return false;
}
int KM()
{
    memarray(linker,-1);
    memarray(ly,0);
    for(int i=1; i<=nx; i++)
    {
        lx[i]=-INF;
        for(int j=1; j<=ny; j++)
            lx[i]=max(lx[i],g[i][j]);
    }
    for(int x=1; x<=nx; x++)
    {
        for(int i=1; i<=ny; i++)
            slack[i]=INF;
        while(true)
        {
            memarray(visx,false);
            memarray(visy,false);
            if(dfs(x))break;
            int d=INF;
            for(int i=1; i<=ny; i++)
                if(!visy[i]&&d>slack[i])
                    d=slack[i];
            for(int i=1; i<=nx; i++)
                if(visx[i])
                    lx[i]-=d;
            for(int i=1; i<=ny; i++)
            {
                if(visy[i])ly[i]+=d;
                else slack[i]-=d;
            }
        }
    }
    int ret=0;
    for(int i=1; i<=ny; i++)
        if(linker[i]!=-1)
            ret+=g[linker[i]][i];
    return ret;
}
void solve()
{
    for(int i=1; i<=n; i++)
        for(int j=1; j<=n; j++)
            g[i][j]=-g[i][j];
    printf("%d\n",-KM());
    for(int i=1; i<=n; i++)
        for(int j=1; j<=n; j++)
            g[i][j]=-g[i][j];
    printf("%d\n",KM());
}
void init()
{
    scanf("%d",&n);
    nx=ny=n;
    for(int i=1; i<=n; i++)
        for(int j=1; j<=n; j++)
            scanf("%d",g[i]+j);
}
```



<span id="19"><h4>19. 最小费用最大流（逐步炒菜）</h4></span>
```cpp
int s,t;

struct Edge
{
    int from,to;
    long long cap,flow,cost;
    Edge(int u,int v,int ca,int f,int co):from(u),to(v),cap(ca),flow(f),cost(co){};
};
int sum=0;
int T[105][105];
struct MCMF {
    int nn,mm,s,t;
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
        this->nn = n;
//        for(int i = 0; i <= n; i++) G[i].clear();
//        edges.clear();
//        r.clear();
//        r.push_back({0,0,0});
    }
    void add_edge(int from,int to, int cap,int cost) {
        edges.push_back((Edge){from,to,cap,0,cost});
        edges.push_back((Edge){to,from,0,0,-cost});
        mm = edges.size();
        G[from].push_back(mm-2);
        G[to].push_back(mm-1);
    }
    int pre[MAXN];
    bool spfa(int s,int t,long long &flow,long long &cost) {
        for(int i = 0; i <= nn; i++) d[i] = INF;
        for(int i = 0; i <= nn; i++) inq[i] = 0;
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
                    /*记录返回*/
                    pre[e.to]=G[u][i];

                    if(!inq[e.to]) {Q.push(e.to); inq[e.to] = 1;}
                }
            }
        }
        if(d[t] == INF) return false;
        flow += a[t];
        cost += d[t] * a[t];
//        r.push_back({flow,cost,d[t]});
        int u = t;
        while(u != s) {
            edges[p[u]].flow += a[t];
            edges[p[u]^1].flow -= a[t];
            u = edges[p[u]].from;
        }
        return true;
    }
    int layer[MAXN];
    long long Mincost(int s,int t) {
        long long cost = 0;
        long long flow=0;
        while(spfa(s,t,flow,cost)){
            /**由于我们跑一次spfa只能找出一次增广路，
            所以我们可以暂时不连不需要的边。
            一开始，我们把所有厨师做倒数第１道菜与所有菜连好，
            然后找一条增广路，这条增广路上一定经过了一个点，
            表示第j个厨师做倒数第１道菜，于是我们添加点
            （第j个厨师做倒数第２道菜），与汇点和所有菜连边，
            以此类推。**/
            int from=edges[pre[t]].from;
            from=(from-n-1)/sum+1;
            layer[from]++;
            if(layer[from]>sum)continue;
            for(int k=1;k<=n;k++){
                add_edge(k,n+(from-1)*sum+layer[from],1,T[k][from]*layer[from]);
            }
        }
        return cost;
    }
} mcmf;
int p[105];
void solve(){
    printf("%lld\n",mcmf.Mincost(s,t));
}
void init(){
    scanf("%d%d",&n,&m);
    for(int i=1;i<=n;i++){
        scanf("%d",p+i);
        sum+=p[i];
    }
    for(int i=1;i<=n;i++){
        for(int j=1;j<=m;j++){
            scanf("%d",T[i]+j);
        }
    }
    s=(m+1)*sum+5;t=s+1;
    mcmf.init(t);
    for(int i=1;i<=n;i++)
        mcmf.add_edge(s,i,p[i],0);
    for(int j=1;j<=m;j++)
        for(int i=1;i<=sum;i++)
            mcmf.add_edge(n+(j-1)*sum+i,t,1,0);
    for(int j=1;j<=m;j++){
        for(int k=1;k<=n;k++){
            mcmf.add_edge(k,n+(j-1)*sum+1,1,T[k][j]*1);
        }
    }
    for(int i=1;i<=m;i++)
        mcmf.layer[i]=1;
}
```
