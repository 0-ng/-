
# 图论
[0.性质](#0)
#### 1.	Johnson 全源最短路
#### 2.	Dijstra最短路
#### 3.	Floyd最短路
#### 4.	LCA
#### 5.	Tarjan找割点
#### 6.	Tarjan找强连通分量个数
#### 7.	Tarjan缩点求路径最大点权和
#### 8.	二分图最大匹配匈牙利算法
#### 9.	拓扑排序
#### 10.	最大流
#### 11.	最小费用最大流
#### 12.	欧拉路径
#### 13.	树链剖分
#### 14.  2-sat输出可行解


<span id="0"><h4>性质</h4></span>
树的重心性质：
1. 树中所有点到某个点的距离和中，到重心的距离和是最小的，如果有两个距离和，他们的距离和一样。
2. 把两棵树通过一条边相连，新的树的重心在原来两棵树重心的连线上。
3. 一棵树添加或者删除一个节点，树的重心最多只移动一条边的位置。
4. 一棵树最多有两个重心，且相邻。


#### 1.	Johnson 全源最短路
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
#### 2.	Dijkstra最短路
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
#### 3.	Floyd最短路
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
#### 4.	LCA
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
#### 5.	Tarjan找割点
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
#### 6.	Tarjan找强连通分量个数
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
#### 7.	Tarjan缩点求路径最大点权和
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
#### 8.	二分图最大匹配匈牙利算法
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
#### 9.	拓扑排序
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
#### 10.	最大流
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
#### 11.	最小费用最大流
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
#### 12.	欧拉路径
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
#### 13.树链剖分
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
#### 15. 2-sat输出可行解
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
