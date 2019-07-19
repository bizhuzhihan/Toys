// Zhihan ZHU
// 201630845294

// 2A (AC) 斐波那契数列(第二章题目)
#include <iostream>
#include <stdio.h>
using namespace std;

int main() {
    int n;
    while (cin >> n) {
        if (n > 0) {
            int N = 1;
            int N_1 = 1;
            int N_f = 1;
            for (int i = 2; i < n; i++) {
                N_f = (N + N_1) % 1997;
                N_1 = N;
                N = N_f;
            }
            cout << N_f << '\n';
        } else {
            break;
        }
    }
    return 0;
}

// 2A (MLE)
#include <iostream>
#include <stdio.h>
using namespace std;

int *fun(int n){
    int a[2];
    a[0] = 1;
    a[1] = 1;
    if (n>2) {
        int *p = fun(n-1);
        a[0] = (*p + *(p+1)) % 1997;
        a[1] = *p % 1997
        ;
    }
    return a;
}

int main() {
    int n;
    while (cin >> n) {
        if (n > 0) {
            cout << *fun(n) << '\n';
        } else {
            break;
        }
    }
    return 0;
}

// 2B (AC) a的b次方模c (第二章题目)
#include <iostream>
#include <stdio.h>
using namespace std;

long long int fun(long long int a,long long int b,long long int c) {
    long long int value = 1;
    if (b>1) {
        if (b>(2*(b/2))) {
            value = fun(a, b/2, c);
            value *= (value * a);
            value = value % c;
        }else{
            value = fun(a, b/2, c);
            value *= value;
            value = value % c;
        }
    }else{
        value = a % c;
    }
    return value;
}

int main() {
    long long int a, b, c;
    while (cin >> a >> b >> c) {
        if ((a>=0) & (b>=0) & (c>1)) {
            long long int r;
            if (b == 0) {
                r = 1;
            }else{
                r = fun(a, b, c);
            }
            cout << r << '\n';
        } else {
            break;
        }
    }
    return 0;
}

// 2C (AC) 输油管道问题(第二章题目)
#include <iostream>
#include <stdio.h>
using namespace std;

int pfun(int Y[], int p, int r){
    int i=p,j=1+r;
    int x=Y[p];
    while (true) {
        while (Y[++i]<x&&i<r);
        while (Y[--j]>x);
        if (i>=j) {
            break;
        }
        swap(Y[i], Y[j]);
    }
    Y[p]=Y[j];
    Y[j]=x;
    return j;
}

void qsorting(int Y[], int p, int r, int n){
    if (p<r) {
        int q=pfun(Y,p,r);
        if (q!=(n/2)) {
            if (q<(n/2)) {
                qsorting(Y, q+1, r, n);
            }
            if (q>(n/2)) {
                qsorting(Y, p, q-1, n);
            }
        }
    }
}

int main() {
    int n;
    cin >> n;
    int Y[n];
    int sum = 0;
    for (int i=0; i<n; i++) {
        int x,y;
        cin >> x >> y;
        if ((x>=(-10000))&(x<=10000)&(y>=(-10000))&(y<=10000)) {
            Y[i] = y;
        } else {
            break;
        }
    }
    qsorting(Y, 0, n-1, n);
    for (int i=0; i<n; i++) {
        if (Y[n/2] > Y[i]) {
            sum += (Y[n/2] - Y[i]);
        }else{
            sum += (Y[i] - Y[n/2]);
        }
    }
    cout<< sum << '\n';
    return 0;
}

// 3A (AC) 回文(第三章题目)
#include <iostream>
#include <stdio.h>
using namespace std;

int LCSr(int n,char x[])
{
    char y[n];
    for (int ii=0; ii<n; ii++) {
        y[ii] = x[n-1-ii];
    }
    int i,j;
    int c[n+1][n+1];
    int b[n+1][n+1];
    for (i=0; i<n+1; i++) {
        c[i][0] = 0;
    }
    for (j=0; j<n+1;j++) {
        c[0][j] = 0;
    }
    for (i=1; i<n+1; i++) {
        for (j=1; j<n+1; j++) {
            if (x[i-1]==y[j-1]) {
                c[i][j] = c[i-1][j-1] + 1;
                b[i][j] = 1;
            }
            else if (c[i-1][j]>=c[i][j-1]) {
                c[i][j] = c[i-1][j];
                b[i][j] = 2;
            }
            else {
                c[i][j] = c[i][j-1];
                b[i][j] = 3;
            }
        }
    }
    int X[n];
    int Y[n];
    for (int i=0; i<n; i++) {
        X[i] = 0;
        Y[i] = 0;
    }
    i = n;
    j = n;
    int counter = 0;
    while (i>0&&j>0) {
        if (b[i][j]==1) {
            counter++;
            X[--i] = counter;
            Y[--j] = counter;
        }else if (b[i][j]==2) {
            i--;
        }else {
            j--;
        }
    }
    int ans=0, ans2=0;
    int p1=0, p2=0, p3=0;
    int pos[2]={n, n+1};
    for (; p1<n; p1++) {
        if (X[p1]>0) {
            for (; p2<=n-1-p1; p2++) {
                if (Y[p2]>0) {
                    ans++;
                    p2++;
                    break;
                }
            }
            for (; p3<n-1-p1; p3++) {
                if (Y[p3]>0) {
                    ans2++;
                    p3++;
                    break;
                }
            }
        }
    }
    int Ans;
    if ((n-2*ans2)<(n+1-2*ans)) {
        Ans = n - 2*ans2;
    }else {
        Ans = n + 1 - 2*ans;
    }
    return Ans;
}

int main() {
    int n;
    while (cin >> n) {
        if (n<=0) {
            break;
        }
        char x[n];
        for (int i=0; i<n; i++) {
            cin >> x[i];
        }
        cout<<LCSr(n, x)<<'\n';
    }
    return 0;
}

// 3B (AC) 石子合并问题(第三章题目)
#include <iostream>
#include <stdio.h>
using namespace std;

void stones(long long int p[], long long int n, long long int *ans) {
    long long int sum[n][n], max[n][n], min[n][n];
    for (long long int i=0; i<n; i++) {
        sum[i][i] = p[i];
        max[i][i] = 0;
        min[i][i] = 0;
    }
    for (long long int l=2; l<=n; l++) {
        for (long long int i=0; i<n; i++) {
            long long int j = (i+l-1) ;
            sum[i][j%n] = sum[(i+1)%n][j%n] + p[i];
            max[i][j%n] = max[(i+1)%n][j%n] + sum[i][j%n];
            min[i][j%n] = min[(i+1)%n][j%n] + sum[i][j%n];
            for (long long int k=(i+1); k<j; k++) {
                long long int M = max[i][k%n] + max[(k+1)%n][j%n] + sum[i][j%n];
                long long int m = min[i][k%n] + min[(k+1)%n][j%n] + sum[i][j%n];
                if (M>max[i][j%n]) {
                    max[i][j%n] = M;
                }
                if (m<min[i][j%n]) {
                    min[i][j%n] = m;
                }
            }
        }
    }
    long long int Min = min[0][n-1];
    long long int Max = max[0][n-1];
    for (long long int i=1; i<n; i++) {
        if (Max<max[i][(i+n-1)%n]) {
            Max = max[i][(i+n-1)%n];
        }
        if (Min>min[i][(i+n-1)%n]) {
            Min = min[i][(i+n-1)%n];
        }
    }
    ans[0] = Min;
    ans[1] = Max;
}

int main() {
    long long int n;
    while (cin>>n) {
        if (n==0) {
            break;
        }
        long long int p[n];
        for (long long int i=0; i<n; i++) {
            cin>>p[i];
        }
        long long int ans[2]={0,0};
        stones(p, n, ans);
        if (n==1) {
            cout<<p[0]<<' '<<p[1]<<'\n';
        }
        cout<<ans[0]<<' '<<ans[1]<<'\n';
    }
    return 0;
}

// 3C (AC) 最长有序子序列(第三章题目)
#include <iostream>
#include <stdio.h>
using namespace std;

int qsorting(int Y[], int p, int r){
    int i = p, j = 1+r;
    int x = Y[p];
    while (true) {
        while (Y[++i]<x&&i<r);
        while (Y[--j]>x);
        if (i>=j) {
            break;
        }
        swap(Y[i], Y[j]);
    }
    Y[p] = Y[j];
    Y[j] = x;
    return j;
}

void Qsorting(int Y[], int p, int r){
    if (p<r) {
        int q = qsorting(Y,p,r);
        Qsorting(Y, q+1, r);
        Qsorting(Y, p, q-1);
    }
}

int LCS(int m,int n,int x[],int y[])
{
    int i,j;
    int c[m+1][n+1];
    for (i=0; i<m+1; i++) {
        c[i][0] = 0;
    }
    for (j=0; j<n+1;j++) {
        c[0][j] = 0;
    }
    for (i=1; i<m+1; i++) {
        for (j=1; j<n+1; j++) {
            if (x[i-1]==y[j-1]) {
                c[i][j] = c[i-1][j-1]+1;
            }
            else if (c[i-1][j]>=c[i][j-1]) {
                c[i][j] = c[i-1][j];
            }
            else {
                c[i][j] = c[i][j-1];
            }
        }
    }
    return c[m][n];
}

int main() {
    int n;
    while (cin>>n) {
        if (n>0) {
            int a0[n], a[n];
            for (int i=0; i<n; i++) {
                cin >> a0[i];
                a[i] =  a0[i];
            }
            Qsorting(a, 0, n-1);
            int i = 0;
            int j = 0;
            while (i<n-1) {
                if (a[i]!=a[i+1]) {
                    a[j] = a[i];
                    j++;
                }
                i++;
            }
            a[j] = a[i];
            cout << LCS(n, j+1, a0, a) <<'\n';
        }else {
            break;
        }
    }
    return 0;
}

// 4A (AC) 财务问题（第四章题目）
#include <iostream>
#include <stdio.h>
using namespace std;

int main() {
    int n,m;
    while (cin>>n>>m) {
        if (n==0&&m==0) {
            break;
        }
        int month[12];
        for (int i=0; i<12; i++) {
            month[i] = 0;
        }
        int max = 0;
        for (int i=4; i>0; i--) {
            if ((i*n-(5-i)*m)<0) {
                max = i;
                break;
            }
        }
        int j = 0;
        for (; j<max; j++) {
            month[j] = 1;
        }
        for (; j<5; j++) {
            month[j] = 0;
        }
        for (; j<12; j++) {
            int sum = n;
            for (int k=j-1; k>=j-4; k--) {
                if (month[k]==0) {
                    sum -= m;
                }else {
                    sum += n;
                }
            }
            if (sum<0) {
                month[j] = 1;
            }else {
                month[j] = 0;
            }
        }
        int profit = 0;
        for (int i=0; i<12; i++) {
            if (month[i]==0) {
                profit -= m;
            }else{
                profit += n;
            }
        }
        if (profit>0) {
            cout << profit << '\n';
        }else{
            cout << "NO" << '\n';
        }
    }
    return 0;
}

// 4B (AC) 游戏预测（第四章题目）
#include <iostream>
#include <stdio.h>
using namespace std;

int pfun(int Y[], int p, int r){
    int i=p,j=1+r;
    int x=Y[p];
    while (true) {
        while (Y[++i]<x&&i<r);
        while (Y[--j]>x);
        if (i>=j) {
            break;
        }
        swap(Y[i], Y[j]);
    }
    Y[p]=Y[j];
    Y[j]=x;
    return j;
}

void qsorting(int Y[], int p, int r){
    if (p<r) {
        int q=pfun(Y,p,r);
        qsorting(Y, q+1, r);
        qsorting(Y, p, q-1);
    }
}

int main() {
    int m,n;
    while (cin>>m>>n) {
        if (n==0&&m==0) {
            break;
        }
        int u[n];
        for (int i=0; i<n; i++) {
            u[i] = 0;
        }
        for (int i=0; i<n; i++) {
            cin >> u[i];
        }
        qsorting(u, 0, n-1);
        int gap[n];
        for (int i=0; i<n; i++) {
            gap[i] = 0;
        }
        gap[n-1] = n*m - u[n-1];
        for (int i=n-2; i>=0; i--) {
            if (gap[i+1]==0) {
                gap[i] = u[i+1] - u[i] + gap[i+1] - 1;
            }else{
                gap[i] = u[i+1] - u[i] + gap[i+1] - 2;
            }
        }
        int time = 0;
        for (int i=0; i<n; i++) {
            time += 1*(gap[i]==0);
        }
        cout << time << '\n';
    }
    return 0;
}

// 4C (AC) 最优合并问题（第四章题目）
#include <iostream>
#include <stdio.h>
using namespace std;

int pfun(int Y[], int p, int r){
    int i=p,j=1+r;
    int x=Y[p];
    while (true) {
        while (Y[++i]<x&&i<r);
        while (Y[--j]>x);
        if (i>=j) {
            break;
        }
        swap(Y[i], Y[j]);
    }
    Y[p]=Y[j];
    Y[j]=x;
    return j;
}

void qsorting(int Y[], int p, int r){
    if (p<r) {
        int q=pfun(Y,p,r);
        qsorting(Y, q+1, r);
        qsorting(Y, p, q-1);
    }
}

int main() {
    int n;
    while (cin>>n) {
        if (n==0) {
            break;
        }
        int m[n],M[n];
        for (int i=0; i<n; i++) {
            int v;
            cin >> v;
            m[i] = v;
            M[i] = v;
        }
        qsorting(m, 0, n-1);
        int min = 0;
        for (int i=1; i<n; i++) {
            m[i] += m[i-1];
            min += m[i] - 1;
            qsorting(m, i, n-1);
        }
        qsorting(M, 0, n-1);
        int max = 0;
        for (int i=n-2; i>=0; i--) {
            M[i] += M[i+1];
            max += M[i] - 1;
        }
        cout<<max<<' '<<min<<'\n';
    }
    return 0;
}

// 5A (AC) 全排列问题（第五章题目）
#include <iostream>
#include <stdio.h>
using namespace std;

class P{
    friend void search(int n,int *p,int *p2);
private:
    bool check(int k);
    void backtrack(int t);
    int n, *x, *y;
};

bool P::check(int t) {
    if (y[t]==1) {
        return false;
    }
    y[t] = 1;
    return true;
}

void P::backtrack(int t)
{
    if (t==n) {
        for (int i=0;i<n-1;i++) {
            cout<<x[i]+1<<' ';
        }
        cout<<x[n-1]+1<<'\n';
    }
    else{
        for (int i=0;i<n;i++) {
            x[t]=i;
            if (check(i)) {
                backtrack(t+1);
                y[i] = 0;
            }
        }
    }
}

void search(int n,int *p,int *p2)
{
    P X;
    X. n=n;
    X.x=p;
    X.y=p2;
    X.backtrack(0);
}

int main() {
    int n;
    while (cin>>n) {
        if (n==0) {
            break;
        }
        int *p2=new int [n];
        for(int i=0;i<n;i++) {
            p2[i]= 0;
        }
        int *p=new int [n];
        for(int i=0;i<n;i++) {
            p[i]= 0;
        }
        search(n,p,p2);
        delete [] p;
        delete [] p2;
    }
    return 0;
}

// 5B (AC) 全组合(第五章题目)
#include <iostream>
#include <stdio.h>
using namespace std;

class C{
    friend void search(int n, int m);
private:
    int n, m, *x;
    void backtrack(int t);
};

void C::backtrack(int t){
    if (t>m) {
        for (int i=1; i<m; i++) {
            cout<<x[i]<<' ';
        }
        cout<<x[m]<<endl;
    }
    else{
        for (int i=x[t-1]+1; i<=n; i++) {
            x[t]=i;
            backtrack(t+1);
        }
    }
}

void search(int n, int m){
    C X;
    X.m = m;
    X.n = n;
    int *p = new int [m+1];
    for(int i=0;i<=m;i++) {
        p[i]= 0;
    }
    X.x=p;
    X.backtrack(1);
    delete [] p;
}

int main() {
    int n, m;
    while (cin>>n>>m) {
        if (((n==0)&&(m==0))||(n<m)) {
            break;
        }
        search(n, m);
    }
    return 0;
}

// 5C (AC) 整数变换问题（第五章题目）
#include <iostream>
#include <stdio.h>
using namespace std;

int f(int i){
    return 3*i;
}

int g(int i){
    return int(i/2);
}

class T{
    friend void search(int n, int m, int num);
private:
    int n, m, num, *x, *y;
    void backtrack(int t);
};

void T::backtrack(int t){
    for (int i=1; i<=2; i++) {
        x[t]=i;
        if (i==1) {
            y[t] = f(y[t-1]);
        }else{
            y[t] = g(y[t-1]);
        }
        if (y[t]==m) {
            num = t;
        }else{
            if (t<num) {
                backtrack(t+1);
            }
        }
    }
}

void search(int n, int m, int num){
    T X;
    X.m = m;
    X.n = n;
    X.num = num;
    int *p = new int [num+1];
    for(int i=0;i<=num;i++) {
        p[i]= 0;
    }
    X.x = p;
    int *p2 = new int [num+1];
    for(int i=0;i<=num;i++) {
        p2[i]= n;
    }
    X.y = p2;
    X.backtrack(1);
    delete [] p;
    cout<< X.num <<'\n';
}

int num_init(int n, int m){
    int num = 0;
    while (n!=m) {
        if (n>m) {
            n = g(n);
        }
        else{
            n = f(n);
        }
        num++;
        if (num >= 25) {
            break;
        }
    }
    return num;
}

int main() {
    int n, m;
    while (cin>>n>>m) {
        if (n==0) {
            break;
        }
        int num = num_init(n, m);
        search(n, m, num);
    }
    return 0;
}

// 6A (AC) 列车问题(第六章题目)
#include <iostream>
#include <stdio.h>
using namespace std;

struct order{
    int begin;
    int end;
    int num;
};

bool toadd(order T, int *M, int k){
    for (int i=T.begin+1; i<=T.end; i++) {
        if ((M[i]-T.num)<0) {
            return false;
        }
    }
    return true;
}

void add(order T, int *M, int k){
    for (int i=T.begin+1; i<=T.end; i++) {
        M[i] -= T.num;
    }
}

void re(order T, int *M, int k){
    for (int i=T.begin+1; i<=T.end; i++) {
        M[i] += T.num;
    }
}

void fun(order *T, int &d, int *M, int &maxprofit, int k, int t, int value){
    if (k<t) {
        for (int i=0; i<=1; i++) {
            if (i==1) {
                if (toadd(T[k], M, k)) {
                    add(T[k], M, k);
                    value += T[k].num*(T[k].end-T[k].begin);
                    maxprofit = max(value, maxprofit);
                    fun(T, d, M, maxprofit, k+1, t, value);
                    value -= T[k].num*(T[k].end-T[k].begin);
                    re(T[k], M, k);
                }else{
                    d -= T[k].num*(T[k].end-T[k].begin);
                    if (d>maxprofit) {
                        fun(T, d, M, maxprofit, k+1, t, value);
                        d += T[k].num*(T[k].end-T[k].begin);
                    }else{
                        d += T[k].num*(T[k].end-T[k].begin);
                    }
                }
            }else{
                d -= T[k].num*(T[k].end-T[k].begin);
                if (d>maxprofit) {
                    fun(T, d, M, maxprofit, k+1, t, value);
                    d += T[k].num*(T[k].end-T[k].begin);
                }else{
                    d += T[k].num*(T[k].end-T[k].begin);
                }
            }
        }
    }
}

int main() {
    int n, m, t;
    while (cin>>n>>m>>t) {
        if (((n==0)&&(m==0)&&(t==0))||(m>=30)||(t>=30)) {
            break;
        }else{
            order *T = new order [t];
            for (int i=0; i<t; i++) {
                cin >> T[i].begin >> T[i].end >> T[i].num;
            }
            int *M = new int [m+1];
            for (int i=0; i<m+1; i++) {
                M[i] = n;
            }
            int d = 0;
            for (int k=0; k<t; k++) {
                d += T[k].num*(T[k].end-T[k].begin);
            }
            int maxprofit = 0;
            fun(T, d, M, maxprofit, 0, t, 0);
            cout<<maxprofit<<endl;
            delete [] T;
            delete [] M;
            d = 0;
        }
    }
    return 0;
}

// 6B (AC) 青蛙跳石头问题（第六章题目）
#include <iostream>
#include <stdio.h>
#include <math.h>
#include<iomanip>

using namespace std;

struct stone{
    double x;
    double y;
};

double distance(stone a, stone b) {
    return sqrt((a.x-b.x)*(a.x-b.x) + (a.y-b.y)*(a.y-b.y));
}

double max_2(double x,double y){
    if (x>y) {
        return x;
    }else{
        return y;
    }
}

void dijkstra(int n, bool *vis, double *dis, double **cost){
    int s = 0;
    int count = n - 1;
    while(count --){
        int k;
        float min_dis = 999999;
        for(int i = 1; i < n; i ++)
            if(!vis[i]){
                if(dis[i] > max(dis[s], cost[i][s]))
                    dis[i] = max(dis[s], cost[i][s]);
                if(min_dis > dis[i]){
                    min_dis = dis[i];
                    k = i;
                }
            }
        if(k == 1)  return;
        s = k;
        vis[k] = true;
    }
}

int main() {
    int n;
    while (cin>>n) {
        if ((n==0)||(n>=100)) {
            break;
        }else if (n<=2){
            stone *stones = new stone[n];
            for (int i=0; i<n; i++) {
                cin>>stones[i].x>>stones[i].y;
            }
            cout<<setiosflags(ios::fixed)<<setprecision(3)<<distance(stones[0], stones[n-1])<<endl;
            delete []stones;
        }else{
            stone *stones = new stone[n];
            for (int i=0; i<n; i++) {
                cin>>stones[i].x>>stones[i].y;
            }

            double **cost = new double *[n];
            for(int i=0;i<n;i++) {
                cost[i]=new double [n];
            }

            double dis[n];
            bool vis[n];
            for (int i=0; i<n; i++) {
                dis[i] = distance(stones[0], stones[i]);
                for (int j=i; j<n; j++) {
                    double t = distance(stones[i], stones[j]);
                    cost[i][j] = t;
                    cost[j][i] = t;
                }
            }
            for(int i=0; i<=n; i++)
            {
                vis[i] = false;
            }
            dijkstra(n, vis, dis, cost);
            cout<<setiosflags(ios::fixed)<<setprecision(3)<<dis[1]<<endl;
            delete []stones;
            for (int i=0; i<n; i++) {
                delete [] cost[i];
            }
        }
    }
    return 0;
}

// 6B (TLE)
#include <iostream>
#include <stdio.h>
#include <math.h>
#include<iomanip>
using namespace std;

struct stone{
    double x;
    double y;
};

double distance(stone a, stone b) {
    return sqrt((a.x-b.x)*(a.x-b.x) + (a.y-b.y)*(a.y-b.y));
}

bool check(bool *p, int i){
    if (p[i]==0) {
        p[i] = 1;
        return true;
    }else{
        return false;
    }
}

void backtrack(stone *stones, int n, double &dismin, double dismax, int i, bool *p){
    for (int j=i+1; j<n; j++) {
        swap(stones[i+1], stones[j]);
        double temp = dismax;
        dismax = max(dismax,distance(stones[i], stones[i+1]));
        if ((dismax<dismin)) {
            if (j!=n-1) {
                backtrack(stones, n, dismin, dismax, i+1, p);
            }else{
                dismin = min(dismin, dismax);
            }
        }
        dismax = temp;
        swap(stones[i+1], stones[j]);
    }
}

int main() {
    int n;
    while (cin>>n) {
        if ((n==0)||(n>=100)) {
            break;
        }else if (n<=2){
            stone *stones = new stone[n];
            for (int i=0; i<n; i++) {
                cin>>stones[i].x>>stones[i].y;
            }
            cout<<setiosflags(ios::fixed)<<setprecision(3)<<distance(stones[0], stones[n-1])<<endl;
            delete []stones;
        }else{
            stone *stones = new stone[n];
            bool *p = new bool[n];
            for (int i=0; i<n; i++) {
                cin>>stones[i].x>>stones[i].y;
                p[i] = 0;
            }
            swap(stones[1], stones[n-1]);
            double dismin = distance(stones[0], stones[n-1]);
            backtrack(stones, n, dismin, 0, 0, p);
            cout<<setiosflags(ios::fixed)<<setprecision(3)<<dismin<<endl;
            delete []stones;
            delete []p;
        }
    }
    return 0;
}

// 6C 王晓东 无优先级运算问题（第六章题目）
#include <iostream>
#include <stdio.h>
using namespace std;

bool found(int k, int m, int *num, int *oper){
    int x = num[0];
    for (int i=0; i<k; i++) {
        switch (oper[i]) {
            case 0:
                x += num[i+1];
                break;
            case 1:
                x -= num[i+1];
                break;
            case 2:
                x *= num[i+1];
                break;
            case 3:
                x /= num[i+1];
                break;
            default:
                break;
        }
    }
    return (x==m);
}

bool search(int dep, int n, int m, int k, int *a, int *num, int *flag, int *oper){
    if (dep>k) {
        if (found(k,m,num,oper)) {
            return true;
        }else{
            return false;
        }
    }
    for (int i=0; i<n; i++) {
        if (flag[i]==0) {
            num[dep] = a[i];
            flag[i] = 1;
            for (int j=0; j<4; j++) {
                oper[dep] = j;
                if (search(dep+1,n,m,k,a,num,flag,oper)) {
                    return true;
                }
            }
            flag[i] = 0;
        }
    }
    return false;
}

int main(){
    int n,m;
    while (cin>>n>>m) {
        if (n==0) {
            break;
        }
        bool panduan = true;
        int a[n];
        int num[n];
        int oper[n];
        int flag[n];
        for (int i=0; i<n; i++) {
            cin>>a[i];
            flag[i] = 0;
        }
        for (int k=0; k<n; k++) {
            if (search(0,n,m,k,a,num,flag,oper)) {
                cout<<k<<endl;
                panduan = false;
                break;
            }
        }
        if (panduan) {cout<<"No Solution!"<<endl;}
    }
    return 0;
}

// 6C (Time Limit Exceed)
#include <iostream>
#include <stdio.h>
#include <queue>
using namespace std;

double fun(int opt, int x, int y){
    switch (opt) {
        case 0:
            return x+y;
        case 1:
            return x*y;
        case 2:
            return y-x;
        case 3:
            return y/x;
        default:
            return 0;
    }
}

int bfs(int n, int m, int now, int *x, int times){
    queue <int> q;
    queue <int> q2;
    for (int s = 0; s<n; s++) {
        int * process = new int [n];
        int ans;
        for (int i=0; i<n; i++) {
            process[i] = 0;
        }
        ans = fun(0, x[s], 0);
        if (ans==m) {
            times = 0;
            break;
        }
        process[s] = 1;
        for (int i=0; i<n; i++) {
            q.push(process[i]);
        }
        q2.push(ans);
        delete [] process;
    }
    while (!q.empty() && times==n) {
        int *node = new int [n];
        for (int i=0; i<n; i++) {
            node[i] = q.front();
            q.pop();
        }
        int a = q2.front();
        q2.pop();
        for (int i=0; i<n; i++) {
            if (node[i]==0) {
                node[i] = 1;
                for (int k=0; k<n; k++) {
                    q.push(node[k]);
                }
                for (int j=0; j<4; j++) {
                    int ans;
                    ans = fun(j, x[i], a);
                    if (ans==m) {
                        int temp = 0;
                        for (int k=0; k<n; k++) {
                            temp += (node[k] != 0);
                        }
                        times = temp - 1;
                        break;
                    }else{
                        q2.push(ans);
                    }

                }
                if (times!=n) {
                    break;
                }
                node[i] = 0;
            }
            if (times != n) {
                break;
            }
        }
        delete [] node;
    }
    while (!q.empty()){
        q.pop();
    }
    while (!q2.empty()){
        q2.pop();
    }
    return times;
}

int main(){
    int n;
    int m;
    while (cin>>n>>m) {
        if (n==0) {
            break;
        }else{
            int *x = new int [n];
            for (int i=0; i<n; i++) {
                cin>>x[i];
            }
            int times = n;
            int now = 0;
            times = bfs(n, m, now, x, times);
            if (times==n) {
                cout<<"No Solution!"<<endl;
            }else{
                cout<<times<<endl;
            }
            delete []x;
        }
    }
    return 0;
}

// 7A (AC) 矩阵乘法检验(第七章题目)
#include <iostream>
#include <stdio.h>
#include <random>

using namespace std;
using std::default_random_engine;

int main() {
    int n;
    while (cin>>n) {
        if ((n==0)||(n>500)) {
            break;
        }else{
            long long int **A = new long long int *[n];
            for(int i=0;i<n;i++) {
                A[i]=new long long  int [n];
            }
            long long int **B = new long long  int *[n];
            for(int i=0;i<n;i++) {
                B[i]=new long long  int [n];
            }
            long long  int **C = new long long  int *[n];
            for(int i=0;i<n;i++) {
                C[i]=new long long  int [n];
            }
            for (int i=0; i<n; i++) {
                for (int j=0; j<n; j++) {
                    cin>>A[i][j];
                }
            }
            for (int i=0; i<n; i++) {
                for (int j=0; j<n; j++) {
                    cin>>B[i][j];
                }
            }
            for (int i=0; i<n; i++) {
                for (int j=0; j<n; j++) {
                    cin>>C[i][j];
                }
            }

            int t=0;

            default_random_engine e;
            bool *bo = new bool[n*n];
            for (int p=0; p<n*n; p++) {
                bo[p] = 0;
            }
            int num = n*n;
            while (num>0) {
                int p = (e())%(n*n);
                if (bo[p]==0) {
                    bo[p] = 1;
                    num--;
                    int i = p/n;
                    int j = p%n;
                    long long int sum = 0;
                    for (int k=0; k<n; k++) {
                        sum += A[i][k]*B[k][j];
                    }
                    if (sum!=C[i][j]) {
                        t++;
                        cout<<"NO!"<<endl;
                        break;
                    }
                }
            }

            for (int i=0; i<n; i++) {
                delete [] A[i];
                delete [] B[i];
                delete [] C[i];
            }
            if (t==0) {
                cout<<"YES!"<<endl;
            }
        }
    }
    return 0;
}

// 7B (AC) 矩阵乘法检验(第七章题目)
#include <iostream>
#include <stdio.h>
#include <random>

using namespace std;
using std::default_random_engine;

void power(unsigned int a, unsigned int p, unsigned int n, unsigned int &result, bool &cp)
{
    unsigned int x;
    if (p==0) {
        result=1;
    }
    else {
        power(a,p/2,n,x,cp);
        result=(x*x)%n;
        if ((result==1)&&(x!=1)&&(x!=n-1))
            cp=true;
        if ((p%2)==1)
            result=(result*a)%n;
    }
}

int main(){
    default_random_engine e;
    int n;
    while (cin>>n) {
        if (n==0) {
            break;
        }else{
            bool cp = false;
            unsigned int result;
            bool b[n-1];
            for (int i=0; i<=n-1; i++) {
                b[i] = 0;
            }
            int num = n-1;
            bool panduan = 1;
            while (num > 0) {
                unsigned int a = (e())%(max(1,n-1)) + 1;
                if (b[a]==0) {
                    b[a] = 1;
                    num--;
                    power(a, n-1, n, result, cp);
                    if (cp||(result!=1)) {
                        panduan = 0;
                        break;
                    }
                }

            }
            if (panduan) {
                cout<<"YES!"<<endl;
            }else{
                cout<<"NO!"<<endl;
            }
        }
    }
    return 0;
}


// 7C (AC) 完全图最大割问题(第七章题目)
#include <iostream>
#include <stdio.h>
using namespace std;

int check(int i,int n,bool *a,int **c){
    int sum = 0;
    for (int j=0; j<n; j++) {
        if (a[i] != a[j]) {
            sum += c[j][i];
        }else{
            sum -= c[j][i];
        }
    }
    return sum;
}

void bt (int n, int temp, int &sum, bool *a, int **c, int k){
    sum = max(sum, temp);
    if (k<n) {
        for (int i=0; i<=1; i++) {
            if (i==0) {
                a[k] = i;
                int p = check(k,n,a,c);
                temp += p;
                bt(n, temp, sum, a, c, k+1);
                temp -= p;
            }else{
                a[k] = i;
                bt(n, temp, sum, a, c, k+1);
            }
        }
    }
}

int main(){
    int n;
    while (cin>>n) {
        if (n==0) {
            break;
        }else{
            int **C = new int *[n];
            for(int i=0;i<n;i++) {
                C[i]=new int [n];
            }
            for (int i=0; i<n; i++) {
                for (int j=0; j<n; j++) {
                    cin>>C[i][j];
                }
            }

            int sum = 0;
            bool *a = new bool [n];
            for (int i=0; i<n; i++) {
                a[i] = true;
            }
            bt(n, sum, sum, a, C, 0);
            cout<<sum<<endl;
            for (int i=0; i<n; i++) {
                delete [] C[i];
            }
            delete [] a;
        }
    }
    return 0;
}
