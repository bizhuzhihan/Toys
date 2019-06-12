//// Zhihan ZHU
//// 201630845294

//// 2A (AC)
//#include <iostream>
//#include <stdio.h>
//using namespace std;
//
//int main() {
//    int n;
//    while (cin >> n) {
//        if (n > 0) {
//            int N = 1;
//            int N_1 = 1;
//            int N_f = 1;
//            for (int i = 2; i < n; i++) {
//                N_f = (N + N_1) % 1997;
//                N_1 = N;
//                N = N_f;
//            }
//            cout << N_f << '\n';
//        } else {
//            break;
//        }
//    }
//    return 0;
//}

//// 2A devide and conquer (Memory Limit Exceed)
//#include <iostream>
//#include <stdio.h>
//using namespace std;
//
//int *fun(int n){
//    int a[2];
//    a[0] = 1;
//    a[1] = 1;
//    if (n>2) {
//        int *p = fun(n-1);
//        a[0] = (*p + *(p+1)) % 1997;
//        a[1] = *p % 1997
//        ;
//    }
//    return a;
//}
//
//int main() {
//    int n;
//    while (cin >> n) {
//        if (n > 0) {
//            cout << *fun(n) << '\n';
//        } else {
//            break;
//        }
//    }
//    return 0;
//}

//// 2B (AC)
//#include <iostream>
//#include <stdio.h>
//using namespace std;
//
//long long int fun(long long int a,long long int b,long long int c) {
//    long long int value = 1;
//    if (b>1) {
//        if (b>(2*(b/2))) {
//            value = fun(a, b/2, c);
//            value *= (value * a);
//            value = value % c;
//        }else{
//            value = fun(a, b/2, c);
//            value *= value;
//            value = value % c;
//        }
//    }else{
//        value = a % c;
//    }
//    return value;
//}
//
//int main() {
//    long long int a, b, c;
//    while (cin >> a >> b >> c) {
//        if ((a>=0) & (b>=0) & (c>1)) {
//            long long int r;
//            if (b == 0) {
//                r = 1;
//            }else{
//                r = fun(a, b, c);
//            }
//            cout << r << '\n';
//        } else {
//            break;
//        }
//    }
//    return 0;
//}

//// 2C (AC)
//#include <iostream>
//#include <stdio.h>
//using namespace std;
//
//int pfun(int Y[], int p, int r){
//    int i=p,j=1+r;
//    int x=Y[p];
//    while (true) {
//        while (Y[++i]<x&&i<r);
//        while (Y[--j]>x);
//        if (i>=j) {
//            break;
//        }
//        swap(Y[i], Y[j]);
//    }
//    Y[p]=Y[j];
//    Y[j]=x;
//    return j;
//}
//
//void qsorting(int Y[], int p, int r, int n){
//    if (p<r) {
//        int q=pfun(Y,p,r);
//        if (q!=(n/2)) {
//            if (q<(n/2)) {
//                qsorting(Y, q+1, r, n);
//            }
//            if (q>(n/2)) {
//                qsorting(Y, p, q-1, n);
//            }
//        }
//    }
//}
//
//int main() {
//    int n;
//    cin >> n;
//    int Y[n];
//    int sum = 0;
//    for (int i=0; i<n; i++) {
//        int x,y;
//        cin >> x >> y;
//        if ((x>=(-10000))&(x<=10000)&(y>=(-10000))&(y<=10000)) {
//            Y[i] = y;
//        } else {
//            break;
//        }
//    }
//    qsorting(Y, 0, n-1, n);
//    for (int i=0; i<n; i++) {
//        if (Y[n/2] > Y[i]) {
//            sum += (Y[n/2] - Y[i]);
//        }else{
//            sum += (Y[i] - Y[n/2]);
//        }
//    }
//    cout<< sum << '\n';
//    return 0;
//}

//// 3A (AC)
//#include <iostream>
//#include <stdio.h>
//using namespace std;
//
//int LCSr(int n,char x[])
//{
//    char y[n];
//    for (int ii=0; ii<n; ii++) {
//        y[ii] = x[n-1-ii];
//    }
//    int i,j;
//    int c[n+1][n+1];
//    int b[n+1][n+1];
//    for (i=0; i<n+1; i++) {
//        c[i][0] = 0;
//    }
//    for (j=0; j<n+1;j++) {
//        c[0][j] = 0;
//    }
//    for (i=1; i<n+1; i++) {
//        for (j=1; j<n+1; j++) {
//            if (x[i-1]==y[j-1]) {
//                c[i][j] = c[i-1][j-1] + 1;
//                b[i][j] = 1;
//            }
//            else if (c[i-1][j]>=c[i][j-1]) {
//                c[i][j] = c[i-1][j];
//                b[i][j] = 2;
//            }
//            else {
//                c[i][j] = c[i][j-1];
//                b[i][j] = 3;
//            }
//        }
//    }
//    int X[n];
//    int Y[n];
//    for (int i=0; i<n; i++) {
//        X[i] = 0;
//        Y[i] = 0;
//    }
//    i = n;
//    j = n;
//    int counter = 0;
//    while (i>0&&j>0) {
//        if (b[i][j]==1) {
//            counter++;
//            X[--i] = counter;
//            Y[--j] = counter;
//        }else if (b[i][j]==2) {
//            i--;
//        }else {
//            j--;
//        }
//    }
//    int ans=0, ans2=0;
//    int p1=0, p2=0, p3=0;
//    int pos[2]={n, n+1};
//    for (; p1<n; p1++) {
//        if (X[p1]>0) {
//            for (; p2<=n-1-p1; p2++) {
//                if (Y[p2]>0) {
//                    ans++;
//                    p2++;
//                    break;
//                }
//            }
//            for (; p3<n-1-p1; p3++) {
//                if (Y[p3]>0) {
//                    ans2++;
//                    p3++;
//                    break;
//                }
//            }
//        }
//    }
//    int Ans;
//    if ((n-2*ans2)<(n+1-2*ans)) {
//        Ans = n - 2*ans2;
//    }else {
//        Ans = n + 1 - 2*ans;
//    }
//    return Ans;
//}
//
//int main() {
//    int n;
//    while (cin >> n) {
//        if (n<=0) {
//            break;
//        }
//        char x[n];
//        for (int i=0; i<n; i++) {
//            cin >> x[i];
//        }
//        cout<<LCSr(n, x)<<'\n';
//    }
//    return 0;
//}

//// 3B (AC)
//#include <iostream>
//#include <stdio.h>
//using namespace std;
//
//void stones(long long int p[], long long int n, long long int *ans) {
//    long long int sum[n][n], max[n][n], min[n][n];
//    for (long long int i=0; i<n; i++) {
//        sum[i][i] = p[i];
//        max[i][i] = 0;
//        min[i][i] = 0;
//    }
//    for (long long int l=2; l<=n; l++) {
//        for (long long int i=0; i<n; i++) {
//            long long int j = (i+l-1) ;
//            sum[i][j%n] = sum[(i+1)%n][j%n] + p[i];
//            max[i][j%n] = max[(i+1)%n][j%n] + sum[i][j%n];
//            min[i][j%n] = min[(i+1)%n][j%n] + sum[i][j%n];
//            for (long long int k=(i+1); k<j; k++) {
//                long long int M = max[i][k%n] + max[(k+1)%n][j%n] + sum[i][j%n];
//                long long int m = min[i][k%n] + min[(k+1)%n][j%n] + sum[i][j%n];
//                if (M>max[i][j%n]) {
//                    max[i][j%n] = M;
//                }
//                if (m<min[i][j%n]) {
//                    min[i][j%n] = m;
//                }
//            }
//        }
//    }
//    long long int Min = min[0][n-1];
//    long long int Max = max[0][n-1];
//    for (long long int i=1; i<n; i++) {
//        if (Max<max[i][(i+n-1)%n]) {
//            Max = max[i][(i+n-1)%n];
//        }
//        if (Min>min[i][(i+n-1)%n]) {
//            Min = min[i][(i+n-1)%n];
//        }
//    }
//    ans[0] = Min;
//    ans[1] = Max;
//}
//
//int main() {
//    long long int n;
//    while (cin>>n) {
//        if (n==0) {
//            break;
//        }
//        long long int p[n];
//        for (long long int i=0; i<n; i++) {
//            cin>>p[i];
//        }
//        long long int ans[2]={0,0};
//        stones(p, n, ans);
//        if (n==1) {
//            cout<<p[0]<<' '<<p[1]<<'\n';
//        }
//        cout<<ans[0]<<' '<<ans[1]<<'\n';
//    }
//    return 0;
//}

//// 3C (AC)
//#include <iostream>
//#include <stdio.h>
//using namespace std;
//
//int qsorting(int Y[], int p, int r){
//    int i = p, j = 1+r;
//    int x = Y[p];
//    while (true) {
//        while (Y[++i]<x&&i<r);
//        while (Y[--j]>x);
//        if (i>=j) {
//            break;
//        }
//        swap(Y[i], Y[j]);
//    }
//    Y[p] = Y[j];
//    Y[j] = x;
//    return j;
//}
//
//void Qsorting(int Y[], int p, int r){
//    if (p<r) {
//        int q = qsorting(Y,p,r);
//        Qsorting(Y, q+1, r);
//        Qsorting(Y, p, q-1);
//    }
//}
//
//int LCS(int m,int n,int x[],int y[])
//{
//    int i,j;
//    int c[m+1][n+1];
//    for (i=0; i<m+1; i++) {
//        c[i][0] = 0;
//    }
//    for (j=0; j<n+1;j++) {
//        c[0][j] = 0;
//    }
//    for (i=1; i<m+1; i++) {
//        for (j=1; j<n+1; j++) {
//            if (x[i-1]==y[j-1]) {
//                c[i][j] = c[i-1][j-1]+1;
//            }
//            else if (c[i-1][j]>=c[i][j-1]) {
//                c[i][j] = c[i-1][j];
//            }
//            else {
//                c[i][j] = c[i][j-1];
//            }
//        }
//    }
//    return c[m][n];
//}
//
//int main() {
//    int n;
//    while (cin>>n) {
//        if (n>0) {
//            int a0[n], a[n];
//            for (int i=0; i<n; i++) {
//                cin >> a0[i];
//                a[i] =  a0[i];
//            }
//            Qsorting(a, 0, n-1);
//            int i = 0;
//            int j = 0;
//            while (i<n-1) {
//                if (a[i]!=a[i+1]) {
//                    a[j] = a[i];
//                    j++;
//                }
//                i++;
//            }
//            a[j] = a[i];
//            cout << LCS(n, j+1, a0, a) <<'\n';
//        }else {
//            break;
//        }
//    }
//    return 0;
//}

//// 3extra
//#include <iostream>
//#include <stdio.h>
//using namespace std;
//
//void LCS(int m,int n,char x[],char y[])
//{
//    int i,j;
//    int c[m+1][n+1];
//    int b[m+1][n+1];
//    for (i=0; i<m+1; i++) {
//        c[i][0] = 0;
//    }
//    for (j=0; j<n+1;j++) {
//        c[0][j] = 0;
//    }
//    for (i=1; i<m+1; i++) {
//        for (j=1; j<n+1; j++) {
//            if (x[i-1]==y[j-1]) {
//                c[i][j] = c[i-1][j-1]+1;
//                b[i][j] = 1;
//            }
//            else if (c[i-1][j]>=c[i][j-1]) {
//                c[i][j] = c[i-1][j];
//                b[i][j] = 2;
//            }
//            else {
//                c[i][j] = c[i][j-1];
//                b[i][j] = 3;
//            }
//        }
//    }
//
//    int x_pos = m;
//    int y_pos = n;
//    while ((x_pos>=0) & (y_pos>=0)) {
//        if (b[x_pos][y_pos]==1) {
//            cout<<x[x_pos-1];
//            x_pos -= 1;
//            y_pos -= 1;
//        }else if (b[x_pos][y_pos]==2) {
//            x_pos -= 1;
//        }else if (b[x_pos][y_pos]==3) {
//            y_pos -= 1;
//        }
//    }
//}
//
//int main() {
//    char a[] = "adfghj";
//    char b[] = "asdfuhjk";
//    LCS(strlen(a), strlen(b), a, b);
//    return 0;
//}


//// 4A (AC)
//#include <iostream>
//#include <stdio.h>
//using namespace std;
//
//int main() {
//    int n,m;
//    while (cin>>n>>m) {
//        if (n==0&&m==0) {
//            break;
//        }
//        int month[12];
//        for (int i=0; i<12; i++) {
//            month[i] = 0;
//        }
//        int max = 0;
//        for (int i=4; i>0; i--) {
//            if ((i*n-(5-i)*m)<0) {
//                max = i;
//                break;
//            }
//        }
////        cout<<'m'<<max<<'\n';
//        int j = 0;
//        for (; j<max; j++) {
//            month[j] = 1;
//        }
//        for (; j<5; j++) {
//            month[j] = 0;
//        }
//        for (; j<12; j++) {
//            int sum = n;
//            for (int k=j-1; k>=j-4; k--) {
//                if (month[k]==0) {
//                    sum -= m;
//                }else {
//                    sum += n;
//                }
//            }
//            if (sum<0) {
//                month[j] = 1;
//            }else {
//                month[j] = 0;
//            }
//        }
//        int profit = 0;
//        for (int i=0; i<12; i++) {
//            if (month[i]==0) {
//                profit -= m;
//            }else{
//                profit += n;
//            }
//        }
//        if (profit>0) {
//            cout << profit << '\n';
//        }else{
//            cout << "NO" << '\n';
//        }
//    }
//    return 0;
//}

//// 4B (AC)
//#include <iostream>
//#include <stdio.h>
//using namespace std;
//
//int pfun(int Y[], int p, int r){
//    int i=p,j=1+r;
//    int x=Y[p];
//    while (true) {
//        while (Y[++i]<x&&i<r);
//        while (Y[--j]>x);
//        if (i>=j) {
//            break;
//        }
//        swap(Y[i], Y[j]);
//    }
//    Y[p]=Y[j];
//    Y[j]=x;
//    return j;
//}
//
//void qsorting(int Y[], int p, int r){
//    if (p<r) {
//        int q=pfun(Y,p,r);
//        qsorting(Y, q+1, r);
//        qsorting(Y, p, q-1);
//    }
//}
//
//int main() {
//    int m,n;
//    while (cin>>m>>n) {
//        if (n==0&&m==0) {
//            break;
//        }
//        int u[n];
//        for (int i=0; i<n; i++) {
//            u[i] = 0;
//        }
//        for (int i=0; i<n; i++) {
//            cin >> u[i];
//        }
//        qsorting(u, 0, n-1);
//        int gap[n];
//        for (int i=0; i<n; i++) {
//            gap[i] = 0;
//        }
//        gap[n-1] = n*m - u[n-1];
//        for (int i=n-2; i>=0; i--) {
//            if (gap[i+1]==0) {
//                gap[i] = u[i+1] - u[i] + gap[i+1] - 1;
//            }else{
//                gap[i] = u[i+1] - u[i] + gap[i+1] - 2;
//            }
//        }
////        for (int i=0; i<n; i++) {
////            cout<<u[i]<<' ';
////        }
////        cout<<'\n';
////        for (int i=0; i<n; i++) {
////            cout<<gap[i]<<' ';
////        }
////        cout<<'\n';
//        int time = 0;
//        for (int i=0; i<n; i++) {
//            time += 1*(gap[i]==0);
//        }
//        cout << time << '\n';
//    }
//    return 0;
//}

//// 4C (AC)
//#include <iostream>
//#include <stdio.h>
//using namespace std;
//
//int pfun(int Y[], int p, int r){
//    int i=p,j=1+r;
//    int x=Y[p];
//    while (true) {
//        while (Y[++i]<x&&i<r);
//        while (Y[--j]>x);
//        if (i>=j) {
//            break;
//        }
//        swap(Y[i], Y[j]);
//    }
//    Y[p]=Y[j];
//    Y[j]=x;
//    return j;
//}
//
//void qsorting(int Y[], int p, int r){
//    if (p<r) {
//        int q=pfun(Y,p,r);
//        qsorting(Y, q+1, r);
//        qsorting(Y, p, q-1);
//    }
//}
//
//int main() {
//    int n;
//    while (cin>>n) {
//        if (n==0) {
//            break;
//        }
//        int m[n],M[n];
//        for (int i=0; i<n; i++) {
//            int v;
//            cin >> v;
//            m[i] = v;
//            M[i] = v;
//        }
//        qsorting(m, 0, n-1);
//        int min = 0;
//        for (int i=1; i<n; i++) {
//            m[i] += m[i-1];
//            min += m[i] - 1;
//            qsorting(m, i, n-1);
//        }
//        qsorting(M, 0, n-1);
//        int max = 0;
//        for (int i=n-2; i>=0; i--) {
//            M[i] += M[i+1];
//            max += M[i] - 1;
//        }
//        cout<<max<<' '<<min<<'\n';
//    }
//    return 0;
//}
