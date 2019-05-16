//// Zhihan ZHU
//// 201630845294

//// 2A
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

//// 2B
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

//// 2C
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

//// 3A
//#include <iostream>
//#include <stdio.h>
//using namespace std;
//
//int LCS(int m,int n,char x[],char y[])
//{
//    int i,j;
//    int c[m+1][n+1];
////    int b[m+1][n+1];
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
////                b[i][j] = 1;
//            }
//            else if (c[i-1][j]>=c[i][j-1]) {
//                c[i][j] = c[i-1][j];
////                b[i][j] = 2;
//            }
//            else {
//                c[i][j] = c[i][j-1];
////                b[i][j] = 3;
//            }
//        }
//    }
//    return c[m][n];
//}
//
//int main() {
//    int n;
//    while (cin >> n) {
//        if (n<=0) {
//            break;
//        }
//        char x[n];
//        int max1 = 0;
//        int max2 = 0;
//        for (int i=0; i<n; i++) {
//            cin >> x[i];
//        }
//        for (int i=0; i<n; i++) {
//            char leftR[i+1];
//            char right[n-i];
//            char right2[n-i-1];
////            cout<<"i="<<i<<'\t';
////            cout<<"left="<<i+1<<'\t';
////            cout<<"right="<<n-i<<'\n';
//            for (int j=0; j<(i+1); j++) {
//                leftR[j] = x[i-j];
//            }
////            cout<<leftR<<'\t';
//            for (int j=0; j<(n-i); j++) {
//                right[j] = x[i+j];
//            }
////            cout<<right<<'\t';
//            for (int j=0; j<(n-i-1); j++) {
//                right2[j] = x[i+j+1];
//            }
////            cout<<right2<<'\n';
////            cout<<LCS(i+1,n-i,leftR,right)<<'\n';
//            int t = LCS(i+1,n-i,leftR,right);
//            if (max1 < t) {
////                max = (n+(n%2))-2*LCS(i+1,n-i,leftR,right);
//                max1 = t;
//            }
//            int t2 = LCS(i+1,n-i-1,leftR,right2);
//            if (max2 < t2) {
//                max2 = t2;
//            }
//        }
//        int lcs1 = n+1 - 2*max1;
//        int lcs2 = n - 2*max2;
//        if (lcs1<lcs2) {
//            cout << lcs1 << '\n';
//        }else{
//            cout << lcs2 << '\n';
//        }
//    }
//    return 0;
//}

//// 3B
//#include <iostream>
//#include <stdio.h>
//using namespace std;
//
//int *stones(int p[], int n) {
//    int sum[n][n], max[n][n], min[n][n];
//    for (int i=0; i<n; i++) {
//        sum[i][i] = p[i];
//        max[i][i] = 0;
//        min[i][i] = 0;
//    }
//    for (int l=2; l<=n; l++) {
//        for (int i=0; i<n; i++) {
//            int j = (i+l-1) % n;
//            sum[i][j] = sum[(i+1)%n][j] + p[i];
//            max[i][j] = max[(i+1)%n][j] + sum[i][j];
//            min[i][j] = min[(i+1)%n][j] + sum[i][j];
//            for (int k=(i+1)%n; k<j; k++) {
//                int M = max[i][k] + max[(k+1)%n][j] + sum[i][j];
//                int m = min[i][k] + min[(k+1)%n][j] + sum[i][j];
//                if (M>max[i][j]) {
//                    max[i][j] = M;
//                }
//                if (m<min[i][j]) {
//                    min[i][j] = m;
//                }
//            }
//        }
//    }
//    int Min = min[0][n-1];
//    int Max = max[0][n-1];
//    for (int i=1; i<n; i++) {
//        if (Max<max[i][(i+n-1)%n]) {
//            Max = max[i][(i+n-1)%n];
//        }
//        if (Min>min[i][(i+n-1)%n]) {
//            Min = min[i][(i+n-1)%n];
//        }
//    }
//    int answer[2] = {Min, Max};
//    return answer;
//}
//
//int main() {
//    int n;
//    while (cin>>n) {
//        if (n==0) {
//            break;
//        }
//        int p[n];
//        for (int i=0; i<n; i++) {
//            cin>>p[i];
//        }
//        cout<<*stones(p, n)<<' '<<*(stones(p, n)+1)<<'\n';
//    }
//    return 0;
//}

//// 3C
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
//            cout << LCS(n, j+1, a0, a) <<'\n';
//        }else {
//            break;
//        }
//    }
//    return 0;
//}
