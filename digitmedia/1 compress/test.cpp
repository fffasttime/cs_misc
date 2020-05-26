/*
Name: AC_template_main
Author: fffasttime 
Date: 
Description: 
*/
//#pragma comment(linker, "/STACK:1024000000,1024000000")
//#pragma GCC optimize("Ofast,no-stack-protector")
//#pragma GCC target("sse,sse2,sse3,ssse3,sse4,popcnt,abm,mmx,avx,tune=native")
//#pragma GCC optimize("unroll-loops")

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <vector>
#include <set>
#include <map>
#include <queue>
#include <complex>
#include <cassert>
#include <algorithm>
using namespace std;
#define USE_ATTR
#ifdef USE_ATTR
#define inc(i,n) for (int i=0;i<n;i++)
#define icc(i,n) for (int i=1;i<=n;i++)
#define dec(i,n) for (int i=n-1;i>=0;i--)
#define dcc(i,n) for (int i=n;i>0;i--)
#define rep(i,a,b) for (int i=a;i<b;i++)
#define rpp(i,a,b) for (int i=a;i<=b;i++)
#define per(i,b,a) for (int i=b-1;i>=a;i--)
#define prr(i,b,a) for (int i=b;i>=a;i--)

#define MS(a,x) memset(a,x,sizeof(a))
//we can use initiallist in c++1x
#define MP make_pair
#define PII pair<int,int>
//std=c++11
#define MT make_tuple
#define TIII tuple<int, int, int>

#endif

typedef long long ll;
typedef double db;
typedef ll T;

//warning: Can't use other input method while using fread.
int rc(){
	//return getchar(); //if use fread, input won't stop until EOF
	static char buf[10000],*p1=buf,*p2=buf;
	return p1==p2&&(p2=(p1=buf)+fread(buf,1,10000,stdin),p1==p2)?EOF:*p1++;
}
int rd(int &x){
	x=0; int f=1,c=rc();
	while (!isdigit(c) && c!=-1) c=='-'?f=-1:0,c=rc();
	if (c==-1) return 0;
	while (isdigit(c)) x=x*10+c-'0', c=rc();
	x*=f; return 1;
}
int rd(char *x){
	int c=rc();
	while (isspace(c) && c!=-1) c=rc();
	if (c==-1) return 0;
	while (!isspace(c) && c!=-1) *(x++)=c,c=rc();
	*x=0; return 1;
}

int qrand(){
	static int seed=2333;
	return seed = (int)((((seed ^ 998244353)+19260817ll)*12345678ll)%1000000007);
}

T gcd(T a, T b){ return b==0?a:gcd(b, a%b);}
int gcd(int a, int b) {return b?gcd(b,a%b):a;}

ll qmul(ll x,ll y,ll p){
	ll t=(x*y-(ll)((long double)x/p*y+0.5)*p);
	return t<0 ? t+p : t;
}
//return a^x
T qpow(T a, int x){
	T ans=1;
	for (;x;a*=a,x>>=1)
		if (x&1) ans*=a;
	return ans;
}
ll qpow(ll a, ll x, ll p){
	ll ans=1;
	for (;x;a=qmul(a,a,p),x>>=1)
		if (x&1) ans=qmul(ans,a,p);
	return ans;
}

const int N=100;
struct Mat{
	ll m[N][N];
	Mat(){memset(m,0,sizeof(m));}
	void I(){for (int i=0;i<N;i++) m[i][i]=1;}
};
Mat mul(const Mat &a, const Mat &b){
	Mat c;
	for (int i=0;i<N;i++)
		for (int j=0;j<N;j++){
			for (int k=0;k<N;k++)
				c.m[i][j]+=a.m[i][k]*b.m[k][j];
			//c.m[i][j]%=p;
		}
	return c;
}
Mat qpow(Mat a, int x){
	Mat ans; ans.I();
	for (;x;a=mul(a,a),x>>=1)
		if (x&1)
			ans=mul(ans,a);
	return ans;
}

//require p is a prime
int inv(int x, int p){
	return qpow(x,p-2,p);
}

namespace NumTheory{

const int maxn=1000010;
int p[maxn],phi[maxn],mu[maxn],pc;
bool isp[maxn];
void genprime(){ //gen prime in range [2,maxn)
	memset(isp,1,sizeof(isp));
	isp[1]=0;
	for (int i=2;i<maxn;i++){
		if (isp[i]) p[pc++]=i;
		for (int j=0;j<pc && i*p[j]<maxn;j++){
			isp[i*p[j]]=0;
			if (i%p[j]==0) break;
		}
	}
}
int d[maxn], c0[maxn]; //d:count of frac,  c0: temp array, i=p1^c0[i] * p2^(...) * ...
void gennumfunc(){
	memset(isp,1,sizeof(isp));
	isp[1]=0;  phi[1]=mu[1]=d[1]=1;
	for (int i=2;i<maxn;i++){
		if (isp[i]) p[pc++]=i,phi[i]=i-1,mu[i]=-1,d[i]=2,c0[i]=1;
		for (int j=0;j<pc && i*p[j]<maxn;j++){
			isp[i*p[j]]=0;
			if (i%p[j]==0){
				phi[i*p[j]]=phi[i]*p[j];
				mu[i*p[j]]=0;
				c0[i*p[j]]=c0[i]+1;
				d[i*p[j]]=d[i]/(c0[i]+1)*(c0[i]+2);
				break;
			}
			phi[i*p[j]]=phi[i]*(p[j]-1); //f(ab)=f(a)f(b), when (a,b)=1
			mu[i*p[j]]=-mu[i];
			c0[i*p[j]]=c0[p[j]];
			d[i*p[j]]=d[i]*2;
		}
	}
}

//sum_i([n/i]) , O(sqrt(n)) 
ll sum_floor(ll n){
	ll ans=0;
	for(int l=1,r;l<=n;l=r){
	    r=n/(n/l)+1; //in section [l,r), floor(n/i) has same value
	    ans+=(r-l)*(n/l); //*(sphi[r]-sphi[l]); //sum([n/i]*phi[i])
	}
	return ans;
}
//another algo of sum_floor: sum(k=1, sqrtint(n), n\k)*2-sqrtint(n)^2

int invt[maxn];
void invTable(int maxl, int p){
	invt[1]=1;
	for (int i=2;i<=maxl;i++)
		invt[i]=invt[p%i]*(p-p/i)%p;
}

//ax+by=gcd(a,b)
ll exgcd(ll a, ll b, ll &x, ll &y){
	if (b==0){
		x=1;y=0;
		return a;
	}
	ll t=exgcd(b,a%b,y,x);
	y-=a/b*x;
	return t;
}
//not require p is prime, but (v,p) should be 1
ll inv_euclid(ll v, ll p){
	ll x,y; exgcd(v,p,x,y); //exgcd(v,p,x,y)==1 required
	return (x+p)%p;
}
//CRT, require (p1,p2,...)=1
//x=a1(mod p1)
//x=a2(mod p2)
//...
//result=sum(bi*Mi*Mi'), MiMi'=1(mod pi)
ll CRT(int n, ll a[], ll p[]){
	ll M=1,x=0;
	for (int i=0;i<n;i++) M*=p[i]; //[!] notice if M>int32, a*b%M may overflow
	for (int i=0;i<n;i++){
		ll w=M/p[i]; //x=pi*k1+a + w*k2
		x=(x+w*inv_euclid(w,p[i])%M*a[i])%M; //get k1, pi*k1=a (Mod w)
		//use inv_euclid() instead qpow() when p[] is not prime 
	}
	return (x+M)%M;
}

//CRT2 equ version
//x=a1(mod m1), x=a2(mod m2)
//result: x(mod p1p2)
int CRT2(int a1, int m1, int a2, int m2){
	int m=m1*m2;    //[!] notice if M>int32
	return (a1*m2%m*inv_euclid(m2,m1)+a2*m1%m*inv_euclid(m1,m2))%m;
}

//EXTCRT
//x=a1(mod p1)
//x=a2(mod p2)
//...
ll EXCRT(int n, ll a[], ll p[]){
	ll n1=p[0],a1=a[0],n2,a2,k1,k2,K,gcd,c,t;
	for (int i=1;i<n;i++){ //merge equs by order
		n2=p[i],a2=a[i]; 
		c=(a2-a1%n2+n2)%n2;
		gcd=exgcd(n1,n2,k1,k2); //n1*k1+n2*k2=gcd(n1,n2)
		if (c%gcd) return -1;
		t=n2/gcd; //n1*K+n2*(c/gcd*k2)=c
		K=c/gcd*k1%t; //K=qmul(c/gcd,k1,t); //if t too large
		a1+=n1*K; n1*=t;
		a1=(a1+n1)%n1;
	}
	return a1; //all answers are a1+LCM(p1,p2,..)
}

//get prime fact, x=mul(pi[i]^pa[i])
int pi[30],pa[30]; //no more 20 prime factor in range 2^64
int getfactor(int x){ //O(sqrt(n)), when x is a prime
	int c=0,m=sqrt(x)+1;
	for (int d=2;x>1 && d<=m;d++) 
	//for (int d=2,pc=0;x>1 && d<m;pc++,d=p[pc]) //faster, O(sqrt(n)/log(n))
		if (x%d==0){
			while (x%d==0) pa[c]++,x/=d;
			pi[c]=d; c++;
		}
	if (x>1) pi[c]=x, pa[c]=1, c++; //x is prime
	return c;
}
//single number phi, O(sqrt(n))
int getphi(int n){
	int ans=n;
	for (int i=2;i*i<=n;i++)
		if (n%i==0){
			ans=ans/i*(i-1);
			while (n%i==0) n/=i;
		}
	if (n>1) ans=ans/n*(n-1);
	return ans;
}

//be sure p=p^k,2p^k
//primitive_root(2)=1, primitive_root(4)=3, special judge
int primitive_root(int p){
	int phi=p-1; //int phi=getphi(p); //when p is not prime
	int pc=getfactor(phi);
	for (int g=2,j;;g++){ //g ~ p^0.25 in average
		for (j=0;j<pc;j++)
			if (qpow(g,phi/pi[j],p)==1)
				break;
		if (j==pc) return g;
	}
	//other solution of g: {g0^k | 0<k<phi, gcd(k, phi)=1}
}
//discrete logarithm
//solve a^x=b(mod p), require gcd(a,p)=1
//if a is primitive root and gcd(b,p)=1, the answer always exists
ll BSGS(ll a, ll b, ll p){
	int m,v,e=1;
	m=(int)sqrt(p+0.5);
	v=inv(qpow(a,m,p),p); //[!] use inv_euclid
	map<int,int> x; //unordered_map -> O(sqrt(N))
	x[1]=0;
	for (int i=1;i<m;i++){
		e=e*a%p;
		if (!x.count(e))
			x[e]=i;
	}
	for (int i=0;i<m;i++){
		if (x.count(b))
			return i*m+x[b];
		b=b*v%p;
	}
	return -1;
}
//find minimum x of a^x=b(mod m)
ll EXBSGS(ll a, ll b, ll m){
	if (b==1 || m==1) return 0;
	ll d,ca=1,k=0;
	while ((d=__gcd(a,m))>1){
		if (b%d) return -1;
		m/=d,b/=d; ca=ca*(a/d)%m; k++;
		if (ca==b) return k; //spj x<k
	}
	ll ret=BSGS(a,b*inv(ca,m)%m,m);
	if (ret==-1) return -1;
	return ret+k;
}

//judge if x^2=n (mod p) has solution
//qpow(n,(p-1)/2,p)=[1 nRp | p-2 n!Rp | 0 n==0]
bool isSquareRemain(int n, int p){
	return qpow(n,(p-1)/2,p)==1;
}

//in uint32, p0={2,7,61} is correct
//p0={2,3,7,61,24251}, only 46856248255981 will wrong
const ll p0[]={2,3,5,7,11,13,17,19,61,2333,24251};
//1. a^(p-1)=1 (mod p) 2. if x^2=1 (mod p) then x={1,p-1}
bool witness(ll a,ll n,ll r,ll s){
	ll x=qpow(a,r,n),pre=x;
	for (int i=0;i<s;i++){
		x=qmul(x,x,n);
		if (x==1 && pre!=1 && pre!=n-1) return 0;
		pre=x;
	}
	return x==1;
}
bool MillerRabin(ll n){
	if (n<=1) return 0;
	if (n==2 || n==3 || n==5 || n==7) return 1;
	if (n%2==0 || n%3==0 || n%5==0 || n%7==0) return 0;
	ll r=n-1,s=0;
	while (!(r&1)) r>>=1,s++;
	for (int i=0;i<10;i++){
		if (p0[i]==n) return 1;
		if (!witness(p0[i],n,r,s)) return 0;
	}
	return 1;
}

//pollard_rho factorization, O(sqrt(sqrt(n))) in expection
ll pol_rho(ll n,ll c){
	if (n%2==0) return 2; if (n%3==0) return 3;
	ll k=2,x=rand()%(n-1)+1,y=x,p=1,val=1;
	for (ll i=1;p==1;i++){
		x=(qmul(x,x,n)+c)%n;
		val=qmul(val,abs(x-y),n); //gcd(ab%n,n)>=gcd(a,n)
		if (!val) return n; //if fail, return before i reach k
		if (!(i&127) || i==k) p=__gcd(val,n);
		if (i==k) y=x,k+=k;
	}
	return p;
}
vector<int> primes;
void fastfactor(ll n){
	if (n==1) return;
	if (MillerRabin(n)) {primes.push_back(n); return;} //n is prime factor
	ll t=n; //if n not always lagre, make a single factor table for int may faster 
	while (t==n) t=pol_rho(n,rand()%(n-1));
	fastfactor(t); fastfactor(n/t);
}

}

struct Q{
	ll p,q;
	Q(ll x){p=x;q=1;}
	void operator=(ll x){p=x;q=1;}
	void simp(){ll t=gcd(p,q); if (t!=1) p/=t,q/=t; if (q<0) p=-p,q=-q;}
	void operator+=(const Q &v){p=p*v.q+v.p*q;q*=v.q;simp();}
	void operator-=(const Q &v){p=p*v.q*v.p*q;q*=v.q;simp();}
	void operator*=(const Q &v){p*=v.p;q*=v.q;simp();}
	void operator/=(const Q &v){p*=v.q;q*=v.p;simp();}
	Q operator+(const Q &y){Q x(*this);x+=y;return x;}
	Q operator-(const Q &y){Q x(*this);x-=y;return x;}
	Q operator*(const Q &y){Q x(*this);x*=y;return x;}
	Q operator/(const Q &y){Q x(*this);x/=y;return x;}
	bool operator<(const Q &v){return p*v.q<v.p*q;}
	double d(){return (double)p/q;}
};
//calc C(n,m), means select m object from n object
namespace LUCAS{
const ll p=10007;
ll fact[p],vfact[p];
ll comb(ll n,ll m){
	if (n<m) return 0;
	return fact[n]*vfact[n-m]%p*vfact[m]%p;
}
ll lucas(ll n,ll m){
	if (m==0) return 1;
	return lucas(n/p,m/p)*comb(n%p,m%p)%p;
}
void pre(){
	fact[0]=1;
	for (int i=1;i<p;i++) fact[i]=fact[i-1]*i%p;
	for (int i=0;i<p;i++) vfact[i]=qpow(fact[i], p-2, p);
}

//----exlucas----
//return C(n,m) mod any number
//time: O(max(p^k)*polylog)
ll calc(ll n, ll x, ll P){ //solve n! mod p^k
	if (!n) return 1;       //x:p ,P:p^k
	ll s=1;
	for (int i=1;i<=P;i++) //main cycle part
		if (i%x) s=s*i%P;
	s=qpow(s,n/P,P);
	for (ll i=1;i<=n%P;i++) //remain part
		if (i%x) s=s*i%P;
	return s*calc(n/x,x,P)%P; //mod p==0 part
}
ll multilucas(ll n, ll m, ll x, ll P) { //solve C(n,m) mod p^k
	ll cnt=0,s1=calc(n,x,P),s2=calc(m,x,P),s3=calc(n-m,x,P);
	for (ll i=n;i;i/=x) cnt+=i/x;
	for (ll i=m;i;i/=x) cnt-=i/x;
	for (ll i=n-m;i;i/=x) cnt-=i/x;
	return qpow(x,cnt,P)%P*s1%P*NumTheory::inv_euclid(s2,P)%P*NumTheory::inv_euclid(s3,P)%P;
}
ll exlucas(ll n, ll m, ll P){
	int cnt=0;
	ll p[20],a[20]; //no more 20 diff prime facter in int64
	for (ll i=2;i*i<=P;i++) //O(sqrt P), use pol_rho when p is large
		if (P%i==0){
			p[cnt]=1;
			while (P%i==0) p[cnt]=p[cnt]*i,P/=i;
			a[cnt]=multilucas(n,m,i,p[cnt]); //solve mod p^k
			cnt++;
		}
	if (P>1) a[cnt]=multilucas(n,m,P,P),p[cnt]=P,++cnt;
	return NumTheory::CRT(cnt,a,p);
}
}

namespace Cantor{
int fac[]={1,1,2,6,24,120,720,5040,40320,362880,3628800};
bool vis[20];
//get rank of permutation
//output range: [0,n!) 
int cantor(int n, int a[]){
	int ret=0;
	for (int i=0;i<n;i++){
		int t=0;
		for (int j=i+1;j<n;j++)
			if (a[i]>a[j]) t++;
		ret+=t*fac[n-i-1];
	}
	return ret;
}
//get kth permutation of {1..n}
//input range: [0,n!) 
void decantor(int n, int k, int ans[]){
	memset(vis,0,sizeof(vis));
	for (int i=0,j;i<n;i++){
		int t=k/fac[n-i-1];
		for (j=1;j<=n;j++)
			if (!vis[j])
				if (t==0) break;
				else t--;
		ans[i]=j;
		vis[j]=1;
		k%=fac[n-i-1];
	}
}
}

namespace Berlekamp_Massey{
const int maxn=2010;

//Berlekamp-Massey alg, O(n^2)
//a[i]=a[i-1]*p[1]+a[i-2]*p[2]+...
//x[]: input,  ps[pn][]: result
db x[maxn],delta[maxn];
vector<db> ps[maxn]; //if memory is limited, use rolling array
int n, cnt, fail[maxn], pn;
void BM(){
	int best=0;
	for(int i=1;i<=n;i++){
		db dt=-x[i];
		for (int j=0;j<ps[pn].size();j++)
			dt+=x[i-j-1]*ps[pn][j];
		delta[i]=dt;
		if (fabs(dt)<=1e-7) continue;
		fail[pn]=i; if (!pn) {ps[++pn].resize(i);continue;}
		vector<db> &ls=ps[best]; db k=-dt/delta[fail[best]]; //get coeff
		vector<db> cur; cur.resize(i-fail[best]-1); //tail 0
		cur.push_back(-k); 
		for(int j=0;j<ls.size();j++) cur.push_back(ls[j]*k); //get k*[last fail]
		if (cur.size()<ps[pn].size()) cur.resize(ps[pn].size());
		for(int j=0;j<ps[pn].size();j++) cur[j]+=ps[pn][j];  //add ps[pn]
        if(i-fail[best]+(int)ps[best].size()>=ps[pn].size()) best=pn;
        ps[++pn]=cur;
	}
}
}

namespace NumericalMethod{
double eps=1e-8;
double f(double x){
	return x;
}
//Three division method, require a convex function
double fargmax(double l, double r){
	while (r-l>eps){
		double d1=(l+l+r)/3,d2=(l+r+r)/3;
		if (f(d1)>f(d2)) r=d2;
		else l=d1;
	}
	return r;
}

//---simpson_intergrate---
double simpson(double l,double r){
	return (r-l)*(f(l)+4*f((l+r)/2)+f(r))/6;	
}
double asr(double l, double r, double ans){
	double mid=(l+r)/2;
	double left=simpson(l,mid),right=simpson(mid,r);
	if (fabs(left+right-ans)<eps) return left+right;
	else return asr(l,mid,left)+asr(mid,r,right);
}
//use simpson method
//warning: avoid odd point
double intergrate(double l, double r){return asr(l,r,simpson(l,r));}

//---SA----
const int maxn=10010;
struct Point{
	double x,y;
}a[maxn];
//Energy function demo
double w[maxn]; int n; 
double E(Point p){
	double ans=0;
	inc(i,n)
		ans+=sqrt((p.x-a[i].x)*(p.x-a[i].x)+(p.y-a[i].y)*(p.y-a[i].y)) * w[i];
	return ans;
}
/*
argmin_SimulateAnneal(): get argmin E(x), while E is not a convex function.
To speed up, the energy function should be nearly smooth and partial convex.
Point: n-d vector, the dimension cannot be too large.
t0:start temprature
tend: end temprature
delta: temprature reduce factor
*/
double rand01() {return (rand()+0.5)/RAND_MAX;}
Point ans; double anse=1e12;
void argmin_SimulateAnneal(double t0=5000, double tend=1e-6, double delta=0.99){
	double t=t0, ne=1e12; //current state energy
	Point p((Point){0,0});
	while (t>=tend){
		Point p1((Point){p.x+(rand01()*t-t/2), p.y+(rand01()*t-t/2)});
		double te=E(p1);
		//cout<<te-ne<<' '<<ne<<' '<<t<<' '<<exp((ne-te)/t)<<' '<<ans.x<<'\n';
		if (te<ne || 0 && exp((ne-te)/t)>rand01()) //disabled jumpout
			p=p1, ne=te; //update
		if (ne<anse)
			ans=p, anse=ne;//cout<<ans.x<<' '<<ans.y<<' '<<anse<<'\n';
		t*=delta;
	}
}
}

namespace Polynomial{
//Attention: n indicates number of coeffs, so deg F[x]=n-1, not n
namespace FFT
{
typedef complex<double> cd;
const int maxl=(1<<20)+1;
const double pi=3.14159265358979;
int rev[maxl];
void get_rev(int bit){
	for (int i=0;i<(1<<bit);i++)
		rev[i]=(rev[i>>1]>>1)|((i&1)<<(bit-1));
}
void fft(cd a[], int n, int dft){
	for(int i=0;i<n;i++) if(i<rev[i]) swap(a[i],a[rev[i]]);
	for (int s=1;s<n;s<<=1){
		cd wn=exp(cd(0,pi*dft/s));
		for (int j=0;j<n;j+=s<<1){
			cd wnk(1,0);
			for (int k=j;k<j+s;k++){
				cd x=a[k],y=wnk*a[k+s];
				a[k]=x+y;
				a[k+s]=x-y;
				wnk*=wn;
			}
		}
	}
	if (dft==-1) for (int i=0;i<n;i++) a[i]/=n;
}
ll G=3,P=998244353;
void ntt(ll *a, int n, int dft){
	for(int i=0;i<n;i++) if(i<rev[i]) swap(a[i],a[rev[i]]);
	for (int s=1;s<n;s<<=1){
		ll wn=qpow(G,dft==1?(P-1)/s/2:P-1-(P-1)/s/2,P);
		for (int j=0;j<n;j+=s<<1){
			ll wnk=1;
			for (int k=j;k<j+s;k++){
				ll x=a[k],y=wnk*a[k+s]%P;
				a[k]=(x+y)%P; //merge
				a[k+s]=(x-y+P)%P;
				wnk=wnk*wn%P;
			}
		}
	}
	if (dft==-1) {
		ll inv=qpow(n,P-2,P);
		for (int i=0;i<n;i++) a[i]=a[i]*inv%P;
	}
}
void conv(cd *fa, cd *fb, int s, cd *ret){
	static cd a[maxl],b[maxl];
	memcpy(a,fa,sizeof(cd)*s); memcpy(b,fb,sizeof(cd)*s);
	fft(a,s,1); fft(b,s,1);
	for (int i=0;i<s;i++) a[i]*=b[i];
	fft(a,s,-1);
	memcpy(ret,a,sizeof(cd)*s);
}
void conv(ll *fa, ll *fb, int s, ll *ret){
	static ll a[maxl],b[maxl];
	memcpy(a,fa,sizeof(ll)*s); memcpy(b,fb,sizeof(ll)*s);
	ntt(a,s,1); ntt(b,s,1);
	for (int i=0;i<s;i++) a[i]*=b[i];
	ntt(a,s,-1);
	memcpy(ret,a,sizeof(ll)*s);
}
int ans[maxl],_;
char s1[100010],s2[100010];
//fast mul
void base10_mul(){
	static cd a[maxl],b[maxl];
	scanf("%s%s",s1,s2);
	int l1=strlen(s1),l2=strlen(s2);
	int s=2,bit=1;
	for (bit=1;(1<<bit)<l1+l2-1;bit++)s<<=1;
	for (int i=0;i<l1;i++) a[i]=s1[l1-i-1]-'0';
	for (int i=0;i<l2;i++) b[i]=s2[l2-i-1]-'0';
	conv(a,b,s,a);
	for (int i=0;i<s;i++) cout<<a[i]<<' '; cout<<'\n';
	for (int i=0;i<s;i++){
		ans[i]+=a[i].real();
		ans[i+1]+=ans[i]/10;
		ans[i]%=10;
	}
	int i;
	for (i=l1+l2;!ans[i]&&i>=0;i--);
	if (i==-1) printf("0");
	for (;i>=0;i--) printf("%d",ans[i]);
	putchar('\n');
}

//C[x]=A[x]*B[x], convolve only
void poly_mul(ll *A, ll *B, int n, ll *C){
	int s=2,bit=1;
    for (bit=1;(1<<bit)<n<<1;bit++)s<<=1; //size n*2
    fill(A+n,A+s,0); fill(B+n,B+s,0);
    get_rev(bit);
    conv(A,B,s,C);
}
void poly_add(ll *A, ll *B, int n, ll *C){
	for (int i=0;i<n;i++) 
		C[i]=(A[i]+B[i])%P;
}
void poly_sub(ll *A, ll *B, int n, ll *C){
	for (int i=0;i<n;i++) 
		C[i]=(A[i]-B[i]+P)%P;
}

#ifdef NO_COMPILE  //3 module ntt version
const ll p=1000000007,P1=1004535809,P2=998244353,P3=469762049;
void poly_mul(ll *A, ll *B, int n, ll *C){
	static ll res[3][maxl];
	int s=2,bit=1;
    for (bit=1;(1<<bit)<n<<1;bit++)s<<=1; //size n*2
    get_rev(bit);
    conv<P1>(A,B,n,res[0]);
    conv<P2>(A,B,n,res[1]);
    conv<P3>(A,B,n,res[2]);
    ll M=P1*P2;
    for (int i=0;i<s;i++){ //merge
    	ll A=(qmul(res[0][i]*P2%M,inv(P2%P1,P1),M)+
			qmul(res[1][i]*P1%M,inv(P1%P2,P2),M))%M;
		ll K=(res[2][i]-A%P3+P3)*inv(M%P3,P3)%P3;
		C[i]=(qmul(K%p,M%p,p)+A)%p;
    }
}
//input A[x], A[x]*B[x]=1 (mod x^n)
void poly_inv(ll *A, ll *B, int n){
    if (n==1) {B[0]=qpow(A[0],p-2,p); return;} //constant element
    poly_inv(A,B,n+1>>1); //divide
    int s=2,bit=1;
    for (bit=1;s<n<<1;bit++)s<<=1; //size n*2
    static ll ta[maxl];//A(x)((2-B(x)A(x))B(x))=1(mod x^2)
    poly_mul(A,B,n,ta);
    for (int i=1;i<s;i++) ta[i]=p-ta[i];
    ta[0]=(2+p-ta[0])%p;
    poly_mul(ta,B,n,B);
    //for (int i=0;i<s;i++)cout<<B[i]<<' '; cout<<'\n';
}
#endif

//input A[x], A[x]*B[x]=1 (mod x^n)
void poly_inv(ll *A, ll *B, int n){
	if (n==1) {B[0]=qpow(A[0],P-2,P); return;} //constant element
	poly_inv(A,B,n+1>>1); //divide
    int s=2,bit=1;
    for (bit=1;(1<<bit)<n<<1;bit++)s<<=1; //size n*2
    get_rev(bit);
    static ll ta[maxl];
    copy(A,A+n,ta); fill(ta+n,ta+s,0);
	ntt(ta,s,1); ntt(B,s,1);
    for (int i=0;i<s;i++) //A(x)(2B(x)-B(x)^2A(x))=1(mod x^2)
    	B[i]=(2-ta[i]*B[i]%P+P)%P*B[i]%P;  
    ntt(B,s,-1); fill(B+n,B+s,0);
    //for (int i=0;i<s;i++)cout<<B[i]<<' '; cout<<'\n';
}

//A[x]=Q[x]*B[x]+R[x], (deg R[x]<deg B[x], deg Q[x]<n-m+1)
void poly_div(ll *A, ll *B, int n, int m, ll *Q, ll *R){
	//conclusion: Qr[x] = Ar[x] * poly_inv(Br[x]) (mod x^(n-m+1))
	static ll X[maxl], Y[maxl], tmp[maxl];
	memset(X,0,sizeof X); memset(Y,0,sizeof Y); memset(tmp,0,sizeof tmp);
	copy(A,A+n,X); copy(B,B+m,Y);
	reverse(X,X+n); reverse(Y,Y+m);
	poly_inv(Y,tmp,n-m+1); 
	poly_mul(X,tmp,n-m+1,Q);
	fill(Q+n-m+1,Q+n,0), reverse(Q,Q+n-m+1);
    copy(A,A+n,X), copy(B,B+m,Y);
    poly_mul(Y,Q,n,tmp), poly_sub(X,tmp,n,R);
	fill(R+m-1,R+n,0);
}

//B[x]=intergrate A[x]dx,  deg B[x] = deg A[x] + 1
void poly_int(ll *A, ll *B, ll n){
	B[0]=0; //constant C
	for (int i=1;i<=n;i++)
		B[i]=A[i-1]*qpow(i,P-2,P)%P;
}
//B[x]=A'[x],  deg B[x] = deg A[x] - 1
void poly_deriv(ll *A, ll *B, ll n){
	B[n-1]=0;
	for (int i=0;i<n-1;i++)
		B[i]=A[i+1]*(i+1)%P;
}

//B[x]=ln A[x]=int(B'[x])=int(A'[x]/A[x])  (mod x^n)
void poly_ln(ll *A, ll *B, ll n){
	static ll X[maxl], Y[maxl], tmp[maxl];
	memset(X,0,sizeof X); memset(Y,0,sizeof Y); memset(tmp,0,sizeof tmp);
	copy(A,A+n,X); poly_deriv(X,Y,n); 
	poly_inv(X,tmp,n); copy(tmp,tmp+n,X); 
	poly_mul(X,Y,n,tmp); poly_int(tmp,B,n);
}

//Fnew[x]=Fb[x](1-ln(Fb[x])+A[x]), O(nlogn)
void poly_exp(ll *A, ll *B, int n) {
    if (n==1) {B[0]=1;return;}
    static ll tmp[maxl]; memset(tmp,0,sizeof tmp);
    poly_exp(A,B,n+1>>1);fill(B+(n+1>>1),B+n,0);
    poly_ln(B,tmp,n); poly_sub(A,tmp,n,tmp);
    poly_mul(tmp,B,n,tmp); poly_add(B,tmp,n,B);
}

//B[x]=sqrt(A[x])=exp(0.5*ln(A[x]))
void poly_sqrt(ll *A, ll *B, int n) {
    static ll tmp[maxl]; memset(tmp,0,sizeof tmp);
	poly_ln(A,tmp,n);
	for (int i=0;i<n;i++) tmp[i]=tmp[i]*qpow(2,P-2,P)%P;
    poly_exp(tmp,B,n);
}

//deleted, please use poly_mul(3 module version) instead
const ll P1=1004535809,P2=998244353,P3=469762049;
ll res[3][maxl];
//conv a sequence with mod p, while p<P1*P2*P3
void conv_mod(){
	static ll a[maxl],b[maxl];
	int l1,l2; ll p;
    rd(l1); rd(l2); l1++; l2++;
    int s=2,bit=1;
    for (bit=1;(1<<bit)<l1+l2-1;bit++)s<<=1;
	get_rev(bit);
    int r; rd(r); p=r;
    for (int i=0;i<l1;i++) rd(r),a[i]=r;
    for (int i=0;i<l2;i++) rd(r),b[i]=r;
    G=3,P=P1; conv(a,b,s,res[0]);
    G=3,P=P2; conv(a,b,s,res[1]);
    G=3,P=P3; conv(a,b,s,res[2]);
    ll M=P1*P2;
    for (int i=0;i<l1+l2-1;i++){ //merge
    	//printf("%lld %lld %lld \n",res[0][i],res[1][i],res[2][i]);
    	ll A=(qmul(res[0][i]*P2%M,inv(P2%P1,P1),M)+
			qmul(res[1][i]*P1%M,inv(P1%P2,P2),M))%M;
		ll K=(res[2][i]-A%P3+P3)*inv(M%P3,P3)%P3;
		//cout<<A<<' '<<K<<' ';
		printf("%lld ", (qmul(K%p,M%p,p)+A)%p);
    }
}
} //namespace FFT

const int maxn=2010;
ll x[maxn],y[maxn];
int n;
//O(n^2) get single point value
//if xi between 1~n, we can optimize it to O(n)
//if xi not in 1~n, we can still preprocess PIj (ki-x[j]+p)%p, 
//  then get multi point value in O(max(n^2,nm))
ll LangrangeInter(ll k, ll p){
    ll sum=0;
    for (int i=0;i<n;i++){
        ll s0=1;
        for (int j=0;j<n;j++)
            if (i!=j) s0=s0*(x[i]-x[j]+p)%p;
        s0=qpow(s0,p-2,p);
        for (int j=0;j<n;j++)
            if (i!=j) s0=s0*(k-x[j]+p)%p;
        sum=(sum+y[i]*s0)%p;
    }
	return sum;
}

namespace FWT{
int N,P,inv2;
const int maxl=1<<18+1;
void fwt_or(int *a,int opt){
    for(int i=1;i<N;i<<=1)
        for(int p=i<<1,j=0;j<N;j+=p)
            for(int k=0;k<i;++k)
                if(opt==1) a[i+j+k]=(a[j+k]+a[i+j+k])%P;
                else a[i+j+k]=(a[i+j+k]+P-a[j+k])%P;
}
void fwt_and(int *a,int opt){
    for(int i=1;i<N;i<<=1)
        for(int p=i<<1,j=0;j<N;j+=p)
            for(int k=0;k<i;++k)
                if(opt==1) a[j+k]=(a[j+k]+a[i+j+k])%P;
                else a[j+k]=(a[j+k]+P-a[i+j+k])%P;
}
void fwt_xor(int *a,int opt){
    for(int i=1;i<N;i<<=1)
        for(int p=i<<1,j=0;j<N;j+=p)
            for(int k=0;k<i;++k){
                int X=a[j+k],Y=a[i+j+k];
                a[j+k]=(X+Y)%P,a[i+j+k]=(X+P-Y)%P;
                if(opt==-1) a[j+k]=(ll)a[j+k]*inv2%P,a[i+j+k]=(ll)a[i+j+k]*inv2%P;
            }
}
void bit_conv(int *fa, int *fb, int *c){
	static int a[maxl],b[maxl];
	memcpy(a,fa,sizeof(int)*N); memcpy(b,fb,sizeof(int)*N);
	fwt_xor(a,1); fwt_xor(b,1);
	for (int i=0;i<N;i++) a[i]=a[i]*b[i]%P;
	fwt_xor(a,-1);
	memcpy(c,a,sizeof(int)*N);
}
}
}

namespace LinAlg{
const int maxn=1010,maxm=1010;
double a[maxn][maxm],b[maxn][maxn],ans[maxn];
int n,m;
const int eps=1e-7;
//require m=n+1 and only one solution
bool gauss_solve(){
	for (int i=0;i<n;i++){
		int maxl=i;
		for (int j=i+1;j<n;j++)
			if (fabs(a[j][i])>fabs(a[maxl][i])) maxl=j;
		if (maxl!=i) swap(a[i],a[maxl]);
		if (fabs(a[i][i])<eps) return 0; //no solution or infinity solution 
		for (int j=i+1;j<n;j++){
			double r=a[j][i]/a[i][i];
			for (int k=i;k<m;k++)
				a[j][k]-=r*a[i][k];
		}
		double r=a[i][i];
		for (int k=i;k<m;k++) a[i][k]/=r;
	}
	for (int i=n-1;i>=0;i--){
		ans[i]=a[i][n];
		for (int j=i+1;j<n;j++)
			ans[i]-=a[i][j]*ans[j];
	}
	return 1;
}
//n*n matrix
bool matinv(){
	memset(b,0,sizeof(b));
	for (int i=0;i<n;i++) b[i][i]=1;	
	for (int i=0;i<n;i++){
		int maxl=i;
		for (int j=i+1;j<n;j++)
			if (fabs(a[j][i])>fabs(a[maxl][i])) maxl=j;
		if (i!=maxl) swap(a[i],a[maxl]),swap(b[i],b[maxl]);
		if (fabs(a[i][i])<eps) return 0; //No trivil
		double r=a[i][i];
		for (int k=0;k<n;k++) a[i][k]/=r,b[i][k]/=r; //k start from 0
		for (int j=i+1;j<n;j++){
			double r=a[j][i]/a[i][i];
			for (int k=0;k<n;k++) //k start from 0
				a[j][k]-=r*a[i][k], b[j][k]-=r*b[i][k];
		}
	}
	return 1;
}
double det(){
	double ans=1;
	for (int i=0;i<n;i++){
		int maxl=i;
		for (int j=i+1;j<n;j++)
			if (fabs(a[j][i])>fabs(a[maxl][i])) maxl=j;
		if (i!=maxl) swap(a[i],a[maxl]),ans=-ans;
		if (fabs(a[i][i])<eps) return 0;
		for (int j=i+1;j<n;j++){
			double r=a[j][i]/a[i][i];
			for (int k=i;k<n;k++)
				a[j][k]-=r*a[i][k];
		}
	}
	for (int i=0;i<n;i++) ans*=a[i][i];
	return ans;
}
int matrank(){
	int l=0; //real line
	for (int i=0;i<m;i++){ //i: col start pos
		int maxl=l;
		for (int j=l+1;j<n;j++)
			if (fabs(a[j][i])>fabs(a[maxl][i])) maxl=j;
		if (l!=maxl) swap(a[l],a[maxl]);
		if (fabs(a[l][i])<eps) continue;
		for (int j=l+1;j<n;j++){
			double r=a[j][i]/a[l][i];
			for (int k=i;k<m;k++)
				a[j][k]-=r*a[i][k];
		}
		l++;
	}
	return l;
}
const ll p=19260817; //const is faster than normal variable
//int det with mod
//used by Matrix-Tree theorem
//M-T theo: a[i][i]=deg i, a[i][j]=-cnt(i->j), n=|vertex|-1
ll detint(ll **a){
	ll ans=1;
	for (int i=0;i<n;i++) for (int j=0;j<n;j++) if (a[i][j]<0) a[i][j]+=p;
	for (int i=0;i<n;i++){
		int maxl=i;
		for (int j=i;j<n;j++)
			if (a[j][i]) {maxl=j;break;}
		if (i!=maxl) swap(a[i],a[maxl]), ans=p-ans;
		if (a[i][i]==0) return 0;
		ans=ans*a[i][i]%p;   //[i] Here update ans before set a[i][i] to 1
		int v=inv(a[i][i],p);
		for (int j=i;j<m;j++) a[i][j]=a[i][j]*v%p;
		for (int j=i+1;j<n;j++){
			ll r1=a[j][i];
			if (!r1) continue;
			for (int k=i;k<n;k++)
				a[j][k]=(a[j][k]-r1*a[i][k]%p+p)%p;
		}
	}
	return ans;
}
//matinv with mod
//require m=2*n, result is a[i][(0..n-1)+n]
bool matinv_int(ll **a){
	m=2*n;
	//a[i][(0..n-1)+n]=0; //set to 0 when necessary
	for (int i=0;i<n;i++) a[i][i+n]=1;	
	for (int i=0;i<n;i++){
		int maxl=i;
		for (int j=i;j<n;j++)
			if (a[j][i]) {maxl=j;break;}
		if (i!=maxl) swap(a[i],a[maxl]);
		if (a[i][i]==0) return 0;
		int v=inv(a[i][i],p);
		for (int j=i;j<m;j++) a[i][j]=a[i][j]*v%p;
		for (int j=0;j<n;j++){
			if (!a[j][i] || j==i) continue;
			ll r1=a[j][i];
			for (int k=i;k<m;k++)
				a[j][k]=(a[j][k]-r1*a[i][k]%p+p)%p;
		}
	}
	return 1;
}
#ifdef NO_COMPILE
int n,m,p;
//solve linear equ in module p
//not require m=n+1, get a possible solution ans[0..m-1)
bool gauss_solve_int(){  //similar as matrank
	int l=0; //real line
	for (int i=0;i<m;i++){ //i: col start pos
		int maxl=l;
		for (int j=l;j<n;j++)
			if (a[j][i]) {maxl=j;break;}
		if (maxl!=l) swap(a[l],a[maxl]);
		if (!a[l][i]) continue; //next col
		int v=inv(a[l][i],p);
		for (int j=i;j<m;j++) a[l][j]=a[l][j]*v%p;
		for (int j=l+1;j<n;j++){
			if (!a[j][i]) continue;
			int r1=a[j][i];
			for (int k=i;k<m;k++)
				a[j][k]=(a[j][k]-r1*a[l][k]%p+p)%p;
		}
		l++;
	} //now l is rank of matrix
	int last=m-1,cur;
	for (int i=l-1;i>=0;i--){
		for (cur=0;cur<m-1 && !a[i][cur];cur++); //first no zero column
		int t=a[i][m-1]; 
		for (int j=last;j<m-1;j++) t=(t-a[i][j]*ans[j]%p+p)%p;
		for (last--;last>cur;last--) //any solution, set to 0
			ans[last]=0; //,t=(t-a[i][j]*ans[last]%p+p)%p; when ans[last] is not 0
		if (cur==m-1 && t) return 0; //no solution
		ans[cur]=t; //a[i][cur]=1, so ans[cur]=t
		last=cur;
	}
	return 1;
}
#endif
}

//-----------------String Algo----------------------

namespace StringHash{
//double module hash
typedef unsigned long long ll;
const ll m1=1000000007;
const int maxn=1000010;
ll h1[maxn],h2[maxn],b1[maxn],b2[maxn];
void pre(){
	b1[0]=b2[0]=1;
	inc(i,maxn-1)
		b1[i+1]=b1[i]*131%m1,
		b2[i+1]=b2[i]*137;
}
void gethash(char *s, int l){
	h1[l]=h2[l]=0;
	dec(i,l){
		h1[i]=(h1[i+1]*131+s[i])%m1;
		h2[i]=h2[i+1]*137+s[i];
		//cout<<h1[i]<<' '<<b1[i]<<'\n';
	}
}
//get substring [l,r) hash value
pair<ll,ll> substr(int l, int r){
	ll r1=h1[l]+m1-h1[r]*b1[r-l]%m1; if (r1>=m1) r1-=m1;
	ll r2=h2[l]-h2[r]*b2[r-l];
	return {r1,r2};
}
}

namespace KMP{
const int maxn=1000010;

int kmp[maxn];
//s is a short partten string
char s[maxn]; int sl; 
void getkmp(){
	int j=0,k=-1;
	kmp[0]=-1;
	while (j<sl)
		if (k==-1 || s[j]==s[k])
			kmp[++j]=++k;
		else
			k=kmp[k];
}
int kmp_idx(char *t, int tl){
	int i=0, j=0;
	while (i<sl && j<tl)
		if (i==-1 || s[i]==t[j])
			i++,j++;
		else
			i=kmp[i];
	if (i==sl) return j-sl;
	else return -1;
}
int kmp_cnt(char *t, int tl){
	int i=0, j=0, cnt=0;
	while (j<tl){
		if (i==-1 || s[i]==t[j])
			i++,j++;
		else
			i=kmp[i];
		if (i==sl) cnt++;
	}
	return cnt;
}
}
//!-- untested
//extend KMP
//extend[i]: LCP lenth between s1[i..l1] and s2
namespace E_KMP{
const int N=100010;
int next[N],extend[N];
void getnext(char *str){
    int i=0,j,po,len=strlen(str);
    next[0]=len;
    while(str[i]==str[i+1] && i+1<len) i++; next[1]=i; //calc next[1]
    po=1;
    for(i=2;i<len;i++)
        if(next[i-po]+i < next[po]+po)
            next[i]=next[i-po];
        else{
            j = next[po]+po-i;
            if(j<0) j=0;
            while(i+j<len && str[j]==str[j+i]) j++; next[i]=j;
            po=i;
        }
}
void exkmp(char *s1,char *s2){
    int i=0,j,po,len=strlen(s1),l2=strlen(s2);
    getnext(s2);
    while(s1[i]==s2[i] && i<l2 && i<len) i++; extend[0]=i;
    po=0;
    for(i=1;i<len;i++)
        if(next[i-po]+i < extend[po]+po)
            extend[i]=next[i-po];
        else //continue try match after e[po]+po
        {
            j = extend[po]+po-i;
            if(j<0) j=0;
            while(i+j<len && j<l2 && s1[j+i]==s2[j]) j++; extend[i]=j;
            po=i; //update po
        }
}
}

namespace ACAM{
const int maxn=100000,alpha=26; //maxn >= sigma len(si)
int ch[maxn][alpha],val[maxn],fail[maxn],lbl[maxn],len[maxn],pc=0;
int cnt[1000]; //str appear times, first element is 1
int strc=0;
void clear(){
	pc=0; strc=0;
	memset(ch,0,sizeof(ch));
	memset(fail,0,sizeof(fail));
	memset(val,0,sizeof(val));
	memset(lbl,0,sizeof(lbl));
}
//Trie construct
void ins(char *s){
	int l=strlen(s), cur=0;
	for (int i=0;i<l;i++){
		int v=s[i]-'a';
		if (!ch[cur][v]) ch[cur][v]=++pc;
		cur=ch[cur][v];
		len[cur]=i+1;
	}
	strc++;
	lbl[cur]=strc;
	val[cur]++;
}
int qu[maxn];
//fail edge add
void build(){
	int qh=0,qt=0;
	for (int i=0;i<alpha;i++)
		if (ch[0][i]) fail[ch[0][i]]=0,qu[qt++]=ch[0][i];
	while (qh<qt){
		int u=qu[qh];qh++;
		for (int i=0;i<alpha;i++)
			if (ch[u][i])
				fail[ch[u][i]]=ch[fail[u]][i],qu[qt++]=ch[u][i];
			else
				ch[u][i]=ch[fail[u]][i]; //opt, move to fail auto. Attention: the multi fail jump will be emitted
	}
}
//count how many mode str appeared in s as substr
int appear(char *s){
	int l=strlen(s),cur=0,ans=0;
	for (int i=0;i<l;i++){
		cur=ch[cur][s[i]-'a']; //the opt trans emitted fail jump chain, do it when necessary
		for (int t=cur;t && ~val[t];t=fail[t]) //the label be sure O(n)
			ans+=val[t],val[t]=-1;
	}
	return ans;
}
//count each mode str in s
//[!] worst O(n^2), a better way is dp on fail tree
void cntall(char *s){
	int l=strlen(s),cur=0;
	memset(cnt,0,sizeof(cnt));
	for (int i=0;i<l;i++){
		cur=ch[cur][s[i]-'a']; //the opt trans emitted fail jump chain, do it when necessary
		for (int t=cur;t;t=fail[t])
			cnt[lbl[t]]++;
	}
}
}

namespace SA{
//sa: pos of ith rk suf, rk: rk of i pos suf, a: s0-'a'
//t1,t2,c: temp array, h: height(LCP of sa[i] and sa[i-1])
int t1[N],t2[N],sa[N],h[N],rk[N],c[N],a[N];
int n,m;
void calcsa(){
	int *x=t1,*y=t2,p=0,f=0;
	icc(i,m) c[i]=0;        
	icc(i,n) c[x[i]=a[i]]++;  //first char sort
	icc(i,m) c[i]+=c[i-1];
	dcc(i,n) sa[c[x[i]]--]=i; //pos of ith first char
	for (int i=1;i<=n && p<=n;i<<=1){
		p=0;
		rpp(j,n-i+1,n) y[++p]=j;  //remain part
		icc(j,n) if (sa[j]>i) y[++p]=sa[j]-i; //main part
		icc(j,m) c[j]=0;
		icc(j,n) c[x[y[j]]]++;
		icc(i,m) c[i]+=c[i-1];
		dcc(j,n) sa[c[x[y[j]]]--]=y[j]; //sort, use dcc because sort should be stable
		swap(x,y);x[sa[1]]=p=1; 
		rpp(j,2,n) x[sa[j]]=y[sa[j]]==y[sa[j-1]]&&y[sa[j]+i]==y[sa[j-1]+i]?p:++p; //refill key
		m=p;
	}
	icc(i,n) rk[sa[i]]=i;
	icc(i,n){ //get height, O(n)
		int j=sa[rk[i]-1];
		if (f) f--; while (a[i+f]==a[j+f]) f++;
		h[rk[i]]=f;
	}
}

char s0[N];
int main(){
	scanf("%s",s0); int l=strlen(s0);
	inc(i,l) a[++n]=s0[i]-' ';
	m=200;
	calcsa();
	icc(i,n) printf("%d ",sa[i]);
	return 0;
}
}

namespace SAM{
const int maxn=100010,alpha=26;

//max state cnt: 2*strlen-1
//max transfer cnt: 3*strlen-4
struct Node{
	int l,num; bool vis;
	Node *p, *tr[alpha];
	//vector<Node *> ch;
	void set(int _l){l=_l;memset(tr,0,sizeof(tr));p=0;num=1;vis=0;/*ch.clear();*/}
}nodes[maxn<<1];
int nodec;
Node *root;
Node *open(int l){
	nodes[nodec++].set(l);
	return nodes+nodec-1;
}
void build(char *s, int l){
	Node *cur;
	cur=root=open(0);
	for (int i=0;i<l;i++){
		int x=s[i]-'a';
		Node *p=cur;
		cur=open(i+1);
		for (;p && !p->tr[x];p=p->p)
			p->tr[x]=cur;
		if (!p) cur->p=root;
		else{
			Node *q=p->tr[x];
			if (p->l+1==q->l) cur->p=q;
			else{
				Node *r=open(-1); r[0]=q[0]; r->l=p->l+1;
				q->p=r; cur->p=r; r->num=0;
				for (;p && p->tr[x]==q;p=p->p) p->tr[x]=r;
			}
		}
	}
}
//get substr last position
int pos(Node *u, char *s){
	if (*s==0) return u->l;
	if (!u->tr[*s-'a']) return -1;
	return pos(u->tr[*s-'a'],s+1);
}

int t[maxn],r[maxn]; //t:temp, r:rank(ith element pos)
//init |right(s)| before cnt
void initnum(int s0l){
	rep(i,0,s0l+1) t[i]=0;
	inc(i,nodec) t[nodes[i].l]++;
	rep(i,1,s0l+1) t[i]+=t[i-1];
	inc(i,nodec) r[--t[nodes[i].l]]=i; //sort by count
	per(i,nodec,1) nodes[r[i]].p->num+=nodes[r[i]].num; //dp
}
//count substr
int cnt(Node *u, char *s){
	if (*s==0) return u->num;
	if (!u->tr[*s-'a']) return 0;
	return cnt(u->tr[*s-'a'],s+1);
}
// longest substring
int lcs(char *x1){
	int lcs=0, ans=0, xl=strlen(x1);
	Node *p=root;
	inc(i,xl){
		int x=x1[i]-'a';
		if (p->tr[x]){
			lcs++;
			p=p->tr[x];
			ans=max(ans,lcs);
			continue;
		}
		for (;p && !p->tr[x];p=p->p);
		if (!p) p=root,lcs=0;
		else{
			lcs=p->l+1;
			p=p->tr[x];
		}
		ans=max(ans,lcs);
	}
	return ans;
}
}

//mulit-str SAM
namespace GSAM{
const int maxn=200010;
struct Node{
	int l, num, las, vis;
	Node *p;
	map<int,Node*> tr; //more time, less memory
	//int tr[26];
}nodes[maxn<<1];
int nodec;
Node *open(int l){
	nodes[nodec++].l=l;
	return nodes+nodec-1;
}
Node *root=open(0);
void build(int *s, int l){
	Node *cur=root;
	for (int i=0;i<l;i++){
		int x=s[i];
		if (cur->tr.count(x)){  //transfer existed
			cur=cur->tr[x];
			continue;
		}
		Node *p=cur;
		cur=open(i+1);
		for (;p && !p->tr.count(x);p=p->p)
			p->tr[x]=cur;
		if (!p) cur->p=root;
		else{
			Node *q=p->tr[x];
			if (p->l+1==q->l) cur->p=q;
			else{
				Node *r=open(-1); r[0]=q[0]; r->l=p->l+1;
				q->p=r; cur->p=r; r->num=0;
				for (;p && p->tr[x]==q;p=p->p) p->tr[x]=r;
			}
		}
	}
}
int len[200010],tot;
vector<int> str[200010];
int ts[200010];
int ans[200010];
//calc how many mode str appear in each query
void upd1(Node *u, int col){ 
	for (;u->las!=col && u!=root;u=u->p)
		u->las=col, u->num++;
}
//calc ecah mode str appear how many query
void upd2(Node *u, int col){ 
	for (;u->las!=col && u!=root;u=u->p)
		u->las=col,ans[col]+=u->vis;
}
//n: ptn str, m: query str
int main(){
	int n,m; cin>>n>>m;
	for (int i=0;i<n;i++){
		++tot;
		scanf("%d", len+tot);
		for (int j=0;j<len[tot];j++){
			scanf("%d",ts+j);
			str[tot].push_back(ts[j]);
		}
		build(ts,len[tot]);
	}
	for (int i=1;i<=tot;i++){
		Node *cur=root;
		for (int j=0;j<len[i];j++)
			cur=cur->tr[str[i][j]],
			upd1(cur,i);
	}
	for (int i=0;i<m;i++){
		int l,x; bool flag=1;
		scanf("%d",&l);
		Node *cur=root;
		for (int j=0;j<l;j++){
			scanf("%d",&x);
			if (flag)
				if (cur->tr.count(x))
					cur=cur->tr[x];
				else //no transfer
					flag=0;
		}
		if (flag)
			printf("%d\n",cur->num),
			cur->vis++;
		else
			printf("0\n");
	}
	for (int i=0;i<nodec;i++) nodes[i].las=0; //!-- clear
	for (int i=1;i<=tot;i++){
		Node *cur=root;
		for (int j=0;j<len[i];j++)
			cur=cur->tr[str[i][j]],
			upd2(cur,i);
	}
	for (int i=1;i<=tot;i++)
		printf("%d ",ans[i]);
	return 0;
}
}

namespace Manacher{
const int maxn=10000000;
//p[]: max palindrome len at pos i
int p[maxn<<1];char s[maxn<<1],s0[maxn];
int sl,s0l;
int manacher(){
    s0l=strlen(s0);
    sl=1; s[0]='$'; s[1]='#';
    inc(i,s0l) s[++sl]=s0[i],s[++sl]='#';
    s[++sl]=0;
    int mx=0,mi=0,ans=0; //mx: max cur pstr right pos, mi: max cur pstr center pos
    rep(i,1,sl){
        p[i]=i<mx?min(p[mi*2-i], mx-i):1;
        while (s[i-p[i]]==s[i+p[i]]) p[i]++;
        if (mx<i+p[i]) mi=i,mx=i+p[i];
        ans=max(ans,p[i]-1);
    }
    return ans;
}
}

namespace PAM{
const int maxn=2000500;

//ch[x]: if cur pstr is a, the ch is xax.  len: len of cur pstr
//fail: longest pstr suffix of cur point.  cnt: count of this pstr.
struct Node{
    int ch[10],fail,len;
	int cnt;
}node[maxn];
int nodec,cur, len[maxn];
char s[maxn];

void pre(){
    node[0].fail=1; node[1].len=-1;
    nodec=2;cur=0;
}
void insert(int p){
    int j, x=s[p]-'0';
    while(s[p-node[cur].len-1]!=s[p])cur=node[cur].fail; //find ch
    if(!node[cur].ch[x]){
        node[nodec].len=node[cur].len+2;
        j=node[cur].fail;
        while(s[p-node[j].len-1]!=s[p])j=node[j].fail; //find fail
        node[nodec].fail=node[j].ch[x];
        node[cur].ch[x]=nodec;
        cur=nodec;
		nodec++;
    }
    else cur=node[cur].ch[x];
    len[p]=node[cur].len;
	node[cur].cnt++;
}
char ts[maxn];
void dfs1(int u, int deep){
	cout<<ts<<' '<<node[u].len<<'\n'; //cur node
	for (int i=0;i<10;i++)
		if (node[u].ch[i]){
			ts[deep]=i+'0';
			dfs1(node[u].ch[i],deep+1);
		}
}
int main(){
	pre();
	scanf("%s",s); int l=strlen(s);
	for (int i=0;i<l;i++) insert(i);
	for (int i=nodec-1;i>0;i--) 
		node[node[i].fail].cnt+=node[i].cnt;
	dfs1(0,0); //even pstr
	dfs1(1,0); //odd pstr
	return 0;
}
}

/*-----------------------data structure------------------------*/

namespace UFSet{
const int maxn=100010;
int fa[maxn];
void clear(){
	for (int i=0;i<maxn;i++) fa[i]=i;
}
int fi(int x){
	if (fa[x]!=x)
		fa[x]=fi(fa[x]);
	return fa[x];
}
void un(int a, int b){
	int ta=fi(a),tb=fi(b);
	if (ta!=tb) fa[ta]=tb;
}
}

namespace Graph{

const int maxn=10010,maxm=100010,inf=0x3f3f3f3f;
int head[maxn],nxt[maxm],to[maxm],co[maxm],ec;
int n;
bool vis[maxn];
int dis[maxn],c[maxn];

void added(int x, int y, int c){
	ec++;
	nxt[ec]=head[x];
	head[x]=ec;
	to[ec]=y;
	co[ec]=c;
}

//wrost O(n^3), but as fast as dijkstra in random data.
int spfa(int s){
	queue<int> q;
	memset(dis,0x3f,sizeof(dis));
	dis[s]=0;
	memset(c,0,sizeof(c)); //judge nagetive loop
	memset(vis,0,sizeof(vis));
	q.push(s); vis[s]=1; c[s]=1;
	while (!q.empty())	{
		int u=q.front(); q.pop();
		vis[u]=0;
		for (int e=head[u];e;e=nxt[e]){
			int v=to[e];
			if (dis[u]+co[e]<dis[v]){
				dis[v]=dis[u]+co[e];
				if (!vis[v]){
					vis[v]=1;
					c[v]++;
					q.push(v);
					if (c[v]>n) return 0; //has nagetive
				}
			}
		}
	}
	return 1;
}
/*
judge negative circle
!-- Discarded, some reasons show it's time complexity is unreliable,
thought it runs fast in random data.
*/
bool spfa_dfsjudge(int u){
	vis[u]=1;
	for (int e=head[u];e;e=nxt[e]){
		int v=to[e];
		if (dis[u]+co[e]<dis[v]){
			dis[v]=dis[u]+co[e];
			if (vis[v] || spfa_dfsjudge(v)) return true;
		}
	}
	vis[u]=0;
	return false;
}

void dijk(int s){
	memset(dis,0x3f,sizeof(dis));
	memset(vis,0,sizeof(vis));
	dis[s]=0;
	priority_queue<pair<int,int> > qu;
	qu.push(make_pair(0,s));
	while (qu.size()){
		int u=qu.top().second, mc=-qu.top().first;
		qu.pop();
		if (vis[u]) continue;
		vis[u]=1;
		for (int e=head[u];e;e=nxt[e])
			if (!vis[to[e]] && mc+co[e]<dis[to[e]]){
				dis[to[e]]=mc+co[e];
				qu.push(make_pair(-dis[to[e]],to[e]));
			}
	}
}
int mp[maxn][maxn];
void dijk_original(int s){
	memset(dis,0x3f,sizeof(dis));
	memset(vis,0,sizeof(vis));
	dis[s]=0;
	for(int i=0;i<n;i++){
		int u=0,md=0x3f3f3f3f;
		for (int j=0;j<n;j++)
			if (!vis[j] && dis[j]<md)
				md=dis[j], u=j;
		vis[u]=1;
		for (int j=0;j<n;j++)
			if (!vis[j])
				dis[j]=min(dis[j],md+mp[u][j]);
	}
}
int d[maxn][maxn];
void floyd(){
	for (int k=0;k<n;k++)
		for (int i=0;i<n;i++)
			for (int j=0;j<n;j++)
				if (d[i][j]>d[i][k]+d[k][j])
					d[i][j]=d[i][k]+d[k][j];
}

//For undirected graph min loop
//!-- In directed graph, use floyd() and d[i][i].
int floyd_minc(){
	int minc=inf;
	for (int k=0;k<n;k++){
		for (int i=0;i<k;i++)
			for (int j=i+1;j<k;j++)
				minc=min(minc,d[i][j]+mp[i][k]+mp[k][j]);
		for (int i=0;i<n;i++)
			for (int j=0;j<n;j++)
				if (d[i][j]>d[i][k]+d[k][j])
					d[i][j]=d[i][k]+d[k][j];
	}
	return minc;
}

vector<int> ed[maxn];
bool ins[maxn];
int st[maxn],stn,dfn[maxn],low[maxn],from[maxn],idx,scn;
vector<int> scc[maxn];

void tarjan(int u){
	st[stn++]=u;
	ins[u]=1;
	dfn[u]=low[u]=++idx;
	for (int i=0;i<ed[u].size();i++)
		if (!dfn[ed[u][i]]){
			tarjan(ed[u][i]);
			low[u]=min(low[u],low[ed[u][i]]);
		}
		else if (ins[ed[u][i]])
			low[u]=min(low[u],low[ed[u][i]]); //either low or dfn are right
	if (low[u]==dfn[u]){
		int v;
		do{
			v=st[--stn];
			scc[scn].push_back(v);
			from[v]=scn;
			ins[v]=0;
		}while (u!=v);
		scn++;
	}
}
int qu[maxn],a[maxn],ind[maxn],sum[maxn],ans[maxn];
vector<int> ed2[maxn];

int circle_dp(){
	int n,m; scanf("%d%d",&n,&m);
	for (int i=0;i<n;i++) scanf("%d",a+i);
	for (int i=0;i<m;i++){
		int a,b; a--; b--;
		scanf("%d%d",&a,&b);
		ed[a].push_back(b);
	}
	inc(i,n)
		if (!dfn[i]) tarjan(i);
	inc(i,n)  //rebuild DAG
		for (int j=0;j<ed[i].size();j++)
			if (from[i]!=from[ed[i][j]])
				ed2[from[i]].push_back(from[ed[i][j]]),
				ind[from[ed[i][j]]]++;
	//for (int i=0;i<scn;i++,cout<<'\n') for (auto u:scc[i]) cout<<u<<' ';
	int head=0,tail=1;
	int tans=-1000000;
	for (int i=0;i<scn;i++) if (!ind[i]) qu[tail++]=i,ans[i]=sum[i],tans=max(tans,ans[i]);
	while (head<tail){
		int u=qu[head];
		for (int i=0;i<ed2[u].size();i++){
			int v=ed2[u][i];
			ind[v]--;
			ans[v]=max(ans[v],ans[u]+sum[v]),tans=max(tans,ans[v]);
			if (ind[v]==0) qu[tail++]=v;
		}
		head++;
	}
	cout<<tans<<'\n';
	return 0;
}

//mark dcc cut point
bool iscut[maxn];
void tarjan_point(int u, int fa){
	int ch=0; //sum of subtree link by this point
	low[u]=dfn[u]=++idx;
	for (int i=0;i<ed[u].size();i++){
		int v=ed[u][i];
		if (!dfn[v]){
			ch++;
			tarjan_point(v,u); //we can use siz[u] to record each size of subtree
			low[u]=min(low[u],low[v]);
			if (low[v]>=dfn[u]) iscut[u]=1;
		}
		else if (v!=fa)
			low[u]=min(low[u],dfn[v]); //dfn
	}
	if (fa<0 && ch==1) iscut[u]=0;
}
//dfs dcc block
int cnt,cucnt;
void dfs_dcc(int u, int dcn){
	vis[u]=dcn; 
	cnt++;
	for (int i=0;i<ed[u].size();i++){
		int v=ed[u][i];
		if (iscut[v] && vis[v]!=dcn) cucnt++,vis[v]=dcn; 
		if (!vis[v] && !iscut[i]) dfs_dcc(v, dcn);
	}
}
//cnt: inner point of current dcc block, cucnt: cut point in cur dcc block
void dcc_caller(){
	//add edge...
    //tarjan_point(1,-1); //mark dcc point
	int dcn=0;
	icc(i,n) if (!vis[i] && !iscut[i]){
		dcn++; cnt=cucnt=0;
		dfs_dcc(i, dcn);
		//cnt, cucnt, ...
	}
}

int ansx[maxn],ansy[maxn],ansc;
void tarjan_ed(int u, int fa){
	low[u]=dfn[u]=++idx;
	for (int i=0;i<ed[u].size();i++){
		int v=ed[u][i];
		if (!dfn[v]){
			tarjan_ed(v,u);
			low[u]=min(low[u],low[v]);
			if (low[v]>dfn[u])
				ansx[ansc]=u,ansy[ansc++]=v;
		}
		else if (v!=fa)
			low[u]=min(low[u],dfn[v]);
	}
}

#ifdef NO_COMPILE
int kruskal(){
	using namespace UFSet;
	int n,m; scanf("%d%d",&n,&m);
	for (int i=1;i<=n;i++) fa[i]=i;
	for (int i=0;i<m;i++) scanf("%d%d%d",&ed[i].x,&ed[i].y,&ed[i].c);
	sort(ed,ed+m);
	int ans=0;
	for (int i=0;i<m;i++){
		int ta=fi(ed[i].x), tb=fi(ed[i].y);
		if (ta!=tb)
			ans+=ed[i].c,
			fa[ta]=tb;
	}
	cout<<ans<<'\n';
	return ans;
}
#endif
//heap opt prim, O((n+m)log(m))
int prim(){
	memset(dis,0x3f,sizeof(dis));
	memset(vis,0,sizeof(vis));
	dis[1]=0;
	priority_queue<pair<int,int> > qu;
	qu.push(make_pair(0,1));
	int ans=0;
	while (qu.size()){
		int u=qu.top().second,c=qu.top().first;
		qu.pop();
		if (vis[u]) continue;
		ans-=c;
		vis[u]=1;
		for (int e=head[u];e;e=nxt[e])
			if (!vis[to[e]]){
				dis[to[e]]=co[e];
				qu.push(make_pair(-co[e],to[e]));
			}
	}
	return ans;
}
//O(n^2)
int prim_original(){
	memset(dis,0x3f,sizeof(dis));
	memset(vis,0,sizeof(vis));
	dis[1]=0; int ans=0;
	for(int i=0;i<n;i++){
		int u=0,md=0x3f3f3f3f;
		for (int j=0;j<n;j++)
			if (!vis[j] && dis[j]<md)
				md=dis[j], u=j;
		vis[u]=1; ans+=md;
		for (int j=0;j<n;j++)
			if (!vis[j])
				dis[j]=min(dis[j],mp[u][j]);
	}
	return ans;
}
#ifdef NO_COMPILE
int deg[maxn],ansp[maxm],al,anse[maxm<<1],el; //ansp has maxm point at most
//O(n+m), euler cycle
//the graph MUST BE only 0 or 2 odd node, and if there is 2 odd node, u must be odd node
//if more than 2 odd point, add virtual edge, and slice result by virtual edge
void euler(int u){
	for (int e=head[u];e;e=ed[e].nxt){
		int v=ed[e].to;
		if (!ed[e].vis && !ed[e^1].vis){
			ed[e].vis=ed[e^1].vis=1;
			euler(v);
			anse[++el]=e;//ed, the road and point will be reversed
		}
	}
	ansp[++al]=u;
}
#endif

#ifdef NO_COMPILE

//min directed tree, chu_liu's alg, O(VE)
const int maxn=1010,inf=0x3f3f3f3f;
struct Edge{
	int u,v,w;
}ed[maxn];
//in[]: min in edge weight, fa[]: min in vertex, id[]:scc id
int in[maxn],fa[maxn],vis[maxn],id[maxn];
int n,m,root;
ll chu_liu(){
	ll ans=0; int cnt,u,v;
	while (1){
		cnt=0; //circnt
		icc(i,n) in[i]=inf, vis[i]=id[i]=0;
		icc(i,m)
			if (ed[i].u!=ed[i].v&&ed[i].w<in[ed[i].v]) //min in edge
				fa[ed[i].v]=ed[i].u,in[ed[i].v]=ed[i].w;
		in[root]=0;
		icc(i,n){
			if (in[i]==inf) return -1;
			ans+=in[i];
			for (u=i;u!=root&&vis[u]!=i&&!id[u];u=fa[u]) //find circle
				vis[u]=i;
			if (u!=root&&!id[u]){
				id[u]=++cnt;
				for (v=fa[u];v!=u;v=fa[v]) id[v]=cnt; //label circle
			}
		}
		if (!cnt) return ans;
		icc(i,n) if(!id[i]) id[i]=++cnt; //single scc
		icc(i,m){ //build new graph on scc
			int laz=in[ed[i].v];
			if ((ed[i].u=id[ed[i].u])!=(ed[i].v=id[ed[i].v])) 
				ed[i].w-=laz; //new weight if exchange out cir edge
		}
		n=cnt; root=id[root];
	}
}
int chu_liu_caller(){
	scanf("%d%d%d",&n,&m,&root);
	icc(i,m) scanf("%d%d%d",&ed[i].u,&ed[i].v,&ed[i].w);
	cout<<chu_liu()<<'\n';
	return 0;
}
#endif
}

namespace BipartiteGraph{

const int maxn=500;
//Info: d[][] 2-d array is O(n^3), forward star is O(ne). 
//d[][]: n->m edge  ;  to[]: set m->n match index
int d[maxn][maxn],to[maxn],n,m;
bool vis[maxn];

//judge whether a graph is BipartiteGraph
bool judge(int u, int col){
	vis[u]=col;
	for (int i=0;i<n;i++)
		if (d[u][i] && (vis[i]==col || !vis[i] && !judge(i,-col)))
			return 0;
	return 1;
}

//index between 0..n-1, so must set to[] -1 before
bool xiong(int u){
	for (int i=0;i<m;i++)
		if (d[u][i] && !vis[i]){
			vis[i]=1;
			if (to[i]==-1 || xiong(to[i])){
				to[i]=u;
				return 1;
			}
		}
	return 0;
}
int match(){
	int ans=0;
	memset(to,-1,sizeof(to));
	for (int i=0;i<n;i++){
		memset(vis,0,sizeof(vis));
		if (xiong(i)) ans++;
	}
	return ans;
}
}
//KM alg, O(n^3), faster than NetworkFlow::MCMF
namespace KM{
const ll inf=0x3f3f3f3f3f3f3f3fll;
const int maxn=410;
ll d[maxn][maxn], to[maxn];
ll ux[maxn], uy[maxn], cx[maxn], cy[maxn], cc[maxn];
int n,m;
bool find(int u){ //Same as bit graph
	ux[u]=1;
	for (int i=1;i<=m;i++){
		if (uy[i]) continue;
		if (cx[u]+cy[i]==d[u][i]){
			uy[i]=1;
			if (!to[i]||find(to[i])) {
				to[i]=u;
				return 1;
			}
		}
		else cc[i]=min(cc[i],cx[u]+cy[i]-d[u][i]);
	}
	return 0;
}
ll km(){
	for (int i=1;i<=n;i++) for (int j=1;j<=m;j++) 
		cx[i]=max(cx[i],d[i][j]);
	for (int i=1;i<=n;i++){
		memset(ux,0,sizeof ux); memset(uy,0,sizeof uy);
		memset(cc,0x3f,sizeof cc);
		if (find(i)) continue;
		ll ms,u;
		while (1){
			ms=inf;
			for (int j=1;j<=n;j++) if (!uy[j]) ms=min(ms,cc[j]);
			for (int j=1;j<=n;j++) if (ux[j]) cx[j]-=ms;
			for (int j=1;j<=n;j++) if (uy[j]) cy[j]+=ms;
				else cc[j]-=ms,cc[j]==0?u=j:0;
			if (!to[u]) break;
			uy[u]=ux[to[u]]=1;
			u=to[u];
			for (int j=1;j<=m;j++) cc[j]=min(cc[j],cx[u]+cy[j]-d[u][j]);
		}
		memset(ux,0,sizeof ux); memset(uy,0,sizeof uy);
		find(i);
	}
	ll ans=0;
	for (int i=1;i<=m;i++) ans+=d[to[i]][i];
	return ans;
}
int fr[maxn];
int main(){
	int n0,m0,k; scanf("%d%d%d",&n0,&m0,&k);
	n=m=max(n0,m0);  //km require n=m
	for (int i=1;i<=k;i++){
		int a,b,x;
		scanf("%d%d%d",&a,&b,&x);
		d[a][b]=x;
	}
	cout<<km()<<'\n';
	for (int i=1;i<=m0;i++) if (d[to[i]][i]) fr[to[i]]=i; //ignore virtual edge
	for (int i=1;i<=n0;i++) printf("%d ",fr[i]); //all match
	return 0;
}
}

namespace MatchOnGraph{
//TreeWithFlower, O(n^3)
const int maxn=1010,maxm=150010;
struct Edge{
	int to,nxt;
}e[maxm<<1];
int head[maxn],ecnt=1;
void added(int a, int b){
	e[ecnt]={b,head[a]};
	head[a]=ecnt++;
}
//match[]:... . id[]:color of a point, -1 uncolored, 0 w, 1 b
//q[]:queue, pre[]: father in bfs tree
int match[maxn],id[maxn],q[maxn],pre[maxn];
int tim,tic[maxn],fa[maxn];
int find(int x){
	if (fa[x]!=x) fa[x]=find(fa[x]);
	return fa[x];
}
int lca(int x, int y){ //find lca without deep, only use pre[] and match[]
	tim++;
	for(tim++;;swap(x,y)) //cross jump up x,y and label road
		if (x){
			x=find(x); //flower point(flower root)
			if (tic[x]==tim) return x;
			else tic[x]=tim, x=find(pre[match[x]]); //jump up
		}
}
int st,ed;
void change(int x, int y, int k){ //circle: x<-->y, x&y as lca k
	while (find(x)!=k){
		pre[x]=y;
		int z=match[x]; id[z]=0; //recolor
		q[ed++]=z; if(ed>=maxn-1) ed=1; //try find new match on each node 
		if (find(z)==z) fa[z]=k; //shrink flower to point
		if (find(x)==x) fa[x]=k; //only find(x)==x is a shrinked point
		y=z;x=pre[y];
	}
}
int n;
bool check(int u){
	for (int i=0;i<=n;i++) fa[i]=i,id[i]=-1,pre[i]=0;
	st=1,ed=2;
	q[st]=u;id[u]=0;
	while (st!=ed){ //bfs argument road
		int x=q[st];
		for (int i=head[x],v=e[i].to;i;i=e[i].nxt,v=e[i].to){
			if (!match[v]&&v!=u){ //get a valid
				pre[v]=x;
				int last,t,now=v;
				while (now){ //cross road upd, same as bit graph
					t=pre[now];
					last=match[t];
					match[t]=now;match[now]=t;
					now=last;
				}
				return 1; //ok
			}
			if (id[v]==-1){ //not visted
				id[v]=1; pre[v]=x;
				id[match[v]]=0;  //cross color
				q[ed++]=match[v];
				if (ed>=maxn-1) ed=1;
			}
			else if(id[v]==0&&find(x)!=find(v)){ //odd circle
				int g=lca(x,v); //find lca in bfs tree
				change(x,v,g);  //shink x and its father
				change(v,x,g);  //shink v and its father
			}
			//even circle is the same as bit graph, so ignored
		}
		st++;
		if (st>=maxn-1) st=1;
	}
	return 0;
}
int main(){ //same as bit graph
	int m; 
	tim=0; memset(match,0,sizeof match);
	scanf("%d%d",&n,&m);
	for (int i=0;i<m;i++){
		int x,y;
		scanf("%d%d",&x,&y);
		added(x,y); added(y,x);
	}
	for (int i=1;i<=n;i++)
		if (!match[i])
			check(i);
	int ans=0;
	for (int i=1;i<=n;i++)
		if (match[i]) ans++;
	printf("%d\n",ans/2);
	for (int i=1;i<=n;i++)
		printf("%d ",match[i]);
	return 0;
}
}

namespace NetFlow{
#define INF 0x3f3f3f3f
const int maxn=1003,maxm=10003<<4;
struct Edge{
	int to,nxt,cap,cost; //assert cap>=0
}ed[maxm];
int head[maxn],ecnt=1,n,m;
void added(int a, int b, int cap){
	ed[++ecnt]=(Edge){b,head[a],cap,0};
	head[a]=ecnt;
	ed[++ecnt]=(Edge){a,head[b],0,0};
	head[b]=ecnt;
}
int s,t,a[maxn],fr[maxn],fp[maxn];
bool vis[maxn];
//deleted O(VE^2)
#ifdef USE_DELETED
int MF_EK(){
	int ans=0;
	while (1){
		memset(vis,0,sizeof(vis));
		memset(a,0,sizeof(a));
		a[s]=INF;
		queue<int> qu;
		qu.push(s);
		vis[s]=1;
		while (qu.size()){
			int u=qu.front(); qu.pop();
			if (u==t) break;
			for (int i=head[u];i;i=ed[i].nxt){
				int v=ed[i].to;
				if (!vis[v] && ed[i].cap){
					vis[v]=1;
					a[v]=min(a[u],ed[i].cap);
					fp[v]=u; fr[v]=i;
					qu.push(v);
				}
			}
		}
		if (!a[t]) break;
		ans+=a[t];
		for (int i=t;i!=s;i=fp[i]){
			ed[fr[i]].cap-=a[t];
			ed[fr[i]^1].cap+=a[t];
		}
	}
	return ans;
}
#endif

//isap, worst O(V^2E), but fast enough on most meaningful graph
//on layered graph is O(VE), and on unit capicity graph is O(E^1.5)
int now[maxn],num[maxn];
int isap_aug(int u, int f){
	if (u==t) return f;
	int ans=0;
	for (int i=now[u],v=ed[i].to;i;i=ed[i].nxt,v=ed[i].to)
		if (a[u]==a[v]+1){
			int w=isap_aug(v,min(f,ed[i].cap));
			ans+=w; f-=w; ed[i].cap-=w; ed[i^1].cap+=w; //aug
			if (!f) return ans;
			now[u]=i;
		}
	if (!(--num[a[u]])) a[s]=n+1; //gap opt
	++num[++a[u]]; now[u]=head[u];
	return ans;
}
int MF_isap(){
	memset(num,0,sizeof(num)); //num: label cnt
	memset(a,0,sizeof(a)); //a: label
	for (int i=1;i<=n;i++) now[i]=head[i];
	static int qu[maxn];
	int ql,qr=1; qu[ql=0]=t;
	++num[a[t]=1];
	while (ql<qr){ //optimize, bfs label at first
		int u=qu[ql++];
		for (int i=head[u],v=ed[i].to;i;i=ed[i].nxt,v=ed[i].to)
			if (!a[v]) ++num[a[v]=a[u]+1],qu[++qr]=v;
	}
	ll ret=isap_aug(s,INF);
	while (a[s]<=n) ret+=isap_aug(s,INF);
	return ret;
}

//deleted, similar as MF_isap
int dinic_dfs(int u, int f){
	int ans=0,w;
	if (u==t) return f;
	for (int i=now[u];i;i=ed[i].nxt){
		int v=ed[i].to;
		if (a[v]==a[u]+1 && ed[i].cap && (w=dinic_dfs(v,min(ed[i].cap,f)))){
			ans+=w;
			ed[i].cap-=w; ed[i^1].cap+=w;
			f-=w; if (f==0) return ans;
			now[u]=i;
		}
	}
	if (!ans) a[u]=-1;
	return ans;
}
int MF_dinic(){
	int ans=0;
	while (1){
		memset(vis,0,sizeof(vis));
		memset(a,0,sizeof(a)); //a: level
		queue<int> qu; qu.push(s); 
		vis[s]=1;
		while (qu.size()){ //BFS
			int u=qu.front(); qu.pop();
			if (u==t) break;
			for (int i=head[u];i;i=ed[i].nxt){
				int v=ed[i].to;
				if (!vis[v] && ed[i].cap){
					qu.push(v);
					a[v]=a[u]+1;
					vis[v]=1;
				}
			}
		}
		if (!vis[t]) break;
		for (int i=1;i<=n;i++) now[i]=head[i];
		ans+=dinic_dfs(s,INF);
	}
	return ans;
}

void added(int a, int b, int cap, int cost){
	ed[++ecnt]=(Edge){b,head[a],cap,cost};
	head[a]=ecnt;
	ed[++ecnt]=(Edge){a,head[b],0,-cost};
	head[b]=ecnt;
}
int dis[maxn];
int MCMF(){
	int ans=0,mf=0;
	while (1){
		memset(vis,0,sizeof(vis));
		memset(dis,0x3f,sizeof(dis));
		queue<int> qu; qu.push(s);
		dis[s]=0; vis[s]=1;
		while (qu.size()){ //spfa
			int u=qu.front(); qu.pop(); vis[u]=0;
			for (int i=head[u];i;i=ed[i].nxt){
				int v=ed[i].to;
				if (ed[i].cap && dis[v]>dis[u]+ed[i].cost){
					dis[v]=dis[u]+ed[i].cost;
					fr[v]=i; fp[v]=u;
					if (!vis[v]) vis[v]=1,qu.push(v);
				}
			}
		}
		if (dis[t]>INF/3) break;
		int mc=INF;
		for (int i=t;i!=s;i=fp[i]) mc=min(mc,ed[fr[i]].cap);
		for (int i=t;i!=s;i=fp[i]){
			ed[fr[i]].cap-=mc;
			ed[fr[i]^1].cap+=mc;
			ans+=mc*ed[fr[i]].cost;
		}
		mf+=mc;
	}
	cout<<mf<<' ';
	return ans;
}
//dijkstra version, more stable. 
//[!] The cost should be positive before first argument, or the
// complexity is not right. Run SPFA on first step is requied. 
int h[maxn];
int MCMF_dijk(){
	memset(h,0,sizeof(h)); //set h[all vertex] to 0
	int ans=0;
	while (1){
		priority_queue<pair<int,int>,vector<pair<int,int>>,greater<pair<int,int>>> qu;
		memset(dis,0x3f,sizeof(dis));
		dis[s]=0; qu.push({0,s});
		while (!qu.empty()){
			int du=qu.top().first, u=qu.top().second;
			qu.pop();
			if (dis[u]<du) continue;
			for (int i=head[u],v=ed[i].to;i;i=ed[i].nxt,v=ed[i].to)
				if (ed[i].cap && dis[v]>dis[u]+ed[i].cost+h[u]-h[v]){
					dis[v]=dis[u]+ed[i].cost+h[u]-h[v];
					fp[v]=u; fr[v]=i;
					qu.push({dis[v],v});
				}
		}
		if (dis[t]>INF/3) break;
		for (int i=0;i<=n;i++) h[i]+=dis[i];
		int mc=INF;
		for (int i=t;i!=s;i=fp[i]) mc=min(mc,ed[fr[i]].cap);
		for (int i=t;i!=s;i=fp[i]){
			ed[fr[i]].cap-=mc;
			ed[fr[i]^1].cap+=mc;
			ans+=mc*ed[fr[i]].cost;
		}
	}
	return ans;
}

//zkw cost flow, faster on wide and small capicity graph
int zkw_ans;
int dfs_aug(int u, int f){
	vis[u]=1;
	if (u==t) return f;
	int w,ad=0; 
	for (int i=head[u];i;i=ed[i].nxt){
		int v=ed[i].to;
		if (!vis[v] && ed[i].cap && dis[u]-ed[i].cost==dis[v] && (w=dfs_aug(v,min(f,ed[i].cap)))){
			zkw_ans+=w*ed[i].cost; ad+=w;
			ed[i].cap-=w; ed[i^1].cap+=w;
			if (f==0) break;
		}
	}
	return ad;
}
int MCMF_zkw(){
	int zkw_mf=0; zkw_ans=0;
	while (1){
		memset(vis,0,sizeof(vis));
		memset(dis,0x3f,sizeof(dis));
		queue<int> qu; qu.push(t);
		dis[t]=0; vis[t]=1;
		while (qu.size()){ //spfa on reverse path
			int u=qu.front(); qu.pop(); vis[u]=0;
			for (int i=head[u];i;i=ed[i].nxt){
				int v=ed[i].to;
				if (ed[i^1].cap && dis[v]>dis[u]-ed[i].cost){
					dis[v]=dis[u]-ed[i].cost;
					if (!vis[v]) vis[v]=1,qu.push(v);
				}
			}
		}
		if (dis[s]>INF/3) break;
		memset(vis,0,sizeof vis);
		zkw_mf+=dfs_aug(s, INF);
	}
	cout<<zkw_mf<<' ';
	return zkw_ans;
}

#undef INF
}

namespace TreeArr{
//WARN: index of tr[] statrs from 1
const int maxn=100010;
int tr[maxn]; int n;
void add(int p, int x){for(;p<=n;p+=p&-p)tr[p]+=x;}
ll sum(int p){ll ret=0;for(;p;p-=p&-p)ret+=tr[p];return ret;}

//section add and section sum version, section [l,r]
template <typename X>    
struct tree_array{    
    struct tree_array_single{    
        X tr[maxn];    
        void add(int p,X x){for(;p<=n;p+=p&-p)tr[p]+=x;}    
        X sum(int p){ll ret=0;for(;p;p-=p&-p)ret+=tr[p];return ret;}    
    }T1,T2;    
    void add(int p,X x){T1.add(p,x);T2.add(p,p*x);}      
    X sum(int p){return (p+1)*T1.sum(p)-T2.sum(p);}
public:
    void update(int l,int r,int x){add(l,x);add(r+1,-x);}  
    X query(int l,int r){return sum(r)-sum(l-1);}    
};
}

namespace SEGT{
const int MAXN=100010;

ll sum[MAXN<<2], tadd[MAXN<<2], tmul[MAXN<<2], a[MAXN];
ll n,m,p;
#define lc u+u+1
#define rc u+u+2
void build(int u, int l, int r){
	tmul[u]=1;
	if (l==r-1){
		sum[u]=a[l];
		return;
	}
	int mid=l+r>>1;
	build(lc,l,mid); build(rc,mid,r);
	sum[u]=(sum[lc]+sum[rc])%p;
}
void upd(int u, int l, int r){
	int mid=l+r>>1;
	sum[lc]*=tmul[u]; sum[lc]+=(mid-l)*tadd[u]; sum[lc]%=p;
	sum[rc]*=tmul[u]; sum[rc]+=(r-mid)*tadd[u]; sum[rc]%=p;
	tadd[lc]*=tmul[u]; tadd[lc]+=tadd[u]; tadd[lc]%=p;
	tmul[lc]*=tmul[u]; tmul[lc]%=p;
	tadd[rc]*=tmul[u]; tadd[rc]+=tadd[u]; tadd[rc]%=p;
	tmul[rc]*=tmul[u]; tmul[rc]%=p;
	tadd[u]=0; tmul[u]=1;
}
void mul(int u, int l, int r, int cl, int cr, ll c){
	if (cl<=l && cr>=r){
		tadd[u]*=c; tadd[u]%=p;
		tmul[u]*=c; tmul[u]%=p;
		sum[u]*=c; sum[u]%=p;
		return;
	}
	if (tadd[u] || tmul[u]!=1) upd(u,l,r);
	int mid=l+r>>1;
	if (cl<mid) mul(lc,l,mid,cl,cr,c);
	if (cr>mid) mul(rc,mid,r,cl,cr,c);
	sum[u]=(sum[lc]+sum[rc])%p;
}
void add(int u, int l, int r, int cl, int cr, ll c){
	if (cl<=l && cr>=r){
		tadd[u]+=c; tadd[u]%=p;
		sum[u]+=c*(r-l)%p; sum[u]%=p;
		return;
	}
	if (tadd[u] || tmul[u]!=1) upd(u,l,r);
	int mid=l+r>>1;
	if (cl<mid) add(lc,l,mid,cl,cr,c);
	if (cr>mid) add(rc,mid,r,cl,cr,c);
	sum[u]=(sum[lc]+sum[rc])%p;
}
ll ask(int u, int l, int r, int cl, int cr){
	if (cl<=l && cr>=r) return sum[u];
	if (tadd[u] || tmul[u]!=1) upd(u,l,r);
	int mid=l+r>>1;
	ll ret=0;
	if (cl<mid) ret+=ask(lc,l,mid,cl,cr);
	if (cr>mid) ret+=ask(rc,mid,r,cl,cr);
	return ret%p;
}
ll pointask(int u, int l, int r, int q){
	if (l==r-1) return sum[u];
	if (tadd[u] || tmul[u]!=1) upd(u,l,r);
	int mid=l+r>>1;
	if (q<mid) return pointask(lc,l,mid,q);
	return pointask(u,mid,r,q);
}

#undef lc
#undef rc

}


namespace FSEGT{
const int maxn=200010;
int sum[maxn<<5], root[maxn], lc[maxn<<5], rc[maxn<<5],trcnt;
int a[maxn],b[maxn];

void build(int &u, int l, int r){
	u=trcnt++;
	if (l==r-1) return;
	int mid=l+r>>1;
	build(lc[u],l,mid); build(rc[u],mid,r);
}
int mod(int u, int l, int r, int c){
	int v=trcnt++;
	lc[v]=lc[u]; rc[v]=rc[u]; sum[v]=sum[u]+1;
	if (l==r-1) return v;
	int mid=l+r>>1;
	if (c<mid) lc[v]=mod(lc[v],l,mid,c);
	else rc[v]=mod(rc[v],mid,r,c);
	return v;
}
int query(int u, int v, int l, int r, int q){
	int mid=l+r>>1, x=sum[lc[v]]-sum[lc[u]];
	if (l==r-1) return l;
	if (x>=q) return query(lc[u],lc[v],l,mid,q);
	return query(rc[u],rc[v],mid,r,q-x);
}
//ask: [l,r) kth number
int main(){
	int n,m;
	cin>>n>>m;
	for (int i=0;i<n;i++)
		scanf("%d", a+i),b[i]=a[i];
	sort(b,b+n);
	int n1=unique(b,b+n)-b;
	build(root[0],0,n1);
	for (int i=0;i<n;i++){
		int q=lower_bound(b,b+n1,a[i])-b;
		root[i+1]=mod(root[i],0,n1,q);
	}
	for (int i=0;i<m;i++){
		int l,r,q;
		scanf("%d%d%d",&l,&r,&q);
		printf("%d\n",b[query(root[l-1],root[r],0,n1,q)]);
	}
	return 0;
}
}

namespace KDT{
const int N=1000010, inf=0x3f3f3f3f;
int n,m,K,rt,ans;

//s[]:tree node  p[2]:2-d coord of leaf node  x[2]:min(LB) coord of a subspace  y[2]:max(RT) coord
struct Node{
	int p[2],x[2],y[2];
	bool operator<(const Node &v)const{
		return p[K]<v.p[K];
	}
}a[N],s[N],q;
int ch[N][2];
#define lc ch[u][0]
#define rc ch[u][1]
void upd(int u){
	inc(i,2){
		if (lc) s[u].x[i]=min(s[u].x[i],s[lc].x[i]),
				s[u].y[i]=max(s[u].y[i],s[lc].y[i]);
		if (rc) s[u].x[i]=min(s[u].x[i],s[rc].x[i]),
				s[u].y[i]=max(s[u].y[i],s[rc].y[i]);
	}
}
void add(int u, Node &t){
	inc(i,2) s[u].x[i]=s[u].y[i]=s[u].p[i]=t.p[i];
}
int disL1Min(int u, Node &t){ //min L1 dis to a Rect of in_tree node
	int ret=0;
	inc(i,2) 
		if (t.p[i]>s[u].y[i]) ret+=t.p[i]-s[u].y[i];
		else if (t.p[i]<s[u].x[i]) ret+=s[u].x[i]-t.p[i];
	return ret;
}
int disL1Max(int u, Node &t){
	int ret=0;
	inc(i,2) ret+=max(abs(t.p[i]-s[u].x[i]),abs(t.p[i]-s[u].y[i]));
	return ret;
}
int sqr(int a){
	return a*a;
}
int disL2Min(int u, Node &t){
	int ret=0;
	inc(i,2) 
		if (t.p[i]>s[u].y[i]) ret+=sqr(t.p[i]-s[u].y[i]);
		else if (t.p[i]<s[u].x[i]) ret+=sqr(t.p[i]-s[u].x[i]);
	return ret;
}
int disL2Max(int u, Node &t){ //max coord dis
	int ret=0;
	inc(i,2) ret+=max(sqr(t.p[i]-s[u].x[i]),sqr(t.p[i]-s[u].y[i]));
	return ret;
}
void build(int &u, int l, int r, int cur){ //O(nlogn)
	u=l+r>>1; K=cur;
	nth_element(a+l,a+u,a+r+1);
	add(u,a[u]);
	if (l<u) build(lc,l,u-1,cur^1);
	if (r>u) build(rc,u+1,r,cur^1);
	upd(u);
}
//Maybe we need to rebuild the tree after unbalanced insert
void ins(int u, int cur){  
	if (q.p[cur]<s[u].p[cur])
		if (lc) ins(lc,cur^1);
		else lc=++n,add(n,q);
	else
		if (rc) ins(rc,cur^1);
		else rc=++n,add(n,q);
	upd(u);
}
void ask(int u){
	ans=min(ans,abs(s[u].p[0]-q.p[0])+abs(s[u].p[1]-q.p[1])); //L1 dis
	int dl=lc?disL1Min(lc,q):inf, dr=rc?disL1Min(rc,q):inf;
	//int dl=lc?disL1Max(lc,q):0, dr=rc?disL1Max(rc,q):0;
	if (dl<dr){ //trim branch, swap > < when search max dis point
		if (dl<ans) ask(lc);
		if (dr<ans) ask(rc);
	}
	else{
		if (dr<ans) ask(rc);
		if (dl<ans) ask(lc);
	}
}
//minDisPoint (L1 dis) with ins operate
//each query asks one nearest point of a giving coord
int main(){
	scanf("%d%d",&n,&m);
	for (int i=1;i<=n;i++) scanf("%d%d",&a[i].p[0],&a[i].p[1]);
	build(rt,1,n,0);
	while (m--){
		int k; scanf("%d%d%d",&k,&q.p[0],&q.p[1]);
		if (k==1) ins(rt,0);
		else{
			ans=inf; ask(rt);
			printf("%d\n",ans);
		}
	}
	return 0;
}
#undef lc
#undef rc
}

namespace Heap{
//Alorithm heap
//Run Fast
const int maxn=1000001;
int heap[maxn+1],hc;
int demo(){
    int n; rd(n);
    for (int i=0;i<n;i++){
        int c; rd(c);
        if (c==1){     //insert
            int x; rd(x);
            heap[hc++]=-x;
            push_heap(heap,heap+hc);
        }
        else if (c==2) //min element
            printf("%d\n",-heap[0]);
        else           //delete
            pop_heap(heap,heap+hc),hc--;
    }
    return 0;
}
}

namespace Treap{
//TT: an ordered struct
typedef int TT;
const int maxn=100001;
struct Node{
	//x: number, s: sum size of cur and subtree, cnt: cnt of cur num
	Node *c[2];
	TT x;
	int s,r,cnt;
	Node(TT _x){c[0]=c[1]=0;s=cnt=1;x=_x;r=rand();}
	Node(){};
}tree[maxn<<1];
#define lc u->c[0]
#define rc u->c[1]
#define lcs (lc?lc->s:0)
#define rcs (rc?rc->s:0)
int trcnt=0;
Node *open(TT x){
	tree[trcnt++]=Node(x);
	return tree+trcnt-1;
}
void upd(Node *u){
	u->s=lcs+rcs+u->cnt;
	//more updates...
}
//rt: set lc to root
void rot(Node* &u, int d){ //0 lt, 1 rt
	Node *t=u->c[d^1]; u->c[d^1]=t->c[d]; t->c[d]=u;
	t->s=u->s; upd(u); u=t;
}
void ins(Node* &u, TT x){
	if (!u){u=open(x);return;}
	if (x==u->x) {u->cnt++;u->s++; return;}
	int d=x>u->x; u->s++;
	ins(u->c[d],x);
	if (u->c[d]->r>u->r) rot(u,d^1);
}
void delx(Node* &u, TT x){
	if (!u) return;
	if (x==u->x){
		if (u->cnt>1) u->cnt--, u->s--;
		else if (!lc || !rc) u=max(lc,rc);
		else{
			rot(u,lc->r>rc->r);
			u->s--,delx(u->c[u->x<x],x);
		}
	}
	else u->s--,delx(u->c[u->x<x],x);
}
int rk(Node *u, TT x){
	if (!u) return 0;
	if (u->x==x) return lcs + 1;
	if (u->x<x) return lcs + u->cnt + rk(rc,x);
	else return rk(lc, x);
}
//get point by element
Node* findx(Node *u, TT x){
	if (!u) return 0;
	if (x==u->x) return u;
	return findx(u->c[u->x<x],x);
}
//get point by rank
//r=(1~tree_size)
Node* findr(Node *u, int r){
	if (!u) return 0;
	if (r<=lcs) return findr(lc,r);
	r-=lcs;
	if (r<=u->cnt) return u;
	r-=u->cnt;
	return findr(rc,r);
}
TT pred(Node *u, TT x){
	if (!u) return -0x3f3f3f3f;
	if (u->x<x) return max(u->x,pred(rc,x));
	else return pred(lc,x);
}
TT succ(Node *u, TT x){
	if (!u) return 0x3f3f3f3f;
	if (x<u->x) return min(u->x,succ(lc,x));
	else return succ(rc,x);
}
void dfs(Node *u, int deep=0){
	if (lc) dfs(lc,deep+1);
	for (int i=0;i<deep;i++) cout<<"   ";
	cout<<u->x<<' '<<u->s<<'\n';
	if (rc) dfs(rc,deep+1);
}
void caller(){
	Node *root=0;
	int T;cin>>T;
	while (T--)	{
		int c,x; scanf("%d%d",&c,&x);
		if (c==1) ins(root,x);
		if (c==2) delx(root,x);
		if (c==3) cout<<rk(root,x)<<'\n';
		if (c==4) cout<<findr(root,x)->x<<'\n';
		if (c==5) cout<<pred(root,x)<<'\n';
		if (c==6) cout<<succ(root,x)<<'\n';
		//dfs(root),cout<<'\n';
	}
}
#undef lc
#undef rc
#undef lcs
#undef rcs
}
namespace Splay{
const int maxn=100010;

int val[maxn],siz[maxn],ch[maxn][2],pa[maxn],cnt[maxn];
bool rev[maxn];
int root,trcnt;
#define lc ch[u][0]
#define rc ch[u][1]
//pushup
void upd(int u){
	siz[u]=cnt[u]+siz[lc]+siz[rc];
}
//lazy tags
void pushdown(int u){
	if (rev[u]){
		rev[lc]^=1;rev[rc]^=1;
		swap(lc,rc);
		rev[u]=0;
	}
}
void rot(int u, int c){
	int p=pa[u];
	ch[p][!c]=ch[u][c];
	pa[ch[u][c]]=p; pa[u]=pa[p];
	if (pa[u]) ch[pa[p]][ch[pa[p]][1]==p]=u;
	ch[u][c]=p; pa[p]=u;
	upd(p); upd(u);
}
//u->under s
void splay(int u, int s){
	while (pa[u]!=s){
		int p=pa[u],pp=pa[p];
		if (pp==s) rot(u,ch[p][0]==u);
		else{
			int c=(ch[pp][0]==p);
			if (ch[p][c]==u) rot(u,!c); else rot(p,c);
			rot(u,c);
		}
	}
	if (s==0) root=u;
}
//rank k->under s
void rk(int k, int s=0){
	int u=root;
	assert(k>=1 && k<=siz[root]);
	while (1){
		pushdown(u);
		if (k<=siz[lc]) u=lc;
		else if (k>siz[lc]+cnt[u]) k-=siz[lc]+cnt[u],u=rc;
		else break;
	}
	splay(u,s);
}
//x->under s
void fi(int x, int s=0){
	int u=root,p;
	while (x!=val[u] && u)
		p=u,u=ch[u][x>val[u]];
	if (u && x==val[u]) splay(u,s);
	else if (!u) splay(p,s);
}

#ifdef NO_COMPILE
//memory restricted open
int avail[maxn],avac;
void open(int &u, int x){
	u=avail[--avac];
	lc=rc=pa[u]=0;
	siz[u]=cnt[u]=1;
	val[u]=x;
}
void pre(){
	for (int i=1;i<maxn;i++) 
		avail[avac++]=maxn-i;
	open(root, -10086); // in section problem, add virtual point is convenient
	int r; open(r, -10086);
	ch[root][1]=r; pa[r]=root;
	upd(root);
}
void release(int u){
	if (!u) return;
	release(lc);
	release(rc);
	avail[avac++]=u;
}
#endif
void open(int &u, int x){
	u=++trcnt;
	lc=rc=pa[u]=0;
	siz[u]=cnt[u]=1;
	val[u]=x;
}
//root, value, parent
void ins(int &u, int x, int p){
	if (!u) open(u,x),pa[u]=p;
	else if (val[u]==x) cnt[u]++,siz[u]++;
	if (!u || val[u]==x){splay(u,0);return;}
	ins(ch[u][val[u]<x],x,u);
	upd(u);
}
//delete root
void del_(){
	int u=root;
	if (rc){
		root=rc; rk(1,0); //right, though it's hard to understand
		ch[root][0]=lc;
		if (ch[root][0]) pa[ch[root][0]]=root;
	}
	else root=lc;
	pa[root]=0;
	upd(root);
}
void del(int x){
	fi(x,0);
	if (val[root]==x)
		if (cnt[root]>1) cnt[root]--,siz[root]--;
		else del_();
}
int pred(int u, int x){
	if (!u) return -0x3f3f3f3f;
	if (val[u]<x) return max(val[u],pred(rc,x));
	return pred(lc,x);
}
int succ(int u, int x){
	if (!u) return 0x3f3f3f3f;
	if (x<val[u]) return min(val[u],succ(lc,x));
	return succ(rc,x);
}
//90 degree rotate debug print
void debug(int u=root, int deep=0){
	if (!u) return;
	debug(rc, deep+1);
	for (int i=0;i<deep;i++) cout<<"  ";
	cout<<val[u]<<' '<<siz[u]<<'\n';
	debug(lc, deep+1);
}
int n,m;
void dfs(int u){
	if (!u) return;
	pushdown(u);
	dfs(lc);
	if (val[u]>0 && val[u]<=n) cout<<val[u]<<' ';
	dfs(rc);
}
void mian(){
	int T; cin>>T;
	while (T--){
		int c,x; scanf("%d%d",&c,&x);
		if (c==1)
			ins(root,x,0);
		else if (c==2)
			del(x);
		else if (c==3){ //get rk of x
			fi(x,0); 
			cout<<siz[ch[root][0]]+1<<'\n';
		}
		else if (c==4){ //get k th
			rk(x,0);
			cout<<val[root]<<'\n';
		}
		else if (c==5){ //pred
			ins(root,x,0);
			rk(siz[ch[root][0]],0);
			cout<<val[root]<<'\n';
			del(x);
		}
		else if (c==6){ //succ
			ins(root,x,0);
			rk(siz[ch[root][0]]+cnt[root]+1,0);
			cout<<val[root]<<'\n';
			del(x);
		}
		//debug(root,0);
	}
}
int main(){ //reverse
	cin>>n>>m;
	for (int i=0;i<=n+1;i++) ins(root,i,0);
	for (int i=0;i<m;i++){
		int l,r; scanf("%d%d",&l,&r);
		rk(l,0); rk(r+2,root);
		rev[ch[ch[root][1]][0]]^=1;
		//debug(root,0);
	}
	dfs(root); putchar('\n');
	return 0;
}
#undef lc
#undef rc
}
namespace NRTreap{
const int maxn=100010;
typedef int TT;
struct Node{
	Node *c[2];
	TT x;
	int s, r;
	bool rev;
}tree[maxn<<1];
typedef pair<Node *,Node *> PD;
int trcnt;
Node *root;
Node *open(int x){
	tree[trcnt++]=(Node){0,0,x,1,rand(),0};
	return tree+trcnt-1;
}
#define lc u->c[0]
#define rc u->c[1]
#define lcs (lc?lc->s:0)
#define rcs (rc?rc->s:0)
void upd(Node *u){
	u->s=lcs+rcs+1;
}
void pushdown(Node *u){
	if (u->rev){
		if (lc) lc->rev^=1;
		if (rc) rc->rev^=1;
		swap(lc,rc);
		u->rev=0;
	}
}
Node *merge(Node *u, Node *v){
	if (!u || !v) return max(u,v);
	pushdown(u); pushdown(v);
	if (u->r<v->r) {rc=merge(rc,v);upd(u);return u;}
	else {v->c[0]=merge(u,v->c[0]);upd(v);return v;} 
}
PD split(Node *u, int k){
	if (!u) return MP((Node *)0,(Node *)0);
	pushdown(u);
	PD t;
	if (k<=lcs){
		t=split(lc,k);
		lc=t.second;
		upd(u);
		t.second=u;
	}
	else{
		t=split(rc,k-lcs-1);
		rc=t.first;
		upd(u);
		t.first=u;
	}
	return t;
}
int rk(Node *u, TT x){
	if (!u) return 0;
	if (u->x<x) return lcs + 1 + rk(rc,x);
	else return rk(lc, x);
}
int findr(Node *u, int r){
	if (!u) return 0;
	if (r<=lcs) return findr(lc,r);	r-=lcs;
	if (r==1) return u->x; r--;
	return findr(rc,r);
}
void ins(TT x){
	int k=rk(root,x);
	PD t=split(root,k);
	Node *u=open(x);
    root=merge(t.first,merge(u,t.second));
}
//t1.second is deleted
void del(TT x){
	int k=rk(root,x);
	PD t1=split(root,k),t2=split(t1.second,1);
	root=merge(t1.first,t2.second);
}
void debug(Node *u, int deep=0){
	if (lc) debug(lc,deep+1);
	for (int i=0;i<deep;i++) cout<<"   ";
	cout<<u->x<<' '<<u->s<<' '<<u->rev<<'\n';
	if (rc) debug(rc,deep+1);
}
int n;
void dfs(Node *u){
	if (!u) return;
	pushdown(u);
	dfs(lc);
	if (u->x>0 && u->x<=n) cout<<u->x<<' ';
	dfs(rc);
}
int mian(){
	int T;cin>>T;
	while (T--)	{
		int c,x; scanf("%d%d",&c,&x);
		if (c==1) ins(x);
		if (c==2) del(x);
		if (c==3) cout<<rk(root,x)+1<<'\n';
		if (c==4) cout<<findr(root,x)<<'\n';
		if (c==5) cout<<findr(root,rk(root,x))<<'\n';
		if (c==6) cout<<findr(root,rk(root,x+1)+1)<<'\n';
		//dfs(root),cout<<'\n';
	}
	return 0;
}
int main(){ //reverse
	int m;cin>>n>>m;
	for (int i=0;i<=n+1;i++) ins(i);
	for (int i=0;i<m;i++){
		int l,r; scanf("%d%d",&l,&r);
		PD x=split(root,l);
		PD y=split(x.second,r-l+1);
		y.first->rev^=1;
		root=merge(x.first,merge(y.first,y.second));
		//dfs(root); putchar('\n');
		//debug(root);
	}
	dfs(root); putchar('\n');
	return 0;
}
#undef lc
#undef rc
}

namespace LCA{
const int maxn=500010;
struct Edge{
	int to,nxt;
	//int w;
}ed[maxn<<1];
int head[maxn],ecnt;
void added(int x, int y){
	ed[++ecnt].to=y;
	ed[ecnt].nxt=head[x];
	head[x]=ecnt;
}
bool vis[maxn];
//online query, sum, minimax
namespace TreeMultiply{
int dep[maxn];
int pa[24][maxn];
//minn[24][maxn], sumv[24][maxn];
int lca(int u, int v){
	//int tu=u,tv=v,sum=0;
	if (dep[u]<dep[v]) swap(u,v);
	for (int k=23;k>=0;k--) 
		if (dep[pa[k][u]]>=dep[v]) u=pa[k][u]; //,sum+=sumv[k][u];
	if (u!=v){
		for (int k=23;k>=0;k--) 
			if (pa[k][u]!=pa[k][v]) u=pa[k][u],v=pa[k][v];//,sum+=sumv[k][u]+sum[k][v];
		u=pa[0][u];//sum+=sumv[0][u]+sumv[0][v];
	}
	//int lenth=dep[tu]+dep[tv]-2*dep[u];
	return u;
}
bool vis[maxn];
void dfs_deep(int u){
	vis[u]=1;
	for (int e=head[u];e;e=ed[e].nxt){
		int v=ed[e].to;
		if (!vis[v]){
			dep[v]=dep[u]+1; pa[0][v]=u;
			for (int k=1;pa[k-1][pa[k-1][v]];k++)
				pa[k][v]=pa[k-1][pa[k-1][v]];
		//	sumv[0][u]=minn[0][u]=maxn[0][u]=ed[e].w;
		//	for (int k=1;pa[k-1][pa[k-1][u]];k++) 
		//		sumv[k][u]=sumv[k-1][u]+sumv[k-1][pa[k-1][u]];
			dfs_deep(v);
		}
	}
}
//no stackoverflow
int qu[maxn];
void bfs_deep(int s){
	int qh=0,qt=1;
	qu[0]=s;dep[s]=1;vis[s]=1;
	while (qh<qt){
		int u=qu[qh];
		for (int e=head[u];e;e=ed[e].nxt){
			int v=ed[e].to;
			if (!vis[v]){
				dep[v]=dep[u]+1;
				vis[v]=1; pa[0][v]=u;
				for (int k=1;pa[k-1][pa[k-1][v]];k++)
					pa[k][v]=pa[k-1][pa[k-1][v]];
			//	sumv[0][u]=minn[0][u]=maxn[0][u]=ed[e].w;
			//	for (int k=1;pa[k-1][pa[k-1][u]];k++) 
			//		sumv[k][u]=sumv[k-1][u]+sumv[k-1][pa[k-1][u]];
				qu[qt++]=v;
			}
		}
		qh++;
	}
}
//n: node, m: query, s: root node
void process(){
	int n,m,s,a,b; cin>>n>>m>>s;
	//for (int i=1;i<=n;i++) rd(w[i]);
	for (int i=0;i<n-1;i++){
		rd(a); rd(b);
		added(a,b);
		added(b,a);
	}
	dep[s]=1;dfs_deep(s);
	//bfs_deep(s);
	for (int i=0;i<m;i++){
		rd(a); rd(b);
		printf("%d\n",lca(a,b));
	}
}
}
//offline lca
namespace tarjan{
using namespace UFSet;
int ans[maxn]; //lca point
struct QEdge{
	int to,nxt,s;
}qc[maxn<<1];
int quh[maxn],qcnt;
void addqu(int i,int a,int b){
	qc[++qcnt]=(QEdge){b,quh[a],i};
	quh[a]=qcnt;
	qc[++qcnt]=(QEdge){a,quh[b],i};
	quh[b]=qcnt;
}
//dfs, notify system stack
void tarjan(int u, int f){
	vis[u]=1;
	for (int e=head[u];e;e=ed[e].nxt){
		int v=ed[e].to;
		if (!vis[v])
			tarjan(v,u);
	}
	for (int i=quh[u];i;i=qc[i].nxt)
		if (vis[qc[i].to])
			ans[qc[i].s]=fi(qc[i].to);
	un(u,f);
}
//n: node, m: query, s: root node
void process(){
	int n,m,s,a,b; cin>>n>>m>>s;
	for (int i=0;i<=n;i++) fa[i]=i;
	for (int i=0;i<n-1;i++){
		rd(a); rd(b);
		added(a,b);
		added(b,a);
	}
	for (int i=0;i<m;i++){
		rd(a); rd(b);
		addqu(i,a,b);
	}
	tarjan(s,0);
	for (int i=0;i<m;i++)
		printf("%d\n",ans[i]);
}
}
}

namespace SplitTree{
const int maxn=100010;
struct Edge{
	int to,nxt;
}ed[maxn<<1];
int head[maxn],ecnt;
void added(int a, int b){
	ed[++ecnt]=(Edge){b,head[a]};
	head[a]=ecnt;
}
int sum[maxn<<2],tadd[maxn<<2],a[maxn],n,P,ucnt;
//son: heavy son, top: chain top, rk: segnode->treenode, id: treenode->segnode
int w[maxn],dep[maxn],fa[maxn],son[maxn],
	siz[maxn],rk[maxn],top[maxn],id[maxn];
#define lc u+u+1
#define rc u+u+2
void build(int u, int l, int r){
	if (l==r-1){
		sum[u]=a[rk[l]];
		return;
	}
	int mid=l+r>>1;
	build(lc,l,mid); build(rc,mid,r);
	sum[u]=(sum[lc]+sum[rc])%P;
}
void upd(int u, int l, int r){
	int mid=l+r>>1;
	sum[lc]+=(mid-l)*tadd[u]; sum[rc]+=(r-mid)*tadd[u];
	tadd[lc]+=tadd[u]; tadd[rc]+=tadd[u];
	sum[lc]%=P; sum[rc]%=P; tadd[lc]%=P; tadd[rc]%=P;
	tadd[u]=0;
}
void add(int u, int l, int r, int cl, int cr, int c){
	if (cl<=l && cr>=r){
		tadd[u]+=c; tadd[u]%=P;
		sum[u]+=c*(r-l)%P; sum[u]%=P;
		return;
	}
	if (tadd[u]) upd(u,l,r);
	int mid=l+r>>1;
	if (cl<mid) add(lc,l,mid,cl,cr,c);
	if (cr>mid) add(rc,mid,r,cl,cr,c);
	sum[u]=(sum[lc]+sum[rc])%P;
}
int ask(int u, int l, int r, int cl, int cr){
	if (cl<=l && cr>=r) return sum[u];
	if (tadd[u]) upd(u,l,r);
	int mid=l+r>>1;
	int ret=0;
	if (cl<mid) ret+=ask(lc,l,mid,cl,cr);
	if (cr>mid) ret+=ask(rc,mid,r,cl,cr);
	return ret%P;
}
void dfs1(int u, int f, int deep){
	fa[u]=f; dep[u]=deep; siz[u]=1;
	for (int i=head[u];i;i=ed[i].nxt){
		int v=ed[i].to; 
		if (v==f) continue;
		dfs1(v,u,deep+1);
		siz[u]+=siz[v];
		if (siz[v]>siz[son[u]]) son[u]=v;
	}
}
void dfs2(int u, int t){
	top[u]=t; id[u]=++ucnt; rk[ucnt]=u; 
	if (son[u]) dfs2(son[u],t);
	for (int i=head[u];i;i=ed[i].nxt){
		int v=ed[i].to;
		if (v!=son[u] && v!=fa[u]) dfs2(v,v);
	}
}
int askt(int x, int y){
	int ans=0;
	int fx=top[x],fy=top[y];
	while (fx!=fy)
		if (dep[fx]>=dep[fy]){
			ans=(ans+ask(0,0,n+1,id[fx],id[x]+1))%P;
			x=fa[fx],fx=top[x];
		}
		else{
			ans=(ans+ask(0,0,n+1,id[fy],id[y]+1))%P;
			y=fa[fy],fy=top[y];
		}
	if (id[x]>id[y]) swap(x,y);
	return (ans+ask(0,0,n+1,id[x],id[y]+1))%P;
}
void addt(int x, int y, int c){
	int fx=top[x],fy=top[y];
	while (fx!=fy)
		if (dep[fx]>=dep[fy]){
			add(0,0,n+1,id[fx],id[x]+1,c);
			x=fa[fx],fx=top[x];
		}
		else{
			add(0,0,n+1,id[fy],id[y]+1,c);
			y=fa[fy],fy=top[y];
		}
	if (id[x]>id[y]) swap(x,y);
	add(0,0,n+1,id[x],id[y]+1,c);
}
//r: root, c: 1-add chain, 2-sum chain, 3-add subtree, 4-sum subtree
void process(){
	int m,r,x,y;
	scanf("%d%d%d%d",&n,&m,&r,&P);
	for (int i=1;i<=n;i++) scanf("%d",a+i);
	for (int i=1;i<n;i++){
		scanf("%d%d",&x,&y);
		added(x,y); added(y,x);
	}
	dfs1(r,0,1);
	dfs2(r,r);
	build(0,0,n+1);
	for (int i=0;i<m;i++){
		int c,z; scanf("%d%d",&c,&x);
		if (c==1){
			scanf("%d%d",&y,&z);
			addt(x,y,z);
		}else if(c==2){
			scanf("%d",&y);
			printf("%d\n",askt(x,y));
		}
		else if (c==3){
			scanf("%d",&y);
			add(0,0,n+1,id[x],id[x]+siz[x],y);
		}
		else
			printf("%d\n",ask(0,0,n+1,id[x],id[x]+siz[x]));
	}
}
}

namespace DividePoint{
const int maxn=20010,maxm=40010;

struct Edge{
	int to,nxt,c;
}e[maxm];
int ec,n,head[maxn];

void added(int a, int b, int c){
	e[++ec]={b,head[a],c};
	head[a]=ec;
	e[++ec]={a,head[b],c};
	head[b]=ec;
}

int query[maxn],q,siz[maxn],ms[maxn];
int MS,root,tn;
bool vis[maxn];

void dfs(int u,int fa, int len){
	;//counter
	for (int i=head[u],v=e[i].to;i;i=e[i].nxt,v=e[i].to)
		if (v!=fa && !vis[v])
			dfs(v,u,len+e[i].c);
}
int calc(int u, int x0){
	dfs(u,u,x0);
	return 0; //return count
}
void getrt(int u, int fa){
	siz[u]=1; ms[u]=0;
	for (int i=head[u],v=e[i].to;i;i=e[i].nxt,v=e[i].to)
		if (v!=fa && !vis[v])
			getrt(v,u),
			siz[u]+=siz[v],ms[u]=max(ms[u],siz[v]);
	ms[u]=max(ms[u],tn-siz[u]);
	if (ms[u]<MS) root=u,MS=ms[u];
}
int ans=0;
void divide(int u){
	vis[u]=1;
	ans+=calc(u,0);
	for (int i=head[u],v=e[i].to;i;i=e[i].nxt,v=e[i].to)
		if (!vis[v]){
			ans-=calc(v,e[i].c); //sub invalid path
			tn=siz[u]; root=0;
			MS=0x3f3f3f3f; getrt(v,u);
			divide(root);
		}
}
int main(){
	scanf("%d",&n);
	for (int i=0;i<n-1;i++){
		int a,b,c; scanf("%d%d%d",&a,&b,&c);
		added(a,b,c);
	}
	tn=n; root=0; MS=0x3f3f3f3f; getrt(1,1); //first point divide
	divide(root);
	cout<<ans<<'\n';
	return 0;
}
}

namespace CalcGeo{

typedef double db;
const db PI=acos(-1);
const db eps=1e-10, inf=1e12;

bool eq(db x){return fabs(x)<eps;}
int sgn(db x){
	if (x<=-eps) return -1;
	return x>=eps;
}

#define Vec const vec &
#define Point const point &
struct vec{
	db x,y;
	vec():x(0),y(0){}
	vec(db x, db y):x(x),y(y){} //[i] init-list is easier to use in c++1x
	vec(db theta):x(cos(theta)),y(sin(theta)){}

	bool operator==(Vec v) const{return eq(x-v.x) && eq(y-v.y);}
	db ang() const{return atan2(y,x);}
	
	vec operator+(Vec v) const{return vec(x+v.x,y+v.y);}
	vec operator-(Vec v) const{return vec(x-v.x,y-v.y);}
	vec operator*(db a) const{return vec(x*a,y*a);}
	vec operator/(db a) const{return vec(x/a,y/a);}
	
	db operator|(Vec v) const{return x*v.x+y*v.y;} //dot
	db operator&(Vec v) const{return x*v.y-y*v.x;} //cross
	db operator!() const{return sqrt(x*x+y*y);}    //len
	
	bool operator<(Vec v) const{return x==v.x?y<v.y:x<v.x;}
	//rotate countclockwise
	vec std() const{return *this/!*this;}
	vec rot(db rad) const{return vec(x*cos(rad)-y*sin(rad), x*sin(rad)+y*cos(rad));}
	vec l90() const{return vec(-y,x);}
	vec r90() const{return vec(y,-x);}
	vec vert() const{return (*this).l90().std();}  //l90 and standard
	
	void rd(){scanf("%lf%lf",&x,&y);}
	void prt(){printf("%f %f\n",x,y);}
	friend ostream& operator<<(ostream &o, Vec v){return o<<v.x<<','<<v.y;}
};
typedef vec point;

//cmp2 sort by angle without atan2
// angle range [-pi,pi)
bool cmp2(Vec a, Vec b){
	int d1=sgn(a.y),d2=sgn(b.y);
	if (d1^d2&&d1&&d2) return d1<d2;
	if (d2==0&&d2==0) return a.x<b.x;
	return (a&b)>0;
}

db angle(Vec a, Vec b){return fabs(atan2(a&b,a|b));}
db cross(Point a, Point b, Point c){return b-a & c-a;}
db dot(Point a, Point b, Point c){return b-a | c-a;}

//cosine theory
db angle(db a, db b, db c){return acos((a*a+b*b-c*c)/(2*a*b));}

//Line: P=P0+t*vp
// Segment: 0<=t<=1
//intersection point of line P and Q
point lineInt(Point p, Vec vp, Point q, Vec vq){
	db t=(vq & p-q)/(vp&vq);
	return p+vp*t;
}
//point projection on line A+tV
point lineProj(Point p, Point s, Vec v){
	return s+v*(v|p-s)/(v|v);
}
//symmetric point of P about line A+tV
point symmetric(Point p, Point s, Vec v){
	return lineProj(p,s,v)*2-p;
}
//distance of p to line A+tV
db lineDis(Point p, Point s, Vec v){
	return fabs(v & p-s)/!v;
}
//distance of p to segment A+tV
db segDis(Point p, Point s, Vec v){
	if (eq(!v)) return !(s-p); //a point
	vec v2=p-s,v3=p-s-v;
	if ((v|v2)<0) return !v2;
	else if ((v|v3)>0) return !v3;
	return fabs(v&v2)/!v;
}
//distance of seg A-B and seg C-D
db segDis(Point a, Point b, Point c, Point d){
	vec u=b-a, v=d-c;
	return min(min(segDis(c,a,u),segDis(d,a,u)),min(segDis(a,c,v),segDis(b,c,v)));
}

//point is on line
bool onLine(Point p, Point a, Point b){return eq(p-a&b-a);}
//point on seg [a,b]
bool onSeg(Point p, Point a, Point b){return onLine(p,a,b) && sgn(a-p|b-p)<=0;}

//fast test before segment cross, 0 indicate the segment are not cross 
bool rectCover(Point a1, Point a2, Point b1, Point b2){return 
	min(a1.x,a2.x)<=max(b1.x,b2.x)+eps &&
	min(b1.x,b2.x)<=max(a1.x,a2.x)+eps &&
	min(a1.y,a2.y)<=max(b1.y,b2.y)+eps &&
	min(b1.y,b2.y)<=max(a1.y,a2.y)+eps;
}
//test if segment A1-A2 B1-B2 is cross
int segCross(Point a1, Point a2, Point b1, Point b2){
	if (!rectCover(a1,a2,b1,b2)) return 0; //not necessary
	db c1=sgn(a2-a1&b1-a1), c2=sgn(a2-a1&b2-a1);
	db c3=sgn(b2-b1&a1-b1), c4=sgn(b2-b1&a2-b1);
	if (c1*c2>0 || c3*c4>0) //no cross
		return 0; 
	if (c1==0 && c2==0||c3==0 && c4==0) //segment on same line
		return -1; 
	if (c1*c2<0 && c3*c4<0) return 1; //normal cross
	return 2; //a point on line
}

//#define const line& Line

struct line{
	point p; vec v;
	line(){}
	line(Point p, db ang):p(p),v(ang){}
	//ax+by+c=0
	line(db a, db b, db c){
		if (eq(b)) p=point(-c/a,0), v=vec(0,1);
		else p=point(0,-c/b),v=vec(1,-a/b);
	}
	line(Point p, Vec v):p(p),v(v){}
	bool operator<(const line &l) const{return v.ang()<l.v.ang();}
};


struct circle{
	point c; db r;
	
	circle(){}
	circle(Point c, db r):c(c),r(r){}
	circle(Point p1, Point p2):c((p1+p2)/2),r(!(p1-p2)/2){}
	//circle passing point P1P2P3
	circle(Point p1, Point p2, Point p3){ //[!] p1,p2,p3 should not on same line
		//c=(p1+lineInt(p2,(p2-p1).l90(),p3,(p3-p1).l90()))/2; //this impl not good
		vec B=p2-p1,C=p3-p1; db D=B&C*2;
		c=vec(C.y*(B|B)-B.y*(C|C),B.x*(C|C)-C.x*(B|B))/D+p1;
		r=!(p1-c);
	}
	//inscribed cricle of triangle P1P2P3
	circle(Point p1, Point p2, Point p3, bool _){
		db x=!(p2-p3),y=!(p3-p1),z=!(p1-p2);
		c=(p1*x+p2*y+p3*z)/(x+y+z);
		r=lineDis(c,p1,p2-p1);
	}
	point angle(db theta){return c+point(theta)*r;}

	bool operator==(const circle &v) const{return c==v.c && eq(r-v.r);}
};


//point in or on circle
bool inCir(Point p, circle c){return sgn(!(c.c-p)-c.r)<=0;}

//return -1,0,1,2, ans[2]
int cirCross(circle A, circle B, point *ans){
	db d=!(A.c-B.c);
	if (eq(d)){
		if (eq(A.r-B.r)) return -1; //same circle
		return 0; //same center
	}
	if (sgn(d-fabs(A.r-B.r))<0) return 0; //inside
	if (sgn(A.r+B.r-d)<0) return 0; //too far
	db a=(B.c-A.c).ang();
	db da=angle(A.r,d,B.r);
	ans[0]=A.angle(a-da),ans[1]=A.angle(a+da);
	if (eq(da) || eq(da-PI)) return 1; //tang
	return 2; //normal inter
}

//get tangent points on circle from point p
//return  ans[2] : tangent point 
int cirTang(Point p, circle c, point *ans){
	db d=!(c.c-p);
	if (sgn(d-c.r)<0) return 0;
	if (eq(d-c.r)){ //[!] notice this time ans[0]-p0 not a line
		ans[0] = p; //ans[0]=(p-c.c).vert()+p; //to get a line
		return 1;
	}
	db base=(p-c.c).ang();
	db ang=acos(c.r/d);
	ans[0]=c.angle(base-ang);
	ans[1]=c.angle(base+ang);
	return 2;
}

//get cir-cir common tangent line
//return  a[4],b[4] : tangent point on circle
//cnt maybe -1(same), 0(in), 1(in tangent), 2(cross), 3(out tangent), 4(out) 
int cirTang(circle A, circle B, point *a, point *b){
	int cnt=0;
	if (A.c==B.c && eq(A.r-B.r)) return -1;
	if (A.r<B.r) swap(A,B),swap(a,b);
	db d=!(A.c-B.c);
	db diff=A.r-B.r, sum=A.r+B.r;
	if (sgn(d-diff)<0) return 0;
	db base=(B.c-A.c).ang();
	if (eq(d-diff)){
		a[0] = A.angle(base);
		b[0] = a[0];
		return 1;
	}
	db ang=acos((A.r-B.r)/d);
	//in common tangent
	a[cnt]=A.angle(base+ang); b[cnt++]=B.angle(base+ang);
	a[cnt]=A.angle(base-ang); b[cnt++]=B.angle(base-ang);
	if (eq(d-sum)){
		a[cnt] = A.angle(base);
		b[cnt] = a[cnt];
		cnt++;
	} else if (sgn(d-sum)>0){ //out common tangent
		ang=acos((A.r+B.r)/d);
		a[cnt]=A.angle(base+ang); b[cnt++]=B.angle(PI+base+ang); 
		a[cnt]=A.angle(base-ang); b[cnt++]=B.angle(PI+base-ang); 
	}
	return cnt;
}

//line A-B cross circle c point
//return  ans[2] : cross or tangent point
int lineInt(Point a, Point b, circle c, point *ans){
	vec v=b-a, u=a-c.c;
	db e=v|v, f=v|u*2, g=(u|u)-c.r*c.r;
	db delta=f*f-4*e*g;
	if (delta<0) return 0;
	if (eq(delta)) return ans[0]=a-v*(f/2/e),1;
	db t1=(-f-sqrt(delta))/(2*e);
	db t2=(-f+sqrt(delta))/(2*e);
	ans[0]=a+v*t1;
	ans[1]=a+v*t2;
	return 2;
}
//seg A-B cross circle c point
//return  ans[2] : cross or tangent point
int segInt(Point a, Point b, circle c, point *ans){
	vec v=b-a, u=a-c.c;
	db e=v|v, f=v|u*2, g=(u|u)-c.r*c.r;
	db delta=f*f-4*e*g;
	if (delta<0) return 0;
	if (eq(delta)){
		db t=f/2/e;
		if (sgn(t)>=0 && sgn(t-1)<=0) return ans[0]=a-v*t,1;
		return 0;
	}
	db t1=(-f-sqrt(delta))/(2*e); //[i] t1 is closer to a because f<0
	db t2=(-f+sqrt(delta))/(2*e);
	point a1=a+v*t1, a2=a+v*t2;
	int cnt=0;
	if (sgn(t1)>=0 && sgn(t1-1)<=0) ans[cnt++]=a1;
	if (sgn(t2)>=0 && sgn(t2-1)<=0) ans[cnt++]=a2;
	return cnt;
}

//insection area of two circle a and b
//!-- test
db cirIntArea(circle a, circle b){
	db d=!(a.c-b.c);
	if (sgn(d-a.r-b.r)>=0) return 0; // too far
	if (sgn(d-fabs(a.r-b.r)<=0)) return PI*min(a.r,b.r)*min(a.r,b.r); //inside
	db hf=(a.r+b.r+d)/2;
	db s=-2*sqrt(hf*(hf-a.r)*(hf-b.r)*(hf-d));
	s+=angle(a.r,d,b.r)*a.r*a.r;
	s+=angle(b.r,d,a.r)*b.r*b.r;
	return s;
}

//UVA12304 get circles with radius r and other conditions
//circle passing A and B with radius r
// return ans[2]: center circle
int getCir(Point a, Point b, db r, point *ans){
	//circle A(a,r),B(b,r); return cirCross(A,B,r); //another implement
	db d=!(a-b)/2;
	if (sgn(d-r)<0) return 0;
	vec v=(b-a)/2;
	if (eq(d-r)) return ans[0]=a+v,1;
	ans[0]=a+v+v.vert()*sqrt(r*r-d*d);
	ans[1]=a+v-v.vert()*sqrt(r*r-d*d);
	return 2;
}
//circle with radius r passing point A and tangent with line P
int getCir(Point a, Point p, vec vp, db r, point *ans){
	if (eq(vp&a-p)){ //special judge point A on line P
		ans[0]=p+vp.vert()*r;
		ans[1]=p-vp.vert()*r;
		return 2;
	}
	//implement by line-cir-intersection 
	//point p1=p+vp.vert()*sgn(vp&a-p)*r; 
	//return lineInt(p1,p1+vp,circle(a,r));
	//independent implement
	point p0=lineProj(a,p,vp); db d=!(a-p0);
	if (sgn(2*r-d)<0) return 0;
	if (eq(2*r-d)) return ans[0]=(a+p0)/2,1;
	point p1=p0+vp.vert()*sgn(vp&a-p)*r;
	d-=r; d=sqrt(r*r-d*d);
	vp=vp.std()*d;
	ans[0]=p1+vp;
	ans[1]=p1-vp;
	return 2;
}
//circle with radius r and tangent with non-parallel line P and Q
int getCir(point p, vec vp, point q, vec vq, db r, point *ans){
	vec mvp=vp.vert()*r; //move dir
	vec mvq=vq.vert()*r;
	ans[0]=lineInt(p-mvp,vp,q-mvq,vq);
	ans[1]=lineInt(p-mvp,vp,q+mvq,vq);
	ans[2]=lineInt(p+mvp,vp,q-mvq,vq);
	ans[3]=lineInt(p+mvp,vp,q+mvq,vq);
	return 4;
}
//circle with radius r and tangent with disjoint circle c1 and c2
int getCir(circle c1, circle c2, db r, point *ans){
	return cirCross(circle(c1.c,c1.r+r),circle(c2.c,c2.r+r),ans);
}

//inverse circle C by  circle P with radius 1
circle inverseCir(Point p, circle c){
	db d=!(c.c-p);
	if (eq(c.r-d)) //get a line
		return circle(vec(inf,inf),inf);
	db d2=1/(d+c.r);
	db d1=1/(d-c.r);
	vec v=c.c-p; v=v/!v;
	return circle(p+v*(d1+d2)/2,(d1-d2)/2);
}
circle inverseCir(Point p, Point s, Vec v){
	point t=lineProj(p,s,v);
	db r=0.5/!(t-p); //radius
	return circle(p+(t-p)/!(t-p)*r,r);
}

//--poly--

//point is in or on polygon
//return  1(in), 0(out), -1(on border)
int inPoly(point p, point *poly, int n){
	int w=0;
	for (int i=0;i<n;i++){
		if (onSeg(p,poly[i],poly[(i+1)%n])) 
			return -1;
		int k=sgn(poly[(i+1)%n]-poly[i] & p-poly[i]);
		int d1=sgn(poly[i].y-p.y);
		int d2=sgn(poly[(i+1)%n].y-p.y);
		if (k>0 && d1<=0 && d2>0) w++;
		if (k<0 && d2<=0 && d1>0) w--;
	}
	return w!=0;
}
//test segment strict in poly, 0 out/border, 1 in
bool inPoly(point p1, point p2, point *poly, int n){
	if (!inPoly(p1,poly,n) || !inPoly(p2,poly,n)) return 0;
	for (int i=0;i<n;i++)
		if (segCross(p1,p2,poly[i],poly[(i+1)%n]))
			return 0;
	return 1;
}
//point at border regard as in poly
const db epr=1e-5; //[!] epr should larger than eps*tan(angle(segment and poly-edge))
bool inPoly2(point p1, point p2, point *poly, int n){
	if (inPoly(p1*epr+p2*(1-epr),poly,n)==0 || inPoly(p2*epr+p1*(1-epr),poly,n)==0) return 0;
	for (int i=0;i<n;i++)
		if (segCross(p1,p2,poly[i],poly[(i+1)%n])==1)
			return 0;
	return 1;
}
bool outPoly(point p1, point p2, point *poly, int n){
	if (inPoly(p1,poly,n) || inPoly(p2,poly,n)) return 0;
	for (int i=0;i<n;i++)
		if (segCross(p1,p2,poly[i],poly[(i+1)%n]))
			return 0;
	return 1;
}
bool outPoly2(point p1, point p2, point *poly, int n){
	if (inPoly(p1*epr+p2*(1-epr),poly,n)==1 || inPoly(p2*epr+p1*(1-epr),poly,n)==1) return 0;
	for (int i=0;i<n;i++)
		if (segCross(p1,p2,poly[i],poly[(i+1)%n])==1)
			return 0;
	return 1;
}

// [!] Require simple polygon, or the result will be strange
db polyArea(point *p, int n){
	db sum=0;
	for (int i=1;i<n-1;i++)
		sum+=cross(p[0],p[i+1],p[i]);
	return fabs(sum)/2;
}

point polyBaryCenter(point *p, int n){
	point ret(0,0);
	for (int i=1;i<n-1;i++)
		ret=ret+(p[0]+p[i]+p[i+1])/3;
	return ret;
}

//convex hull, Andrew algo
// return  ans[m]
int convex(point *p, int n, point *ans){
	sort(p,p+n);
	n=unique(p,p+n)-p;
	int m=0;
	for (int i=0;i<n;i++){
		while (m>1 && cross(ans[m-2],ans[m-1],p[i])<=0) m--;
		ans[m++]=p[i];
	}
	int k=m;
	for (int i=n-2;i>=0;i--){
		while (m>k && cross(ans[m-2],ans[m-1],p[i])<=0) m--;
		ans[m++]=p[i];
	}
	if (n>1) m--; //p[0]==p[m]
	return m;
}

//test point P strictly in convex polygon, o(nlogn)
bool inConvex(point *p, int n, point q){ //require countclockwise convex hull
	if (sgn(cross(p[0],q,p[1]))>=0 || sgn(cross(p[0],p[n-1],q))>=0) return 0;
	int s=lower_bound(p+1,p+n,q,[&](Point a, Point b){return sgn(cross(p[0],a,b))>0;})-p;
	return sgn(cross(p[s-1],p[s],q))>0;
}

//cut convex polygon by line A, return right side of remain poly
int convexCut(point *p, int n, point a, vec v, point *ans){
	int c=0;
	for (int i=0;i<n;i++){
		int d1=sgn(v&p[i]-a);
		int d2=sgn(v&p[(i+1)%n]-a);
		if (d1>=0) ans[c++]=p[i];
		if (d1*d2<0) ans[c++]=lineInt(a,v,p[i],p[(i+1)%n]-p[i]); //cut
	}
	return c;
}
//Minkowski sum
int msum(point *p, int n, point *q, int m, point *ans){ 
	int i,j,res=0; //res=m+n
	#if 1          //flip q by (0,0) when calc moving q convex clip p
	inc(i,m) q[i]=(vec){0,0}-q[i];
	rotate(q,min_element(q,q+m),q+m); //be sure 0 leftest
	#endif
	p[n]=p[0]; q[m]=q[0];
	for (i=0,j=0;i<n && j<m;)
		if ((p[i+1]-p[i]&q[j+1]-q[j])>0) ans[res++]=p[++i]+q[j];
		else ans[res++]=p[i]+q[++j];
	for (;i<n;) ans[res++]=p[++i]+q[j];
	for (;j<m;) ans[res++]=p[i]+q[++j];
	return res; //[!] notice coliner point in result
}
//weather point p is lefter than line l
bool onleft(point p, const line &l){return sgn(l.v&p-l.p)>0;}
const int maxp=1001;
line Q[maxp<<1]; //deque of lines
point T[maxp<<1]; //deque of points(result)
//intersection of left side half plane, return countclockwise polygon point
//[!] The result area can't be unlimited.
int halfplaneInt(line *l, int n, point *ans){
	sort(l,l+n); //[!] This operation changed input
	int head=0,tail=0; //rangeof Q:[head,tail] ; range of T: [head, tail)
	Q[0]=l[0];
	for (int i=1;i<n;i++){
		while (head<tail && !onleft(T[tail-1],l[i])) tail--;
		while (head<tail && !onleft(T[head],l[i])) head++;
		Q[++tail]=l[i];
		if (eq(Q[tail].v&Q[tail-1].v)){ //same direction
			--tail;
			if (onleft(l[i].p,l[i])) Q[tail]=l[i]; //replace righter line
		}
		if (head<tail) //get point
			T[tail-1]=lineInt(Q[tail-1].p,Q[tail-1].v,Q[tail].p,Q[tail].v);		
	}
	while (head<tail && !onleft(T[tail-1],Q[head])) tail--; 
	if (head>=tail-1) return 0;  //m<3, no available area
	T[tail]=lineInt(Q[head].p,Q[head].v,Q[tail].p,Q[tail].v); //head cross tail
	int m=0;
	for (int i=head;i<=tail;i++) ans[m++]=T[i];
	return m;
}
//half plane intersection with unlimted space judge
int halfplaneInt_(line *l, int n, point *ans){ //[!] array l should have 4 extra space
	l[n]=line(point(-inf,-inf),vec(1,0));
	l[n+1]=line(point(inf,-inf),vec(0,1));
	l[n+2]=line(point(inf,inf),vec(-1,0));
	l[n+3]=line(point(-inf,inf),vec(0,-1));
	int ret=halfplaneInt(l,n+4,ans);
	for (int i=0;i<ret;i++)
		if (fabs(ans[i].x)>inf/2 || fabs(ans[i].y)>inf/2)
			return -1; //unlimited
	return ret;
}

//--rotating stuck--

const int maxn=100010;
//max dis point pair on poly
// (farthest point pair on plane)
db polyDiam(point *p0, int n0){
	static point p[maxn];
	int n=convex(p0,n0,p); //[!] p0 changed
	p[n]=p[0];
	int opp=1; db ans=!(p[0]-p[1]);
	for (int i=0;i<n;i++){
		while (cross(p[i],p[i+1],p[opp+1])>cross(p[i],p[i+1],p[opp])) opp=(opp+1)%n;
		ans=max(ans, max(!(p[opp]-p[i]),!(p[opp]-p[i+1])));
	}
	return ans;
}
//min dis between parallel lines clip polygon
db polyWidth(point *p0, int n0){
	static point p[maxn];
	int n=convex(p0,n0,p); //[!] p0 changed
	p[n]=p[0];
	int opp=1; db ans=1e10;
	for (int i=0;i<n;i++){
		while (cross(p[i],p[i+1],p[opp+1])>cross(p[i],p[i+1],p[opp])) opp=(opp+1)%n;
		ans=min(ans, lineDis(p[opp],p[i],p[i+1]-p[i]));
	}
	return ans;
}
//min rectangle area cover polygon
db minRectCover(point *p0, int n0){
	static point p[maxn];
	int n=convex(p0,n0,p); //[!] p0 changed
	if (n<3) return 0;
	p[n]=p[0];
	db ans=-1;
	int h=1,r=1,l;
	for (int i=0;i<n;i++){
		while (cross(p[i],p[i+1],p[h+1])-cross(p[i],p[i+1],p[h])>=0) h=(h+1)%n; //farest
		while (dot(p[i],p[i+1],p[r+1])-dot(p[i],p[i+1],p[r])>=0) r=(r+1)%n; //rightest
		if (i==0) l=h;
		while (dot(p[i],p[i+1],p[l+1])-dot(p[i],p[i+1],p[l])<=0) l=(l+1)%n; //leftest
		db t=p[i+1]-p[i]|p[i+1]-p[i];
		db s=cross(p[i],p[i+1],p[h])*(dot(p[i],p[i+1],p[r])-dot(p[i],p[i+1],p[l]))/t; //rect area
		//min circumference of rectangle
		//db c=2*(cross(p[i],p[i+1],p[h])+dot(p[i],p[i+1],p[r])-dot(p[i],p[i+1],p[l]))/!(p[i+1]-p[i]);
		if (ans<0 || ans>s) ans=s;
	}
	return ans;
}
//minimum convex hull distanse (actually what support vector machine do on plane)
//[!] require non-cross countclockwise convex hull
db minConvexDis(point *p, int n, point *q, int m){ 
	p[n]=p[0]; q[n]=q[0];
	db ans=inf; int r=0;
	for (int i=1;i<m;i++)
		if (cross(p[0],p[1],q[i])>cross(p[0],p[1],q[r]))
			r=i;
	for (int i=0;i<n;i++){
		while (cross(p[i],p[i+1],q[r+1])>cross(p[i],p[i+1],q[r])){
			r=(r+1)%m;
			ans=min(ans,segDis(p[i],p[i+1],q[r],q[r+1]));
		}
	}
	return ans;
}
//inner common tangent line of two convex hull, O(n+m)
//return one of tangent line (postion in input array)
//[!] require non-corss countclockwise convex hull)
pair<int,int> convexInnerTang(point *p, int n, point *q, int m){ 
	p[n]=p[0]; q[n]=q[0];
	int r=0;
	for (int i=1;i<m;i++)
		if (cross(p[0],p[1],q[i])>cross(p[0],p[1],q[r]))
			r=i;
	for (int i=0;i<n;i++){
		while (cross(p[i],p[i+1],q[r+1])>cross(p[i],p[i+1],q[r])){
			r=(r+1)%m;
			if (cross(p[(i+n-1)%n],p[i],q[r])>=0 && cross(p[i],p[i+1],q[r])<0 &&
				cross(q[(r+m-1)%m],q[r],p[i])>=0 && cross(q[r],q[r+1],p[i])<0) //change here to get another tangent line
				return {i,r};
		}
	}
	throw;
}

//---complex---

//sector a~b of radius r
db secArea(point a, point b, db r){return r*r*angle(a,b)/2;}
db triArea(point a, point b){return fabs(a&b)/2;}
//intersection area of circle C and triangle P1-P2-C
db tri_cirArea(point p1, point p2, circle c){
	db r=c.r;
	p1=p1-c.c; p2=p2-c.c;
	c.c=vec(0,0);
	point p[2];
	if (sgn(!p1-r)<0){ //p1 in circle
		if (sgn(!p2-r)<0) return triArea(p1,p2);
		segInt(p1,p2,c,p);
		return triArea(p1,p[0])+secArea(p[0],p2,r);
	}
	if (sgn(!p2-r)<0){ //p2 in circle
		segInt(p1,p2,c,p);
		return secArea(p1,p[0],r)+triArea(p[0],p2);
	}
	int pc=segInt(p1,p2,c,p);
	if (pc==2) return secArea(p1,p[0],r)+triArea(p[0],p[1])+secArea(p[1],p2,r);
	return secArea(p1,p2,r);	
}
//intersection area of polygon P and circle C
db poly_cirArea(point *p, int n, circle c){
	db ans=0;
	for (int i=0;i<n;i++){
		db d=sgn(cross(c.c,p[i],p[(i+1)%n]));
		ans+=d*tri_cirArea(p[i],p[(i+1)%n],c);
	}
	return fabs(ans);
}

//min circle corver point set p
//average O(n)
circle mincirCover(point *p, int n){
    random_shuffle(p,p+n); //[!] This operation changed input
    circle c(p[0],0);
    for (int i=1;i<n;i++)
        if (sgn(!(c.c-p[i])-c.r)>0){
            c=circle(p[i],0);
            for (int j=0;j<i;j++)
                if (sgn(!(c.c-p[j])-c.r)>0){
                    c=circle(p[i],p[j]);
                    for (int k=0;k<j;k++)
                        if (sgn(!(c.c-p[k])-c.r)>0)
                            c=circle(p[i],p[j],p[k]);
                }
        }
    return c;
}


//union area of circles
namespace Circles{
	const int N=1010; //O(n^2log(n))
	circle c[N];
	db ans[N],pre[N];
	int n;
	//remove inside or same circles
	void init(){ //[!] c[N] changed
		sort(c,c+n,[](const circle &a, const circle &b){return a.c==b.c?a.r<b.r:a.c<b.c;});
		n=unique(c,c+n)-c; //use circle::operator==
	}
	db arcarea(db rad, db r){return 0.5*r*r*(rad-sin(rad));}
	//union area of circles
	// ans[1] is union area
	// ans[i]-ans[i+1] is k times intersection area; [!] the circles should be unique
	db areaunion(){
		memset(ans,0,sizeof ans);
		vector<pair<db,int>> v; //int 1: start of section | -1: end of section
		init(); //delete inside; [!] should NOT init() when get k-times intersection
		for (int i=0;i<n;i++){
			v.clear();
			v.emplace_back(-PI,1); //default [-PI,PI] full circle
			v.emplace_back(PI,-1); 
			for (int j=0;j<n;j++) //label arc secions
				if (i^j){
					point q=c[j].c-c[i].c;
					db d=!q, x=c[i].r, y=c[j].r;
					if (sgn(d+x-y)<=0){ //cover by circle[j]
						v.emplace_back(-PI,1);
						v.emplace_back(PI,-1);
						continue;
					}
					if (sgn(d-x+y)<=0) continue; //cover circle[j]
					if (sgn(d-x-y)>0) continue; //too far
					db base=q.ang(), ang=angle(x,d,y);
					db a0=base-ang;	if (sgn(a0+PI)<0) a0+=2*PI;
					db a1=base+ang; if (sgn(a1-PI)>0) a1-=2*PI;
					v.emplace_back(a0,1);
					if (sgn(a0-a1)>0){ //arc across 180 degree
						v.emplace_back(PI,-1);
						v.emplace_back(-PI,1);
					}
					v.emplace_back(a1,-1);
				}
			sort(v.begin(),v.end());
			int cur=0;
			for (auto &a:v){ //point
				if (cur && sgn(a.first-pre[cur])){
					ans[cur]+=arcarea(a.first-pre[cur],c[i].r); //arcarea
					ans[cur]+=(c[i].angle(pre[cur])&c[i].angle(a.first))/2; //piece of center polygon area(signed)
				}
				cur+=a.second;
				pre[cur]=a.first;
			}
		}
		//for (int i=1;i<n;i++)
		//	ans[i]-=ans[i+1];
		return ans[1];
	}
}

void test(){
	vec a(1.2,2.5);
	vec b(1.4,1.3);
	vec c(1,2),vc(0,1);
	vec d(3,1),vd(-3,1.5);
	vec ep(eps/2,-eps/2);
	cout<<a+b<<" expect 2.6 3.8\n";
	cout<<a-b<<" expect -0.2 1.2\n";
	cout<<a*2<<" expect 2.4 5\n";
	cout<<b/2<<" expect 0.7 0.65\n";
	cout<<(a|b)<<" expect 4.93\n";
	cout<<(a&b)<<" expect -1.94\n";
	cout<<(b&a)<<" expect 1.94\n";
	cout<<(a==b)<<" expect 0\n";
	cout<<(a==a+ep)<<" expect 1\n";
	cout<<!a<<" expect 2.77308\n";
	cout<<(a|a)<<" expect 7.69\n";
	cout<<(c.ang())<<" expect 1.10715\n";
	cout<<(c.rot(PI/2))<<" expect -2 1\n";
	cout<<(c.rot(-PI/2))<<" expect 2 -1\n";
	cout<<c.vert()<<" expect -0.8944 0.4472\n";
	cout<<angle(c,d)<<" expect "<<c.ang()-d.ang()<<'\n';
	cout<<lineInt(c,vc,d,vd)<<" expect 1 2\n";
	cout<<lineInt(d,vd,c,vc)<<" expect 1 2\n";
	cout<<lineDis(point(0,0),d,vec(0,2.5)-d)<<" expect 2.23607\n";
	cout<<segDis(point(0,0),d,vec(0,2.5)-d)<<" expect 2.23607\n";
	cout<<segDis(point(0,5),d,vec(0,2.5)-d)<<" expect 2.5\n";
	cout<<lineProj(point(0,0),d,vec(4,0)-d)<<" expect 2 2\n";
	
	cout<<onLine(point(2,2),d,vec(4,0))<<" expect 1\n";
	cout<<onSeg(point(2,2),d,vec(4,0))<<" expect 0\n";
	cout<<onSeg(point(3.5,0.5),d,vec(4,0))<<" expect 1\n";
	cout<<onSeg(point(4,0),d,vec(4,0))<<" expect 1\n";
	
	cout<<segCross(point(2,2),point(0,0),d,vec(0,4))<<" expect 2\n";
	cout<<segCross(point(3,3),point(0,0),d,vec(0,4))<<" expect 1\n";
	cout<<segCross(point(0,4),point(0,0),d,vec(0,4))<<" expect 2\n";
	cout<<segCross(point(1,1),point(0,0),d,vec(0,4))<<" expect 0\n";
	cout<<segCross(point(2,2),point(-1,5),d,vec(0,4))<<" expect -1\n";
	cout<<segCross(point(0,4),point(-1,5),d,vec(0,4))<<" expect -1\n";
	
	point ans[2];
	circle c1(point(0,1),1),c2(point(0,0),1);
	cout<<cirCross(c1,c2,ans)<<" expect 2\n";
	cout<<ans[0]<<' '<<ans[1]<<" expect -0.866 0.5 0.866 0.5\n";
	
	c1=circle(point(0,1),1),c2=c1;
	cout<<cirCross(c1,c2,ans)<<" expect -1\n";
	
	c1=circle(point(0,1),1),c2=circle(point(4,4),1);
	cout<<cirCross(c1,c2,ans)<<" expect 0\n";
	
	c1=circle(point(0,1),1),c2=circle(point(0,0),2);
	cout<<cirCross(c1,c2,ans)<<" expect 1\n";
	cout<<ans[0]<<" expect 0 2\n";
	
	cout<<cirTang(vec(0,0),c1,ans)<<" expect 1\n";
	cout<<ans[0]<<" expect 0 0\n";
	
	cout<<cirTang(vec(1,0),c1,ans)<<" expect 2\n";
	cout<<ans[0]<<' '<<ans[1]<<" expect 1 1 0 0 or 0 0 1 1\n";
	
	c1=circle(point(0,0),2);
	cout<<cirTang(vec(-4,0),c1,ans)<<" expect 2\n";
	cout<<ans[0]<<' '<<ans[1]<<" expect -1 1.73205 -1 -1.73205\n";
	
	cout<<lineInt(vec(-4,4),vec(4,-4),c1, ans)<<" expect 2\n";
	cout<<ans[0]<<' '<<ans[1]<<" expect -1.414 1.414 1.414 -1.414\n";
	
	//cout<<segInt(vec(0,0),vec(4,0),c1)<<" expect 2 0\n";
	//cout<<segInt(vec(4,0),vec(0,0),c1)<<" expect 2 0\n";
	
	c2=circle(point(0,-1),1);
	point xa[4],xb[4];
	cout<<cirTang(c1,c2,xa,xb)<<" expect 1\n";
	cout<<xa[0]<<' '<<xb[0]<<" expect 0 -2 0 -2\n";
	
	c2=circle(point(2,2),2);
	cout<<cirTang(c1,c2,xa,xb)<<" expect 2\n";
	cout<<xa[0]<<' '<<xb[0]<<' '<<xa[1]<<' '<<xb[1]<<" expect -1.414 1.414 0.586 3.414 1.414 -1.414 3.414 0.586\n";
	
	c2=circle(point(4,0),2);
	cout<<cirTang(c1,c2,xa,xb)<<" expect 3\n";
	cout<<xa[0]<<' '<<xb[0]<<' '<<xa[1]<<' '<<xb[1]<<' '<<xa[2]<<' '<<xb[2]<<
		" expect 0 2 4 2 0 -2 4 -2 2 0\n";
	
	c1=circle(point(-2,0),sqrt(2));c2=circle(point(2,0),sqrt(2));
	cout<<cirTang(c1,c2,xa,xb)<<" expect 4\n";
	cout<<xa[2]<<' '<<xb[2]<<' '<<xa[3]<<' '<<xb[3]<<" expect -1 1 1 -1 -1 -1 1 1\n";
	
	a=vec(PI*0.75);
	cout<<a<<" expect -0.707 0.707\n";
	
	c1=circle(point(0,0),point(1,2));
	cout<<c1.c<<' '<<c1.r<<" expect 0.5 1 1.118\n";

	c1=circle(point(0,2),point(0,0),point(1,1));
	cout<<c1.c<<' '<<c1.r<<" expect 0 1 1\n";
	c1=circle(point(0,2),point(1,sqrt(3)),point(-sqrt(3),-1));
	cout<<c1.c<<' '<<c1.r<<" expect 0 0 2\n";
	
	point poly[4]={{-1,0},{2,1},{1,0},{2,-1}};
	cout<<inPoly({0,0},poly,4)<<' '<<inPoly({-2,0},poly,4)<<' '<<inPoly({1,0},poly,4)<<" expect 1 0 -1\n";
	cout<<inPoly({0,-0.5},poly,4)<<' '<<inPoly({1,0.2},poly,4)<<' '<<inPoly({1.5,0.2},poly,4)<<" expect 0 1 0\n";
	cout<<inPoly({1.5,0.5},poly,4)<<' '<<polyArea(poly,4)<<" expect -1 2\n";
	
	point aa[4];
	point polyt[4]={{-1,0},{2,1},{1,0},{2,-1}};
	cout<<convex(polyt,4,aa)<<" expect 3\n";
	cout<<inConvex(aa,3,{0,0})<<" expect 1\n";
	
	cout<<mincirCover(polyt,4).c<<" expect "<<circle(poly[0],poly[1],poly[3]).c<<'\n';
	cout<<mincirCover(polyt,4).r<<" expect "<<circle(poly[0],poly[1],poly[3]).r<<'\n';
	
	cout<<poly_cirArea(poly, 4, {{0,0},1})<<" expect ???\n";

	{
	using namespace Circles;
	n=2;
	Circles::c[0]=circle(vec(0,0),1);
	Circles::c[1]=circle(vec(0,1),1);
	areaunion();
	cout<<Circles::ans[1]<<" expect 5.048156\n";
	}
}

//cdq func for minDisPoint
point tp[200010],use[200010];
db cdq(point *p,int l, int r){
	if (l==r-1) return 1e12;
	if (l==r-2) {
		if (p[l].y>p[l+1].y) swap(p[l],p[l+1]);
		return !(p[l]-p[l+1]);
	}
	int mid=l+r+1>>1;
	int uc=0; Point pmid=p[mid];
	db d=min(cdq(p,l,mid),cdq(p,mid,r));
	for (int cl=l,cr=mid,cc=l;cc<r;cc++){
		if (cr>=r || cl<mid && p[cl].y<=p[cr].y)
			tp[cc]=p[cl++];
		else
			tp[cc]=p[cr++];
		if (fabs(tp[cc].x-pmid.x)<=d+eps)
			use[uc++]=tp[cc];
	}
	inc(i,uc)
		rep(j,i+1,uc){
			if (use[j].y>use[i].y+d+eps) break;
			d=min(!(use[i]-use[j]),d);
		}
	rep(i,l,r) p[i]=tp[i];
	return d;
}
db minDisPoint(point *p, int n){
	sort(p,p+n);
	return cdq(p,0,n);
}

}

namespace Geo3D{

typedef double db;
const db PI=acos(-1);
const db eps=1e-10, inf=1e12;

bool eq(db x){return fabs(x)<eps;}
int sgn(db x){
	if (x<=-eps) return -1;
	return x>=eps;
}

#define Vec const vec &
#define Point const point &
struct vec{
	db x,y,z;
	vec():x(0),y(0){}
	vec(db x, db y, db z):x(x),y(y),z(z){}
	vec(db theta):x(cos(theta)),y(sin(theta)){}

	bool operator==(Vec v) const{return eq(x-v.x) && eq(y-v.y) && eq(z-v.z);}
	
	vec operator+(Vec v) const{return vec(x+v.x,y+v.y,z+v.z);}
	vec operator-(Vec v) const{return vec(x-v.x,y-v.y,z-v.z);}
	vec operator*(db a) const{return vec(x*a,y*a,z*a);}
	vec operator/(db a) const{return vec(x/a,y/a,z/a);}
	
	db operator|(Vec v) const{return x*v.x+y*v.y+z*v.z;} //dot
	vec operator&(Vec v) const{return vec(y*v.z-z*v.y,z*v.x-x*v.z,x*v.y-y*v.x);} //cross
	db operator!() const{return sqrt(x*x+y*y+z*z);} //len
	
	friend ostream& operator<<(ostream &o, Vec v){
		o<<v.x<<','<<v.y<<','<<v.z;
		return o;
	}
};
typedef vec point;

db angle(Vec a, Vec b){return atan2(!(a&b),a|b);}
vec cross(Point a, Point b, Point c){return b-a & c-a;}
db dot(Point a, Point b, Point c){return b-a | c-a;}

//mixtured product; 6-times directed volume
db vol6(Point a, Point b, Point c, Point d){
	return b-a&c-a|d-a;
}

//point projection on line S+tV
point lineProj(Point p, Point s, Vec v){
	return s+v*(v|p-s)/(v|v);
}
//symmetric point of P about line S+tV
point symmetric(Point p, Point s, Vec v){
	return lineProj(p,s,v)*2-p;
}
//distance of p to line S+tV
db lineDis(Point p, Point s, Vec v){
	return !(v & p-s)/!v;
}
//distance of p to segment S+tV
db segDis(Point p, Point s, Vec v){
	if (eq(!v)) return !(s-p); //single point
	vec v2=p-s,v3=p-s-v;
	if ((v|v2)<0) return !v2;
	else if ((v|v3)>0) return !v3;
	return !(v&v2)/!v;
}
//distance of seg A-B and seg C-D
db segDis(Point a, Point b, Point c, Point d){
	vec u=b-a, v=d-c;
	return min(min(segDis(c,a,u),segDis(d,a,u)),min(segDis(a,c,v),segDis(b,c,v)));
}
//point is on line
bool onLine(Point p, Point a, Point b){return eq(!(p-a&b-a));}
//point on seg [a,b]
bool onSeg(Point p, Point a, Point b){return eq(!(p-a&b-a)) && sgn(a-p|b-p)<=0;}

//rot point P by line S+tV ang rads clockwise(see from s to t)
point rot(Point p, Point s, Vec v, db ang){
	if (eq(!(v&p-s))) return p;
	point f1=v&p-s;
	point f2=f1&v;
	f1=f1/!v; 
	f2=f2/!f2*!f1;
	return p-f2+f1*sin(ang)+f2*cos(ang);
}

struct plane{
	point p;
	vec v; //normal vector
	plane(){}
	plane(Point p, Vec v):p(p),v(v){}
	plane(Point a, Point b, Point c):p(a),v(cross(a,b,c)){}
	//ax+by+cz+d=0
	plane(db a, db b, db c, db d){
		v=vec(a,b,c);
		if (sgn(a)) p=point((-d-c-b)/a,1,1);
		else if (sgn(b)) p=point(1,(-d-c-a)/b,1);
		else p=point(1,1,(-d-a-b)/c);
	}
};
//point is on plane
bool onPlane(Point p, plane f){
	return eq(dot(f.p,p,f.v));
}
//line s cross plane f
int lineInt(point s, vec v, plane f, point &ans){
	db d=v|f.v;
	if (eq(d)) return 0; //parallel
	ans=s+v/d*(f.p-s|f.v);
	return 1;
}
//porjection of point p on plane f
point planeProj(point p, plane f){
	db d=f.v|f.v;
	return p+f.v/d*(f.p-p|f.v);
}
//plane a cross plane b, get a line
int planeInt(plane a, plane b, point &s, point &v){
	v=a.v&b.v;
	if (eq(!v)) return 0; //parallel
	point t=a.v&v;
	s=a.p+t/fabs(b.v|t)*(b.v|b.p-a.p); //s is cent pos
	return 1;
}

//area of triangle on unit sphere
db angle3d_sptri(Point x, Point y, Point z){
	vec a=x&y,b=y&z,c=x&z;
	return angle(a,b)+angle(b,c)+angle(a,c)-PI;
}
//triangle projection on unit sphere
db angle3d_tri(Point x, Point y, Point z){
	db a=angle(x,y),b=angle(y,z),c=angle(x,z);
	db s=(a+b+c)/2;
	return 4*atan(sqrt(tan(s/2)*tan(s/2-a/2)*tan(s/2-b/2)*tan(s/2-c/2)));
}

struct sphere{
	point o; db r;
	sphere(){}
	sphere(Point o, db r):o(o),r(r){}
	sphere(Point a, Point b):o((a+b)/2),r(!(a-b)/2){}
	//min sphere passing point A,B,C
	//[!] a,b,c should not on same line
	sphere(Point a, Point b, Point c){ 
		vec h1=b-a,h2=c-a,h3=b&c; //three plane intersection
		vec g=vec(h1|h1,h2|h2,0)/2;   //ax+by+cz=g
		vec g1=vec(h1.x,h2.x,h3.x); //transfer
		vec g2=vec(h1.y,h2.y,h3.y);
		vec g3=vec(h1.z,h2.z,h3.z);
		db s=g1&g2|g3;              //Cramer's Rule
		o=vec(g&g2|g3,g1&g|g3,g1&g2|g)/s + a; 
		r=!(a-o);
	}
	 //[!] a,b,c,d should not collinear or coplanear
	sphere(Point a, Point b, Point c, Point d){
		vec h1=b-a,h2=c-a,h3=d-a; //three plane intersection
		vec g=vec(h1|h1,h2|h2,h3|h3)/2;   //ax+by+cz=g
		vec g1=vec(h1.x,h2.x,h3.x); //transfer
		vec g2=vec(h1.y,h2.y,h3.y);
		vec g3=vec(h1.z,h2.z,h3.z);
		db s=g1&g2|g3;              //Cramer's Rule
		o=vec(g&g2|g3,g1&g|g3,g1&g2|g)/s + a; 
		r=!(a-o);
	}
};

//convex hull 3D
namespace CH3D{

const int N=1010; //O(n^2)
point p[N];
struct face{
	int v[3]; //index on p
};
vector<face> ans;
bool vis[N][N];
void convex(int n){ //[i] as result, cross(p[v[0]],p[v[1]],p[v[2]]) towards outside of poly
	vector<face> nxt;
	//make first face not collineration; [!] point p changed
	for (int i=2;i<n;i++) if (sgn(!cross(p[0],p[1],p[i]))){swap(p[2],p[i]);break;}
	for (int i=3;i<n;i++) if (sgn(vol6(p[0],p[1],p[2],p[i]))) {swap(p[3],p[i]);break;}
	if (eq(vol6(p[0],p[1],p[2],p[3]))) return; //all on same line
	ans.push_back((face){{1,2,0}}); 
	ans.push_back((face){{2,1,0}}); //another direction. algo will select one auto.
	for (int i=3;i<n;i++){ //adding points
		nxt.clear();
		for (auto &f:ans){ //remove visable face
			bool see=sgn(vol6(p[f.v[0]],p[f.v[1]],p[f.v[2]],p[i]))>=0; //assume coplanear face visable, so previous coplanear point will be deleted
			if (!see) nxt.push_back(f);
			for (int k=0;k<3;k++) vis[f.v[k]][f.v[(k+1)%3]]=see; //label edges
		}
		if (nxt.size()==ans.size()) continue;
		for (auto &f:ans)
			for (int k=0;k<3;k++){
				int a=f.v[k],b=f.v[(k+1)%3];
				if (!vis[b][a] && vis[a][b])
					nxt.push_back((face){{a,b,i}}),vis[a][b]=1;
			}
		ans.swap(nxt);//update to ans
	}
}

//--polyhedron--

db volume(){ //[!] the input face should towards same side
	db sum=0;
	for (auto &f:ans)
		sum+=vol6(p[0],p[f.v[0]],p[f.v[1]],p[f.v[2]]);
	return fabs(sum/6);
}
point barycenter(){ //[!] the input face should towards same side
	point s(0,0,0);
	db sum=0;
	for (auto &f:ans){
		db v=vol6(p[0],p[f.v[0]],p[f.v[1]],p[f.v[2]]);
		sum+=v;
		s=s+(p[0]+p[f.v[0]]+p[f.v[1]]+p[f.v[2]])/4*v;
	}
	return s/sum;
}

//point s is in or on polygon
//return  1(in), 0(out), -1(on border)
/*
int inPoly(Point s){ //[!] the input face should towards outside
	auto rdi=[](){return rand()%10+1;};
start:
	vec v=vec(rdi(),rdi(),rdi());
	int w=0;
	for (auto &f: ans){
		point cp;
		int ret=lineInt(s,v,plane(p[f.v[0]],p[f.v[1]],p[f.v[2]]),cp);
		if (!ret || sgn(cp-s|v)<0) continue;
		int s1=sgn(!cross(p[f.v[0]],p[f.v[1]],cp)); //TODO : bug here
		int s2=sgn(!cross(p[f.v[0]],cp,p[f.v[2]])); //how to test point in 3d triangle?
		int s3=sgn(!cross(p[f.v[1]],cp,p[f.v[0]]));
		int s4=sgn(!cross(p[f.v[1]],p[f.v[2]],cp));
		if (s1==1 && s2==1 && s3==1 && s4==1){
			if (sgn(cp-s|v)==0) return -1;
			w+=sgn(vol6(p[f.v[0]],p[f.v[1]],p[f.v[2]],s));
		}
		if ((s1||s2||s3||s4) == 0){
			if (sgn(cp-s|v)==0) return -1;
			goto start;
		}
	}
	return w!=0;
}
*/
//another impl, return sum of angle3d
db inPoly2(Point s){
	db w=0;
	for (auto &f: ans) 
		w+=sgn(vol6(s,p[f.v[0]],p[f.v[1]],p[f.v[2]]))*angle3d_tri(p[f.v[0]]-s,p[f.v[1]]-s,p[f.v[2]]-s);
	return w;
}

}// namespace CH3D

void test(){
	point p1(0,0,0),p2(1,0,0),p3(0,1,0),p4(0,0,1);
	printf("%f expected 1\n",vol6(p1,p2,p3,p4));
	cout<<rot(p3,p1,p2,PI/4)<<" expect 0,0.707,0.707\n";
	cout<<rot(p3,p1,p2,-PI/4)<<" expect 0,0.707,-0.707\n";
	cout<<rot(p3,p1,p2,PI/2)<<" expect 0,0,1\n";
	cout<<rot(p3,p1,p2,PI)<<" expect 0,-1,0\n";
	plane f(p1,vec(1,1,1));
	vec ans;
	lineInt(p4,p1-p3-p2,f,ans);
	cout<<ans<<" expect -0.5,-0.5,1\n";
	cout<<planeProj(p4,f)<<" expect -0.3333,-0.3333,0.6667\n";
	point pv;
	planeInt(f,plane(p2,p2),ans,pv);
	cout<<ans<<' '<<pv<<" expect 1,0.5,0.5 0,k,-k\n";
	sphere ball({1,1,0},p2,p3);
	cout<<ball.o<<' '<<ball.r<<" expect 0.5,0.5,0 0.707\n";
	ball=sphere({1,1,1},{1,0,1},{0,1,1},{1,1,0});
	cout<<ball.o<<' '<<ball.r<<" expect 0.5,0.5,0.5 0.866\n";
	
	{
		using namespace CH3D;
		inc(i,2) inc(j,2) inc(k,2) p[i*4+j*2+k]=vec(i,j,k);
		convex(8);
		for (auto &f:CH3D::ans)
			cout<<p[f.v[0]]<<' '<<p[f.v[1]]<<' '<<p[f.v[2]]<<'\n';

		cout<<volume()<<" expect 1\n";
		cout<<barycenter()<<" expect 0.5,0.5,0.5\n";

		//cout<<inPoly({0.5,0.5,0.5})<<" expect 1\n";
		//cout<<inPoly({0,0.5,0.5})<<" expect -1\n";
		//cout<<inPoly({1.5,0.5,0.5})<<" expect 0\n";
		
		cout<<inPoly2({0.5,0.5,0.5})<<" expect 4*pi\n";
	}
}

} //end namespace Geo3D

//-------------------MISC------------------------
namespace DateTime{

int gettime(int h, int m, int s){
	return h*3600+m*60+s;
}

bool isleapyear(int y){
	//if (y<0) return isleapyear(-y-1);
	//if (y%3200==0) return y%172800==0; 
	return y%4==0 && y%100 || y%400==0;
}

int mm[13]={0,31,28,31,30,31,30,31,31,30,31,30,31};
//get day diff from 0000/01/01 (BC 0001y), but require y>0
int getday(int y, int m, int d){
	if (m<3) y--,m+=12;
	return (d+30*m+3*(m+1)/5+y*365+y/4-y/100+y/400)-33;
}
//inverse function of getday()
void getdate(int d0, int &y, int &m, int &d){
	int y1=(d0)/146097;
	int y2=(d0-1-y1*146097)/36524;
	int y3=(d0-1-y1*146097-y2*36524)/1461;
	y=y1*400+y2*100+y3*4+(d0-1-y1*146097-y2*36524-y3*1461)/365;
	d=d0-y*365-(y-1)/4+(y-1)/100-(y-1)/400; m=1;
	if (y%4==0&&y%100||y%400==0) mm[2]++;
	while (d-mm[m]>0) d-=mm[m++];
	if (y%4==0&&y%100||y%400==0) mm[2]--;
}

//get week by date,1 for Monday
//[!] Because 1582/10/05~1582/10/14 is not existed
// the formula is correct after that day
int getweek(int y, int m, int d){
	if (m<3) y--,m+=12;
	return (d+2*m+3*(m+1)/5+y+y/4-y/100+y/400+1)%7;
}

}

namespace scannerLine{
const int maxn=100010;
struct Line{
	double l,r,h; int c;
	bool operator<(const Line &v) const{
		return h<v.h;
	}
}li[maxn];
int lic;

#define lc u+u+1
#define rc u+u+2
double len[maxn<<2]; int cnt[maxn<<2];
double x[maxn<<1]; int xc;

void calc(int u, int l, int r){
	if (cnt[u])
		len[u]=x[r]-x[l];
	else if (l==r-1)
		len[u]=0;
	else
		len[u]=len[lc]+len[rc];
}
void add(int u, int l, int r, int cl, int cr, int c){
	if (cl<=l && cr>=r){
		cnt[u]+=c;
		calc(u,l,r);
		return;
	}
	int mid=l+r>>1;
	if (cl<mid) add(lc,l,mid,cl,cr,c);
	if (cr>mid) add(rc,mid,r,cl,cr,c);
	calc(u,l,r);
}

double x0[maxn],y0[maxn],x1[maxn],y1[maxn];
void rectInt(int n){
	xc=lic=0;
	memset(len,0,sizeof(len));
	memset(cnt,0,sizeof(cnt));
	for (int i=0;i<n;i++){
		double x1,y1,x2,y2;
		scanf("%lf%lf%lf%lf",&x1,&y1,&x2,&y2);
		x[xc++]=x1; x[xc++]=x2;
		li[lic++]=(Line){x1,x2,y1,1};
		li[lic++]=(Line){x1,x2,y2,-1};
	}
	sort(li,li+lic);
	sort(x,x+xc);
	double ans=0,last=0;
	for (int i=0;i<lic;i++){
		int l=lower_bound(x,x+xc,li[i].l)-x;
		int r=lower_bound(x,x+xc,li[i].r)-x;
		add(0,0,xc,l,r,li[i].c);
		ans+=len[0]*(li[i+1].h-li[i].h);
		//sum of lenth on sx
		//ans+=fabs(len[0]-last); last=len[0];
	}
	printf("%.2f\n",ans);
}
#undef lc
#undef rc
}

//O(nlogn)
namespace HFMT{
const int maxn=30;

//finally tn is the root of tree
int a[maxn],ch[maxn<<1][2],n,tn; //idx from 1..n
int sum[maxn<<1]; //not very necessary
//input  a[1..n]:frequency of each character, n: |character|
//result ch[maxn<<1][2], the path to leaf node is the encoding of input char
void build(){
	priority_queue<PII> qu;
	for (int i=1;i<=n;i++) qu.emplace(-a[i],i),sum[i]=a[i];
	tn=n;
	while (qu.size()>1){
		int x1=qu.top().first, p1=qu.top().second;
		qu.pop();
		int x2=qu.top().first, p2=qu.top().second;
		qu.pop();
		ch[++tn][0]=p1; ch[tn][1]=p2;
		sum[tn]=-x1-x2;
		qu.emplace(x1+x2, tn);
	}
}
int len[maxn];
//dfs: debug, and label lenth of encode after
void dfs(int u=tn, int deep=0){
	if (!u) return;
	if (u<=n) len[u]=deep;
	dfs(ch[u][1],deep+1);
	//for (int i=0;i<deep;i++) printf("  "); printf("%d\n",sum[u]);
	dfs(ch[u][0],deep+1);
}
};

namespace SquareTransform{
const int N=100;
typedef int Arr[N][N];int n;
void cp(Arr to,Arr from){inc(i,n)inc(j,n) to[i][j]=from[i][j];}
Arr _t;
//clockwise 90 deg
void rot(Arr a){
	inc(i,n)inc(j,n) _t[j][n-i-1]=a[i][j];
	cp(a,_t);
}
//LR flip
void flip(Arr a){
	inc(i,n)inc(j,n) _t[i][n-j-1]=a[i][j];
	cp(a,_t);
}
bool same(Arr a, Arr b){
	inc(i,n) inc(j,n) if (a[i][j]!=b[i][j]) return 0;
	return 1;
}
}

namespace Expr{
//Easy experission, calc +-*/^()
#define CP cin.peek()
#define CG cin.get()
#define CS while (CP==' ') CG;

int S();
int V(){CS
	int ans=0;
	if (CP=='('){
		CG;
		ans=S();
		CS;CG;
	}
	else cin>>ans;
	return ans;
}
int U(){
	int ans=V(); CS;
	while (CP=='^'){
		CG;
		int v=V(),d=ans;
		if (v==0) ans=1;
		for (int i=1;i<v;i++)
			ans*=d;
	}
	return ans;
}
int T(){
	int ans=U(); CS;
	while (CP=='*' || CP=='/'){
		if (CG=='*') ans*=U();
		else ans/=U();
	}
	return ans;
}
int S(){
	int ans=0; CS;
	if (CP=='-'){
		CG; ans=-T();
	}
	else ans=T();
	CS;
	while (CP=='+' || CP=='-'){
		if (CG=='+') ans+=T();
		else ans-=T();
	}
	return ans;
}

#undef CG
#undef CP
#undef CS
}

int main(){
	return 0;
}

