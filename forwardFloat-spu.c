#include <spu_intrinsics.h>
#include <spu_mfcio.h>
#include <simdmath.h>
//#include <simdmath/fmaf4.h>
#include <simdmath/sqrtf4.h>
#include <simdmath/rsqrtf4.h>
#include <simdmath/fmaxf4.h>
//#include <simdmath/fabsf4.h>
#include <simdmath/cosf4.h>
#include <simdmath/sinf4.h>
#include <simdmath/powf4.h>
//#include <simdmath/fminf4.h>
#include <math.h>   // smallpt, a Path Tracer by Kevin Beason, 2008
#include <stdlib.h> // Make : g++ -O3 -fopenmp smallpt.cpp -o smallpt
#include <stdio.h>  //        Remove "-fopenmp" for g++ version < 4.2
#include "spu-rng.h"

volatile unsigned int status_termination __attribute__((aligned(16)));

volatile unsigned int status[128/sizeof(unsigned int)] __attribute__((aligned(128)));

typedef struct _contexto {
	unsigned long long retorno, status_termination, status; //endereço de retorno
	unsigned long long seed;
	int samps, w, h, dummy;
} contexto;

volatile contexto ctxt __attribute__((aligned(16)));

typedef vector float Vec;

Vec vzero={0.f,0.f,0.f,0.f};
Vec vones={1.f,1.f,1.f,1.f};
Vec vtwos={2.f,2.f,2.f,2.f};
Vec vdotone={.1f,.1f,.1f,.1f};
Vec vhalf={.5f,.5f,.5f,.5f};
Vec vquarter={.25f,.25f,.25f,.25f};
Vec vfactor={140.f,140.f,140.f,140.f};
Vec veps={1e-4, 1e-4, 1e-4, 1e-4};
Vec vinf={1e20,1e20,1e20,1e20};

vector unsigned char vpat1={0x04,0x05,0x06,0x07,
							0x00,0x01,0x02,0x03,
							0x0c,0x0d,0x0e,0x0f,
							0x08,0x09,0x0a,0x0b};
vector unsigned char vpat2={0x08,0x09,0x0a,0x0b,
							0x0c,0x0d,0x0e,0x0f,
							0x00,0x01,0x02,0x03,
							0x04,0x05,0x06,0x07};
vector unsigned char vpat3={0x04,0x05,0x06,0x07,
							0x08,0x09,0x0a,0x0b,
							0x00,0x01,0x02,0x03,
							0x00,0x01,0x02,0x03,};
vector unsigned char vpat4={0x04,0x05,0x06,0x07,
							0x08,0x09,0x0a,0x0b,
							0x0c,0x0d,0x0e,0x0f,
							0x00,0x01,0x02,0x03,};
Vec VLoad(float x, float y, float z) {
	vector float result=(vector float){x, y, z, 0.0};
	return  result;
}
Vec VPromote(float x) {
	return VLoad(x,x,x);
}

float VGetX(Vec v)  {
	return spu_extract(v, 0);
}
float VGetY(Vec v)  {
	return spu_extract(v, 1);
}
float VGetZ(Vec v)  {
	return spu_extract(v, 2);
}

int maxdepth;

const Vec vec(float x, float y, float z) {
	return VLoad(x,y,z);
}

const Vec vecmult(const Vec a, const Vec b) {
	return spu_mul(a,b);
}

const Vec vecdbl(const Vec a, float b) {
	return spu_mul(a, VLoad(b, b, b));
}

const Vec vecsum(const Vec a, const Vec b) {
	return spu_add(a,b);
}

const Vec vecsub(const Vec a, const Vec b) {
	return spu_sub(a,b);
}

const Vec vec_cross(const Vec a, const Vec b) {
	Vec x1=VLoad(VGetY(a),VGetZ(a),VGetX(a));
	Vec x2=VLoad(VGetZ(b),VGetX(b),VGetY(b));
	Vec x3=VLoad(VGetZ(a),VGetX(a),VGetY(a));
	Vec x4=VLoad(VGetY(b),VGetZ(b),VGetX(b));
	Vec result=vecmult(x1,x2);
	return spu_nmsub(x3,x4,result);
}

float dot(const Vec a, const Vec b) {
	Vec d=vecmult(a,b);
	float x=VGetX(d);
	float y=VGetY(d);
	float z=VGetZ(d);
	return x+y+z;
}

const Vec vecnorm(const Vec a) {
	vector float rroot=_rsqrtf4(spu_splats(dot(a,a)));
	return vecmult(a,rroot);
}

typedef struct _Ray {
	Vec o, d;
} Ray;

Ray criaRay(Vec o, Vec d) {
	Ray r;
	r.o=o;
	r.d=d;
	return r;

}

typedef enum _Refl_t { DIFF, SPEC, REFR } Refl_t;  // material types, used in radiance()

typedef struct _Sphere {
	float rad;       // radius
	Vec p, e, c;      // position, emission, color
	Refl_t refl;      // reflection type (DIFFuse, SPECular, REFRactive)
} Sphere;

const Sphere sphere(float rad_, Vec p_, Vec e_, Vec c_, Refl_t refl_) {
	Sphere s;
	s.rad=rad_;
	s.p=p_;
	s.e=e_;
	s.c=c_;
	s.refl=refl_;
	return s;
}
#define NUM_SPHERES (12)
Sphere spheres[NUM_SPHERES];
void initSpheres() {
	spheres[0]=sphere(1e4, vec( 1e4+1,40.8,81.6), vec(0,0,0),vec(.75,.25,.25),DIFF);//Left
	spheres[1]=sphere(1e4, vec(-1e4+99,40.8,81.6),vec(0,0,0),vec(.25,.25,.75),DIFF);//Rght
	spheres[2]=sphere(1e4, vec(50,40.8, 1e4),     vec(0,0,0),vec(.75,.75,.75),DIFF);//Back
	spheres[3]=sphere(1e4, vec(50,40.8,-1e4+170), vec(0,0,0),vec(0,0,0),           DIFF);//Frnt
	spheres[4]=sphere(1e4, vec(50, 1e4, 81.6),    vec(0,0,0),vec(.75,.75,.75),DIFF);//Botm
	spheres[5]=sphere(1e4, vec(50,-1e4+81.6,81.6),vec(0,0,0),vec(.75,.75,.75),DIFF);//Top
	spheres[6]=sphere(16.5,vec(27,16.5,47),       vec(0,0,0),vecdbl(vec(1,1,1),.999), SPEC);//Mirr
	spheres[7]=sphere(16.5,vec(73,16.5,78),       vec(0,0,0),vecdbl(vec(1,1,1),.999), REFR);//Glas
	spheres[8]=sphere(8, vec(50,81.6-16.5,81.6),vec(12,12,12),  vec(0,0,0), DIFF); //Lite
	spheres[9]=sphere(6, vec(50,681.6-.27,81.6),vec(12,12,12),  vec(0,0,0), DIFF); //Dummy
	spheres[10]=sphere(6, vec(50,681.6-.27,81.6),vec(12,12,12),  vec(0,0,0), DIFF); //Dummy
	spheres[11]=sphere(6, vec(50,681.6-.27,81.6),vec(12,12,12),  vec(0,0,0), DIFF); //Dummy
}

float clamp(float x) {
	return x<0 ? 0 : x>1 ? 1 : x;
}
int toInt(float x) {
	return (long)(pow(clamp(x),1/2.2)*255+.5);
}
/*vector float vecClamp(vector float x){
	return spu_sel(spu_sel(x,vones,spu_cmpgt(x,vones)),vzero,spu_cmpgt(vzero,x));
}
vector int vecToInt(vector float x){
	float aux=1/2.2;
	return spu_convts(spu_add(spu_mul(_powf4(x,spu_splats(aux)),spu_splats(255.f)),vhalf),0);
}*/
typedef struct Vecsoa_ {
	vector float x;
	vector float y;
	vector float z;
} Vecsoa;
Vecsoa vecsoa(Vec a, Vec b, Vec c, Vec d){
	Vecsoa r = {{VGetX(a),VGetX(b), VGetX(c), VGetX(d)},
				{VGetY(a),VGetY(b), VGetY(c), VGetY(d)},
				{VGetZ(a),VGetZ(b), VGetZ(c), VGetZ(d)}};
	return r;
}
typedef struct Spheresoa_ {
	vector float rad;
	Vecsoa p, e, c;
	vector signed int refl;
	vector signed int id;
} Spheresoa;
#define NUM_SPHERESOA (3)
Spheresoa sphsoa[]={
		{
	/*spheres[0]=sphere(1e4, vec( 1e4+1,40.8,81.6), vec(0,0,0),vec(.75,.25,.25),DIFF);//Left
	spheres[1]=sphere(1e4, vec(-1e4+99,40.8,81.6),vec(0,0,0),vec(.25,.25,.75),DIFF);//Rght
	spheres[2]=sphere(1e4, vec(50,40.8, 1e4),     vec(0,0,0),vec(.75,.75,.75),DIFF);//Back
	spheres[3]=sphere(1e4, vec(50,40.8,-1e4+170), vec(0,0,0),vec(0,0,0),           DIFF);//Frnt*/
			{1e4, 1e4, 1e4, 1e4},
			{
				{1e4+1, 1e4+99, 50, 50},
				{40.8, 40.8, 40.8, 40.8},
				{81.6, 81.6, 1e4, -1e4+170}
			},
			{
				{0,0,0,0},
				{0,0,0,0},
				{0,0,0,0}
			},
			{
				{.75, .25, .75, 0},
				{.25, .25, .75, 0},
				{.25, .75, .75, 0}
			},
			{DIFF, DIFF, DIFF, DIFF},
			{0,1,2,3}
		},
		/*spheres[4]=sphere(1e4, vec(50, 1e4, 81.6),    vec(0,0,0),vec(.75,.75,.75),DIFF);//Botm
	spheres[5]=sphere(1e4, vec(50,-1e4+81.6,81.6),vec(0,0,0),vec(.75,.75,.75),DIFF);//Top
	spheres[6]=sphere(16.5,vec(27,16.5,47),       vec(0,0,0),vecdbl(vec(1,1,1),.999), SPEC);//Mirr
	spheres[7]=sphere(16.5,vec(73,16.5,78),       vec(0,0,0),vecdbl(vec(1,1,1),.999), REFR);//Glas*/
		{
			{1e4, 1e4, 16.5, 16.5},
			{
				{50,	50,		27,		73},
				{1e4,	1e4+81,	16.5,	16.5},
				{81.6,	81.6,	47,		78}
			},
			{
				{0,0,0,0},
				{0,0,0,0},
				{0,0,0,0}
			},
			{
				{.75f, .75f, .999f, .999f},
				{.75f, .75f, .999f, .999f},
				{.75f, .75f, .999f, .999f}
			},
			{DIFF, DIFF, SPEC, REFR},
			{4,5,6,7}
		},
		/*spheres[8]=sphere(8, vec(50,81.6-16.5,81.6),vec(12,12,12),  vec(0,0,0), DIFF); //Lite
	spheres[9]=sphere(6, vec(50,681.6-.27,81.6),vec(12,12,12),  vec(0,0,0), DIFF); //Dummy
	spheres[10]=sphere(6, vec(50,681.6-.27,81.6),vec(12,12,12),  vec(0,0,0), DIFF); //Dummy
	spheres[11]=sphere(6, vec(50,681.6-.27,81.6),vec(12,12,12),  vec(0,0,0), DIFF); //Dummy*/
		{
			{8, 6, 6, 6},
			{
				{50,		50,			50,			50},
				{81.6-16.5,	681.6-.27,	681.6-.27,	681.6-.27},
				{81.6,		81.6,		81.6,		81.6}
			},
			{
				{12,12,12,12},
				{12,12,12,12},
				{12,12,12,12}
			},
			{
				{0,0,0,0},
				{0,0,0,0},
				{0,0,0,0}
			},
			{DIFF, DIFF, DIFF, DIFF},
			{8,9,10,11}
		}
	};
typedef struct _Raysoa {
	Vecsoa o, d;
} Raysoa;
Raysoa criaRaysoa(Ray r1, Ray r2, Ray r3, Ray r4){
	Raysoa r;
	r.o=vecsoa(r1.o, r2.o, r3.o, r4.o);
	r.d=vecsoa(r1.d, r2.d, r3.d, r4.d);
	return r;
}
Raysoa criaRaysoa2(Vecsoa o, Vecsoa d){
	Raysoa r;
	r.o=o;
	r.d=d;
	return r;
}
Spheresoa spheresoa(Sphere s1, Sphere s2, Sphere s3, Sphere s4){
	Spheresoa sphsoa;
	sphsoa.p=vecsoa(s1.p,s2.p,s3.p,s4.p);
	sphsoa.e=vecsoa(s1.e,s2.e,s3.e,s4.e);
	sphsoa.c=vecsoa(s1.c,s2.c,s3.c,s4.c);
	
	sphsoa.rad=spu_insert(s1.rad, sphsoa.rad, 0);
	sphsoa.rad=spu_insert(s2.rad, sphsoa.rad, 1);
	sphsoa.rad=spu_insert(s3.rad, sphsoa.rad, 2);
	sphsoa.rad=spu_insert(s4.rad, sphsoa.rad, 3);
	
	sphsoa.refl=spu_insert(s1.refl,sphsoa.refl,0);
	sphsoa.refl=spu_insert(s2.refl,sphsoa.refl,1);
	sphsoa.refl=spu_insert(s3.refl,sphsoa.refl,2);
	sphsoa.refl=spu_insert(s4.refl,sphsoa.refl,3);
	return sphsoa;
}
void loadSpheresSOA(){
	int i, j=0;
	for(i=0;i<NUM_SPHERESOA;i++){
		sphsoa[i].p=vecsoa(spheres[j].p,spheres[j+1].p,spheres[j+2].p,spheres[j+3].p);
		sphsoa[i].e=vecsoa(spheres[j].e,spheres[j+1].e,spheres[j+2].e,spheres[j+3].e);
		sphsoa[i].c=vecsoa(spheres[j].c,spheres[j+1].c,spheres[j+2].c,spheres[j+3].c);
		
		sphsoa[i].rad=spu_insert(spheres[j].rad, sphsoa[i].rad, 0);
		sphsoa[i].rad=spu_insert(spheres[j+1].rad, sphsoa[i].rad, 1);
		sphsoa[i].rad=spu_insert(spheres[j+2].rad, sphsoa[i].rad, 2);
		sphsoa[i].rad=spu_insert(spheres[j+3].rad, sphsoa[i].rad, 3);
		
		sphsoa[i].refl=spu_insert(spheres[j].refl,sphsoa[i].refl,0);
		sphsoa[i].refl=spu_insert(spheres[j+1].refl,sphsoa[i].refl,1);
		sphsoa[i].refl=spu_insert(spheres[j+2].refl,sphsoa[i].refl,2);
		sphsoa[i].refl=spu_insert(spheres[j+3].refl,sphsoa[i].refl,3);
		
		sphsoa[i].id=spu_insert(j++,sphsoa[i].id,0);
		sphsoa[i].id=spu_insert(j++,sphsoa[i].id,1);
		sphsoa[i].id=spu_insert(j++,sphsoa[i].id,2);
		sphsoa[i].id=spu_insert(j++,sphsoa[i].id,3);
	}
}
Vecsoa SOAPromote(Vec a){
	Vecsoa r={	{VGetX(a),VGetX(a),VGetX(a),VGetX(a)},
				{VGetY(a),VGetY(a),VGetY(a),VGetY(a)},
				{VGetZ(a),VGetZ(a),VGetZ(a),VGetZ(a)}};
	return r;
}
Vecsoa soa_vecsub(Vecsoa a, Vecsoa b){
	Vecsoa tmp;
	tmp.x=spu_sub(a.x,b.x);
	tmp.y=spu_sub(a.y,b.y);
	tmp.z=spu_sub(a.z,b.z);
	return tmp;
}
Vecsoa soa_vecsum(Vecsoa a, Vecsoa b){
	Vecsoa tmp;
	tmp.x=spu_add(a.x, b.x);
	tmp.y=spu_add(a.y, b.y);
	tmp.z=spu_add(a.z, b.z);
	return tmp;
}
Vecsoa soa_vecmul(Vecsoa a, Vecsoa b){
	Vecsoa tmp;
	tmp.x=spu_mul(a.x, b.x);
	tmp.y=spu_mul(a.y, b.y);
	tmp.z=spu_mul(a.z, b.z);
	return tmp;
}
Vecsoa soa_vecdbl(Vecsoa a, Vec b){
	Vecsoa tmp;
	tmp.x=spu_mul(a.x,b);
	tmp.y=spu_mul(a.y,b);
	tmp.z=spu_mul(a.z,b);
	return tmp;
	//return soa_vecmul(a,tmp);
}
Vec soa_dot(Vecsoa a, Vecsoa b){
	vector float x = spu_madd(a.x,b.x,vzero);
	x=spu_madd(a.y,b.y,x);
	x=spu_madd(a.z,b.z,x);
	return x;
}
Vecsoa soa_vecnorm(Vecsoa a) {
	vector float length=_rsqrtf4(soa_dot(a,a));
	Vecsoa rroot;
	rroot.x=spu_mul(a.x,length);
	rroot.y=spu_mul(a.y,length);
	rroot.z=spu_mul(a.z,length);
	return rroot;
}
Vecsoa soa_veccross(Vecsoa a, Vecsoa b){
	Vecsoa r;
	r.x=vecsub(vecmult(a.y,b.z),vecmult(a.z,b.y));
	r.y=vecsub(vecmult(a.z,b.x),vecmult(a.x,b.z));
	r.z=vecsub(vecmult(a.x,b.y),vecmult(a.y,b.x));
	return r;
}
Vecsoa soa_select(Vecsoa a, Vecsoa b, vector unsigned int mask){
	Vecsoa r;
	r.x=spu_sel(a.x,b.x,mask);
	r.y=spu_sel(a.y,b.y,mask);
	r.z=spu_sel(a.z,b.z,mask);
	return r;
}
Vec soa_intersect(Vecsoa v, Raysoa *r, vector float raio, vector signed int *id){
	Vecsoa v1, v2, v3;
	v1.x=spu_shuffle(v.x,vzero,vpat4);
	v1.y=spu_shuffle(v.y,vzero,vpat4);
	v1.z=spu_shuffle(v.z,vzero,vpat4);
	v2.x=spu_shuffle(v1.x,vzero,vpat4);
	v2.y=spu_shuffle(v1.y,vzero,vpat4);
	v2.z=spu_shuffle(v1.z,vzero,vpat4);
	v3.x=spu_shuffle(v2.x,vzero,vpat4);
	v3.y=spu_shuffle(v2.y,vzero,vpat4);
	v3.z=spu_shuffle(v2.z,vzero,vpat4);
	vector float 	raio1=spu_shuffle(raio,vzero,vpat4),
					raio2=spu_shuffle(raio1,vzero,vpat4),
					raio3=spu_shuffle(raio2,vzero,vpat4);
	Vecsoa op1 = soa_vecsub(v, r->o),
			op2 = soa_vecsub(v1, r->o),
			op3 = soa_vecsub(v2, r->o),
			op4 = soa_vecsub(v3, r->o);
	vector float	t11,t12,t21,t22,t31,t32,t41,t42,
	        b1=soa_dot(op1, r->d),
	        b2=soa_dot(op2, r->d),
	        b3=soa_dot(op3, r->d),
	        b4=soa_dot(op4, r->d),
	        det1=spu_add(spu_sub(spu_mul(b1,b1),soa_dot(op1, op1)),spu_mul(raio,raio)),
	        det2=spu_add(spu_sub(spu_mul(b2,b2),soa_dot(op2, op2)),spu_mul(raio1,raio1)),
	        det3=spu_add(spu_sub(spu_mul(b3,b3),soa_dot(op3, op3)),spu_mul(raio2,raio2)),
	        det4=spu_add(spu_sub(spu_mul(b4,b4),soa_dot(op4, op4)),spu_mul(raio3,raio3)),
			rdet1,rdet2,rdet3,rdet4;
	rdet1=_sqrtf4(det1);
	rdet2=_sqrtf4(det2);
	rdet3=_sqrtf4(det3);
	rdet4=_sqrtf4(det4);
	t11=spu_sub(b1,rdet1);
	t12=spu_add(b1,rdet1);
	t21=spu_sub(b2,rdet2);
	t22=spu_add(b2,rdet2);
	t31=spu_sub(b3,rdet3);
	t32=spu_add(b3,rdet3);
	t41=spu_sub(b4,rdet4);
	t42=spu_add(b4,rdet4);
	vector unsigned int mask0=spu_cmpgt(vzero,det1);
	vector unsigned int mask1=spu_cmpgt(t11,veps);
	vector unsigned int mask2=spu_cmpgt(t12,veps);
	
	vector float result1 = spu_sel(vzero,t12,mask2);
	result1 = spu_sel(result1, t11, mask1);
	result1 = spu_sel(result1,vzero,mask0);
	
	mask0=spu_cmpgt(vzero,det2);
	mask1=spu_cmpgt(t21,veps);
	mask2=spu_cmpgt(t22,veps);
	
	vector float result2 = spu_sel(vzero,t22,mask2);
	result2 = spu_sel(result2, t21, mask1);
	result2 = spu_sel(result2,vzero,mask0);
	
	mask0=spu_cmpgt(vzero,det3);
	mask1=spu_cmpgt(t31,veps);
	mask2=spu_cmpgt(t32,veps);
	
	vector float result3 = spu_sel(vzero,t32,mask2);
	result3 = spu_sel(result3, t31, mask1);
	result3 = spu_sel(result3,vzero,mask0);
	
	mask0=spu_cmpgt(vzero,det4);
	mask1=spu_cmpgt(t41,veps);
	mask2=spu_cmpgt(t42,veps);
	
	vector float result4 = spu_sel(vzero,t42,mask2);
	result4 = spu_sel(result4, t41, mask1);
	result4 = spu_sel(result4,vzero,mask0);
	vector signed int id1=spu_shuffle(*id,*id,vpat4);
	vector signed int id2=spu_shuffle(id1,*id,vpat4);
	vector signed int id3=spu_shuffle(id2,*id,vpat4);
	
	vector unsigned int mask6=spu_cmpabsgt(result1,vzero);
	vector unsigned int mask5=spu_cmpgt(vinf,result1);
	vector signed int negone=spu_splats(-1);
	
	vector float distancia1 = spu_sel(vinf,result1,spu_and(mask5,mask6));
	vector signed int rid1=spu_sel(negone,*id,spu_and(mask5,mask6));
	
	vector unsigned int mask8=spu_cmpabsgt(result2,vzero);
	vector unsigned int mask7=spu_cmpgt(distancia1,result2);
	vector float distancia2 = spu_sel(distancia1, result2, spu_and(mask7,mask8));
	vector signed int rid2=spu_sel(rid1,id1,spu_and(mask7,mask8));
	
	vector unsigned int mask9=spu_cmpabsgt(result3,vzero);
	vector unsigned int mask10=spu_cmpgt(distancia2,result3);
	vector float distancia3=spu_sel(distancia2,result3,spu_and(mask9,mask10));
	vector signed int rid3=spu_sel(rid2,id2,spu_and(mask9,mask10));
	
	vector unsigned int mask11=spu_cmpabsgt(result4,vzero);
	vector unsigned int mask12=spu_cmpgt(distancia3,result4);
	vector float distancia = spu_sel(distancia3, result4, spu_and(mask11,mask12));
	*id=spu_sel(rid3,id3,spu_and(mask11,mask12));
	
	return distancia;
}
vector unsigned int vecsoa_intersect(Raysoa *r, vector float *t, vector signed int *id ){
	vector float d, vt={1e20,1e20,1e20,1e20};
	vector signed int aux,vid=spu_splats(0);
	int i;
	for(i=0; i<NUM_SPHERESOA; i++){
		aux=sphsoa[i].id;
		d=soa_intersect(sphsoa[i].p,r,sphsoa[i].rad,&aux);
		vector unsigned int mask1=spu_cmpgt(vt,d);
		vt=spu_sel(vt,d,mask1);
		vid=spu_sel(vid,aux,mask1);
	}
	*t=vt;
	*id=spu_sel(spu_splats(0),vid,(vector unsigned int)spu_cmpgt(vinf,vt));
	vector unsigned int retorno=spu_cmpgt(vinf,vt);
	return retorno;
}
/* Sample the light sphere uniformly over the solid angle it subtends. */
Vecsoa direct_light(Vecsoa x, Vecsoa nl, Vecsoa throughput,
                    SpuRng *rng, vector unsigned int active) {
	Vecsoa center=SOAPromote(spheres[8].p);
	Vecsoa to_light=soa_vecsub(center,x);
	Vec dist2=soa_dot(to_light,to_light);
	Vec cos_max=_sqrtf4(_fmaxf4(vzero,vones-spu_splats(64.f)*spu_re(dist2)));
	Vec cos_theta=vones-spu_rng_uniform(rng)*(vones-cos_max);
	Vec sin_theta=_sqrtf4(_fmaxf4(vzero,vones-cos_theta*cos_theta));
	Vec phi=spu_rng_uniform(rng)*spu_splats((float)(2*M_PI));
	Vecsoa w=soa_vecnorm(to_light);
	Vecsoa unit_x=SOAPromote(vec(1,0,0));
	Vecsoa unit_y=SOAPromote(vec(0,1,0));
	vector unsigned int use_y=spu_cmpabsgt(w.x,vdotone);
	Vecsoa u=soa_vecnorm(soa_veccross(soa_select(unit_x,unit_y,use_y),w));
	Vecsoa v=soa_veccross(w,u);
	Vecsoa direction=soa_vecnorm(soa_vecsum(
		soa_vecsum(soa_vecdbl(u,_cosf4(phi)*sin_theta),
		           soa_vecdbl(v,_sinf4(phi)*sin_theta)),
		soa_vecdbl(w,cos_theta)));
	Vec cosine=_fmaxf4(vzero,soa_dot(nl,direction));
	Raysoa shadow=criaRaysoa2(soa_vecsum(x,soa_vecdbl(nl,spu_splats(1e-3f))),direction);
	Vec distance;
	vector signed int hit_id;
	vecsoa_intersect(&shadow,&distance,&hit_id);
	vector unsigned int visible=spu_and(active,spu_cmpeq(hit_id,spu_splats(8)));
	visible=spu_and(visible,spu_cmpgt(cosine,vzero));
	Vec weight=spu_splats(2.f)*(vones-cos_max)*cosine;
	Vecsoa contribution=soa_vecdbl(
		soa_vecmul(throughput,SOAPromote(spheres[8].e)),weight);
	return soa_select(SOAPromote(vzero),contribution,visible);
}
Vec soa_extract(Vecsoa a, int index){
	return vec(((float*)&a.x)[index],((float*)&a.y)[index],((float*)&a.z)[index]);
}
Vecsoa soa_insert(Vecsoa a, Vec b, int index){
	Vecsoa r=a;
	((float*)&r.x)[index]=VGetX(b);
	((float*)&r.y)[index]=VGetY(b);
	((float*)&r.z)[index]=VGetZ(b);
	return r;
}
Vecsoa radiance(const Raysoa r_, SpuRng *rng){
  vector float t;                               // distance to intersection
  vector signed int id;                        // id of intersected object
  vector signed int depth={0,0,0,0};
  Vecsoa cl=SOAPromote(vec(0,0,0));   // accumulated color
  Vecsoa cf=SOAPromote(vec(1,1,1));  // accumulated reflectance
  vector unsigned int mask3={0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF};
  vector unsigned int previous_diffuse=spu_splats(0u);
  Raysoa r=r_;
  while (1){
	vector unsigned int aux1;
	
	aux1=vecsoa_intersect(&r,&t,&id);
	mask3=spu_and(mask3,aux1);
	if(spu_extract(spu_gather(mask3),0)==0x00) return cl;
	unsigned int * eid=(unsigned int*)&id;

	Spheresoa obj=spheresoa(spheres[*eid],spheres[eid[1]],spheres[eid[2]],spheres[eid[3]]);
	
	Vecsoa 	x=soa_vecsum(r.o, soa_vecdbl(r.d,t)),
			n=soa_vecnorm(soa_vecsub(x,obj.p)),
			nl,
			f=obj.c;
	Vec teste1=soa_dot(n,r.d);
	vector unsigned int mask1=spu_cmpgt(vzero,teste1);
	nl.x=spu_sel(vecdbl(n.x,-1.f),n.x,mask1);
	nl.y=spu_sel(vecdbl(n.y,-1.f),n.y,mask1);
	nl.z=spu_sel(vecdbl(n.z,-1.f),n.z,mask1);
	
	Vec p = _fmaxf4(f.x,f.y);
	p=_fmaxf4(p,f.z);
	
	/* Direct sampling already accounts for light hits after a diffuse bounce. */
	Vecsoa cl1 = soa_vecsum(cl,soa_vecmul(cf,obj.e));
	cl=soa_select(cl,cl1,spu_and(mask3,spu_xor(previous_diffuse,0xffffffffu)));
	mask3=spu_and(mask3,spu_cmpgt(p,vzero));
	if(spu_extract(spu_gather(mask3),0)==0) return cl;
	
	depth=spu_add(depth,spu_splats(1));
	vector unsigned int mask4=spu_cmpgt(depth,spu_splats(5));
	vector unsigned int mask5=spu_splats(0xffffffffu);
	if(spu_extract(spu_gather(spu_and(mask3,mask4)),0))
		mask5=spu_cmpgt(p,spu_rng_uniform(rng));
	f=soa_select(f,soa_vecdbl(f,spu_re(p)),spu_and(mask4,mask5));
	vector unsigned int teste2=spu_and(mask4,spu_xor(mask5,0xffffffff));
	mask3=spu_and(mask3,spu_xor(teste2,0xffffffff));
	if(spu_extract(spu_gather(mask3),0)==0) return cl;
	cf=soa_select(cf,soa_vecmul(cf,f),mask3);
	
	//Calculating REFR only when an active lane needs it.
	vector unsigned int refr_mask=spu_and(mask3,spu_cmpeq(obj.refl,REFR));
	if(spu_extract(spu_gather(refr_mask),0)) {
		Raysoa reflRay=criaRaysoa2(x, soa_vecsub(r.d,soa_vecdbl(n,vecdbl(soa_dot(n,r.d),2))));
		vector unsigned int into = spu_cmpgt(soa_dot(n,nl),vzero);
		Vec nc=spu_splats(1.f), nt=spu_splats(1.5f), nnt, ddn=soa_dot(r.d,nl), cos2t;
		nnt=spu_sel(vecmult(nt,spu_re(nc)),vecmult(nc,spu_re(nt)), into);
		cos2t=vecsub(vones,nnt*nnt*(vecsub(vones,ddn*ddn)));
		Vecsoa tdir=soa_vecnorm(soa_vecsub(soa_vecdbl(r.d,nnt), soa_vecdbl(n, spu_sel(spu_splats(-1.f),vones,into)*vecsum(ddn*nnt,_sqrtf4(cos2t)))));
		Vec a=vecsub(nt,nc),
			b=vecsum(nt,nc),
			R0=a*a*spu_re(b*b),
			c=vecsub(vones, spu_sel(soa_dot(tdir,n),spu_splats(-1.f)*ddn,into));
		Vec Re=R0+(vones-R0)*c*c*c*c*c,
			Tr=vecsub(vones,Re),
			P=vecsum(vquarter,vhalf*Re),
			RP=Re*spu_re(P),
			TP=Tr*spu_re(vecsub(vones,P));
		vector unsigned int mask6=spu_cmpgt(P, spu_rng_uniform(rng));
		vector unsigned int mask7=spu_cmpgt(vzero,cos2t);
		//mask7 e mask8 verdadeiros mantem o cf constante e muda o raio para reflRay
		
		Vecsoa cf1=soa_select(soa_select(soa_vecdbl(cf,TP),soa_vecdbl(cf,RP),mask6),cf,mask7);
		Raysoa r1=criaRaysoa2(x,tdir);
		Raysoa result=criaRaysoa2(soa_select(r1.o,reflRay.o,spu_or(mask6,mask7)), soa_select(r1.d,reflRay.d,spu_or(mask6,mask7)));
		r=criaRaysoa2(soa_select(r.o,result.o,refr_mask),soa_select(r.d,result.d,refr_mask));
		cf=soa_select(cf,cf1,refr_mask);
	}
	
	//Calculating DIFF only when an active lane needs it.
	vector unsigned int diff_mask=spu_and(mask3,spu_cmpeq(obj.refl,DIFF));
	if(spu_extract(spu_gather(diff_mask),0)) {
		Vecsoa safe_x=soa_select(SOAPromote(vec(50,40.8,81.6)),x,diff_mask);
		cl=soa_vecsum(cl,direct_light(safe_x,nl,cf,rng,diff_mask));
		Vec r1=vecdbl(spu_rng_uniform(rng),2*M_PI), r2=spu_rng_uniform(rng), r2s=_sqrtf4(r2);
		Vecsoa w=nl,
				u,
				v;
		vector unsigned int mask8=spu_cmpabsgt(w.x,vdotone);
		Vecsoa aux100=SOAPromote(vec(1,0,0)),
				aux010=SOAPromote(vec(0,1,0));
		u=soa_vecnorm(soa_veccross(soa_select(aux100,aux010,mask8),w));
		v=soa_veccross(w,u);
		Vecsoa d=soa_vecnorm(soa_vecsum(soa_vecsum(soa_vecdbl(u,_cosf4(r1)*r2s), soa_vecdbl(v, _sinf4(r1)*r2s)), soa_vecdbl(w,_sqrtf4(vecsub(vones,r2)))));
		Raysoa raio1=criaRaysoa2(x,d);
		r.o=soa_select(r.o,raio1.o,diff_mask);
		r.d=soa_select(r.d,raio1.d,diff_mask);
	}
	previous_diffuse=diff_mask;
	
	//Calculating SPEC
	{
		Raysoa r1=criaRaysoa2(x,soa_vecsub(r.d,soa_vecdbl(n,vecdbl(soa_dot(n,r.d),2))));
		vector unsigned int mask10=spu_cmpeq(obj.refl,SPEC);
		r.o = soa_select(r.o, r1.o, spu_and(mask3,mask10));
		r.d = soa_select(r.d, r1.d, spu_and(mask3,mask10));
	}
	}
}

typedef vector signed int Vecint;

Vecint vecint(int x, int y, int z){
	Vecint r={x,y,z,0};
	return r;
}
vector float vsx={0,0,1,1};
vector float vsy={0,1,0,1};

int main(unsigned long long spe_id, unsigned long long data_ea, unsigned long long engp) {
	Ray cam=criaRay(vec(50,52,295.6), vecnorm(vec(0,-0.042612,-1))); // cam pos, dir
	int samps, w, h, x, y, sx, sy, idx=0;
	Vec cx, cy;
	Vecint *c_int[2];

	initSpheres();
	loadSpheresSOA();
	
	mfc_get(&ctxt, data_ea, sizeof(contexto), 0,0,0);
	spu_writech(MFC_WrTagMask, 1);
	spu_mfcstat(MFC_TAG_UPDATE_ALL);

	samps=ctxt.samps;
	w=ctxt.w;
	h=ctxt.h;
	
	cx=vec(w*.5135/h,0,0);
	cy=vecdbl(vecnorm(vec_cross(cx,cam.d)),.5135);
	
	unsigned int statushi = (unsigned int)(ctxt.status>>32);
	unsigned int statuslo = (unsigned int)ctxt.status;

	c_int[0]=(Vecint*)memalign(128, ctxt.w*sizeof(Vecint));
	c_int[1]=(Vecint*)memalign(128, ctxt.w*sizeof(Vecint));
	if (!c_int[0] || !c_int[1]) {
		fprintf(stderr, "SPU row allocation failed\n");
		return 1;
	}
	
	Vec avgsamples=spu_re(spu_splats((float)samps));
	
	do {
		unsigned int mfcstatus;

		do {
			spu_mfcdma64(status, statushi, statuslo, 128, 0, MFC_GETLLAR_CMD);
			spu_readch(MFC_RdAtomicStat);

			y=(status[0]++);

			spu_mfcdma64(status, statushi, statuslo, 128, 0, MFC_PUTLLC_CMD);
			mfcstatus = spu_readch(MFC_RdAtomicStat) & MFC_PUTLLC_STATUS;
		} while(mfcstatus);
		
		if(y>=h){
			status_termination=1;
			spu_mfcdma64(&status_termination, (unsigned int)(ctxt.status_termination>>32),
		             (unsigned int)ctxt.status_termination, sizeof(unsigned int), 1, MFC_PUTB_CMD);
			break;
		}

		spu_writech(MFC_WrTagMask, 1<<idx);
		spu_mfcstat(MFC_TAG_UPDATE_ALL);
		
		fprintf(stderr,"\rRendering (%d spp) %5.2f%%",samps*4,100.*y/(h-1));
		int linha=(h-y-1)*w;
		for (x=0; x<w; x++) {  // Loop cols
			SpuRng rng;
			spu_rng_seed(&rng,ctxt.seed,(uint64_t)y*w+x);
			Vec c=vzero;
			int s;
			Vec r=vzero;
			Vecsoa rsoa=SOAPromote(vzero);
			for (s=0; s<samps; s++) {
				vector float vrand1=spu_rng_uniform(&rng);
				vector float vrand2=spu_rng_uniform(&rng);
				vector float vr1=vrand1*vtwos;
				vector float vdx=spu_sel(spu_sub(vones,_sqrtf4(spu_sub(vtwos,vr1))),spu_sub(_sqrtf4(vr1),vones),spu_cmpgt(vones,vr1));
				vector float vr2=vrand2*vtwos;
				vector float vdy=spu_sel(spu_sub(vones,_sqrtf4(spu_sub(vtwos,vr2))),spu_sub(_sqrtf4(vr2),vones),spu_cmpgt(vones,vr2));
				vector float auxx=spu_sub(spu_mul(spu_add(spu_mul(spu_add(spu_add(vsx,vhalf),vdx),vhalf),spu_splats((float)x)),spu_re(spu_splats((float)w))),vhalf);
				vector float auxy=spu_sub(spu_mul(spu_add(spu_mul(spu_add(spu_add(vsy,vhalf),vdy),vhalf),spu_splats((float)y)),spu_re(spu_splats((float)h))),vhalf);
				Vecsoa soad=soa_vecsum(soa_vecsum(soa_vecdbl(SOAPromote(cx),auxx),SOAPromote(cam.d)),
										soa_vecdbl(SOAPromote(cy),auxy));
				Raysoa raysoa;
				raysoa.o=soa_vecsum(SOAPromote(cam.o), soa_vecdbl(soad,vfactor));
				raysoa.d=soa_vecnorm(soad);
				rsoa=soa_vecsum(rsoa, soa_vecdbl(radiance(raysoa, &rng),avgsamples));
				
			} // Camera rays are pushed ^^^^^ forward to start in interior
			int jj=0;
			for(; jj<4; jj++){
				Vec xxx=soa_extract(rsoa,jj);
				c=vecsum(c, vecmult(vec(clamp(VGetX(xxx)),clamp(VGetY(xxx)),clamp(VGetZ(xxx))),vquarter));
			}
			
			c_int[idx][x]=vecint(toInt(VGetX(c)),toInt(VGetY(c)),toInt(VGetZ(c)));
		}
		unsigned long long ea_retorno=ctxt.retorno+(linha*sizeof(Vecint));
		spu_mfcdma64(c_int[idx], (unsigned int)(ea_retorno>>32),
			(unsigned int)ea_retorno, ctxt.w*sizeof(Vecint), idx, MFC_PUT_CMD);
		idx^=1;
	} while(1);
	spu_writech(MFC_WrTagMask, 0x03);
	spu_mfcstat(MFC_TAG_UPDATE_ALL);
	free(c_int[0]);
	free(c_int[1]);
	return 0;
}
