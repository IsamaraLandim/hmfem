#include <stdio.h>

#define N 64
#define L 32

struct add_scale_params {
	double *x;
	double *y;
	double s;
	unsigned int m;
};
__global__ void __device_add_scale(struct add_scale_params parms)
{
	double s0 = parms.s;
	//for(unsigned int i = N * blockIdx.x + threadIdx.x; i < parms.m; i += N * gridDim.x)
	unsigned int i = N * blockIdx.x + threadIdx.x;
	if(i < parms.m)
	{
		parms.x[i] += parms.y[i] * s0;
	}
}
extern "C" {
void _device_add_scale(double *x, double *y, double s, unsigned int m)
{
	struct add_scale_params parms;
	parms.x = x;
	parms.y = y;
	parms.s = s;
	parms.m = m;
	//__device_add_scale<<< M, N >>>(parms);
	__device_add_scale<<< (m + N - 1)/N, N >>>(parms);
}
}
///////////////////////////////

struct sub_scale_params {
	double *x;
	double *y;
	double s;
	unsigned int m;
};
__global__ void __device_sub_scale(struct sub_scale_params parms)
{
	double s0 = parms.s;
	//for(unsigned int i = N * blockIdx.x + threadIdx.x; i < parms.m; i += N * gridDim.x)
	unsigned int i = N * blockIdx.x + threadIdx.x;
	if(i < parms.m)
	{
		parms.x[i] -= parms.y[i] * s0;
	}
}
extern "C" {
void _device_sub_scale(double *x, double *y, double s, unsigned int m)
{
	struct sub_scale_params parms;
	parms.x = x;
	parms.y = y;
	parms.s = s;
	parms.m = m;
	//__device_sub_scale<<< M, N >>>(parms);
	__device_sub_scale<<< (m + N - 1)/N, N >>>(parms);
}
}
///////////////////////////////

struct scale_add_params {
	double *x;
	double *y;
	double s;
	unsigned int m;
};
__global__ void __device_scale_add(struct scale_add_params parms)
{
	double s0 = parms.s;
	//for(unsigned int i = N * blockIdx.x + threadIdx.x; i < parms.m; i += N * gridDim.x)
	unsigned int i = N * blockIdx.x + threadIdx.x;
	if(i < parms.m)
	{
		parms.x[i] = parms.x[i] * s0 + parms.y[i];
	}
}
extern "C" {
void _device_scale_add(double *x, double *y, double s, unsigned int m)
{
	struct scale_add_params parms;
	parms.x = x;
	parms.y = y;
	parms.s = s;
	parms.m = m;
	//__device_scale_add<<< M, N >>>(parms);
	__device_scale_add<<< (m + N - 1)/N, N >>>(parms);
}
}
///////////////////////////////

struct vector_product_params {
    const double *u;
    const double *v;
    double *w;
    unsigned int n;
};
__global__ void __device_vector_product(struct vector_product_params parms) 
{
	__shared__ double _acc[N];        
	unsigned int tid = threadIdx.x;

    double s = 0.0;
    for (unsigned int j = N * blockIdx.x + threadIdx.x; j < parms.n; j += N * gridDim.x)
	{
        s += parms.u[j] * parms.v[j];
    }
    _acc[tid] = s;

	for (unsigned int i = N >> 1; i > 0; i >>= 1)
	{
        __syncthreads();
        if (tid < i) _acc[tid] += _acc[tid + i];
    }
    if (tid == 0) parms.w[blockIdx.x] = _acc[tid];
}

extern "C" {
void _device_vector_product(double *u, double *v, double *w, unsigned int n, unsigned int m)
{
	struct vector_product_params parms;
	parms.u = u;
	parms.v = v;
	parms.w = w;
	parms.n = n;
	//__device_vector_product<<<M, N>>>(parms);
	__device_vector_product<<<m, N>>>(parms);
}
}

///////////////////////////////
struct linear_operator_params {
	const int *cnt;
	const int *dsp;
	const int *col;
	const double *ele;
	const double *u;
	double *v;
	int n;
};

texture<int2, 1> tex_u1;
texture<int2, 1> tex_u2;
texture<int2, 1> tex_e1;
texture<int2, 1> tex_e2;
__global__ void __device_linear_operator(struct linear_operator_params parms)
{
#if 1
	//for(unsigned int j = N * blockIdx.x + threadIdx.x; j < parms.n; j += N * gridDim.x)
	unsigned int j = N * blockIdx.x + threadIdx.x;
	if(j < parms.n)
	{
		unsigned int blkStart = parms.dsp[j];
		unsigned int blkStop = blkStart + L * parms.cnt[j];

		double s = 0.0;
		for(unsigned int i = blkStart; i < blkStop; i += L)
		{
			unsigned int q = parms.col[i];
			double a = parms.ele[i];
			//double b = parms.u[q];
			int2 c = tex1Dfetch(tex_u1, q);
			double b = __hiloint2double(c.y, c.x);
			s += a * b;
		}
		parms.v[j] = s;
	}
#else
	//for(unsigned int j = N * blockIdx.x + threadIdx.x; j < parms.n; j += N * gridDim.x)
	unsigned int j = N * blockIdx.x + threadIdx.x;
	if(j < parms.n)
	{
		unsigned int blkStart = parms.dsp[j];
		unsigned int blkStop = blkStart + parms.cnt[j];

		double s = 0.0;
		for(unsigned int i = blkStart; i < blkStop; i++)
		{
			unsigned int q = parms.col[i];
			//double a = parms.ele[i];
			int2 e = tex1Dfetch(tex_e1, i);
			double a = __hiloint2double(e.y, e.x);
			//double b = parms.u[q];
			int2 c = tex1Dfetch(tex_u1, q);
			double b = __hiloint2double(c.y, c.x);
			s += a * b;
		}
		parms.v[j] = s;
	}
/*
	unsigned int p = blockIdx.x * blockDim.x + threadIdx.x;
	if(p < n)
	{
		unsigned int l = dsp[p];
		unsigned int r = l + cnt[p];

		double s = 0.0;
		for(unsigned int j = l; j < r; j++)
		{
			unsigned int q = col[j];
			double a = ele[j];
			double b = tex1Dfetch(tex_u, q);
			s += a * b;
		}
		v[p] = s;
	}
*/
#endif
}
extern "C" {
void _device_linear_operator(int *cnt, int *dsp, int *col, double *ele, int l, int m, int n, double *u, double *v)
{
	cudaBindTexture(0, tex_u1, (int2*)u, sizeof(double) * m);
	//cudaBindTexture(0, tex_e1, (int2*)ele, sizeof(double) * l);

	struct linear_operator_params parms;
	parms.cnt = cnt;
	parms.dsp = dsp;
	parms.col = col;
	parms.ele = ele;
	parms.u = u;
	parms.v = v;
	parms.n = n;
	__device_linear_operator<<< (n + N - 1)/N, N >>>(parms);
	//__device_linear_operator<<< M, N >>>(parms);

	//cudaUnbindTexture(tex_e1);
	cudaUnbindTexture(tex_u1);
}
}


///////////////////////////////
struct jacobi_params {
	const int *cnt;
	const int *dsp;
	const int *col;
	const double *ele;
	const double *f;
	const double *u;
	double *v;
	int n;
};

__global__ void __device_jacobi(struct jacobi_params parms)
{
#if 1
	//for(unsigned int j = N * blockIdx.x + threadIdx.x; j < parms.n; j += N * gridDim.x)
	unsigned int j = N * blockIdx.x + threadIdx.x;
	if(j < parms.n)
	{
		unsigned int blkStart = parms.dsp[j];
		unsigned int blkStop = blkStart + L * parms.cnt[j];

		double s = parms.f[j];
		for(unsigned int i = blkStart; i < blkStop; i += L)
		{
			unsigned int q = parms.col[i];
			double a = parms.ele[i];
			//double b = parms.u[q];
			int2 c = tex1Dfetch(tex_u2, q);
			double b = __hiloint2double(c.y, c.x);
			s -= a * b;
		}
		parms.v[j] = s;
	}
#else
	//for(unsigned int j = N * blockIdx.x + threadIdx.x; j < parms.n; j += N * gridDim.x)
	unsigned int j = N * blockIdx.x + threadIdx.x;
	if(j < parms.n)
	{
		unsigned int blkStart = parms.dsp[j];
		unsigned int blkStop = blkStart + parms.cnt[j];

		double s = parms.f[j];
		for(unsigned int i = blkStart; i < blkStop; i++)
		{
			unsigned int q = parms.col[i];
			//double a = parms.ele[i];
			int2 e = tex1Dfetch(tex_e2, i);
			double a = __hiloint2double(e.y, e.x);
			//double b = parms.u[q];
			int2 c = tex1Dfetch(tex_u2, q);
			double b = __hiloint2double(c.y, c.x);
			s -= a * b;
		}
		parms.v[j] = s;
	}
/*
void operator()(const vector<S> &_f, vector<S> &_u)
{
	const T *acnt = &_acnt[0], *acol = &_acol[0];
	const S *aele = &_aele[0], *adia = &_adia[0];
	const S *f = &_f[0];
	S *u = &_u[0];
	S *v = &_v[0];

	S omega = _omega;
	int dsize = (int)_acnt.size();
	for(int i = 0; i < dsize; i++)
	{
		S s = f[i];
		int csize = acnt[i];
		for(int j = 0; j < csize; j++)
		{
			T c = *acol++;
			S t = *aele++;
			s -= t * u[c];
		}
		v[i] = s;
	}

	_com.accumulate(_v);
	for(int i = 0; i < dsize; i++)
	{
		u[i] += omega * v[i] * adia[i];
	}
}
*/
#endif
}
extern "C" {
void _device_jacobi(int *cnt, int *dsp, int *col, double *ele, int l, int n, double *f, double *u, double *v)
{
	cudaBindTexture(0, tex_u2, (int2*)u, sizeof(double) * n);
	//cudaBindTexture(0, tex_e2, (int2*)ele, sizeof(double) * l);

	struct jacobi_params parms;
	parms.cnt = cnt;
	parms.dsp = dsp;
	parms.col = col;
	parms.ele = ele;
	parms.f = f;
	parms.u = u;
	parms.v = v;
	parms.n = n;
	__device_jacobi<<< (n + N - 1)/N, N >>>(parms);
	//__device_jacobi<<< M, N >>>(parms);

	//cudaUnbindTexture(tex_e2);
	cudaUnbindTexture(tex_u2);
}
}

struct set_mul_params {
	double *x;
	double *y;
	double *z;
	double w;
	unsigned int m;
};
__global__ void __device_set_mul(struct set_mul_params parms)
{
	//for(unsigned int i = N * blockIdx.x + threadIdx.x; i < parms.m; i += N * gridDim.x)
	unsigned int i = N * blockIdx.x + threadIdx.x;
	if(i < parms.m)
	{
		parms.x[i] = parms.w * parms.y[i] * parms.z[i];
	}
}
extern "C" {
void _device_set_mul(double *x, double *y, double *z, double w, unsigned int m)
{
	struct set_mul_params parms;
	parms.x = x;
	parms.y = y;
	parms.z = z;
	parms.w = w;
	parms.m = m;
	//__device_set_mul<<< M, N >>>(parms);
	__device_set_mul<<< (m + N - 1)/N, N >>>(parms);
}
}

struct add_mul_params {
	double *x;
	double *y;
	double *z;
	double w;
	unsigned int m;
};
__global__ void __device_add_mul(struct add_mul_params parms)
{
	//for(unsigned int i = N * blockIdx.x + threadIdx.x; i < parms.m; i += N * gridDim.x)
	unsigned int i = N * blockIdx.x + threadIdx.x;
	if(i < parms.m)
	{
		parms.x[i] += parms.w * parms.y[i] * parms.z[i];
	}
}
extern "C" {
void _device_add_mul(double *x, double *y, double *z, double w, unsigned int m)
{
	struct add_mul_params parms;
	parms.x = x;
	parms.y = y;
	parms.z = z;
	parms.w = w;
	parms.m = m;
	//__device_add_mul<<< M, N >>>(parms);
	__device_add_mul<<< (m + N - 1)/N, N >>>(parms);
}
}
////
struct diagonal_solver_params {
	double *x;
	double *y;
	double *z;
	unsigned int m;
};
__global__ void __device_diagonal_solver(struct diagonal_solver_params parms)
{
	//for(unsigned int i = N * blockIdx.x + threadIdx.x; i < parms.m; i += N * gridDim.x)
	unsigned int i = N * blockIdx.x + threadIdx.x;
	if(i < parms.m)
	{
		parms.x[i] = parms.y[i] * parms.z[i];
	}
}
extern "C" {
void _device_diagonal_solver(double *x, double *y, double *z, unsigned int m)
{
	struct diagonal_solver_params parms;
	parms.x = x;
	parms.y = y;
	parms.z = z;
	parms.m = m;
	//__device_diagonal_solver<<< M, N >>>(parms);
	__device_diagonal_solver<<< (m + N - 1)/N, N >>>(parms);
}
}
////
struct com_snd_params {
	double *x;
	double *svec;
	int *rcom;
	unsigned int n;
};
__global__ void __device_com_snd(struct com_snd_params parms)
{
	//for(int i = 0; i < (int)_rcom.size(); i++) _svec[i] = _x[_rcom[i]];
	//for(unsigned int i = N * blockIdx.x + threadIdx.x; i < parms.n; i += N * gridDim.x)
	unsigned int i = N * blockIdx.x + threadIdx.x;
	if(i < parms.n)
	{
		parms.svec[i] = parms.x[parms.rcom[i]];
	}
}
extern "C" {
void _device_com_snd(double *x, double *svec, int *rcom, unsigned int n)
{
	struct com_snd_params parms;
	parms.x = x;
	parms.svec = svec;
	parms.rcom = rcom;
	parms.n = n;
	//__device_com_snd<<< M, N >>>(parms);
	__device_com_snd<<< (n + N - 1)/N, N >>>(parms);
}
}

struct com_rcv_params {
	double *x;
	double *rvec;
	int *rcom;
	unsigned int n;
};
__device__ inline void atomic_add(double *p, double d)
{
	unsigned long long int a, b, t;
	a = *(unsigned long long int*)p;
	do
	{
		t = a;
		b = __double_as_longlong(__longlong_as_double(a) + d);
		a = atomicCAS((unsigned long long int*)p, a, b);
	}
	while(a != t);
}
__global__ void __device_com_rcv(struct com_rcv_params parms)
{
	//for(int i = 0; i < (int)_rcom.size(); i++) _x[_rcom[i]] += _rvec[i];
	//for(unsigned int i = N * blockIdx.x + threadIdx.x; i < parms.n; i += N * gridDim.x)
	unsigned int i = N * blockIdx.x + threadIdx.x;
	if(i < parms.n)
	{
		//parms.x[parms.rcom[i]] += parms.rvec[i];
		atomic_add(parms.x + parms.rcom[i], parms.rvec[i]);
	}
}
extern "C" {
void _device_com_rcv(double *x, double *rvec, int *rcom, unsigned int n)
{
	struct com_rcv_params parms;
	parms.x = x;
	parms.rvec = rvec;
	parms.rcom = rcom;
	parms.n = n;
	//__device_com_rcv<<< M, N >>>(parms);
	__device_com_rcv<<< (n + N - 1)/N, N >>>(parms);
}
}

struct sub_params {
	double *x;
	double *y;
	unsigned int m;
};
__global__ void __device_sub(struct sub_params parms)
{
	//for(unsigned int i = N * blockIdx.x + threadIdx.x; i < parms.m; i += N * gridDim.x)
	unsigned int i = N * blockIdx.x + threadIdx.x;
	if(i < parms.m)
	{
		parms.y[i] -= parms.x[i];
	}
}
extern "C" {
void _device_sub(double *x, double *y, unsigned int m)
{
	struct sub_params parms;
	parms.x = x;
	parms.y = y;
	parms.m = m;
	//__device_sub<<< M, N >>>(parms);
	__device_sub<<< (m + N - 1)/N, N >>>(parms);
}
}

struct add_params {
	double *x;
	double *y;
	unsigned int m;
};
__global__ void __device_add(struct add_params parms)
{
	//for(unsigned int i = N * blockIdx.x + threadIdx.x; i < parms.m; i += N * gridDim.x)
	unsigned int i = N * blockIdx.x + threadIdx.x;
	if(i < parms.m)
	{
		parms.y[i] += parms.x[i];
	}
}
extern "C" {
void _device_add(double *x, double *y, unsigned int m)
{
	struct add_params parms;
	parms.x = x;
	parms.y = y;
	parms.m = m;
	//__device_add<<< M, N >>>(parms);
	__device_add<<< (m + N - 1)/N, N >>>(parms);
}
}

struct rsb_params {
	double *x;
	double *y;
	unsigned int m;
};
__global__ void __device_rsb(struct rsb_params parms)
{
	//for(unsigned int i = N * blockIdx.x + threadIdx.x; i < parms.m; i += N * gridDim.x)
	unsigned int i = N * blockIdx.x + threadIdx.x;
	if(i < parms.m)
	{
		parms.y[i] = parms.x[i] - parms.y[i];
	}
}
extern "C" {
void _device_rsb(double *x, double *y, unsigned int m)
{
	struct rsb_params parms;
	parms.x = x;
	parms.y = y;
	parms.m = m;
	//__device_rsb<<< M, N >>>(parms);
	__device_rsb<<< (m + N - 1)/N, N >>>(parms);
}
}

__global__ void __device_matrix_conversion(int *_acnt, int *_adsp, int *_acol, double *_aele, int *_tdsp, int *_tcol, double *_tele, int _n)
{
	//for(int j = N * blockIdx.x + threadIdx.x; j < _n; j += N * gridDim.x)
	int j = N * blockIdx.x + threadIdx.x;
	if(j < _n)
	{
		int tdsp = _tdsp[j];
		int acnt = _acnt[j];
		int adsp = _adsp[j];
		for(int k = adsp, i = tdsp; i < tdsp + acnt; i++, k += L)
		{
			_acol[k] = _tcol[i];
			_aele[k] = _tele[i];
		}
	}
}
extern "C" {
void _device_matrix_conversion(int *_acnt, int *_adsp, int *_acol, double *_aele, int *_tdsp, int *_tcol, double *_tele, int _n)
{
	//__device_matrix_conversion<<< M, N >>>(_acnt, _adsp, _acol, _aele, _tdsp, _tcol, _tele, _n);
	__device_matrix_conversion<<< (_n + N - 1)/N, N >>>(_acnt, _adsp, _acol, _aele, _tdsp, _tcol, _tele, _n);
}
}


































