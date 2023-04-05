template<class T, class S> class device_generic_operator
{
public:
	virtual void operator()(const device_vector<S> &_u, device_vector<S> &_v) const = 0;
};

extern "C" {
//void _device_linear_operator(int *cnt, int *dsp, int *col, float *ele, int l, int m, int n, float *u, float *v);
void _device_linear_operator(int *cnt, int *dsp, int *col, double *ele, int l, int m, int n, double *u, double *v);
}
template<class T, class S> class device_linear_operator : public device_generic_operator<T, S>
{
private:
	const device_vector<T> &_acnt;
	const device_vector<T> &_adsp;
	const device_vector<T> &_acol;
	const device_vector<S> &_aele;
public:
	device_linear_operator(const device_vector<T> &_acnt, const device_vector<T> &_adsp, const device_vector<T> &_acol, const device_vector<S> &_aele)
		: _acnt(_acnt), _adsp(_adsp), _acol(_acol), _aele(_aele)
	{
	}
	void operator()(const device_vector<S> &_u, device_vector<S> &_v) const
	{
/*
		device_timer timer;
		timer.start();
*/
		//_device_linear_operator(&_acnt[0], &_adsp[0], &_acol[0], &_aele[0], _aele.size(), _u.size(), _v.size(), &_u[0], &_v[0]);
		_device_linear_operator(_acnt.data(), _adsp.data(), _acol.data(), _aele.data(), _aele.size(), _u.size(), _v.size(), _u.data(), _v.data());
/*
		cudaThreadSynchronize();
 		timer.stop();
		cudaThreadSynchronize();

		S time = timer.time();
		S flops = 1000.0 * 2 * _acol.size() / (time * 1000.0 * 1000.0 * 1000.0);
		S gbs = 1000.0 * (3 * sizeof(S) * _acol.size() + sizeof(S) * _acnt.size()) / (time * 1024.0 * 1024.0 * 1024.0);
		cout << "Processing time: " << time << " (ms)" << endl;
		cout << "Performance: " << flops << " (GFLOPS)" << endl;
		cout << "Performance: " << gbs << " (GB/s)" << endl;
*/
	}
};

#if 0
template<class T, class S> class simple_prolongation : public generic_operator<T, S>
{
private:
	const vector<T> &_rcnt;
	const vector<T> &_rcol;
public:
	simple_prolongation(const vector<T> &_rcnt, const vector<T> &_rcol)
		: _rcnt(_rcnt), _rcol(_rcol)
	{
	}
	void operator()(const vector<S> &_u, vector<S> &_s) const
	{
		const S *u = &_u[0];
		S *s = &_s[0];
		int size = (int)_rcnt.size();
		const T *rcnt = &_rcnt[0];
		const T *rcol = &_rcol[0];
		S t;
		for(int i = 0; i < size; i++)
		{
			_mm_prefetch((char*)&rcol[64],_MM_HINT_NTA);
			int c = rcnt[i];
			switch(c)
			{
			case 1:
				s[i] = u[*rcol++];
				break;
			case 2:
				t = u[*rcol++];
				t += u[*rcol++];
				s[i] = t / S(2.0);
				break;
			case 3:
				t = u[*rcol++];
				t += u[*rcol++];
				t += u[*rcol++];
				s[i] = t / S(3.0);
				break;
			case 4:
				t = u[*rcol++];
				t += u[*rcol++];
				t += u[*rcol++];
				t += u[*rcol++];
				s[i] = t / S(4.0);
				break;
			default:
				S t = S(0.0);
				for(int j = 0; j < c; j++)
				{
					t += u[*rcol++];
				}
				s[i] = t / S(c);
				break;
			}
		}
	}
};

template<class T, class S> class simple_restriction : public generic_operator<T, S>
{
private:
	const vector<T> &_rcnt;
	const vector<T> &_rcol;
public:
	simple_restriction(const vector<T> &_rcnt, const vector<T> &_rcol)
		: _rcnt(_rcnt), _rcol(_rcol)
	{
	}
	void operator()(const vector<S> &_u, vector<S> &_s) const
	{
		const S *u = &_u[0];
		S *s = &_s[0];
		int dsize = (int)_s.size();
		for(int i = 0; i < dsize; i++) s[i] = S(0.0);
		int size = (int)_rcnt.size();
		const T *rcnt = &_rcnt[0];
		const T *rcol = &_rcol[0];
		S t;
		for(int i = 0; i < size; i++)
		{
			_mm_prefetch((char*)&rcol[64],_MM_HINT_NTA);
			int c = rcnt[i];
			switch(c)
			{
			case 1:
				s[*rcol++] += u[i];
				break;
			case 2:
				t = u[i] / S(2.0);
				s[*rcol++] += t;
				s[*rcol++] += t;
				break;
			case 3:
				t = u[i] / S(3.0);
				s[*rcol++] += t;
				s[*rcol++] += t;
				s[*rcol++] += t;
				break;
			case 4:
				t = u[i] / S(4.0);
				s[*rcol++] += t;
				s[*rcol++] += t;
				s[*rcol++] += t;
				s[*rcol++] += t;
				break;
			default:
				t = u[i] / S(c);
				for(int j = 0; j < c; j++)
				{
					s[*rcol++] += t;
				}
				break;
			}
		}
	}
};

template<class T, class S> class prolongation : public generic_operator<T, S>
{
private:
	const vector<T> &_rcnt;
	const vector<T> &_rcol;
	const vector<S> &_rele;
public:
	prolongation(const vector<T> &_rcnt, const vector<T> &_rcol, const vector<S> &_rele)
		: _rcnt(_rcnt), _rcol(_rcol), _rele(_rele)
	{
	}
	void operator()(const vector<S> &_u, vector<S> &_s) const
	{
		const T *rcol = &_rcol[0];
		const S *rele = &_rele[0];
		for(int i = 0; i < (int)_rcnt.size(); i++)
		{
			int c = _rcnt[i];
			S s = 0.0;
			for(int j = 0; j < c; j++)
			{
				s += *rele++ * _u[*rcol++];
			}
			_s[i] = s;
		}
	}
};

template<class T, class S> class restriction : public generic_operator<T, S>
{
private:
	const vector<T> &_rcnt;
	const vector<T> &_rcol;
	const vector<S> &_rele;
public:
	restriction(const vector<T> &_rcnt, const vector<T> &_rcol, const vector<S> &_rele)
		: _rcnt(_rcnt), _rcol(_rcol), _rele(_rele)
	{
	}
	void operator()(const vector<S> &_u, vector<S> &_s) const
	{
		memset(&_s[0], 0, sizeof(S) * _s.size());
		const T *rcol = &_rcol[0];
		const S *rele = &_rele[0];
		for(int i = 0; i < (int)_rcnt.size(); i++)
		{
			int c = _rcnt[i];
			S s = _u[i];
			for(int j = 0; j < c; j++)
			{
				_s[*rcol++] += *rele++ * s;
			}
		}
	}
};

template<class T, class S> class diagonal_operator : public generic_operator<T, S>
{
private:
	const vector<S> &_adia;
	const S _asca;
public:
	diagonal_operator(const vector<S> &_adia, const S _asca)
		: _adia(_adia), _asca(_asca)
	{
	}
	void operator()(const vector<S> &_u, vector<S> &_v) const
	{
		S *adia = &_adia[0];
		S *u = &_u[0];
		S *v = &_v[0];
		int m = (int)_adia.size();

		for(int i = 0; i < m; i++)
		{
			S val = *adia++;
			v[i] = _asca * val * u[i];
		}
	}
};

template<class T, class S> class linear_operator : public generic_operator<T, S>
{
private:
	const vector<T> &_acnt;
	const vector<T> &_acol;
	const vector<S> &_aele;
public:
	linear_operator(const vector<T> &_acnt, const vector<T> &_acol, const vector<S> &_aele)
		: _acnt(_acnt), _acol(_acol), _aele(_aele)
	{
	}
	void operator()(const vector<S> &_u, vector<S> &_v) const
	{
		const T *__restrict acnt = &_acnt[0], *__restrict acol = &_acol[0];
		const S *__restrict aele = &_aele[0];
		const S *__restrict u = &_u[0];
		S *__restrict v = &_v[0];

		int size = (int)_acnt.size();
		for(int i = 0; i < size; i++)
		{
			S s = S(0.0);
			int csize = acnt[i];
			for(int j = 0; j < csize; j++)
			{
				_mm_prefetch((char*)&acol[64],_MM_HINT_NTA);
				_mm_prefetch((char*)&aele[64],_MM_HINT_NTA);
				T q = *acol++;
				S a = *aele++;
				s += a * u[q];
			}
			v[i] = s;
		}

	}
};

template<class T, class S> class linear_operator2 : public generic_operator<T, S>
{
private:
	const vector<T> &_acnt;
	const vector<T> &_acol;
	const vector<S> &_aele;
public:
	linear_operator2(const vector<T> &_acnt, const vector<T> &_acol, const vector<S> &_aele)
		: _acnt(_acnt), _acol(_acol), _aele(_aele)
	{
	}
	void operator()(const vector<S> &_u, vector<S> &_v) const
	{
		const T *__restrict acnt = &_acnt[0], *__restrict acol = &_acol[0];
		const S *__restrict aele = &_aele[0];
		const S *__restrict u = &_u[0];
		S *__restrict v = &_v[0];

		int size = (int)_acnt.size();
		for(int i = 0; i < size; i++)
		{
			S s0 = S(0.0);
			S s1 = S(0.0);
			int csize = acnt[i];
			for(int j = 0; j < csize; j++)
			{
				_mm_prefetch((char*)&acol[64],_MM_HINT_NTA);
				_mm_prefetch((char*)&aele[64],_MM_HINT_NTA);
				T q = *acol++;
				S a = *aele++;
				s0 += a * u[2*q+0];
				s1 += a * u[2*q+1];
			}
			v[2*i+0] = s0;
			v[2*i+1] = s1;
		}

	}
};

template<class T, class S> class linear_operator4 : public generic_operator<T, S>
{
private:
	const vector<T> &_acnt;
	const vector<T> &_acol;
	const vector<S> &_aele;
public:
	linear_operator4(const vector<T> &_acnt, const vector<T> &_acol, const vector<S> &_aele)
		: _acnt(_acnt), _acol(_acol), _aele(_aele)
	{
	}
	void operator()(const vector<S> &_u, vector<S> &_v) const
	{
		const T *__restrict acnt = &_acnt[0], *__restrict acol = &_acol[0];
		const S *__restrict aele = &_aele[0];
		const S *__restrict u = &_u[0];
		S *__restrict v = &_v[0];

		int size = (int)_acnt.size();
		for(int i = 0; i < size; i++)
		{
			S s0 = S(0.0);
			S s1 = S(0.0);
			S s2 = S(0.0);
			S s3 = S(0.0);
			int csize = acnt[i];
			for(int j = 0; j < csize; j++)
			{
				_mm_prefetch((char*)&acol[64],_MM_HINT_NTA);
				_mm_prefetch((char*)&aele[64],_MM_HINT_NTA);
				T q = *acol++;
				S a = *aele++;
				s0 += a * u[4*q+0];
				s1 += a * u[4*q+1];
				s2 += a * u[4*q+2];
				s3 += a * u[4*q+3];
			}
			v[4*i+0] = s0;
			v[4*i+1] = s1;
			v[4*i+2] = s2;
			v[4*i+3] = s3;
		}

	}
};

template<class T, class S> class linear_operator8 : public generic_operator<T, S>
{
private:
	const vector<T> &_acnt;
	const vector<T> &_acol;
	const vector<S> &_aele;
public:
	linear_operator8(const vector<T> &_acnt, const vector<T> &_acol, const vector<S> &_aele)
		: _acnt(_acnt), _acol(_acol), _aele(_aele)
	{
	}
	void operator()(const vector<S> &_u, vector<S> &_v) const
	{
		const T *__restrict acnt = &_acnt[0], *__restrict acol = &_acol[0];
		const S *__restrict aele = &_aele[0];
		const S *__restrict u = &_u[0];
		S *__restrict v = &_v[0];

		int size = (int)_acnt.size();
		for(int i = 0; i < size; i++)
		{
			S s0 = S(0.0);
			S s1 = S(0.0);
			S s2 = S(0.0);
			S s3 = S(0.0);
			S s4 = S(0.0);
			S s5 = S(0.0);
			S s6 = S(0.0);
			S s7 = S(0.0);
			int csize = acnt[i];
			for(int j = 0; j < csize; j++)
			{
				_mm_prefetch((char*)&acol[64],_MM_HINT_NTA);
				_mm_prefetch((char*)&aele[64],_MM_HINT_NTA);
				T q = *acol++;
				S a = *aele++;
				s0 += a * u[8*q+0];
				s1 += a * u[8*q+1];
				s2 += a * u[8*q+2];
				s3 += a * u[8*q+3];
				s4 += a * u[8*q+4];
				s5 += a * u[8*q+5];
				s6 += a * u[8*q+6];
				s7 += a * u[8*q+7];
			}
			v[8*i+0] = s0;
			v[8*i+1] = s1;
			v[8*i+2] = s2;
			v[8*i+3] = s3;
			v[8*i+4] = s4;
			v[8*i+5] = s5;
			v[8*i+6] = s6;
			v[8*i+7] = s7;
		}

	}
};

template<class T, class S> class symmetric_operator : public generic_operator<T, S>
{
private:
	const vector<T> &_acnt;
	const vector<T> &_acol;
	const vector<S> &_aele;
	const vector<S> &_adia;
public:
	symmetric_operator(const vector<T> &_acnt, const vector<T> &_acol, const vector<S> &_aele, const vector<S> &_adia)
		: _acnt(_acnt), _acol(_acol), _aele(_aele), _adia(_adia)
	{
	}
	void operator()(const vector<S> &_u, vector<S> &_v) const
	{
		const T *__restrict acnt = &_acnt[0], *__restrict acol = &_acol[0];
		const S *__restrict adia = &_adia[0], *__restrict aele = &_aele[0];
		const S *__restrict u = &_u[0];
		S *__restrict v = &_v[0];

		int dsize = (int)_adia.size();
		for(int i = 0; i < dsize; i++) *v++ = *adia++ * *u++;
		v -= dsize;
		u -= dsize;
		for(int i = 0; i < dsize; i++)
		{
			S s = v[i];
			S r = u[i];
			int csize = acnt[i];
			for(int j = 0; j < csize; j++)
			{
				_mm_prefetch((char*)&acol[64],_MM_HINT_NTA);
				_mm_prefetch((char*)&aele[64],_MM_HINT_NTA);
				T q = *acol++;
				S a = *aele++;
				S t = v[q];
				t += a * r;
				s += a * u[q];
				v[q] = t;
			}
			v[i] = s;
		}
	}
};

template<class T, class S> class symmetric_operator2 : public generic_operator<T, S>
{
private:
	const vector<T> &_acnt;
	const vector<T> &_acol;
	const vector<S> &_aele;
	const vector<S> &_adia;
public:
	symmetric_operator2(const vector<T> &_acnt, const vector<T> &_acol, const vector<S> &_aele, const vector<S> &_adia)
		: _acnt(_acnt), _acol(_acol), _aele(_aele), _adia(_adia)
	{
	}
	void operator()(const vector<S> &_u, vector<S> &_v) const
	{
		const T *__restrict acnt = &_acnt[0], *__restrict acol = &_acol[0];
		const S *__restrict adia = &_adia[0], *__restrict aele = &_aele[0];
		const S *__restrict u = &_u[0];
		S *__restrict v = &_v[0];

		int dsize = (int)_adia.size();
		for(int i = 0; i < dsize; i++) *v++ = *adia++ * *u++;
		v -= dsize;
		u -= dsize;
		for(int i = 0; i < dsize; i++)
		{
			S s0 = v[2*i+0];
			S s1 = v[2*i+1];
			S r0 = u[2*i+0];
			S r1 = u[2*i+1];
			int csize = acnt[i];
			for(int j = 0; j < csize; j++)
			{
				_mm_prefetch((char*)&acol[64],_MM_HINT_NTA);
				_mm_prefetch((char*)&aele[64],_MM_HINT_NTA);
				T q = *acol++;
				S a = *aele++;
				S t0 = v[2*q+0];
				S t1 = v[2*q+1];
				t0 += a * r0;
				t1 += a * r1;
				s0 += a * u[2*q+0];
				s1 += a * u[2*q+1];
				v[2*q+0] = t0;
				v[2*q+1] = t1;
			}
			v[2*i+0] = s0;
			v[2*i+1] = s1;
		}
	}
};

template<class T, class S> class symmetric_operator4 : public generic_operator<T, S>
{
private:
	const vector<T> &_acnt;
	const vector<T> &_acol;
	const vector<S> &_aele;
	const vector<S> &_adia;
public:
	symmetric_operator4(const vector<T> &_acnt, const vector<T> &_acol, const vector<S> &_aele, const vector<S> &_adia)
		: _acnt(_acnt), _acol(_acol), _aele(_aele), _adia(_adia)
	{
	}
	void operator()(const vector<S> &_u, vector<S> &_v) const
	{
		const T *__restrict acnt = &_acnt[0], *__restrict acol = &_acol[0];
		const S *__restrict adia = &_adia[0], *__restrict aele = &_aele[0];
		const S *__restrict u = &_u[0];
		S *__restrict v = &_v[0];

		int dsize = (int)_adia.size();
		for(int i = 0; i < dsize; i++) *v++ = *adia++ * *u++;
		v -= dsize;
		u -= dsize;
		for(int i = 0; i < dsize; i++)
		{
			S s0 = v[4*i+0];
			S s1 = v[4*i+1];
			S s2 = v[4*i+2];
			S s3 = v[4*i+3];
			S r0 = u[4*i+0];
			S r1 = u[4*i+1];
			S r2 = u[4*i+2];
			S r3 = u[4*i+3];
			int csize = acnt[i];
			for(int j = 0; j < csize; j++)
			{
				_mm_prefetch((char*)&acol[64],_MM_HINT_NTA);
				_mm_prefetch((char*)&aele[64],_MM_HINT_NTA);
				T q = *acol++;
				S a = *aele++;
				S t0 = v[4*q+0];
				S t1 = v[4*q+1];
				S t2 = v[4*q+2];
				S t3 = v[4*q+3];
				t0 += a * r0;
				t1 += a * r1;
				t2 += a * r2;
				t3 += a * r3;
				s0 += a * u[4*q+0];
				s1 += a * u[4*q+1];
				s2 += a * u[4*q+2];
				s3 += a * u[4*q+3];
				v[4*q+0] = t0;
				v[4*q+1] = t1;
				v[4*q+2] = t2;
				v[4*q+3] = t3;
			}
			v[4*i+0] = s0;
			v[4*i+1] = s1;
			v[4*i+2] = s2;
			v[4*i+3] = s3;
		}
	}
};

template<class T, class S> class hilbert_operator : public generic_operator<T, S>
{
private:
	const vector<T> &_arow;
	const vector<T> &_acol;
	const vector<S> &_aele;
public:
	hilbert_operator(const vector<T> &_arow, const vector<T> &_acol, const vector<S> &_aele)
		: _arow(_arow), _acol(_acol), _aele(_aele)
	{
	}
	void operator()(const vector<S> &_u, vector<S> &_v) const
	{
#if 1
		const T *__restrict arow = &_arow[0], *__restrict acol = &_acol[0];
		const S *__restrict aele = &_aele[0];
		const S *__restrict u = &_u[0];
		S *__restrict v = &_v[0];

		int dsize = (int)_v.size();
		for(int i = 0; i < dsize; i++) *v++ = S(0.0);
		v -= dsize;
		int size = (int)_aele.size();
		for(int i = 0; i < size; i++)
		{
			_mm_prefetch((char*)&arow[64],_MM_HINT_NTA);
			_mm_prefetch((char*)&acol[64],_MM_HINT_NTA);
			_mm_prefetch((char*)&aele[64],_MM_HINT_NTA);
			T p = *arow++;
			T q = *acol++;
			S a = *aele++;
			S s = v[p];
			s += a * u[q];
			v[p] = s;
		}
#else
		int dsize = (int)_v.size();
		for(int i = 0; i < dsize; i++)
		{
			_v[i] = S(0.0);
		}
		int size = (int)_aele.size();
		for(int i = 0; i < size; i++)
		{
			T p = _arow[i];
			S s = _v[p];
			T q = _acol[i];
			S a = _aele[i];
			s += a * _u[q];
			_v[p] = s;
		}
#endif
	}
};

template<class T, class S> class hilbert_symmetric_operator : public generic_operator<T, S>
{
private:
	const vector<T> &_arow;
	const vector<T> &_acol;
	const vector<S> &_aele;
	const vector<S> &_adia;
public:
	hilbert_symmetric_operator(const vector<T> &_arow, const vector<T> &_acol, const vector<S> &_aele, const vector<S> &_adia)
		: _arow(_arow), _acol(_acol), _aele(_aele), _adia(_adia)
	{
	}
	void operator()(const vector<S> &_u, vector<S> &_v) const
	{
#if 1
		const T *__restrict arow = &_arow[0], *__restrict acol = &_acol[0];
		const S *__restrict aele = &_aele[0], *__restrict adia = &_adia[0];
		const S *__restrict u = &_u[0];
		S *__restrict v = &_v[0];

		int dsize = (int)_adia.size();
		for(int i = 0; i < dsize; i++) *v++ = *adia++ * *u++;
		v -= dsize;
		u -= dsize;
		int size = (int)_aele.size();
		for(int i = 0; i < size; i++)
		{
			_mm_prefetch((char*)&arow[64],_MM_HINT_NTA);
			_mm_prefetch((char*)&acol[64],_MM_HINT_NTA);
			_mm_prefetch((char*)&aele[64],_MM_HINT_NTA);
			T p = *arow++;
			T q = *acol++;
			S a = *aele++;
			S s = v[p];
			S t = v[q];
			s += a * u[q];
			t += a * u[p];
			v[p] = s;
			v[q] = t;
		}
#else
		int dsize = (int)_adia.size();
		for(int i = 0; i < dsize; i++) _v[i] = _adia[i] * _u[i];
		int size = (int)_aele.size();
		for(int i = 0; i < size; i++)
		{
			T p = _arow[i];
			T q = _acol[i];
			S a = _aele[i];
			S s = _v[p];
			S t = _v[q];
			s += a * _u[q];
			t += a * _u[p];
			_v[p] = s;
			_v[q] = t;
		}
#endif
	}
};

template<class T, class S> class hilbert_struct_operator : public generic_operator<T, S>
{
private:
	const vector<T> &_aent;
	const vector<S> &_aele;
public:
	hilbert_struct_operator(const vector<T> &_aent, const vector<S> &_aele)
		: _aent(_aent), _aele(_aele)
	{
	}
	void operator()(const vector<S> &_u, vector<S> &_v) const
	{
#if 1
		const T *__restrict aent = &_aent[0];
		const S *__restrict aele = &_aele[0];
		const S *__restrict u = &_u[0];
		S *__restrict v = &_v[0];

		int dsize = (int)_v.size();
		for(int i = 0; i < dsize; i++) *v++ = S(0.0);
		v -= dsize;
		int size = (int)_aele.size();
		for(int i = 0; i < size; i++)
		{
			_mm_prefetch((char*)&aent[64],_MM_HINT_NTA);
			_mm_prefetch((char*)&aele[64],_MM_HINT_NTA);
			T p = aent[0];
			T q = aent[1];
			aent += 2;
			S a = *aele++;
			S s = v[p];
			s += a * u[q];
			v[p] = s;
		}
#else
		int dsize = (int)_v.size();
		for(int i = 0; i < dsize; i++)
		{
			_v[i] = S(0.0);
		}
		int size = (int)_aele.size();
		for(int i = 0; i < size; i++)
		{
			T p = _aent[2*i+0];
			T q = _aent[2*i+1];
			S a = _aele[i];
			S s = _v[p];
			s += a * _u[q];
			_v[p] = s;
		}
#endif
	}
};

template<class T, class S> class hilbert_struct_symmetric_operator : public generic_operator<T, S>
{
private:
	const vector<T> &_aent;
	const vector<S> &_aele;
	const vector<S> &_adia;
public:
	hilbert_struct_symmetric_operator(const vector<T> &_aent, const vector<S> &_aele, const vector<S> &_adia)
		: _aent(_aent), _aele(_aele), _adia(_adia)
	{
	}
	void operator()(const vector<S> &_u, vector<S> &_v) const
	{
#if 1
		const T *__restrict aent = &_aent[0];
		const S *__restrict aele = &_aele[0], *__restrict adia = &_adia[0];
		const S *__restrict u = &_u[0];
		S *__restrict v = &_v[0];

		int dsize = (int)_adia.size();
		for(int i = 0; i < dsize; i++) *v++ = *adia++ * *u++;
		v -= dsize;
		u -= dsize;
		int size = (int)_aele.size();
		for(int i = 0; i < size; i++)
		{
			_mm_prefetch((char*)&aent[64],_MM_HINT_NTA);
			_mm_prefetch((char*)&aele[64],_MM_HINT_NTA);
			T p = aent[0];
			T q = aent[1];
			aent += 2;
			S a = *aele++;
			S s = v[p];
			S t = v[q];
			s += a * u[q];
			t += a * u[p];
			v[p] = s;
			v[q] = t;
		}
#else
		int dsize = (int)_adia.size();
		for(int i = 0; i < dsize; i++) _v[i] = _adia[i] * _u[i];
		int size = (int)_aele.size();
		for(int i = 0; i < size; i++)
		{
			T p = _aent[2*i+0];
			T q = _aent[2*i+1];
			S a = _aele[i];
			S s = _v[p];
			S t = _v[q];
			s += a * _u[q];
			t += a * _u[p];
			_v[p] = s;
			_v[q] = t;
			//_mm_prefetch((char*)&_aent[2*i]+256,_MM_HINT_NTA);
			//_mm_prefetch((char*)&_aele[i]+256,_MM_HINT_NTA);
		}
#endif
	}
};


template<class T, class S> class hilbert_struct_operator2 : public generic_operator<T, S>
{
private:
	const vector<T> &_aent;
	const vector<S> &_aele;
public:
	hilbert_struct_operator2(const vector<T> &_aent, const vector<S> &_aele)
		: _aent(_aent), _aele(_aele)
	{
	}
	void operator()(const vector<S> &_u, vector<S> &_v) const
	{
		const T *__restrict aent = &_aent[0];
		const S *__restrict aele = &_aele[0];
		const S *__restrict u = &_u[0];
		S *__restrict v = &_v[0];

		int dsize = (int)_v.size();
		for(int i = 0; i < dsize; i++) *v++ = S(0.0);
		v -= dsize;
		int size = (int)_aele.size();
		for(int i = 0; i < size; i++)
		{
			_mm_prefetch((char*)&aent[64],_MM_HINT_NTA);
			_mm_prefetch((char*)&aele[64],_MM_HINT_NTA);
			T p = 2 * aent[0];
			T q = 2 * aent[1];
			aent += 2;
			S a = *aele++;
			S s0 = v[p+0];
			S s1 = v[p+1];
			s0 += a * u[q+0];
			s1 += a * u[q+1];
			v[p+0] = s0;
			v[p+1] = s1;
		}
	}
};

template<class T, class S> class hilbert_struct_symmetric_operator2 : public generic_operator<T, S>
{
private:
	const vector<T> &_aent;
	const vector<S> &_aele;
	const vector<S> &_adia;
public:
	hilbert_struct_symmetric_operator2(const vector<T> &_aent, const vector<S> &_aele, const vector<S> &_adia)
		: _aent(_aent), _aele(_aele), _adia(_adia)
	{
	}
	void operator()(const vector<S> &_u, vector<S> &_v) const
	{
		const T *__restrict aent = &_aent[0];
		const S *__restrict aele = &_aele[0], *__restrict adia = &_adia[0];
		const S *__restrict u = &_u[0];
		S *__restrict v = &_v[0];

		int dsize = (int)_adia.size();
		for(int i = 0; i < dsize; i++)
		{
			S a = *adia++;
			v[0] = a * u[0];
			v[1] = a * u[1];
			u+=2; v+=2;
		}
		v -= 2 * dsize;
		u -= 2 * dsize;
		int size = (int)_aele.size();
		for(int i = 0; i < size; i++)
		{
			_mm_prefetch((char*)&aent[64],_MM_HINT_NTA);
			_mm_prefetch((char*)&aele[64],_MM_HINT_NTA);
			T p = 2 * aent[0];
			T q = 2 * aent[1];
			aent += 2;
			S a = *aele++;
			S s0 = v[p+0];
			S s1 = v[p+1];
			S t0 = v[q+0];
			S t1 = v[q+1];
			s0 += a * u[q+0];
			s1 += a * u[q+1];
			t0 += a * u[p+0];
			t1 += a * u[p+1];
			v[p+0] = s0;
			v[p+1] = s1;
			v[q+0] = t0;
			v[q+1] = t1;
		}
	}
};


template<class T, class S> class hilbert_struct_operator4 : public generic_operator<T, S>
{
private:
	const vector<T> &_aent;
	const vector<S> &_aele;
public:
	hilbert_struct_operator4(const vector<T> &_aent, const vector<S> &_aele)
		: _aent(_aent), _aele(_aele)
	{
	}
	void operator()(const vector<S> &_u, vector<S> &_v) const
	{
		const T *__restrict aent = &_aent[0];
		const S *__restrict aele = &_aele[0];
		const S *__restrict u = &_u[0];
		S *__restrict v = &_v[0];

		int dsize = (int)_v.size();
		for(int i = 0; i < dsize; i++) *v++ = S(0.0);
		v -= dsize;
		int size = (int)_aele.size();
		for(int i = 0; i < size; i++)
		{
			_mm_prefetch((char*)&aent[64],_MM_HINT_NTA);
			_mm_prefetch((char*)&aele[64],_MM_HINT_NTA);
			T p = 4 * aent[0];
			T q = 4 * aent[1];
			aent += 2;
			S a = *aele++;
			S s0 = v[p+0];
			S s1 = v[p+1];
			S s2 = v[p+2];
			S s3 = v[p+3];
			s0 += a * u[q+0];
			s1 += a * u[q+1];
			s2 += a * u[q+2];
			s3 += a * u[q+3];
			v[p+0] = s0;
			v[p+1] = s1;
			v[p+2] = s2;
			v[p+3] = s3;
		}
	}
};
#if 0
template<class T, class S> class hilbert_struct_symmetric_operator4 : public generic_operator<T, S>
{
private:
	const vector<T> &_aent;
	const vector<S> &_aele;
	const vector<S> &_adia;
public:
	hilbert_struct_symmetric_operator4(const vector<T> &_aent, const vector<S> &_aele, const vector<S> &_adia)
		: _aent(_aent), _aele(_aele), _adia(_adia)
	{
	}
	void operator()(const vector<S> &_u, vector<S> &_v) const
	{
		const T *__restrict aent = &_aent[0];
		const S *__restrict aele = &_aele[0], *__restrict adia = &_adia[0];
		const S *__restrict u = &_u[0];
		S *__restrict v = &_v[0];

		if(sizeof(S) == 4)
		{
			int dsize = (int)_adia.size();
			for(int i = 0; i < dsize; i++)
			{
				__m128 a = _mm_load1_ps((float*)adia);
				adia++;
				__m128 s0 = _mm_load_ps((float*)u);
				s0 = _mm_mul_ps(a, s0);
				//_mm_store_ps((float*)v, s0);
				_mm_stream_ps((float*)v, s0);
				u+=4; v+=4;

				//*v++ = *adia++ * *u++;
			}
			v -= 4 * dsize;
			u -= 4 * dsize;
			int size = (int)_aele.size();
			for(int i = 0; i < size; i++)
			{
				_mm_prefetch((char*)&aent[64],_MM_HINT_NTA);
				_mm_prefetch((char*)&aele[64],_MM_HINT_NTA);
				T p = 4 * aent[0];
				T q = 4 * aent[1];
				aent += 2;
				__m128 a = _mm_load1_ps((float*)aele);
				aele++;
				__m128 e = _mm_load_ps((float*)&u[p]);
				__m128 f = _mm_load_ps((float*)&u[q]);
				__m128 s = _mm_load_ps((float*)&v[p]);
				__m128 t = _mm_load_ps((float*)&v[q]);
				f = _mm_mul_ps(a, f);
				e = _mm_mul_ps(a, e);
				s = _mm_add_ps(s, f);
				t = _mm_add_ps(t, e);
				_mm_store_ps((float*)&v[p], s);
				_mm_store_ps((float*)&v[q], t);
			}		
		}
		else
		{
			int dsize = (int)_adia.size();
			for(int i = 0; i < dsize; i++)
			{
				__m128d a = _mm_load1_pd((double*)adia);
				adia++;
				__m128d s0 = _mm_load_pd((double*)u);
				s0 = _mm_mul_pd(a, s0);
				//_mm_store_pd((double*)v, s0);
				_mm_stream_pd((double*)v, s0);
				u+=2; v+=2;
				__m128d s1 = _mm_load_pd((double*)u);
				s1 = _mm_mul_pd(a, s1);
				//_mm_store_pd((double*)v, s1);
				_mm_stream_pd((double*)v, s1);
				u+=2; v+=2;

				//*v++ = *adia++ * *u++;
			}
			v -= 4 * dsize;
			u -= 4 * dsize;
			int size = (int)_aele.size();
			for(int i = 0; i < size; i++)
			{
				_mm_prefetch((char*)&aent[64],_MM_HINT_NTA);
				_mm_prefetch((char*)&aele[64],_MM_HINT_NTA);
				T p = 4 * aent[0];
				T q = 4 * aent[1];
				aent += 2;
				__m128d a = _mm_load1_pd((double*)aele);
				aele++;
				__m128d s0 = _mm_load_pd((double*)&v[p]);
				__m128d t0 = _mm_load_pd((double*)&v[q]);
				__m128d e0 = _mm_load_pd((double*)&u[p]);
				__m128d f0 = _mm_load_pd((double*)&u[q]);
				f0 = _mm_mul_pd(a, f0);
				e0 = _mm_mul_pd(a, e0);
				s0 = _mm_add_pd(s0, f0);
				t0 = _mm_add_pd(t0, e0);
				_mm_store_pd((double*)&v[p], s0);
				_mm_store_pd((double*)&v[q], t0);
				p+=2;
				q+=2;
				__m128d s1 = _mm_load_pd((double*)&v[p]);
				__m128d t1 = _mm_load_pd((double*)&v[q]);
				__m128d e1 = _mm_load_pd((double*)&u[p]);
				__m128d f1 = _mm_load_pd((double*)&u[q]);
				f1 = _mm_mul_pd(a, f1);
				e1 = _mm_mul_pd(a, e1);
				s1 = _mm_add_pd(s1, f1);
				t1 = _mm_add_pd(t1, e1);
				_mm_store_pd((double*)&v[p], s1);
				_mm_store_pd((double*)&v[q], t1);
			}
		}
	}
};
#else
template<class T, class S> class hilbert_struct_symmetric_operator4 : public generic_operator<T, S>
{
private:
	const vector<T> &_aent;
	const vector<S> &_aele;
	const vector<S> &_adia;
public:
	hilbert_struct_symmetric_operator4(const vector<T> &_aent, const vector<S> &_aele, const vector<S> &_adia)
		: _aent(_aent), _aele(_aele), _adia(_adia)
	{
	}
	void operator()(const vector<S> &_u, vector<S> &_v) const
	{
		const T *__restrict aent = &_aent[0];
		const S *__restrict aele = &_aele[0], *__restrict adia = &_adia[0];
		const S *__restrict u = &_u[0];
		S *__restrict v = &_v[0];

		int dsize = (int)_adia.size();
		for(int i = 0; i < dsize; i++)
		{
			S a = *adia++;
			v[0] = a * u[0];
			v[1] = a * u[1];
			v[2] = a * u[2];
			v[3] = a * u[3];
			u+=4; v+=4;
		}
		v -= 4 * dsize;
		u -= 4 * dsize;
		int size = (int)_aele.size();
		for(int i = 0; i < size; i++)
		{
			_mm_prefetch((char*)&aent[64],_MM_HINT_NTA);
			_mm_prefetch((char*)&aele[64],_MM_HINT_NTA);
			T p = 4 * aent[0];
			T q = 4 * aent[1];
			aent += 2;
			S a = *aele++;
			S s0 = v[p+0];
			S s1 = v[p+1];
			S s2 = v[p+2];
			S s3 = v[p+3];
			S t0 = v[q+0];
			S t1 = v[q+1];
			S t2 = v[q+2];
			S t3 = v[q+3];
			s0 += a * u[q+0];
			s1 += a * u[q+1];
			s2 += a * u[q+2];
			s3 += a * u[q+3];
			t0 += a * u[p+0];
			t1 += a * u[p+1];
			t2 += a * u[p+2];
			t3 += a * u[p+3];
			v[p+0] = s0;
			v[p+1] = s1;
			v[p+2] = s2;
			v[p+3] = s3;
			v[q+0] = t0;
			v[q+1] = t1;
			v[q+2] = t2;
			v[q+3] = t3;
		}
	}
};
#endif
#endif
