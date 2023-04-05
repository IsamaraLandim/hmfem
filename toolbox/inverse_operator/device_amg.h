extern "C" {
void _device_jacobi(int *cnt, int *dsp, int *col, double *ele, int l, int n, double *f, double *u, double *v);
void _device_set_mul(double *x, double *y, double *z, double omega, unsigned int m);
void _device_add_mul(double *x, double *y, double *z, double omega, unsigned int m);
void _device_diagonal_solver(double *x, double *y, double *z, unsigned int m);
void _device_sub(double *x, double *y, unsigned int m);
void _device_add(double *x, double *y, unsigned int m);
void _device_rsb(double *x, double *y, unsigned int m);
void _device_matrix_conversion(int *_acnt, int *_adsp, int *_acol, double *_aele, int *_tdsp, int *_tcol, double *_tele, int _n);
}
template<class T, class S> class device_jacobi
{
private:
	const device_vector<T> &_acnt;
	const device_vector<T> &_adsp;
	const device_vector<T> &_acol;
	const device_vector<S> &_aele;
	const device_vector<S> &_adia;
	device_vector<S> &_v;
	const S _omega;
	device_communicator<T, S> &_com;
public:
	device_jacobi(const device_vector<T> &_acnt, const device_vector<T> &_adsp, const device_vector<T> &_acol, const device_vector<S> &_aele, const device_vector<S> &_adia, device_vector<S> &_v, const S _omega, device_communicator<T, S> &_com)
		: _acnt(_acnt), _adsp(_adsp), _acol(_acol), _aele(_aele), _adia(_adia), _v(_v), _omega(_omega), _com(_com)
	{
	}
	void operator()(const device_vector<S> &_f, device_vector<S> &_u) const
	{
		//_com.accumulate_begin(_v);
		_device_jacobi(_acnt.data(), _adsp.data(), _acol.data(), _aele.data(), _aele.size(), _acnt.size(), _f.data(), _u.data(), _v.data());
		//_com.accumulate_end(_v);
		_com.accumulate(_v);
		_device_add_mul(_u.data(), _v.data(), _adia.data(), _omega, _acnt.size());
	}
    void operator()(const device_vector<S> &_f, device_vector<S> &_u, device_vector<S> &_r)
	{
		_device_set_mul(_u.data(), _f.data(), _adia.data(), _omega, _acnt.size());
		_com.accumulate(_u);
		//_com.accumulate_begin(_r);
		_device_jacobi(_acnt.data(), _adsp.data(), _acol.data(), _aele.data(), _aele.size(), _acnt.size(), _f.data(), _u.data(), _r.data());
		//_com.accumulate_end(_r);
	}
};
#if 0
template<class T, class S> class device_matrix_conversion
{
public:
	void operator()(const vector<T> &_cnt, const vector<T> &_col, const vector<S> &_ele, device_vector<T> &dev_cnt, device_vector<T> &dev_dsp, device_vector<T> &dev_col, device_vector<S> &dev_ele)
	{
		vector<T> _dsp(_cnt.size(), 0);
		bucket_sort_offset(_cnt, _dsp);

		dev_cnt.assign(_cnt.begin(), _cnt.end());
		dev_dsp.assign(_dsp.begin(), _dsp.end());
		dev_col.assign(_col.begin(), _col.end());
		dev_ele.assign(_ele.begin(), _ele.end());
	}
};
#else
template<class T, class S> class device_matrix_conversion
{
public:
	void operator()(const vector<T> &_cnt, const vector<T> &_col, const vector<S> &_ele, device_vector<T> &dev_cnt, device_vector<T> &dev_dsp, device_vector<T> &dev_col, device_vector<S> &dev_ele)
	{
#if 1
		int size = (int)_cnt.size();
		vector<T> _dsp(size);
		bucket_sort_offset(_cnt, _dsp);
		device_vector<T> tmp_dsp(_dsp);
		device_vector<T> tmp_col(_col);
		device_vector<S> tmp_ele(_ele);
		int warp = L;
		int m = 0, max = 0;
		for(int i = 0; i < size; i++)
		{
			T cnt = _cnt[i];
			if(i % warp == 0) m += max * warp, max = 0;
			_dsp[i] = m++;
			if(max < cnt) max = cnt;
		}
		m += max * warp;
		dev_cnt.assign(_cnt.begin(), _cnt.end());
		dev_dsp.assign(_dsp.begin(), _dsp.end());
		dev_col.assign(m);
		dev_ele.assign(m);
		_device_matrix_conversion(dev_cnt.data(), dev_dsp.data(), dev_col.data(), dev_ele.data(), tmp_dsp.data(), tmp_col.data(), tmp_ele.data(), size);
		//cout << "Device: " << (m * (sizeof(T) + sizeof(S)) + 2 * sizeof(T) * size) / 1024.0 / 1024.0 << " MB" << endl;
#else
		int size = (int)_cnt.size();
		vector<T> _dsp(size);
		bucket_sort_offset(_cnt, _dsp);

		vector<T> glo_dsp(size);
		int k = 0, l = 0;
		for(int i = 0; i < size; i++)
		{
			if(i % L == 0)
			{
				k += l * L;
				l = 0;
			}
			glo_dsp[i] = k;
			k++;
			if(l < _cnt[i]) l = _cnt[i];
		}
		k += l * L;
		vector<T> glo_col(k, 0);
		vector<S> glo_ele(k, 0.0f);
		for(int i = 0; i < size; i++)
		{
			for(int j = 0; j < _cnt[i]; j++)
			{
				glo_col[glo_dsp[i] + L * j] = _col[_dsp[i] + j];
				glo_ele[glo_dsp[i] + L * j] = _ele[_dsp[i] + j];
			}
		}
		//cout << "Efficiency: " << 1.0f * glo_ele.size() / _ele.size() << endl;
		//cout << "Overhead: " << 100.0f * glo_ele.size() / _ele.size() - 100.0f << "%" << endl;

		dev_cnt.assign(_cnt.begin(), _cnt.end());
		dev_dsp.assign(glo_dsp.begin(), glo_dsp.end());
		dev_col.assign(glo_col.begin(), glo_col.end());
		dev_ele.assign(glo_ele.begin(), glo_ele.end());
#endif
	}
};
#endif
template<class T, class S> class device_diagonal_solver : public device_inverse_operator<T, S>
{
private:
	const device_vector<S> &_adia;
	device_communicator<T, S> &_com;
public:
	device_diagonal_solver(const device_vector<S> &_adia, device_communicator<T, S> &_com)
		: _adia(_adia), _com(_com)
	{
	}
	void operator()(const device_vector<S> &_u, device_vector<S> &_v)
	{
		_device_diagonal_solver(_v.data(), _u.data(), _adia.data(), _adia.size());
		_com.accumulate(_v);
	}
};

template<class T, class S> class device_amg_solver : public device_inverse_operator<T, S>
{
private:
	vector<vector<char> > _rnod;
	vector<vector<T> > _rcnt;
	vector<vector<T> > _rcol;
	vector<vector<S> > _rele;
	vector<vector<T> > _anod;
	vector<vector<T> > _acnt;
	vector<vector<T> > _acol;
	vector<vector<S> > _aele;
	vector<vector<S> > _adia;
	vector<vector<char> > _abou;

	vector<vector<S> > _u;
	vector<vector<S> > _f;
	vector<vector<S> > _s;
	vector<vector<S> > _r;
	vector<vector<S> > _v;

	vector<communicator<T, S> > _com;
	vector<S> _norm;
	int _level;
	const S _epsilon;
	const S _omega;
	const bool _gauss_seidel;
	const bool _cascadic;

#ifdef SUPERLU
    superlu_solver<T, S> *_superlu;
#else
#ifdef MUMPS
    mumps_solver<T, S> *_mumps;
#endif
#endif

	vector<device_communicator<T, S> > dev_com;
	vector<device_vector<T> > dev_rcnt;
	vector<device_vector<T> > dev_rdsp;
	vector<device_vector<T> > dev_rcol;
	vector<device_vector<S> > dev_rele;
	vector<device_vector<T> > dev_scnt;
	vector<device_vector<T> > dev_sdsp;
	vector<device_vector<T> > dev_scol;
	vector<device_vector<S> > dev_sele;
	vector<device_vector<T> > dev_acnt;
	vector<device_vector<T> > dev_adsp;
	vector<device_vector<T> > dev_acol;
	vector<device_vector<S> > dev_aele;
	vector<device_vector<S> > dev_adia;

	vector<device_vector<S> > dev_u;
	vector<device_vector<S> > dev_f;
	vector<device_vector<S> > dev_s;
	vector<device_vector<S> > dev_r;
	vector<device_vector<S> > dev_v;

public:
	device_amg_solver(const vector<T> &_anod0, const vector<T> &_acnt0, const vector<T> &_acol0, const vector<S> &_aele0, 
		const int max_level, const int max_size, const S _epsilon, const S _omega, const bool _gauss_seidel = false, const bool _cascadic = false)
		: _epsilon(_epsilon), _omega(_omega), _gauss_seidel(_gauss_seidel), _cascadic(_cascadic)
	{
		_rnod.resize(max_level);
		_rcnt.resize(max_level);
		_rcol.resize(max_level);
		_rele.resize(max_level);
		_anod.resize(max_level);
		_acnt.resize(max_level);
		_acol.resize(max_level);
		_aele.resize(max_level);
		_adia.resize(max_level);
		_abou.resize(max_level);

		_u.resize(max_level);
		_f.resize(max_level);
		_s.resize(max_level);
		_r.resize(max_level);
		_v.resize(max_level);
		_com.resize(max_level);
		_norm.resize(max_level);

		dev_com.resize(max_level);
		dev_rcnt.resize(max_level);
		dev_rdsp.resize(max_level);
		dev_rcol.resize(max_level);
		dev_rele.resize(max_level);
		dev_scnt.resize(max_level);
		dev_sdsp.resize(max_level);
		dev_scol.resize(max_level);
		dev_sele.resize(max_level);
		dev_acnt.resize(max_level);
		dev_adsp.resize(max_level);
		dev_acol.resize(max_level);
		dev_aele.resize(max_level);
		dev_adia.resize(max_level);

		dev_u.resize(max_level);
		dev_f.resize(max_level);
		dev_s.resize(max_level);
		dev_r.resize(max_level);
		dev_v.resize(max_level);

#if 1
		_anod[0].assign(_anod0.begin(), _anod0.end());
		_acnt[0].assign(_acnt0.begin(), _acnt0.end());
		_acol[0].assign(_acol0.begin(), _acol0.end());
		_aele[0].assign(_aele0.begin(), _aele0.end());

		device_matrix_conversion<T, S> conversion;
		conversion(_acnt[0], _acol[0], _aele[0], dev_acnt[0], dev_adsp[0], dev_acol[0], dev_aele[0]);
#else
		_anod[0] = _anod0;
		_acnt[0] = _acnt0;
		_acol[0] = _acol0;
		_aele[0] = _aele0;
#endif
		int rank;
		network::rank(rank);

		amg_restriction<T, S> amg_restriction;
		int size = (int)_acnt[0].size();
		_u[0].resize(size, 0.0);
		_f[0].resize(size, 0.0);
		_s[0].resize(size, 0.0);
		_r[0].resize(size, 0.0);
		_v[0].resize(size, 0.0);

		dev_u[0].assign(_u[0].begin(), _u[0].end());
		dev_f[0].assign(_f[0].begin(), _f[0].end());
		dev_s[0].assign(_s[0].begin(), _s[0].end());
		dev_r[0].assign(_r[0].begin(), _r[0].end());
		dev_v[0].assign(_v[0].begin(), _v[0].end());

		int l = 0;
		int gsize;
		network::size(gsize);
		vector<int> _gcnt(gsize, 0);
		network::allgather(&size, sizeof(T), &_gcnt[0], sizeof(T));
		int sum = 0;
		for(int i = 0; i < (int)_gcnt.size(); i++) sum += _gcnt[i];
		int gcnt = sum / gsize;
		double rel = 0.0;
		while(l < max_level - 1 && gcnt > max_size && rel < 0.95)
		//while(l < max_level - 1 && gcnt > max_size)
		{
			_com[l].assign(_anod[l]);
			dev_com[l].assign(_anod[l]);
			amg_restriction(_epsilon, _anod[l], _acnt[l], _acol[l], _aele[l], _adia[l], _abou[l], _com[l], _rnod[l], _rcnt[l], _rcol[l], _rele[l], _anod[l+1], _acnt[l+1], _acol[l+1], _aele[l+1]);
			//rel = (float)_acnt[l+1].size()/(float)_acnt[l].size();

			conversion(_acnt[l+1], _acol[l+1], _aele[l+1], dev_acnt[l+1], dev_adsp[l+1], dev_acol[l+1], dev_aele[l+1]);
			conversion(_rcnt[l], _rcol[l], _rele[l], dev_rcnt[l], dev_rdsp[l], dev_rcol[l], dev_rele[l]);
			////
			vector<T> _row(_rcol[l]), _col(_rcol[l]);
			vector<S> _ele(_rele[l]);
			for(int i = 0, k = 0; i < (int)_rcnt[l].size(); i++) for(int j = 0; j < _rcnt[l][i]; j++, k++) _col[k] = i;
			binary_sort_sort_copy(_row, _col, _ele);
			vector<T> _cnt(_acnt[l+1].size(), 0);
			for(int i = 0; i < (int)_row.size(); i++) _cnt[_row[i]]++;
			////
			conversion(_cnt, _col, _ele, dev_scnt[l], dev_sdsp[l], dev_scol[l], dev_sele[l]);

			for(int i = 0; i < (int)_adia[l].size(); i++) _adia[l][i] = 1.0 / _adia[l][i];
			//cout << "Rank: " << rank << " Size: (" << _acol[l].size() << " | " << _acnt[l].size() << ")" << endl;
			dev_adia[l].assign(_adia[l].begin(), _adia[l].end());
			l++;
			size = (int)_acnt[l].size();
			_u[l].resize(size, 0.0);
			_f[l].resize(size, 0.0);
			_s[l].resize(size, 0.0);
			_r[l].resize(size, 0.0);
			_v[l].resize(size, 0.0);

			dev_u[l].assign(_u[l].begin(), _u[l].end());
			dev_f[l].assign(_f[l].begin(), _f[l].end());
			dev_s[l].assign(_s[l].begin(), _s[l].end());
			dev_r[l].assign(_r[l].begin(), _r[l].end());
			dev_v[l].assign(_v[l].begin(), _v[l].end());

			network::allgather(&size, sizeof(T), &_gcnt[0], sizeof(T));
			int sum = 0;
			for(int i = 0; i < (int)_gcnt.size(); i++) sum += _gcnt[i];
			int gold = gcnt;
			gcnt = sum / gsize;
			rel = (double)gcnt / (double)gold;
			if(rank == 0) cout << "Rank: " << rank << " Relative: " << rel << " Nodes: " << sum << " Level: " << l << endl;
		}
		_com[l].assign(_anod[l]);
		dev_com[l].assign(_anod[l]);
		_adia[l].resize(_acnt[l].size());
		for(int t = 0, i = 0; i < (int)_acnt[l].size(); i++)
		{
			for(int j = 0; j < _acnt[l][i]; j++, t++)
			{
				if(_acol[l][t] == i) _adia[l][i] = _aele[l][t];
			}
		}
		_com[l].accumulate(_adia[l]);
		for(int i = 0; i < (int)_adia[l].size(); i++) _adia[l][i] = 1.0 / _adia[l][i];
		//cout << "Rank: " << rank << " Size: (" << _acol[l].size() << " | " << _acnt[l].size() << ")" << endl;
		dev_adia[l].assign(_adia[l].begin(), _adia[l].end());
		_level = l;
		//cout << "Rank: " << rank << " Level: " << _level << endl;

#ifdef SUPERLU
        _superlu = new superlu_solver<T, S>(_anod[l], _acnt[l], _acol[l], _aele[l]);
#else
#ifdef MUMPS
        _mumps = new mumps_solver<T, S>(_anod[l], _acnt[l], _acol[l], _aele[l]);
#endif
#endif

	}
	void sub(const device_vector<S> &_a, device_vector<S> &_b)
	{
		_device_sub(_a.data(), _b.data(), _a.size());
	}
	void rsb(const device_vector<S> &_a, device_vector<S> &_b)
	{
		_device_rsb(_a.data(), _b.data(), _a.size());
	}
	void add(const device_vector<S> &_a, device_vector<S> &_b)
	{
		_device_add(_a.data(), _b.data(), _a.size());
	}
	void operator()(const device_vector<S> &_f0, device_vector<S> &_u0)
	{
		unsigned int size = _f0.size();
		CUDA_SAFE_CALL(cudaMemcpy(dev_f[0].data(), _f0.data(), size * sizeof(S), cudaMemcpyDeviceToDevice) );
		multigrid(0);
		CUDA_SAFE_CALL(cudaMemcpy(_u0.data(), dev_u[0].data(), size * sizeof(S), cudaMemcpyDeviceToDevice) );
	}
	void multigrid(int _k)
	{
		int l = _k++;
		S delta = _omega;

		if(l < _level)
		{
			device_linear_operator<T, S> _K(dev_acnt[l], dev_adsp[l], dev_acol[l], dev_aele[l]);
			device_jacobi<T, S> _J(dev_acnt[l], dev_adsp[l], dev_acol[l], dev_aele[l], dev_adia[l], dev_v[l], delta, dev_com[l]);
			device_linear_operator<T, S> _P(dev_rcnt[l], dev_rdsp[l], dev_rcol[l], dev_rele[l]);
			device_linear_operator<T, S> _R(dev_scnt[l], dev_sdsp[l], dev_scol[l], dev_sele[l]);
			if(_cascadic)
			{
				//if(_gauss_seidel) cout << "NOT IMPLEMENTED ON GPUS!" << endl;
				_R(dev_f[l], dev_f[l + 1]);
				multigrid(_k);
				_P(dev_u[l + 1], dev_u[l]);
				_J(dev_f[l], dev_u[l]);
				_J(dev_f[l], dev_u[l]);
			}
			else
			{
				//if(_gauss_seidel) cout << "NOT IMPLEMENTED ON GPUS!" << endl;
#if 1
				_J(dev_f[l], dev_u[l], dev_r[l]);
#else
				CUDA_SAFE_CALL(cudaMemset(dev_u[l].data(), 0, dev_u[l].size() * sizeof(S)) );
				_J(dev_f[l], dev_u[l]);
				_K(dev_u[l], dev_r[l]);
				rsb(dev_f[l], dev_r[l]);
#endif
				_R(dev_r[l], dev_f[l + 1]);
				multigrid(_k);
				_P(dev_u[l + 1], dev_s[l]);
				add(dev_s[l], dev_u[l]);
				_J(dev_f[l], dev_u[l]);
			}

		}
		else
		{
#ifdef SUPERLU
            (*_superlu)(dev_f[l], dev_u[l]);
#else
#ifdef MUMPS
            (*_mumps)(dev_f[l], dev_u[l]);
#else
			device_diagonal_solver<T, S> _D(dev_adia[l], dev_com[l]);
			_D(dev_f[l], dev_u[l]);
#endif
#endif
		}


	}
};

