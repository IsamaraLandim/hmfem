#ifdef SUPERLUMT
#include <pdsp_defs.h>
#else
#include <slu_ddefs.h>
#endif
template<class T, class S> class superlu_solver : public inverse_operator<T, S>
{
private:
	fact_t _fact;
    trans_t _trans;
    SuperMatrix _A, _P, _L, _U, _B;
    int _info, _nrhs, _n;
#ifdef SUPERLUMT
    superlumt_options_t _options;
	Gstat_t _gstat;
#else
    SuperLUStat_t _stat;
    superlu_options_t _options;
#endif
    vector<T> _perm_c;
    vector<T> _perm_r;
    vector<T> _etree;

    vector<T> _enod;
    vector<T> _emap;
    vector<T> _ecnt;
    vector<T> _edsp;
    vector<S> _e;
    vector<S> _r;
	vector<S> _s;

    vector<T> _sptr;
    vector<T> _srow;
    vector<T> _scol;
    vector<S> _sele;

    vector<T> _tptr;
    vector<T> _trow;
    vector<T> _tcol;
    vector<S> _tele;

public:
	superlu_solver(const vector<T> &_anod, const vector<T> &_acnt, const vector<T> &_acol, const vector<S> &_aele)
    {
		int rank, size;
		network::rank(rank);
		network::size(size);

#ifdef NOMPI

        _info = 0;
        _nrhs = 1;
        _n = _acnt.size();
        _scol.assign(_acol.begin(), _acol.end());
        _sele.assign(_aele.begin(), _aele.end());
		_sptr.resize(_n + 1);
		_sptr[0] = 0;
		for(int t = 0, i = 0; i < _n; i++) t += _acnt[i], _sptr[i + 1] = t;
        dCreate_CompCol_Matrix(&_A, _n, _n, _sele.size(), _sele.data(), _scol.data(), _sptr.data(), SLU_NC, SLU_D, SLU_GE);
		_r.resize(_anod.size());
#else
		_s.resize(_anod.size());
        vector<T> _snod(_anod);
        int ssize = _snod.size();
        _ecnt.resize(size, 0);
        _edsp.resize(size);
        MPI_Gather(&ssize, sizeof(int), MPI_BYTE, _ecnt.data(), sizeof(int), MPI_BYTE, 0, MPI_COMM_WORLD);
        bucket_sort_offset(_ecnt, _edsp);
        int esize = bucket_sort_size(_ecnt);
        _emap.resize(esize);

        for(int i = 0; i < size; i++) _ecnt[i] *= sizeof(T), _edsp[i] *= sizeof(T);
        MPI_Gatherv(_snod.data(), ssize * sizeof(T), MPI_BYTE, _emap.data(), _ecnt.data(), _edsp.data(), MPI_BYTE, 0, MPI_COMM_WORLD);
        for(int i = 0; i < size; i++) _ecnt[i] /= sizeof(T), _edsp[i] /= sizeof(T);

        if(rank == 0)
        {
            _enod.assign(_emap.begin(), _emap.end());

            vector<T> _eloc(esize);
            for(int i = 0; i < esize; i++) _eloc[i] = i;
            binary_sort_copy(_enod, _eloc);
            for(int k = 0, i = 0; i < esize; k++)
            {
                T nod = _enod[i];
                while(i < esize && _enod[i] == nod) _emap[_eloc[i]] = k, i++;
            }
            unique_resize(_enod);
            _r.resize(_enod.size());
            _e.resize(_emap.size());
        }

        for(int i = 0; i < size; i++) _ecnt[i] *= sizeof(T), _edsp[i] *= sizeof(T);
        MPI_Scatterv(_emap.data(), _ecnt.data(), _edsp.data(), MPI_BYTE, _snod.data(), ssize * sizeof(T), MPI_BYTE, 0, MPI_COMM_WORLD);
        for(int i = 0; i < size; i++) _ecnt[i] /= sizeof(T), _edsp[i] /= sizeof(T);
        for(int i = 0; i < size; i++) _ecnt[i] *= sizeof(S), _edsp[i] *= sizeof(S);

        int asize = _aele.size();
        _srow.resize(asize);
        for(int k = 0, i = 0; i < _acnt.size(); i++) for(int j = 0; j < _acnt[i]; j++, k++) _srow[k] = _snod[i];
        _scol.resize(asize);
        vector<T> _tnod(_acol);
        vector<T> _tloc(asize);
        for(int i = 0; i < asize; i++) _tloc[i] = i;
        binary_sort_copy(_tnod, _tloc);
        for(int k = 0, i = 0; i < asize; k++)
        {
            T nod = _tnod[i];
            while(i < asize && _tnod[i] == nod) _scol[_tloc[i]] = _snod[k], i++;
        }

        vector<T> _tcnt(size, 0);
        vector<T> _tdsp(size);
        MPI_Gather(&asize, sizeof(int), MPI_BYTE, _tcnt.data(), sizeof(int), MPI_BYTE, 0, MPI_COMM_WORLD);
        bucket_sort_offset(_tcnt, _tdsp);
        int tsize = bucket_sort_size(_tcnt);

        if(rank == 0)
        {
            _trow.resize(tsize);
            _tcol.resize(tsize);
            _tele.resize(tsize);
        }
        for(int i = 0; i < size; i++) _tcnt[i] *= sizeof(T), _tdsp[i] *= sizeof(T);
        MPI_Gatherv(_srow.data(), asize * sizeof(T), MPI_BYTE, _trow.data(), _tcnt.data(), _tdsp.data(), MPI_BYTE, 0, MPI_COMM_WORLD);
        MPI_Gatherv(_scol.data(), asize * sizeof(T), MPI_BYTE, _tcol.data(), _tcnt.data(), _tdsp.data(), MPI_BYTE, 0, MPI_COMM_WORLD);
        for(int i = 0; i < size; i++) _tcnt[i] /= sizeof(T), _tdsp[i] /= sizeof(T);
        for(int i = 0; i < size; i++) _tcnt[i] *= sizeof(S), _tdsp[i] *= sizeof(S);
        MPI_Gatherv(_aele.data(), asize * sizeof(S), MPI_BYTE, _tele.data(), _tcnt.data(), _tdsp.data(), MPI_BYTE, 0, MPI_COMM_WORLD);
        for(int i = 0; i < size; i++) _tcnt[i] /= sizeof(S), _tdsp[i] /= sizeof(S);

        if(rank == 0)
        {
            binary_sort_sort_copy(_trow, _tcol, _tele);
            unique_accumulate(_trow, _tcol, _tele);

            _info = 0;
            _nrhs = 1;
            _n = _enod.size();
			int psize = _tele.size();
			vector<T> _pcnt(psize, 0);
			for(int i = 0; i < psize; i++) _pcnt[_trow[i]]++;
            _tptr.resize(_n + 1);
            _tptr[0] = 0;
            for(int t = 0, i = 0; i < _n; i++) t += _pcnt[i], _tptr[i + 1] = t;
            dCreate_CompCol_Matrix(&_A, _n, _n, _tele.size(), _tele.data(), _tcol.data(), _tptr.data(), SLU_NC, SLU_D, SLU_GE);
		}
#endif

		if(rank == 0)
		{
            dCreate_Dense_Matrix(&_B, _n, _nrhs, _r.data(), _n, SLU_DN, SLU_D, SLU_GE);

            _perm_c.resize(_n);
            _perm_r.resize(_n);
            _etree.resize(_n);

			int _panel_size = sp_ienv(1);
			int _relax = sp_ienv(2);
#ifdef SUPERLUMT
			int _nprocs = 2;
			_fact               = EQUILIBRATE;
			yes_no_t _refact             = NO;
			_trans              = NOTRANS;
			_info = 0;

			get_perm_c(3, &_A, _perm_c.data());
//			sp_preorder(&_options, &_A, _perm_c.data(), _etree.data(), &_P);

			StatAlloc(_n, _nprocs, _panel_size, _relax, &_gstat);
			StatInit(_n, _nprocs, &_gstat);

			pdgstrf_init(_nprocs, _fact, _trans, _refact, _panel_size, _relax, 1.0, NO, 0.0, _perm_c.data(), _perm_r.data(), NULL, 0, &_A, &_P, &_options, &_gstat);
			pdgstrf(&_options, &_P, _perm_r.data(), &_L, &_U, &_gstat, &_info);
			//dgstrs(_trans, &_L, &_U, _perm_r.data(), _perm_c.data(), &_B, &_gstat, &_info);

			//pxgstrf_finalize(&_options, &_P);
			PrintStat(&_gstat);
			StatFree(&_gstat);

			//pdgstrf(&_options, &_P, 0.0, _relax, _panel_size, _etree.data(), NULL, 0, _perm_c.data(), _perm_r.data(), &_L, &_U, &_stat, &_info);
#else
			StatInit(&_stat);
            set_default_options(&_options);
            _options.SymmetricMode = YES;
			get_perm_c(3, &_A, _perm_c.data());
			sp_preorder(&_options, &_A, _perm_c.data(), _etree.data(), &_P);
			dgstrf(&_options, &_P, 0.0, _relax, _panel_size, _etree.data(), NULL, 0, _perm_c.data(), _perm_r.data(), &_L, &_U, &_stat, &_info);
			StatFree(&_stat);
#endif
			//Destroy_SuperMatrix_Store(&_A);
			//Destroy_CompCol_Permuted(&_P);
            cout << "Rank: " << rank << " Coarse System: " << _n << endl;
		}
	}

	~superlu_solver()
	{
		int rank;
		network::rank(rank);
		if(rank == 0)
		{
//		Destroy_CompCol_Matrix(&_A);
//		Destroy_SuperMatrix_Store(&_B);
//		Destroy_SuperNode_Matrix(&_L);
//		Destroy_CompCol_Matrix(&_U);
		}
#ifdef SUPERLU
#endif
	}
	void operator()(const vector<S> &_f, vector<S> &_u)
	{
		_trans = NOTRANS;
#ifdef NOMPI
        for(int i = 0; i < _f.size(); i++) _r[i] = _f[i];
#ifdef SUPERLUMT
		dgstrs(_trans, &_L, &_U, _perm_c.data(), _perm_r.data(), &_B, &_gstat, &_info);
#else
		StatInit(&_stat);
        dgstrs(_trans, &_L, &_U, _perm_c.data(), _perm_r.data(), &_B, &_stat, &_info);
		StatFree(&_stat);
#endif
        for(int i = 0; i < _u.size(); i++) _u[i] = _r[i];
#else
        int rank;
        network::rank(rank);
        MPI_Gatherv((S*)_f.data(), _f.size() * sizeof(S), MPI_BYTE, _e.data(), _ecnt.data(), _edsp.data(), MPI_BYTE, 0, MPI_COMM_WORLD);
        if(rank == 0)
        {
			for(int i = 0; i < _r.size(); i++) _r[i] = 0.0;
			for(int i = 0; i < _e.size(); i++) _r[_emap[i]] += _e[i];
#ifdef SUPERLUMT
            dgstrs(_trans, &_L, &_U, _perm_c.data(), _perm_r.data(), &_B, &_gstat, &_info);
#else
			StatInit(&_stat);
            dgstrs(_trans, &_L, &_U, _perm_c.data(), _perm_r.data(), &_B, &_stat, &_info);
			StatFree(&_stat);
#endif
			for(int i = 0; i < _e.size(); i++) _e[i] = _r[_emap[i]];
        }
        MPI_Scatterv(_e.data(), _ecnt.data(), _edsp.data(), MPI_BYTE, _u.data(), _u.size() * sizeof(S), MPI_BYTE, 0, MPI_COMM_WORLD);
#endif
	}
#ifdef CUDA
	void operator()(const device_vector<S> &_f, device_vector<S> &_u)
	{
		_trans = NOTRANS;
#ifdef NOMPI
		_f.download(_r.data(), _r.size());
#ifdef SUPERLUMT
		dgstrs(_trans, &_L, &_U, _perm_c.data(), _perm_r.data(), &_B, &_gstat, &_info);
#else
		StatInit(&_stat);
        dgstrs(_trans, &_L, &_U, _perm_c.data(), _perm_r.data(), &_B, &_stat, &_info);
		StatFree(&_stat);
#endif
		_u.upload(_r.data(), _r.size());
#else
        int rank;
        network::rank(rank);
		_f.download(_s.data(), _s.size());
        MPI_Gatherv((S*)_s.data(), _s.size() * sizeof(S), MPI_BYTE, _e.data(), _ecnt.data(), _edsp.data(), MPI_BYTE, 0, MPI_COMM_WORLD);
        if(rank == 0)
        {
			for(int i = 0; i < _r.size(); i++) _r[i] = 0.0;
			for(int i = 0; i < _e.size(); i++) _r[_emap[i]] += _e[i];
#ifdef SUPERLUMT
            dgstrs(_trans, &_L, &_U, _perm_c.data(), _perm_r.data(), &_B, &_gstat, &_info);
#else
			StatInit(&_stat);
            dgstrs(_trans, &_L, &_U, _perm_c.data(), _perm_r.data(), &_B, &_stat, &_info);
			StatFree(&_stat);
#endif
			for(int i = 0; i < _e.size(); i++) _e[i] = _r[_emap[i]];
        }
        MPI_Scatterv(_e.data(), _ecnt.data(), _edsp.data(), MPI_BYTE, _s.data(), _s.size() * sizeof(S), MPI_BYTE, 0, MPI_COMM_WORLD);
		_u.upload(_s.data(), _s.size());
#endif
	}
#endif
};
