#include <dmumps_c.h>
template<class T, class S> class mumps_solver : public inverse_operator<T, S>
{
private:
	DMUMPS_STRUC_C _id;

	vector<T> _enod;
	vector<T> _emap;
	vector<T> _ecnt;
	vector<T> _edsp;
	vector<S> _e;
	vector<S> _r;
	vector<S> _s;

	vector<T> _srow;
	vector<T> _scol;
	vector<S> _sele;

	vector<T> _trow;
	vector<T> _tcol;
	vector<S> _tele;

public:
	mumps_solver(const vector<T> &_anod, const vector<T> &_acnt, const vector<T> &_acol, const vector<S> &_aele)
	{
		int rank, size;
		network::rank(rank);
		network::size(size);

		_id.comm_fortran = -987654; //Magic number
		_id.par = 1; //Host participates in computation
		_id.sym = 0; //0,1,2
		_id.job = -1;
		dmumps_c(&_id);
		_id.icntl[1 - 1] = -1;
		_id.icntl[2 - 1] = -1;
		_id.icntl[3 - 1] = -1;
		_id.icntl[4 - 1] = 0;

#ifdef NOMPI
		_srow.resize(_acol.size());
		_scol.assign(_acol.begin(), _acol.end());
		_sele.assign(_aele.begin(), _aele.end());
		for(int k = 0, i = 0; i < _acnt.size(); i++)
			for(int j = 0; j < _acnt[i]; j++, k++) _srow[k] = i + 1, _scol[k] = _acol[k] + 1;
		_id.n = _acnt.size();
		_id.nz = _sele.size();
		_id.irn = _srow.data();
		_id.jcn = _scol.data();
		_id.a = _sele.data();
		_id.job = 4; //Analysis and factorization
		dmumps_c(&_id);
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

		for(int i = 0; i < _snod.size(); i++) _snod[i] += 1;//Fortran!

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
#if 0
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
			_id.n = _enod.size();
			_id.nz = _tele.size();
			_id.irn = _trow.data();
			_id.jcn = _tcol.data();
			_id.a = _tele.data();
		}
		_id.job = 4; //Analysis and factorization
		dmumps_c(&_id);

#else
		_sele.assign(_aele.begin(), _aele.end());

		_id.n = _enod.size();
		_id.nz_loc = _sele.size();
		_id.irn_loc = _srow.data();
		_id.jcn_loc = _scol.data();
		_id.a_loc = _sele.data();
		_id.icntl[18 - 1] = 3;//Distributed matrix
		_id.job = 4; //Analysis and factorization
		dmumps_c(&_id);
#endif
#endif
		if(rank == 0) cout << "Rank: " << rank << " Coarse System: " << _id.n << endl;
	}
	~mumps_solver()
	{
		_id.job = -2;
		dmumps_c(&_id);
	}
	void operator()(const vector<S> &_f, vector<S> &_u)
	{
#ifdef NOMPI
		_u.assign(_f.begin(), _f.end());
		_id.rhs = _u.data();
		_id.job = 3; //Solve
		dmumps_c(&_id);
#else
		int rank;
		network::rank(rank);
		MPI_Gatherv((S*)_f.data(), _f.size() * sizeof(S), MPI_BYTE, _e.data(), _ecnt.data(), _edsp.data(), MPI_BYTE, 0, MPI_COMM_WORLD);
		if(rank == 0)
		{
			for(int i = 0; i < _r.size(); i++) _r[i] = 0.0;
			for(int i = 0; i < _e.size(); i++) _r[_emap[i]] += _e[i];
			_id.rhs = _r.data();
		}
		_id.job = 3; //Solve
		dmumps_c(&_id);
		if(rank == 0) for(int i = 0; i < _e.size(); i++) _e[i] = _r[_emap[i]];
		MPI_Scatterv(_e.data(), _ecnt.data(), _edsp.data(), MPI_BYTE, _u.data(), _u.size() * sizeof(S), MPI_BYTE, 0, MPI_COMM_WORLD);
#endif
	}
#ifdef CUDA
	void operator()(const device_vector<S> &_f, device_vector<S> &_u)
	{
#ifdef NOMPI
		_f.download(_r.data(), _r.size());
		_id.rhs = _r.data();
		_id.job = 3; //Solve
		dmumps_c(&_id);
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
			_id.rhs = _r.data();
		}
		_id.job = 3; //Solve
		dmumps_c(&_id);
		if(rank == 0) for(int i = 0; i < _e.size(); i++) _e[i] = _r[_emap[i]];
		MPI_Scatterv(_e.data(), _ecnt.data(), _edsp.data(), MPI_BYTE, _s.data(), _s.size() * sizeof(S), MPI_BYTE, 0, MPI_COMM_WORLD);
		_u.upload(_s.data(), _s.size());
#endif
	}
#endif
};
