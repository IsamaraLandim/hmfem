template<class T> class node_selection
{
public:
	void operator()(const char _sel, const vector<T> &_acnt, const vector<T> &_adsp, const vector<T> &_acol, const vector<char> &_acox, vector<char> &_asel, vector<char> &_abou, const bool _bou)
	{
		int anodes = (int)_acnt.size();
		for(int i = 0; i < anodes; i++)
		{
			if((_bou || _abou[i] != 0) && _asel[i] == _sel)
			{
				if(_sel == 0) _asel[i] = 1; // master node
				T cnt = _acnt[i], dsp = _adsp[i];
				for(int j = dsp; j < dsp + cnt; j++)
				{
					T col = _acol[j];
					if(_asel[col] == 0)
					{
						//
						T acnt = _acnt[col], adsp = _adsp[col];
						for(int k = adsp; k < adsp + acnt; k++)
						{
							if(_acol[k] == i && _acox[k] != 0) _asel[col] = 2, k = adsp + acnt; // slave node
						}
						//
						//if(_acox[j] != 0) _asel[col] = 2; // slave node
					}
				}

			}
		}
	}
};

template<class T, class S> class parallel_coarsening
{
public:
	void operator()(const vector<T> &_anod, const vector<T> &_acnt, const vector<T> &_acol, const vector<S> &_aele, vector<S> &_adia, vector<char> &_abou, vector<char> &_asel, vector<char> &_acox, communicator<T, S> &_acom, const S _epsilon)
	{
		int size, rank;
		network::size(size);
		network::rank(rank);

		int anodes = (int)_anod.size();
		vector<T> _adsp(anodes);
		bucket_sort_offset(_acnt, _adsp);
		_abou.resize(anodes, 0);
		_asel.resize(anodes, 0);
		_acom.boundary(_abou);

		vector<T> _rcbn, _rloc, _rcnt, _rdsp;
		_acom.boundary_nodes(_rcbn, _rloc, _rcnt, _rdsp);
		vector<T> _rmap(anodes, -1);

		vector<T> _scnt(size, 0);
		vector<T> _sdsp(size);
		for(int i = 0; i < size; i++)
		{
			T rcnt = _rcnt[i], rdsp = _rdsp[i];
			for(int j = rdsp; j < rdsp + rcnt; j++) _rmap[_rcbn[j]] = _rloc[j];
			for(int j = rdsp; j < rdsp + rcnt; j++)
			{
				T l = _rcbn[j];
				T acnt = _acnt[l], adsp = _adsp[l];
				for(int k = adsp; k < adsp + acnt; k++)
				{
					if(_rmap[_acol[k]] != -1) _scnt[i]++;
				}
			}
			for(int j = rdsp; j < rdsp + rcnt; j++) _rmap[_rcbn[j]] = -1;
		}
		bucket_sort_offset(_scnt, _sdsp);
		int ssize = bucket_sort_size(_scnt);
		vector<T> _scol(ssize);
		vector<T> _srow(ssize);
		vector<S> _sele(ssize);
		for(int m = 0, i = 0; i < size; i++)
		{
			T rcnt = _rcnt[i], rdsp = _rdsp[i];
			for(int j = rdsp; j < rdsp + rcnt; j++) _rmap[_rcbn[j]] = _rloc[j];
			for(int j = rdsp; j < rdsp + rcnt; j++)
			{
				T l = _rcbn[j];
				T acnt = _acnt[l], adsp = _adsp[l];
				for(int k = adsp; k < adsp + acnt; k++)
				{
					if(_rmap[_acol[k]] != -1)
					{
						_scol[m] = _rmap[_acol[k]];
						_srow[m] = _rmap[l];
						_sele[m] = _aele[k];
						m++;
					}
				}
			}
			for(int j = rdsp; j < rdsp + rcnt; j++) _rmap[_rcbn[j]] = -1;
		}
		vector<T> _qcnt(size, 0);
		vector<T> _qdsp(size);
		network::alltoall(&_scnt[0], sizeof(T), &_qcnt[0], sizeof(T));
		bucket_sort_offset(_qcnt, _qdsp);
		int qsize = bucket_sort_size(_qcnt);
		vector<T> _brow(qsize);
		vector<T> _bcol(qsize);
		vector<S> _bele(qsize);
		for(int i = 0; i < size; i++) _scnt[i] *= sizeof(T), _sdsp[i] *= sizeof(T);
		for(int i = 0; i < size; i++) _qcnt[i] *= sizeof(T), _qdsp[i] *= sizeof(T);
		network::alltoallv(_scol.data(), _scnt.data(), _sdsp.data(), _bcol.data(), _qcnt.data(), _qdsp.data());
		network::alltoallv(_srow.data(), _scnt.data(), _sdsp.data(), _brow.data(), _qcnt.data(), _qdsp.data());
		for(int i = 0; i < size; i++) _scnt[i] /= sizeof(T), _sdsp[i] /= sizeof(T);
		for(int i = 0; i < size; i++) _qcnt[i] /= sizeof(T), _qdsp[i] /= sizeof(T);
		for(int i = 0; i < size; i++) _scnt[i] *= sizeof(S), _sdsp[i] *= sizeof(S);
		for(int i = 0; i < size; i++) _qcnt[i] *= sizeof(S), _qdsp[i] *= sizeof(S);
		network::alltoallv(_sele.data(), _scnt.data(), _sdsp.data(), _bele.data(), _qcnt.data(), _qdsp.data());
		for(int i = 0; i < size; i++) _scnt[i] /= sizeof(S), _sdsp[i] /= sizeof(S);
		for(int i = 0; i < size; i++) _qcnt[i] /= sizeof(S), _qdsp[i] /= sizeof(S);
		for(int i = 0; i < size; i++)
		{
			T cnt = _qcnt[i], dsp = _qdsp[i];
			for(int j = dsp; j < dsp + cnt; j++) _bcol[j] = (_bcol[j] << 8) | i;
		}
		binary_sort_sort_copy(_brow, _bcol, _bele);
		for(int i = 0; i < qsize; i++) _bcol[i] = _bcol[i] >> 8;
		unique_accumulate(_brow, _bcol, _bele);

		int bsize = (int)_bele.size();
		vector<T> _bcnt(anodes, 0);
		vector<T> _bdsp(anodes);
		for(int i = 0; i < bsize; i++) _bcnt[_brow[i]]++;
		bucket_sort_offset(_bcnt, _bdsp);

// Phase (a)
		node_selection<T> selection;

		_adia.resize(anodes, S(0.0));
		for(int i = 0; i < anodes; i++)
		{
			if(_abou[i] != 0)
			{
				T cnt = _bcnt[i], dsp = _bdsp[i];
				for(int j = dsp; j < dsp + cnt; j++) if(i == _bcol[j]) _adia[i] = _bele[j];
			}
			else
			{
				T cnt = _acnt[i], dsp = _adsp[i];
				for(int j = dsp; j < dsp + cnt; j++) if(i == _acol[j]) _adia[i] = _aele[j];
			}
		}
		vector<S> _asda(anodes);
		for(int i = 0; i < anodes; i++) _asda[i] = fabs(_epsilon * _adia[i]);
//
		vector<char> _bcox(bsize, 0);
		for(int i = 0; i < anodes; i++)
		{
			T cnt = _bcnt[i], dsp = _bdsp[i];
			for(int j = dsp; j < dsp + cnt; j++)
			{
				//if(fabs(_bele[j]) >= _asda[_bcol[j]]) _bcox[j] = 1; // strong connection
				if(fabs(_bele[j]) >= _asda[i]) _bcox[j] = 1; // strong connection
			}
		}
		int asize = (int)_aele.size();
		_acox.resize(asize, 0);
		for(int i = 0; i < anodes; i++)
		{
			T cnt = _acnt[i], dsp = _adsp[i];
			for(int j = dsp; j < dsp + cnt; j++)
			{
				//if(fabs(_aele[j]) >= _asda[_acol[j]]) _acox[j] = 1; // strong connection
				if(fabs(_aele[j]) >= _asda[i]) _acox[j] = 1; // strong connection
			}
		}
		for(int i = 0; i < anodes; i++)
		{
			if(_abou[i] != 0)
			{
				T cnt = _acnt[i], dsp = _adsp[i];
				for(int j = dsp; j < dsp + cnt; j++)
				{
					T bcnt = _bcnt[i], bdsp = _bdsp[i];
					for(int k = bdsp; k < bdsp + bcnt; k++)
					{
						if(_acol[j] == _bcol[k]) _acox[j] = _bcox[k], k = bdsp + bcnt;
					}
				}
			}
		}
//
		vector<char> _rsel(_rcbn.size(), 0);
		for(int n = 0; n < size; n++)
		{
			if(n == rank)
			{
				selection(1, _bcnt, _bdsp, _bcol, _bcox, _asel, _abou, false);
				selection(0, _bcnt, _bdsp, _bcol, _bcox, _asel, _abou, false);

				for(int i = n + 1; i < size; i++)
				{
					T rcnt = _rcnt[i], rdsp = _rdsp[i];
					for(int j = rdsp; j < rdsp + rcnt; j++)
					{
						_rsel[j] = _asel[_rcbn[j]];
					}
				}
			}
			_rcnt[n] = 0;
			T vsize = 0;
			network::scatter(_rcnt.data(), sizeof(T), &vsize, sizeof(T), n);
			vector<char> _vsel(vsize);
			network::scatterv(_rsel.data(), _rcnt.data(), _rdsp.data(), _vsel.data(), vsize, n);
			if(n < rank)
			{
				for(int i = 0; i < vsize; i++)
				{
					if(_vsel[i] != 0) _asel[_rcbn[i + _rdsp[n]]] = _vsel[i];
				}
			}
		}

// Phase (b)
		selection(1, _acnt, _adsp, _acol, _acox, _asel, _abou, true);

// Phase (c)
		selection(0, _acnt, _adsp, _acol, _acox, _asel, _abou, true);

	}
};

template<class T, class S> class parallel_interpolation
{
public:
	void operator()(const vector<T> &_anod, const vector<T> &_acnt, const vector<T> &_acol, const vector<S> &_aele, const vector<S> &_adia, const vector<char> &_abou, const vector<char> &_asel, const vector<char> &_acox, communicator<T, S> &_acom, vector<T> &_mnod, vector<T> &_mcnt, vector<T> &_mcol, vector<S> &_mele)
	{
		int size, rank;
		network::size(size);
		network::rank(rank);

		int anodes = (int)_anod.size();
		vector<T> _adsp(anodes);
		bucket_sort_offset(_acnt, _adsp);

		vector<T> _rcbn, _rloc, _rcnt, _rdsp;
		_acom.boundary_nodes(_rcbn, _rloc, _rcnt, _rdsp);
		vector<T> _rmap(anodes, -1);

		vector<T> _scnt(size, 0);
		vector<T> _sdsp(size);
		for(int i = 0; i < size; i++)
		{
			T rcnt = _rcnt[i], rdsp = _rdsp[i];
			//for(int j = rdsp; j < rdsp + rcnt; j++) _rmap[_rcbn[j]] = _rloc[j];
			for(int j = rdsp; j < rdsp + rcnt; j++)
			{
				T l = _rcbn[j];
				T acnt = _acnt[l], adsp = _adsp[l];
				for(int k = adsp; k < adsp + acnt; k++)
				{
					if((_asel[_acol[k]] == 1) && (_acox[k] != 0)) _scnt[i]++; // strong coarse node
				}
			}
			//for(int j = rdsp; j < rdsp + rcnt; j++) _rmap[_rcbn[j]] = -1;
		}
		bucket_sort_offset(_scnt, _sdsp);
		int ssize = bucket_sort_size(_scnt);
		vector<T> _scol(ssize);
		vector<T> _srow(ssize);
		//vector<S> _sele(ssize);
		for(int m = 0, i = 0; i < size; i++)
		{
			T rcnt = _rcnt[i], rdsp = _rdsp[i];
			for(int j = rdsp; j < rdsp + rcnt; j++) _rmap[_rcbn[j]] = _rloc[j];
			for(int j = rdsp; j < rdsp + rcnt; j++)
			{
			T l = _rcbn[j];
				T acnt = _acnt[l], adsp = _adsp[l];
				for(int k = adsp; k < adsp + acnt; k++)
				{
					if((_asel[_acol[k]] == 1) && (_acox[k] != 0))
					{
						_scol[m] = _anod[_acol[k]];
						_srow[m] = _rmap[l];
						//_sele[m] = _aele[k];
						m++;
					}
				}
			}
			for(int j = rdsp; j < rdsp + rcnt; j++) _rmap[_rcbn[j]] = -1;
		}
		vector<T> _qcnt(size, 0);
		vector<T> _qdsp(size);
		network::alltoall(&_scnt[0], sizeof(T), &_qcnt[0], sizeof(T));
		bucket_sort_offset(_qcnt, _qdsp);
		int qsize = bucket_sort_size(_qcnt);
		vector<T> _ccol(qsize);
		vector<T> _crow(qsize);
		//vector<S> _cele(qsize);
		for(int i = 0; i < size; i++) _scnt[i] *= sizeof(T), _sdsp[i] *= sizeof(T);
		for(int i = 0; i < size; i++) _qcnt[i] *= sizeof(T), _qdsp[i] *= sizeof(T);
		network::alltoallv(_scol.data(), _scnt.data(), _sdsp.data(), _ccol.data(), _qcnt.data(), _qdsp.data());
		network::alltoallv(_srow.data(), _scnt.data(), _sdsp.data(), _crow.data(), _qcnt.data(), _qdsp.data());
		for(int i = 0; i < size; i++) _scnt[i] /= sizeof(T), _sdsp[i] /= sizeof(T);
		for(int i = 0; i < size; i++) _qcnt[i] /= sizeof(T), _qdsp[i] /= sizeof(T);
/*
		for(int i = 0; i < size; i++) _scnt[i] *= sizeof(S), _sdsp[i] *= sizeof(S);
		for(int i = 0; i < size; i++) _qcnt[i] *= sizeof(S), _qdsp[i] *= sizeof(S);
		network::alltoallv(_sele.data(), _scnt.data(), _sdsp.data(), _cele.data(), _qcnt.data(), _qdsp.data());
		for(int i = 0; i < size; i++) _scnt[i] /= sizeof(S), _sdsp[i] /= sizeof(S);
		for(int i = 0; i < size; i++) _qcnt[i] /= sizeof(S), _qdsp[i] /= sizeof(S);
		for(int i = 0; i < size; i++)
		{
			T cnt = _qcnt[i], dsp = _qdsp[i];
			for(int j = dsp; j < dsp + cnt; j++) _ccol[j] = (_ccol[j] << 8) | i;
		}
		binary_sort_sort_copy(_crow, _ccol, _cele);
		for(int i = 0; i < qsize; i++) _ccol[i] = _ccol[i] >> 8;
		unique_accumulate(_crow, _ccol, _cele);
*/
		binary_sort_sort(_crow, _ccol);//
		unique_resize(_crow, _ccol);//

		int csize = (int)_crow.size();
		vector<T> _ccnt(anodes, 0);
		vector<T> _cdsp(anodes);
		for(int i = 0; i < csize; i++) _ccnt[_crow[i]]++;
		bucket_sort_offset(_ccnt, _cdsp);
////
		_mcnt.resize(anodes, 0);
		vector<T> _mdsp(anodes);
		for(int i = 0; i < anodes; i++)
		{
			if(_asel[i] == 1) _mcnt[i] = 1;
			else
			{
				if(_abou[i] == 0)
				{
					T acnt = _acnt[i], adsp = _adsp[i];
					for(int j = adsp; j < adsp + acnt; j++)
					{
						if((_asel[_acol[j]] == 1) && (_acox[j] != 0)) _mcnt[i]++;
					}
				}
				else
				{
					_mcnt[i] = _ccnt[i];
				}
			}
		}
		bucket_sort_offset(_mcnt, _mdsp);
		int msize = bucket_sort_size(_mcnt);
		_mcol.resize(msize);
		_mele.resize(msize);
		for(int i = 0; i < anodes; i++)
		{
			if(_asel[i] == 1)
			{
				T mdsp = _mdsp[i];
				_mcol[mdsp] = _anod[i];
				_mele[mdsp] = S(1.0);
			}
			else
			{
				if(_abou[i] == 0)
				{
					T mcnt = _mcnt[i], mdsp = _mdsp[i];
					T acnt = _acnt[i], adsp = _adsp[i];
					for(int j = adsp; j < adsp + acnt; j++)
					{
						if((_asel[_acol[j]] == 1) && (_acox[j] != 0))
						{
							_mcol[mdsp] = _anod[_acol[j]];
							_mele[mdsp] = S(1.0) / mcnt;
							mdsp++;
						}
					}
				}
				else
				{
					T mcnt = _mcnt[i], mdsp = _mdsp[i];
					T ccnt = _ccnt[i], cdsp = _cdsp[i];
					for(int j = cdsp; j < cdsp + ccnt; j++)
					{
						_mcol[mdsp] = _ccol[j];
						_mele[mdsp] = S(1.0) / mcnt;
						mdsp++;
					}
				}
			}
		}

		for(int i = 0; i < anodes; i++)
		{
			if(_mcnt[i] == 0) cout << "Rank: " << rank << " ! " << i << " " << (int)_abou[i] << " " << (int)_asel[i] << endl;
		}

		_mnod.assign(_mcol.begin(), _mcol.end());
		vector<T> _mloc(msize);
		for(int i = 0; i < msize; i++) _mloc[i] = i;
		binary_sort_copy(_mnod, _mloc);
		for(int k = 0, i = 0; i < msize; k++)
		{
			T nod = _mnod[i];
			while(i < msize && _mnod[i] == nod) _mcol[_mloc[i]] = k, i++;
		}
		unique_resize(_mnod);

	}
};

template<class T, class S> class sparse_multiply_transpose
{
public:
	void operator()(const vector<T> &_acnt, const vector<T> &_acol, const vector<S> &_aele, const vector<T> &_bcnt, const vector<T> &_bcol, const vector<S> &_bele, const vector<T> &_cnod, vector<T> &_ccnt, vector<T> &_ccol, vector<S> &_cele)
	{
#if 0
		int arows = (int)_acnt.size();
		vector<T> _adsp(arows);
		bucket_sort_offset(_acnt, _adsp);
		int brows = (int)_bcnt.size();
		vector<T> _bdsp(brows);
		bucket_sort_offset(_bcnt, _bdsp);
		int crows = (int)_cnod.size();
		_ccnt.assign(crows, 0);

		vector<T> _clst(crows, -1);
		for(int i = 0; i < arows; i++)
		{
			T acnt = _acnt[i], adsp = _adsp[i];
			for(int j = adsp; j < adsp + acnt; j++)
			{
				_mm_prefetch((char*)(_acol.data()+j+64),_MM_HINT_NTA);
				T c = _acol[j];
				T bcnt = _bcnt[c], bdsp = _bdsp[c];
				for(int k = bdsp; k < bdsp + bcnt; k++)
				{
					T d = _bcol[k];
					if(_clst[d] != i)
					{
						_clst[d] = i;
						_ccnt[d]++;
					}
				}
			}
		}
		int csize = bucket_sort_size(_ccnt);
		_ccol.resize(csize);
		_cele.resize(csize);

		_clst.assign(crows, -1);
		vector<T> _cdsp(crows);
		bucket_sort_offset(_ccnt, _cdsp);

		for(int i = 0; i < arows; i++)
		{
			T acnt = _acnt[i], adsp = _adsp[i];
			for(int j = adsp; j < adsp + acnt; j++)
			{
				_mm_prefetch((char*)(_acol.data()+j+64),_MM_HINT_NTA);
				_mm_prefetch((char*)(_aele.data()+j+64),_MM_HINT_NTA);
				T c = _acol[j];
				S s = _aele[j];
				T bcnt = _bcnt[c], bdsp = _bdsp[c];
				for(int k = bdsp; k < bdsp + bcnt; k++)
				{
					T d = _bcol[k];
					T e = _cdsp[d];
					S t = s * _bele[k];
					if(_clst[d] != i)
					{
						_clst[d] = i;
						_ccol[e] = i;
						_cele[e] = t;
						_cdsp[d] = e + 1;
					}
					else
					{
						_cele[e - 1] += t;
					}
				}
			}
		}
#else
		int arows = (int)_acnt.size();
		const T *acnt = _acnt.data();
		const T *acol = _acol.data();
		const S *aele = _aele.data();
		const T *bcnt = _bcnt.data();
		const T *bcol = _bcol.data();
		const S *bele = _bele.data();

		int brows = (int)_bcnt.size();
		vector<T> _bdsp(brows);
		T *bdsp = _bdsp.data();
		for(int t = 0, i = 0; i < brows; i++)
		{
			bdsp[i] = t;
			t += bcnt[i];
		}
		int bcols = (int)_cnod.size();
		//int bcols = *max_element(_bcol.begin(), _bcol.end()) + 1;

		_ccnt.assign(bcols, 0);
		vector<T> _clst(bcols, -1);
		T *ccnt = _ccnt.data();
		T *clst = _clst.data();
		int csize = 0;
		for(int m = 0, i = 0; i < arows; i++)
		{
			for(int j = 0; j < acnt[i]; j++)
			{
				T c = acol[m];
				m++;
				for(int k = bdsp[c]; k < bdsp[c] + bcnt[c]; k++)
				{
					T d = bcol[k];
					if(clst[d] != i)
					{
						clst[d] = i;
						ccnt[d]++;
						csize++;
					}
				}
			}
		}

		_clst.assign(bcols, -1);
		clst = _clst.data();
		vector<T> _cdsp(bcols);
		T *cdsp = _cdsp.data();
		for(int t = 0, i = 0; i < bcols; i++)
		{
			cdsp[i] = t;
			t += ccnt[i];
		}
		_ccol.resize(csize);
		_cele.resize(csize);
		T *ccol = _ccol.data();
		S *cele = _cele.data();
		for(int m = 0, i = 0; i < arows; i++)
		{
			for(int j = 0; j < acnt[i]; j++)
			{
				T c = acol[m];
				S s = aele[m];
				m++;
				for(int k = bdsp[c]; k < bdsp[c] + bcnt[c]; k++)
				{
					T d = bcol[k];
					T e = cdsp[d];
					S t = s * bele[k];
					if(clst[d] != i)
					{
						clst[d] = i;
						ccol[e] = i;
						cele[e] = t;
						cdsp[d] = e + 1;
					}
					else
					{
						cele[e - 1] += t;
					}
				}
			}
		}
#endif
	}
};

template<class T, class S> class new_triple_product
{
public:
	void operator()(const vector<T> &_cnt, const vector<T> &_col, const vector<S> &_ele, const vector<T> &_mcnt, const vector<T> &_mcol, const vector<S> &_mele, const vector<T> &_anod, vector<T> &_acnt, vector<T> &_acol, vector<S> &_aele)
	{
		int rank;
		network::rank(rank);
		vector<T> _tcnt;
		vector<T> _tcol;
		vector<S> _tele;
		sparse_multiply_transpose<T, S> multiply;

		multiply(_cnt, _col, _ele, _mcnt, _mcol, _mele, _anod, _tcnt, _tcol, _tele);
		multiply(_tcnt, _tcol, _tele, _mcnt, _mcol, _mele, _anod, _acnt, _acol, _aele);
	}
};

template<class T, class S> class matrix_reduction
{
public:
	void operator()(vector<T> &_acnt, vector<T> &_acol, vector<S> &_aele, const S _alpha)
	{
		int anumnod = (int)_acnt.size();
		int anumele = (int)_aele.size();
		vector<T> _adsp(_acnt.size());
		bucket_sort_offset(_acnt, _adsp);

		//diagonal
		vector<S> _adia(anumnod, 0.0);
		for(int i = 0; i < anumnod; i++)
		{
			for(int c = _adsp[i]; c < _adsp[i] + _acnt[i]; c++)
			{
				T j = _acol[c];
				S a = _aele[c];
				if(i == j) _adia[i] = a;
			}
		}

		vector<T> _cnt(anumnod, 0);
		vector<T> _col(anumele);
		vector<S> _ele(anumele);

		int k = 0;
		for(int i = 0; i < anumnod; i++)
		{
			for(int c = _adsp[i]; c < _adsp[i] + _acnt[i]; c++)
			{
				T j = _acol[c];
				S a = _aele[c];
				S s = fabs(_adia[i]) * _alpha;
				S t = fabs(_adia[j]) * _alpha;
				if((i == j) || ((fabs(a) >= s) && (fabs(a) >= t))) { _cnt[i] += 1; _col[k] = j; _ele[k] = a; k++; }
			}
		}
		_col.resize(k);
		_ele.resize(k);
		_acnt.assign(_cnt.begin(), _cnt.end());
		_acol.assign(_col.begin(), _col.end());
		_aele.assign(_ele.begin(), _ele.end());
	}
};

template<class T, class S> class amg_restriction
{
public:
	void operator()(const S _epsilon, const vector<T> &_anod, const vector<T> &_acnt, const vector<T> &_acol, const vector<S> &_aele, vector<S> &_adia, vector<char> &_abou, communicator<T, S> &_acom, vector<char> &_mnod, vector<T> &_mcnt, vector<T> &_mcol, vector<S> &_mele, vector<T> &_cnod, vector<T> &_ccnt, vector<T> &_ccol, vector<S> &_cele)
	{
		static int once = 0;
		int oncen = 4;
		vector<char> _asel, _acox;

		parallel_coarsening<T, S> coarsening;
		coarsening(_anod, _acnt, _acol, _aele, _adia, _abou, _asel, _acox, _acom, _epsilon);
/*
		if(once == 0){
		binary_write(_anod.begin(), _anod.end(), "./anod.bin");
		binary_write(_acnt.begin(), _acnt.end(), "./acnt.bin");
		binary_write(_acol.begin(), _acol.end(), "./acol.bin");
		binary_write(_aele.begin(), _aele.end(), "./aele.bin");
		binary_write(_adia.begin(), _adia.end(), "./adia.bin");
		binary_write(_abou.begin(), _abou.end(), "./abou.bin");
		binary_write(_asel.begin(), _asel.end(), "./asel.bin");
		binary_write(_acox.begin(), _acox.end(), "./acox.bin");
		}
*/
		parallel_interpolation<T, S> interpolation;
		interpolation(_anod, _acnt, _acol, _aele, _adia, _abou, _asel, _acox, _acom, _cnod, _mcnt, _mcol, _mele);
/*
		if(once == 4){
        binary_write(_mcnt.begin(), _mcnt.end(), "./mcnt.bin");
        binary_write(_mcol.begin(), _mcol.end(), "./mcol.bin");
        binary_write(_mele.begin(), _mele.end(), "./mele.bin");
		}
*/
		new_triple_product<T, S> triple_product;
		triple_product(_acnt, _acol, _aele, _mcnt, _mcol, _mele, _cnod, _ccnt, _ccol, _cele);
/*
		if(once == 3){
		binary_write(_cnod.begin(), _cnod.end(), "./cnod.bin");
        binary_write(_ccnt.begin(), _ccnt.end(), "./ccnt.bin");
        binary_write(_ccol.begin(), _ccol.end(), "./ccol.bin");
        binary_write(_cele.begin(), _cele.end(), "./cele.bin");
		}
		once++;
*/
	}
};

////

template<class T, class S> class gauss_seidel_f
{
private:
	const vector<T> &_acnt;
	const vector<T> &_acol;
	const vector<S> &_aele;
	const vector<S> &_adia;
	const vector<char> &_abou;
	vector<S> &_v;
	const S _omega;
	communicator<T, S> &_com;
public:
	gauss_seidel_f(const vector<T> &_acnt, const vector<T> &_acol, const vector<S> &_aele, const vector<S> &_adia, const vector<char> &_abou, vector<S> &_v, const S _omega, communicator<T, S> &_com)
		: _acnt(_acnt), _acol(_acol), _aele(_aele), _adia(_adia), _abou(_abou), _v(_v), _omega(_omega), _com(_com)
	{
	}
	void operator()(const vector<S> &_f, vector<S> &_u)
	{
		const T *acnt = &_acnt[0], *acol = &_acol[0];
		const S *aele = &_aele[0], *adia = &_adia[0];
		const S *f = &_f[0];
		const char *abou = &_abou[0];
		S *u = &_u[0];
		S *v = &_v[0];

		S omega = _omega;
		int dsize = (int)_acnt.size();
		for(int i = 0; i < dsize; i++)
		{
			int csize = acnt[i];
			if(abou[i] != 0)
			{
				S s = f[i];
				for(int j = 0; j < csize; j++)
				{
					//_mm_prefetch((char*)&acol[64],_MM_HINT_NTA);
					//_mm_prefetch((char*)&aele[64],_MM_HINT_NTA);
					T q = *acol++;
					S a = *aele++;
					s -= a * u[q];
				}
				v[i] = s * omega;
			}
			else
			{
				acol += csize;
				aele += csize;
			}
		}
		_com.accumulate(_v);
		for(int i = 0; i < dsize; i++)
		{
			if(abou[i] != 0) u[i] += v[i] * adia[i];
		}

		acol = &_acol[0]; aele = &_aele[0];
		for(int i = 0; i < dsize; i++)
		{
			int csize = acnt[i];
			if(abou[i] == 0)
			{
				S s = f[i];
				for(int j = 0; j < csize; j++)
				{
					_mm_prefetch((char*)&acol[64],_MM_HINT_NTA);
					_mm_prefetch((char*)&aele[64],_MM_HINT_NTA);
					T q = *acol++;
					S a = *aele++;
					s -= a * u[q];
				}
				u[i] += s * adia[i];
				//u[i] += s * omega * adia[i];
			}
			else
			{
				acol += csize;
				aele += csize;

			}
		}
	}
	void operator()(const vector<S> &_f, vector<S> &_u, vector<S> &_r)
	{
		const T *acnt = &_acnt[0], *acol = &_acol[0];
		const S *aele = &_aele[0], *adia = &_adia[0];
		const S *f = &_f[0];
		const char *abou = &_abou[0];
		S *u = &_u[0];
		S *r = &_r[0];

		S omega = _omega;
		int dsize = (int)_acnt.size();
		for(int i = 0; i < dsize; i++) r[i] = f[i];
		for(int i = 0; i < dsize; i++)
		{
			int csize = acnt[i];
			S s = adia[i] * r[i];
			u[i] = s;
			for(int j = 0; j < csize; j++)
			{
				_mm_prefetch((char*)&acol[64],_MM_HINT_NTA);
				_mm_prefetch((char*)&aele[64],_MM_HINT_NTA);
				T c = *acol++;
				S t = *aele++;
				r[c] -= t * s;
			}
		}

/*
		for(int i = 0; i < dsize; i++)
		{
			if(abou[i] != 0) u[i] = omega * f[i] * adia[i];
		}
		_com.accumulate(_u);
		for(int i = 0; i < dsize; i++)
		{
			S s = f[i];
			int csize = acnt[i];
			if(abou[i] != 0)
			{
				for(int j = 0; j < csize; j++)
				{
					_mm_prefetch((char*)&acol[64],_MM_HINT_NTA);
					_mm_prefetch((char*)&aele[64],_MM_HINT_NTA);
					T c = *acol++;
					if(abou[c] != 0)
					{
						S t = *aele++;
						s -= t * u[c];
					}
				}
			}
			else
			{
				acol += csize;
				aele += csize;
			}
			r[i] = s;
		}

		//for(int i = 0; i < dsize; i++) r[i] = f[i];
		acol = &_acol[0]; aele = &_aele[0];
		for(int i = 0; i < dsize; i++)
		{
			int csize = acnt[i];
			if(abou[i] == 0)
			{
				S s = adia[i] * r[i];
				u[i] = s;
				for(int j = 0; j < csize; j++)
				{
					_mm_prefetch((char*)&acol[64],_MM_HINT_NTA);
					_mm_prefetch((char*)&aele[64],_MM_HINT_NTA);
					T c = *acol++;
					S t = *aele++;
					r[c] -= t * s;
				}
			}
			else
			{
				acol += csize;
				aele += csize;
			}
		}
*/
	}
};

template<class T, class S> class gauss_seidel_b
{
private:
	const vector<T> &_acnt;
	const vector<T> &_acol;
	const vector<S> &_aele;
	const vector<S> &_adia;
	const vector<char> &_abou;
	vector<S> &_v;
	const S _omega;
	communicator<T, S> &_com;
public:
	gauss_seidel_b(const vector<T> &_acnt, const vector<T> &_acol, const vector<S> &_aele, const vector<S> &_adia, const vector<char> &_abou, vector<S> &_v, const S _omega, communicator<T, S> &_com)
		: _acnt(_acnt), _acol(_acol), _aele(_aele), _adia(_adia), _abou(_abou), _v(_v), _omega(_omega), _com(_com)
	{
	}
	void operator()(const vector<S> &_f, vector<S> &_u)
	{
		const T *acnt = &_acnt[0], *acol = &_acol[0];
		const S *aele = &_aele[0], *adia = &_adia[0];
		const S *f = &_f[0];
		const char *abou = &_abou[0];
		S *u = &_u[0];
		S *v = &_v[0];

		int size = (int)_aele.size();
		acol += size;
		aele += size;

		S omega = _omega;
		int dsize = (int)_acnt.size();
		for(int i = dsize - 1; i >= 0; i--)
		{
			int csize = acnt[i];
			if(abou[i] == 0)
			{
				S s = f[i];
				for(int j = 0; j < csize; j++)
				{
					_mm_prefetch((char*)&acol[-64],_MM_HINT_NTA);
					_mm_prefetch((char*)&aele[-64],_MM_HINT_NTA);
					T q = *--acol;
					S a = *--aele;
					s -= a * u[q];
				}
				u[i] += s * adia[i];
				//u[i] += s * omega * adia[i];
			}
			else
			{
				acol -= csize;
				aele -= csize;
			}
		}

		acol = &_acol[0]; aele = &_aele[0];
		for(int i = 0; i < dsize; i++)
		{
			int csize = acnt[i];
			if(abou[i] != 0)
			{
				S s = f[i];
				for(int j = 0; j < csize; j++)
				{
					//_mm_prefetch((char*)&acol[64],_MM_HINT_NTA);
					//_mm_prefetch((char*)&aele[64],_MM_HINT_NTA);
					T q = *acol++;
					S a = *aele++;
					s -= a * u[q];
				}
				v[i] = s * omega;
			}
			else
			{
				acol += csize;
				aele += csize;
			}
		}
		_com.accumulate(_v);
		for(int i = 0; i < dsize; i++)
		{
			if(abou[i] != 0) u[i] += v[i] * adia[i];
		}
	}
	void operator()(const vector<S> &_f, vector<S> &_u, vector<S> &_s)
	{
		const T *acnt = &_acnt[0], *acol = &_acol[0];
		const S *aele = &_aele[0], *adia = &_adia[0];
		const S *f = &_f[0];
		const char *abou = &_abou[0];
		S *u = &_u[0];
		S *s = &_s[0];

		int size = (int)_aele.size();
		acol += size;
		aele += size;

		S omega = _omega;
		int dsize = (int)_acnt.size();
		for(int i = 0; i < dsize; i++) u[i] += s[i];
		for(int i = dsize - 1; i >= 0; i--)
		{
			int csize = acnt[i];
			S s = f[i];
			for(int j = 0; j < csize; j++)
			{
				_mm_prefetch((char*)&acol[-64],_MM_HINT_NTA);
				_mm_prefetch((char*)&aele[-64],_MM_HINT_NTA);
				T q = *--acol;
				S a = *--aele;
				s -= a * u[q];
			}
			u[i] += s * adia[i];
		}

	}
};

template<class T, class S> class jacobi
{
private:
	const vector<T> &_acnt;
	const vector<T> &_acol;
	const vector<S> &_aele;
	const vector<S> &_adia;
	vector<S> &_v;
	const S _omega;
	communicator<T, S> &_com;
public:
	jacobi(const vector<T> &_acnt, const vector<T> &_acol, const vector<S> &_aele, const vector<S> &_adia, vector<S> &_v, const S _omega, communicator<T, S> &_com)
		: _acnt(_acnt), _acol(_acol), _aele(_aele), _adia(_adia), _v(_v), _omega(_omega), _com(_com)
	{
	}
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
				_mm_prefetch((char*)&acol[64],_MM_HINT_NTA);
				_mm_prefetch((char*)&aele[64],_MM_HINT_NTA);
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
	void operator()(const vector<S> &_f, vector<S> &_u, vector<S> &_r)
	{
		const T *acnt = &_acnt[0], *acol = &_acol[0];
		const S *aele = &_aele[0], *adia = &_adia[0];
		const S *f = &_f[0];
		S *u = &_u[0];
		S *r = &_r[0];

		S omega = _omega;
		int dsize = (int)_acnt.size();
		for(int i = 0; i < dsize; i++)
		{
			u[i] = omega * f[i] * adia[i];
		}
		_com.accumulate(_u);

		for(int i = 0; i < dsize; i++)
		{
			S s = f[i];
			int csize = acnt[i];
			for(int j = 0; j < csize; j++)
			{
				_mm_prefetch((char*)&acol[64],_MM_HINT_NTA);
				_mm_prefetch((char*)&aele[64],_MM_HINT_NTA);
				T c = *acol++;
				S t = *aele++;
				s -= t * u[c];
			}
			r[i] = s;
		}
	}
};

template<class T, class S> class pseudo_jacobi
{
private:
	const vector<S> &_adia;
	const S _omega;
	communicator<T, S> &_com;
public:
	pseudo_jacobi(const vector<S> &_adia, const S _omega, communicator<T, S> &_com)
		: _adia(_adia), _omega(_omega), _com(_com)
	{
	}
	void operator()(const vector<S> &_f, vector<S> &_u)
	{
		const S *adia = &_adia[0];
		const S *f = &_f[0];
		S *u = &_u[0];

		_com.distribute(_u);
		S omega = _omega;
		int dsize = (int)_adia.size();
		for(int i = 0; i < dsize; i++)
		{
			*u++ += omega * *f++ * *adia++;
		}
		_com.accumulate(_u);
	}
};



////

template<class T, class S> class diagonal_solver : public inverse_operator<T, S>
{
private:
	const vector<S> &_adia;
	communicator<T, S> &_com;
public:
	diagonal_solver(const vector<S> &_adia, communicator<T, S> &_com)
		: _adia(_adia), _com(_com)
	{
	}
	void operator()(const vector<S> &_u, vector<S> &_v)
	{
		const S *adia = _adia.data();
		const S *u = _u.data();
		S *v = _v.data();

		int asize = (int)_adia.size();
		for(int i = 0; i < asize; i++)
		{
			v[i] = adia[i] * u[i];
		}
		_com.accumulate(_v);
	}
};

template<class T, class S> class amg_solver : public inverse_operator<T, S>
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
public:
	amg_solver(const vector<T> &_anod0, const vector<T> &_acnt0, const vector<T> &_acol0, const vector<S> &_aele0, 
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

#if 1
		_anod[0].assign(_anod0.begin(), _anod0.end());
		_acnt[0].assign(_acnt0.begin(), _acnt0.end());
		_acol[0].assign(_acol0.begin(), _acol0.end());
		_aele[0].assign(_aele0.begin(), _aele0.end());
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
		int l = 0;
		int gsize;
		network::size(gsize);
		vector<int> _gcnt(gsize, 0);
		network::allgather(&size, sizeof(T), &_gcnt[0], sizeof(T));
		int sum = 0;
		for(int i = 0; i < (int)_gcnt.size(); i++) sum += _gcnt[i];
		int gcnt = sum / gsize;
		double rel = 0.0;
		while(l < max_level - 1 && gcnt > max_size && rel < 1.0)
		//while(l < max_level - 1 && gcnt > max_size)
		{
			_com[l].assign(_anod[l]);
			amg_restriction(_epsilon, _anod[l], _acnt[l], _acol[l], _aele[l], _adia[l], _abou[l], _com[l], _rnod[l], _rcnt[l], _rcol[l], _rele[l], _anod[l+1], _acnt[l+1], _acol[l+1], _aele[l+1]);
			//rel = (float)_acnt[l+1].size()/(float)_acnt[l].size();

			for(int i = 0; i < (int)_adia[l].size(); i++) _adia[l][i] = 1.0 / _adia[l][i];
			//cout << "Rank: " << rank << " Size: (" << _acol[l].size() << " | " << _acnt[l].size() << ")" << endl;
			l++;

			size = (int)_acnt[l].size();
			_u[l].resize(size, 0.0);
			_f[l].resize(size, 0.0);
			_s[l].resize(size, 0.0);
			_r[l].resize(size, 0.0);
			_v[l].resize(size, 0.0);

			network::allgather(&size, sizeof(T), &_gcnt[0], sizeof(T));
			int sum = 0;
			for(int i = 0; i < (int)_gcnt.size(); i++) sum += _gcnt[i];
			int gold = gcnt;
			gcnt = sum / gsize;
			rel = (double)gcnt / (double)gold;
		    //Paola Commented on 20-03-2019	
            //if(rank == 0) cout << "Rank: " << rank << " Relative: " << rel << " Nodes: " << sum << " Level: " << l << endl;

		}
		_com[l].assign(_anod[l]);
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
	void sub(const vector<S> &_a, vector<S> &_b)
	{
		for(int i = 0; i < _a.size(); i++)
		{
			_b[i] -= _a[i];
		}
	}
	void rsb(const vector<S> &_a, vector<S> &_b)
	{
		for(int i = 0; i < (int)_a.size(); i++)
		{
			_b[i] = _a[i] - _b[i];
		}
	}
	void add(const vector<S> &_a, vector<S> &_b)
	{
		for(int i = 0; i < (int)_a.size(); i++)
		{
			_b[i] += _a[i];
		}
	}
	void operator()(const vector<S> &_f0, vector<S> &_u0)
	{
		int size = (int)_f0.size();
		//for(int i = 0; i < size; i++) _f[0][i] = _f0[i];
		memcpy(_f[0].data(), _f0.data(), _f0.size() * sizeof(S));
		multigrid(0);
		memcpy(_u0.data(), _u[0].data(), _u[0].size() * sizeof(S));
		//for(int i = 0; i < size; i++) _u0[i] = _u[0][i];
	}
	void multigrid(int _k)
	{
		int l = _k++;
		S delta = _omega;

		if(l < _level)
		{
			linear_operator<T, S> _K(_acnt[l], _acol[l], _aele[l]);
			gauss_seidel_f<T, S> _S(_acnt[l], _acol[l], _aele[l], _adia[l], _abou[l], _v[l], delta, _com[l]);
			gauss_seidel_b<T, S> _T(_acnt[l], _acol[l], _aele[l], _adia[l], _abou[l], _v[l], delta, _com[l]);
			jacobi<T, S> _J(_acnt[l], _acol[l], _aele[l], _adia[l], _v[l], delta, _com[l]);
			pseudo_jacobi<T, S> _I(_adia[l], delta, _com[l]);
			simple_prolongation<T, S> _P(_rcnt[l], _rcol[l]);
			simple_restriction<T, S> _R(_rcnt[l], _rcol[l]);
			//prolongation<T, S> _P(_rcnt[l], _rcol[l], _rele[l]);
			//restriction<T, S> _R(_rcnt[l], _rcol[l], _rele[l]);
			if(_cascadic)
			{
				_R(_f[l], _f[l + 1]);
				multigrid(_k);
				_P(_u[l + 1], _u[l]);

				if(l < 0)
				{
					_I(_f[l], _u[l]);
				}
				else
				{
					//for(int i = 0; i < l * l + 1; i++)
					{
						if(_gauss_seidel)
						{
							_S(_f[l], _u[l]);
							_T(_f[l], _u[l]);
						}
						else
						{
							_J(_f[l], _u[l]);
							_J(_f[l], _u[l]);
						}
					}
				}
			}
			else
			{
				if(_gauss_seidel)
				{
#ifdef NOMPI
					_S(_f[l], _u[l], _r[l]);
					_R(_r[l], _f[l + 1]);
					multigrid(_k);
					_P(_u[l + 1], _s[l]);
					_T(_f[l], _u[l], _s[l]);
#else
					memset(&_u[l][0], 0, _u[l].size() * sizeof(S));
					_S(_f[l], _u[l]);
					_K(_u[l], _r[l]);
					rsb(_f[l], _r[l]);
					_R(_r[l], _f[l + 1]);
					multigrid(_k);
					_P(_u[l + 1], _s[l]);
					add(_s[l], _u[l]);
					_T(_f[l], _u[l]);
#endif
				}
				else
				{
#if 1
					_J(_f[l], _u[l], _r[l]);
#else
					memset(&_u[l][0], 0, _u[l].size() * sizeof(S));
					_J(_f[l], _u[l]);
					_K(_u[l], _r[l]);
					rsb(_f[l], _r[l]);
#endif
					_R(_r[l], _f[l + 1]);
					multigrid(_k);
					_P(_u[l + 1], _s[l]);
					add(_s[l], _u[l]);
					_J(_f[l], _u[l]);
				}
			}
		}
		else
		{
#ifdef SUPERLU
			(*_superlu)(_f[l], _u[l]);
#else
#ifdef MUMPS
			(*_mumps)(_f[l], _u[l]);
#else
			diagonal_solver<T, S> _D(_adia[l], _com[l]);
			_D(_f[l], _u[l]);
#endif
#endif
		}

	}
};

template<class T, class S> class fast_amg_solver : public inverse_operator<T, S>
{
private:
	amg_solver<T, float> &_amg;
	vector<float> _f;
	vector<float> _u;
public:
	fast_amg_solver(amg_solver<T, float> &_amg)
		: _amg(_amg)
	{
	}
	void operator()(const vector<S> &_f0, vector<S> &_u0)
	{
		int size = (int)_f0.size();
		_f.resize(size);
		_u.resize(size);
		for(int i = 0; i < size; i++) _f[i] = (float)_f0[i];
		_amg(_f, _u);
		for(int i = 0; i < size; i++) _u0[i] = (S)_u[i];
	}
};
