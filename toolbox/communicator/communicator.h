/**
 * The communicator class provides all communication functions for the parallel toolbox.
 * \author Manfred Liebmann
 */
template <class T, class S> class communicator
{
private:
	vector<T> _rpos;	/**< Vector of ordered local node numbers for the distribute procedure */
	vector<S> _rsca;	/**< Vector of the reciprocal multiplicity for the distribute procedure */
	vector<T> _rcom;	/**< Vector of all shared local nodes ordered by processor for the accumulate procedure */
	vector<T> _rpro;	/**< Vector of processor numbers */
	vector<T> _rloc;	/**< Vector of all shared foreign nodes ordered by processor */
	vector<int> _rcnt;	/**< Vector of shared node counts for the accumulate procedure */
	vector<int> _rdsp;	/**< Vector of shared node displacements for the accumulate procedure. */
	vector<S> _rvec;	/**< Receive buffer for the accumulate procedure */
	vector<S> _svec;	/**< Send buffer for the accumulate procedure */
	vector<S> _rcol;	/**< Receive buffer for the collect procedure */
	//vector<T> _rglo;
	vector<int> _pcnt;
	vector<int> _pdsp;
	vector<T> _ploc;
	vector<T> _ppro;
	vector<T> _ppos;
public:
	/**
	* The collect procedure calculates the global sum of all input values.
	* \param _t Input: Numerical value, Output: Sum of all input values.
	*/
	void collect(S &_t)
	{
		network::allgather(&_t, sizeof(S), _rcol.data(), sizeof(S));
		S s = 0.0;
		for(int i = 0; i < (int)_rcol.size(); i++) s += _rcol[i];
		_t = s;
	}
	/**
	* The collect procedure calculates the global sum of all input vectors.
	* \param _t Input: Numerical vector, Output: Sum of all input vectors.
	*/
	void collect(vector<S> &_t)
	{
		int m = (int)_t.size();
		int n = (int)_rcnt.size();
		_rcol.resize(n * m);
		network::allgather(_t.data(), m * sizeof(S), _rcol.data(), m * sizeof(S));

		S *r = &_rcol[0];
		for(int i = 0; i < m; i++)
		{
			S s = 0.0;
			for(int j = 0; j < n; j++) s += r[m * j + i];
			_t[i] = s;
		}

	}
	/**
	* The distribute procedure converts an accumulated vector to a distributed vector.
	* \param _x Input: Accumulated vector, Output: Distributed vector.
	*/
	void distribute(vector<S> &_x)
	{
		if(_rpos.size() == 0) return;
		S *x = _x.data();
		for(int i = 0; i < (int)_rpos.size(); i++)
		{
			x[_rpos[i]] *= _rsca[i];
		}
	}
	/**
	* The accumulate procedure converts a distributed vector to an accumulated vector.
	* \param _x Input: Distributed vector, Output: Accumulated vector.
	*/
	void accumulate(vector<S> &_x)
	{
		if(_rcom.size() == 0) return;
		for(int i = 0; i < (int)_rcom.size(); i++) _svec[i] = _x[_rcom[i]];
		for(int i = 0; i < (int)_rcnt.size(); i++) _rcnt[i] *= sizeof(S), _rdsp[i] *= sizeof(S);
		network::alltoallv(_svec.data(), _rcnt.data(), _rdsp.data(), _rvec.data(), _rcnt.data(), _rdsp.data());
		for(int i = 0; i < (int)_rcnt.size(); i++) _rcnt[i] /= sizeof(S), _rdsp[i] /= sizeof(S);
		for(int i = 0; i < (int)_rcom.size(); i++) _x[_rcom[i]] += _rvec[i];
	}
	/**
	* The assign procedure calculates the communication data structures from a vector of global node numbers.
	* \param _gnod Input: Vector of global node numbers.
	*/
	void assign(const vector<T> &_gnod)
	{
		//_rglo.assign(_gnod.begin(), _gnod.end());
		int size, rank;
		network::size(size);
		network::rank(rank);

		int snodes = (int)_gnod.size();
		vector<T> _snod(snodes);
		vector<T> _spro(snodes);
		vector<T> _sloc(snodes);
		vector<int> _scnt(size, 0);
		vector<int> _sdsp(size);
		for(int i = 0; i < snodes; i++) 
		{
			//T pro = _gnod[i] & (size - 1);
			T pro = _gnod[i] % size;
			_scnt[pro]++;
			_spro[i] = pro;
		}
		bucket_sort_offset(_scnt, _sdsp);
		for(int i = 0; i < snodes; i++)
		{
			T pro = _spro[i];
			T nod = _gnod[i];
			int dsp = _sdsp[pro];
			_snod[dsp] = nod;
			_sloc[dsp] = i;
			_sdsp[pro]++;
		}
		bucket_sort_offset(_scnt, _sdsp);

		vector<int> _qcnt(size, 0);
		vector<int> _qdsp(size);
		network::alltoall(_scnt.data(), sizeof(int), _qcnt.data(), sizeof(int));
		bucket_sort_offset(_qcnt, _qdsp);
		int qnodes = bucket_sort_size(_qcnt);
		vector<T> _qnod(qnodes);
		vector<T> _qloc(qnodes);

		for(int i = 0; i < size; i++) _scnt[i] *= sizeof(T), _sdsp[i] *= sizeof(T);
		for(int i = 0; i < size; i++) _qcnt[i] *= sizeof(T), _qdsp[i] *= sizeof(T);
		network::alltoallv(_sloc.data(), _scnt.data(), _sdsp.data(), _qloc.data(), _qcnt.data(), _qdsp.data());
		network::alltoallv(_snod.data(), _scnt.data(), _sdsp.data(), _qnod.data(), _qcnt.data(), _qdsp.data());
		for(int i = 0; i < size; i++) _scnt[i] /= sizeof(T), _sdsp[i] /= sizeof(T);
		for(int i = 0; i < size; i++) _qcnt[i] /= sizeof(T), _qdsp[i] /= sizeof(T);

		vector<T> _qpro(qnodes);
		for(int i = 0; i < size; i++)
		{
			int dsp = _qdsp[i];
			int cnt = _qcnt[i];
			for(int j = dsp; j < dsp + cnt; j++) _qpro[j] = i;
		}
		binary_sort_sort_copy(_qnod, _qpro, _qloc);
		vector<int> _qacc(qnodes, 1);
		unique_accumulate(_qnod, _qacc);

		int qnodes2 = (int)_qacc.size();
		vector<int> _tcnt(size, 0);
		vector<int> _tdsp(size);
		for(int k = 0, i = 0; i < qnodes2; i++)
		{
			int acc = _qacc[i];
			for(int j = 0; j < acc; j++)
			{
				_tcnt[_qpro[k + j]] += acc - 1;
			}
			k += acc;
		}
		bucket_sort_offset(_tcnt, _tdsp);

		int tnodes = bucket_sort_size(_tcnt);
		vector<T> _tloc(tnodes);
		vector<T> _tpro(tnodes);
		for(int k = 0, i = 0; i < qnodes2; i++)
		{
			int acc = _qacc[i];
			for(int j = 0; j < acc; j++)
			{
				for(int l = 0; l < acc; l++)
				{
					if(l != j)
					{
						T pro = _qpro[k + j];
						int dsp = _tdsp[pro];
						_tloc[dsp] = _qloc[k + j];
						_tpro[dsp] = _qpro[k + l];
						_tdsp[pro]++;
					}
				}
			}
			k += acc;
		}
		bucket_sort_offset(_tcnt, _tdsp);

		_pcnt.resize(size, 0);
		_pdsp.resize(size);
		network::alltoall(_tcnt.data(), sizeof(int), _pcnt.data(), sizeof(int));
		bucket_sort_offset(_pcnt, _pdsp);
		int pnodes = bucket_sort_size(_pcnt);
		_ploc.resize(pnodes);
		_ppro.resize(pnodes);

		for(int i = 0; i < size; i++) _tcnt[i] *= sizeof(T), _tdsp[i] *= sizeof(T);
		for(int i = 0; i < size; i++) _pcnt[i] *= sizeof(T), _pdsp[i] *= sizeof(T);
		network::alltoallv(_tloc.data(), _tcnt.data(), _tdsp.data(), _ploc.data(), _pcnt.data(), _pdsp.data());
		network::alltoallv(_tpro.data(), _tcnt.data(), _tdsp.data(), _ppro.data(), _pcnt.data(), _pdsp.data());
		for(int i = 0; i < size; i++) _tcnt[i] /= sizeof(T), _tdsp[i] /= sizeof(T);
		for(int i = 0; i < size; i++) _pcnt[i] /= sizeof(T), _pdsp[i] /= sizeof(T);
		binary_sort_copy(_ploc, _ppro);

		_rcom.assign(_ploc.begin(), _ploc.end());
		_rpro.assign(_ppro.begin(), _ppro.end());
		binary_sort_sort(_rpro, _rcom);

		_rcnt.resize(size, 0);
		_rdsp.resize(size);
		for(int i = 0; i < pnodes; i++) _rcnt[_rpro[i]]++;
		bucket_sort_offset(_rcnt, _rdsp);

		_ppos.resize(pnodes);
		for(int i = 0; i < pnodes; i++)
		{
			T pro = _ppro[i];
			int dsp = _rdsp[pro];
			_rdsp[pro]++;
			_ppos[i] = dsp;
		}
		bucket_sort_offset(_rcnt, _rdsp);
////
		int rsize = (int)_rcom.size();
		_rloc.resize(rsize);

		for(int i = 0; i < size; i++) _rcnt[i] *= sizeof(T), _rdsp[i] *= sizeof(T);
		network::alltoallv(_rcom.data(), _rcnt.data(), _rdsp.data(), _rloc.data(), _rcnt.data(), _rdsp.data());
		for(int i = 0; i < size; i++) _rcnt[i] /= sizeof(T), _rdsp[i] /= sizeof(T);
////
		//for(int i = 0; i < size; i++) _pcnt[i] = _rcnt[i] * sizeof(S);
		//bucket_sort_offset(_pcnt, _pdsp);

		_rvec.resize(pnodes);
		_svec.resize(pnodes);
		_rcol.resize(size);

		_rpos.assign(_ploc.begin(), _ploc.end());
		_rsca.resize(pnodes, 1.0);
		unique_accumulate(_rpos, _rsca);
		int pnodes2 = (int)_rsca.size();
		for(int i = 0; i < pnodes2; i++)
		{
			_rsca[i] = S(1.0)/(S(1.0) + _rsca[i]);
		}

	}
	/**
	* The constructor calculates the communication data structures from a vector of global node numbers.
	* \param _gnod Input: Vector of global node numbers.
	*/
	communicator(const vector<T> &_gnod)
	{
		this->assign(_gnod);
	}
	/**
	* The constructor calculates the communication data structures from a vector of global node numbers.
	* \param _gnod Input: Vector of global node numbers.
	*/
	communicator()
	{
	}

	void boundary(vector<char> &_abou)
	{
		//int nodes = (int)_rcom.size();
		//for(int i = 0; i < nodes; i++) _abou[_rcom[i]] = 1;
		int nodes = (int)_rpos.size();
		for(int i = 0; i < nodes; i++) _abou[_rpos[i]] = 1;
	}
	void boundary_nodes(vector<T> &_acbn, vector<T> &_aloc, vector<T> &_acnt, vector<T> &_adsp)
	{
		int size, rank;
		network::size(size);
		network::rank(rank);

		int rnodes = _rpos.size();
		_acnt.assign(_rcnt.begin(), _rcnt.end());
		_acnt[rank] = rnodes;
		_adsp.resize(size);
		bucket_sort_offset(_acnt, _adsp);
		int asize = bucket_sort_size(_acnt);
		_acbn.resize(asize);
		_aloc.resize(asize);
		for(int k = 0, i = 0; i < size; i++)
		{
			if(i == rank)
			{
				for(int j = 0; j < rnodes; j++, k++) _acbn[k] = _aloc[k] = _rpos[j];
			}
			else
			{
				T rcnt = _rcnt[i], rdsp = _rdsp[i];
				for(int j = rdsp; j < rdsp + rcnt; j++, k++)
				{
					_acbn[k] = _rcom[j];
					_aloc[k] = _rloc[j];
				}
			}
		}
/*
		_acbn.assign(_rcom.begin(), _rcom.end());
		_aloc.assign(_rloc.begin(), _rloc.end());
		_acnt.assign(_rcnt.begin(), _rcnt.end());
		_adsp.assign(_rdsp.begin(), _rdsp.end());
*/
	}
};
