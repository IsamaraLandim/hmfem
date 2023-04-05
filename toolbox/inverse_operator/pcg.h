/**
* The conjugate_gradient function implements the preconditioned conjugate gradient algorithm.
* \param _K Input: System matrix.
* \param _C Input: Preconditioner matrix.
* \param _S Input: Scalar product.
* \param _u Input: Solution vector in packed vector format.
* \param _f Input: Right hand side in packed vector format.
* \param _eps Input: Desired accuracy for residual norm.
* \param _max Input: Maximum number of iterations.
* \return Output: Number of iterations.
*/
#if 1
template<class T, class S> class conjugate_gradient : public inverse_operator<T, S>
{
private:
	const generic_operator<T, S> &_K;
	inverse_operator<T, S> &_C;
	binary_operator<T, S> &_S;
	const S _eps;
	const int _max;
	const bool _show;
	int _loop;
public:
	conjugate_gradient(const generic_operator<T, S> &_K, inverse_operator<T, S> &_C, binary_operator<T, S> &_S, const S _eps, const int _max, const bool _show = true)
		: _K(_K), _C(_C), _S(_S), _eps(_eps), _max(_max), _show(_show)
	{
		_loop = 0;
	}
	int iterations()
	{
		return _loop;
	}
	void operator()(const vector<S> &_f, vector<S> &_u)
	{
		int rank;
		network::rank(rank);
		int n = 1;
		vector<S> alpha(n, 1.0), beta(n, 0.0), sigma(n, 0.0), sigma_first(n, 0.0), sigma_last(n, 0.0), tau(n, 0.0);
		vector<S> _r(_f);
		vector<S> _s(_u);
		vector<S> _v(_r.size());

		_K(_s, _v);
		sub_scale(_r, _v, alpha);
		_C(_r, _v);
		_S(_v, _r, sigma);
		if(sigma[0] < 0.0) { sigma[0] = fabs(sigma[0]); if(_show) if(rank == 0) cout << "[!!!]"; }
		scale_add(_s, _v, beta);
		memcpy(&sigma_last[0], &sigma[0], sizeof(S) * sigma.size());
		memcpy(&sigma_first[0], &sigma[0], sizeof(S) * sigma.size());
		for(int i = 0; i < (int)alpha.size(); i++)
		{
			if(sigma[i] == 0.0) sigma_first[i] = 1.0;
			tau[i] = sigma[i] / sigma_first[i];
		}
		_loop = 0;
		while((++_loop < _max) && (*max_element(tau.begin(), tau.end()) > _eps * _eps))
		{
			if(_show) if(rank == 0) cout << "**Rank: " << rank << " Loop: " << _loop << " Pre Norm2: " << sqrt(sigma[0]) << " " << sqrt(sigma[0] / sigma_first[0]) << endl;
			_K(_s, _v);
			_S(_s, _v, tau);
			if(tau[0] < 0.0) { tau[0] = fabs(tau[0]); if(_show) if(rank == 0) cout << "[!]"; }
			for(int i = 0; i < (int)alpha.size(); i++) alpha[i] = sigma[i] / tau[i];
			add_scale(_u, _s, alpha);
			sub_scale(_r, _v, alpha);
			_C(_r, _v);
			_S(_v, _r, sigma);
			if(sigma[0] < 0.0) { sigma[0] = fabs(sigma[0]); if(_show) if(rank == 0) cout << "[!!]"; }
			for(int i = 0; i < (int)beta.size(); i++) beta[i] = sigma[i] / sigma_last[i];
			scale_add(_s, _v, beta);
			memcpy(&sigma_last[0], &sigma[0], sizeof(S) * sigma.size());
			for(int i = 0; i < (int)alpha.size(); i++) tau[i] = sigma[i] / sigma_first[i];
		}
		if(_show) if(rank == 0) cout << "**Rank: " << rank << " Loop: " << _loop << " Pre Norm2: " << sqrt(sigma[0]) << " " << sqrt(sigma[0] / sigma_first[0]) << endl;
	}
};
#else
template<class T, class S> class conjugate_gradient : public inverse_operator<T, S>
{
private:
	const generic_operator<T, S> &_K;
	inverse_operator<T, S> &_C;
	binary_operator<T, S> &_S;
	const S _eps;
	const int _max;
	const bool _show;
	int _loop;
public:
	conjugate_gradient(const generic_operator<T, S> &_K, inverse_operator<T, S> &_C, binary_operator<T, S> &_S, const S _eps, const int _max, const bool _show = true)
		: _K(_K), _C(_C), _S(_S), _eps(_eps), _max(_max), _show(_show)
	{
		_loop = 0;
	}
	int iterations()
	{
		return _loop;
	}
	void operator()(const vector<S> &_f, vector<S> &_u)
	{
		int rank;
		network::rank(rank);
		int n = 1;
		vector<S> alpha(n, 1.0), beta(n, 0.0), sigma(n, 0.0), sigma_first(n, 0.0), sigma_last(n, 0.0), tau(n, 0.0);
		vector<S> _r(_f);
		vector<S> _s(_u);
		vector<S> _v(_r.size());
		vector<S> gamma(n, 0.0);

		_K(_s, _v);
		sub_scale(_r, _v, alpha);
		_S(_r, _r, gamma);//BUG SEQ
		_C(_r, _v);
		_S(_v, _r, sigma);
		if(sigma[0] < 0.0) { sigma[0] = fabs(sigma[0]); if(_show) if(rank == 0) cout << "[!!!]"; }
		scale_add(_s, _v, beta);
		memcpy(&sigma_last[0], &sigma[0], sizeof(S) * sigma.size());
		memcpy(&sigma_first[0], &sigma[0], sizeof(S) * sigma.size());
		for(int i = 0; i < (int)alpha.size(); i++) tau[i] = sigma[i] / sigma_first[i];

		_loop = 0;
		//while((++_loop < _max) && (*max_element(tau.begin(), tau.end()) > _eps * _eps))
		while((++_loop < _max) && (*max_element(gamma.begin(), gamma.end()) > _eps * _eps))
		{
			if(_show) if(rank == 0) cout << "**Rank: " << rank << " Loop: " << _loop << " Pre Norm2: " << sqrt(gamma[0]) << endl;
			_K(_s, _v);
			_S(_s, _v, tau);
			if(tau[0] < 0.0) { tau[0] = fabs(tau[0]); if(_show) if(rank == 0) cout << "[!]"; }
			for(int i = 0; i < (int)alpha.size(); i++) alpha[i] = sigma[i] / tau[i];
			add_scale(_u, _s, alpha);
			sub_scale(_r, _v, alpha);
			_S(_r, _r, gamma);//BUG SEQ
			_C(_r, _v);
			_S(_v, _r, sigma);
			if(sigma[0] < 0.0) { sigma[0] = fabs(sigma[0]); if(_show) if(rank == 0) cout << "[!!]"; }
			for(int i = 0; i < (int)beta.size(); i++) beta[i] = sigma[i] / sigma_last[i];
			scale_add(_s, _v, beta);
			memcpy(&sigma_last[0], &sigma[0], sizeof(S) * sigma.size());
			for(int i = 0; i < (int)alpha.size(); i++) tau[i] = sigma[i] / sigma_first[i];
		}
		if(_show) if(rank == 0) cout << "**Rank: " << rank << " Loop: " << _loop << " Pre Norm2: " << sqrt(gamma[0]) << endl;
	}
};
#endif
