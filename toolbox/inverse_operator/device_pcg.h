extern "C" {
//void _device_add_scale(float *x, float *y, float s, unsigned int m);
//void _device_sub_scale(float *x, float *y, float s, unsigned int m);
//void _device_scale_add(float *x, float *y, float s, unsigned int m);
void _device_add_scale(double *x, double *y, double s, unsigned int m);
void _device_sub_scale(double *x, double *y, double s, unsigned int m);
void _device_scale_add(double *x, double *y, double s, unsigned int m);
}

/**
* The add_scale procedure calculates the sum of two packed vectors, scaling the second packed vector. x += y * s
* \param _x Input: Numerical vector, first argument. Output: Scaled vector sum. x += y * s
* \param _y Input: Numerical vector, second argument.
* \param _s Input: Vector of scaling factors.
*/
template<class S>
void device_add_scale(device_vector<S> &_x, const device_vector<S> &_y, const vector<S> &_s)
{
	//_device_add_scale(&_x[0], &_y[0], _s[0], _x.size());
	_device_add_scale(_x.data(), _y.data(), _s[0], _x.size());
}
/**
* The sub_scale procedure calculates the difference of two packed vectors, scaling the second packed vector. x -= y * s
* \param _x Input: Numerical vector, first argument. Output: Scaled vector sum. x -= y * s
* \param _y Input: Numerical vector, second argument.
* \param _s Input: Vector of scaling factors.
*/
template<class S>
void device_sub_scale(device_vector<S> &_x, const device_vector<S> &_y, const vector<S> &_s)
{
	//_device_sub_scale(&_x[0], &_y[0], _s[0], _x.size());
	_device_sub_scale(_x.data(), _y.data(), _s[0], _x.size());
}

/**
* The scale_add procedure calculates the sum of two vectors, scaling the first vector. x = x * s + y
* \param _x Input: Numerical vector, first argument. Output: Scaled vector sum. x = x * s + y
* \param _y Input: Numerical vector, second argument.
* \param _s Input: Scaling factor.
*/
template<class S>
void device_scale_add(device_vector<S> &_x, const device_vector<S> &_y, const vector<S> &_s)
{
	//_device_scale_add(&_x[0], &_y[0], _s[0], _x.size());
	_device_scale_add(_x.data(), _y.data(), _s[0], _x.size());
}
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
template<class T, class S> class device_conjugate_gradient : public device_inverse_operator<T, S>
{
private:
	const device_generic_operator<T, S> &_K;
	device_inverse_operator<T, S> &_C;
	device_binary_operator<T, S> &_S;
	const S _eps;
	const int _max;
	const bool _show;
	int _loop;
public:
	device_conjugate_gradient(const device_generic_operator<T, S> &_K, device_inverse_operator<T, S> &_C, device_binary_operator<T, S> &_S, const S _eps, const int _max, const bool _show = true)
		: _K(_K), _C(_C), _S(_S), _eps(_eps), _max(_max), _show(_show)
	{
		_loop = 0;
	}
	int iterations()
	{
		return _loop;
	}
	void operator()(const device_vector<S> &_f, device_vector<S> &_u)
	{
		int rank;
		network::rank(rank);
		int n = 1;
		vector<S> alpha(n, 1.0), beta(n, 0.0), sigma(n, 0.0), sigma_first(n, 0.0), sigma_last(n, 0.0), tau(n, 0.0);
		device_vector<S> _r(_f);
		device_vector<S> _s(_u);
		device_vector<S> _v(_r.size());

		_K(_s, _v);
		device_sub_scale(_r, _v, alpha);
		_C(_r, _v);
		_S(_v, _r, sigma);
		if(sigma[0] < 0.0) { sigma[0] = fabs(sigma[0]); if(_show) if(rank == 0) cout << "[!!!]"; }
		device_scale_add(_s, _v, beta);
		memcpy(&sigma_last[0], &sigma[0], sizeof(S) * sigma.size());
		memcpy(&sigma_first[0], &sigma[0], sizeof(S) * sigma.size());
		for(int i = 0; i < (int)alpha.size(); i++) tau[i] = sigma[i] / sigma_first[i];

		_loop = 0;
		while((++_loop < _max) && (*max_element(tau.begin(), tau.end()) > _eps * _eps))
		{
			if(_show) if(rank == 0) cout << "**Rank: " << rank << " Loop: " << _loop << " Pre Norm2: " << sqrt(sigma[0]) << " " << sqrt(sigma[0] / sigma_first[0]) << endl;
			_K(_s, _v);
			_S(_s, _v, tau);
			if(tau[0] < 0.0) { tau[0] = fabs(tau[0]); if(_show) if(rank == 0) cout << "[!]"; }
			for(int i = 0; i < (int)alpha.size(); i++) alpha[i] = sigma[i] / tau[i];
			device_add_scale(_u, _s, alpha);
			device_sub_scale(_r, _v, alpha);
			_C(_r, _v);
			_S(_v, _r, sigma);
			if(sigma[0] < 0.0) { sigma[0] = fabs(sigma[0]); if(_show) if(rank == 0) cout << "[!!]"; }
			for(int i = 0; i < (int)beta.size(); i++) beta[i] = sigma[i] / sigma_last[i];
			device_scale_add(_s, _v, beta);
			memcpy(&sigma_last[0], &sigma[0], sizeof(S) * sigma.size());
			for(int i = 0; i < (int)alpha.size(); i++) tau[i] = sigma[i] / sigma_first[i];
		}
		if(_show) if(rank == 0) cout << "**Rank: " << rank << " Loop: " << _loop << " Pre Norm2: " << sqrt(sigma[0]) << " " << sqrt(sigma[0] / sigma_first[0]) << endl;
	}
};
#else
template<class T, class S> class device_conjugate_gradient : public device_inverse_operator<T, S>
{
private:
	const device_generic_operator<T, S> &_K;
	device_inverse_operator<T, S> &_C;
	device_binary_operator<T, S> &_S;
	const S _eps;
	const int _max;
	const bool _show;
	int _loop;
public:
	device_conjugate_gradient(const device_generic_operator<T, S> &_K, device_inverse_operator<T, S> &_C, device_binary_operator<T, S> &_S, const S _eps, const int _max, const bool _show = true)
		: _K(_K), _C(_C), _S(_S), _eps(_eps), _max(_max), _show(_show)
	{
		_loop = 0;
	}
	int iterations()
	{
		return _loop;
	}
	void operator()(const device_vector<S> &_f, device_vector<S> &_u)
	{
		int rank;
		network::rank(rank);
		int n = 1;
		vector<S> alpha(n, 1.0), beta(n, 0.0), gamma(n, 0.0), sigma(n, 0.0), sigma_first(n, 0.0), sigma_last(n, 0.0), tau(n, 0.0);
		//cout << "**" << endl;
		device_vector<S> _r(_f);
		device_vector<S> _s(_u);
		//cout << "**" << endl;
		device_vector<S> _v(_r.size());

		_K(_s, _v);
		device_sub_scale(_r, _v, alpha);
		_S(_r, _r, gamma);//BUG SEQ
		_C(_r, _v);
		_S(_v, _r, sigma);
		if(sigma[0] < 0.0) { sigma[0] = fabs(sigma[0]); if(_show) if(rank == 0) cout << "[!!!]"; }
		device_scale_add(_s, _v, beta);
		
		//memcpy(&sigma_last[0], &sigma[0], sizeof(S) * sigma.size());
		//memcpy(&sigma_first[0], &sigma[0], sizeof(S) * sigma.size());
		for(int i = 0; i < (int)sigma.size(); i++) sigma_last[i] = sigma[i];
		for(int i = 0; i < (int)sigma.size(); i++) sigma_first[i] = sigma[i];
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
			device_add_scale(_u, _s, alpha);
			device_sub_scale(_r, _v, alpha);
			_S(_r, _r, gamma);//BUG SEQ
			_C(_r, _v);
			_S(_v, _r, sigma);
			if(sigma[0] < 0.0) { sigma[0] = fabs(sigma[0]); if(_show) if(rank == 0) cout << "[!!]"; }
			for(int i = 0; i < (int)beta.size(); i++) beta[i] = sigma[i] / sigma_last[i];
			device_scale_add(_s, _v, beta);
			//memcpy(&sigma_last[0], &sigma[0], sizeof(S) * sigma.size());
			for(int i = 0; i < (int)sigma.size(); i++) sigma_last[i] = sigma[i];
			for(int i = 0; i < (int)alpha.size(); i++) tau[i] = sigma[i] / sigma_first[i];
		}
		if(_show) if(rank == 0) cout << "**Rank: " << rank << " Loop: " << _loop << " Pre Norm2: " << sqrt(gamma[0]) << endl;
	}
};
#endif
