template<class T, class S> class device_binary_operator
{
public:
	virtual void operator()(const device_vector<S> &_u, const device_vector<S> &_v, vector<S> &_w) const = 0;
};
extern "C" {
//void _device_vector_product(float *u, float *v, float *w, int n);
void _device_vector_product(double *u, double *v, double *w, unsigned int n, unsigned int m);
}
template<class T, class S> class device_vector_product : public device_binary_operator<T, S>
{
private:
#ifdef NOZEROCOPY
	device_vector<S> _b;
	vector<S> _a;
#else
	host_vector<S> _a;
	S* dev_a;
#endif
	device_communicator<T, S> &_com;
public:
	//device_vector_product(device_vector<S> &_b, vector<S> &_a, device_communicator<T, S> &_com)
	//	: _b(_b), _a(_a), _com(_com)
	//{
	//}
	device_vector_product(device_communicator<T, S> &_com)
		: _com(_com)
	{
		_a.assign(M);
#ifdef NOZEROCOPY
		_b.assign(M);
#else
		CUDA_SAFE_CALL(cudaHostGetDevicePointer((void**)&dev_a, _a.data(), 0));
#endif
	}
	void operator()(const device_vector<S> &_u, const device_vector<S> &_v, vector<S> &_w) const
	{
#if 1
		unsigned int size = (int)_a.size();
#ifdef NOZEROCOPY
		_device_vector_product(_u.data(), _v.data(), _b.data(), _u.size(), _a.size());
		_b.download(&_a[0], size);
#else
		_device_vector_product(_u.data(), _v.data(), dev_a, _u.size(), _a.size());
		cudaThreadSynchronize();
#endif

		S acc = 0.0;
		for (unsigned int i = 0; i < size; i++) acc += _a[i];
		_w[0] = acc;
		_com.collect(_w);

#else
		_w[0] = 1.0;
#endif
	}
};
