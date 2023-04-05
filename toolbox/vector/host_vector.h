#ifndef NOZEROCOPY
template<class T> class host_vector
{
public:
	host_vector()
	{
		_data = 0;
		_capacity = 0;
		_size = 0;
	}
	host_vector(int size)
	{
		//CUDA_SAFE_CALL(cudaMallocHost((void**) &_data, size * sizeof(T)));
		if(size != 0) CUDA_SAFE_CALL(cudaHostAlloc((void**) &_data, size * sizeof(T), cudaHostAllocMapped));
		else _data = 0;
		_capacity = size;
		_size = size;
		//for(int i = 0; i < size; i++) _data[i] = 0;
	}
	host_vector(int size, T val)
	{
		//CUDA_SAFE_CALL(cudaMallocHost((void**) &_data, size * sizeof(T)));
		if(size != 0) CUDA_SAFE_CALL(cudaHostAlloc((void**) &_data, size * sizeof(T), cudaHostAllocMapped));
		else _data = 0;
		_capacity = size;
		_size = size;
		for(int i = 0; i < size; i++) _data[i] = val;
	}
	host_vector(const vector<T> &vec)
	{
		int size = vec.size();
		T* data = vec.data();
		//CUDA_SAFE_CALL(cudaMallocHost((void**) &_data, size * sizeof(T)));
		if(size != 0) CUDA_SAFE_CALL(cudaHostAlloc((void**) &_data, size * sizeof(T), cudaHostAllocMapped));
		else _data = 0;
		_capacity = size;
		_size = size;
		CUDA_SAFE_CALL(cudaMemcpy(_data, data, size * sizeof(T), cudaMemcpyHostToHost) );
	}
	host_vector(const device_vector<T> &vec)
	{
		int size = vec.size();
		T* data = vec.data();
		//CUDA_SAFE_CALL(cudaMallocHost((void**) &_data, size * sizeof(T)));
		if(size != 0) CUDA_SAFE_CALL(cudaHostAlloc((void**) &_data, size * sizeof(T), cudaHostAllocMapped));
		else _data = 0;
		_capacity = size;
		_size = size;
		CUDA_SAFE_CALL(cudaMemcpy(_data, data, size * sizeof(T), cudaMemcpyDeviceToHost) );
	}
	virtual ~host_vector()
	{
		if(_data != 0) CUDA_SAFE_CALL(cudaFreeHost(_data));
	}
	T& operator[](int i) const
	{
		return _data[i];
	}
	T* data() const
	{
		return _data;
	}
	int size() const
	{
		return _size;
	}
	int capacity() const
	{
		return _capacity;
	}
	T* begin() const
	{
		return _data;
	}
	T* end() const
	{
		return _data + _size;
	}
	void assign(T* s, T* e)
	{
		int size = (int)(e - s);
		if(_capacity < size)
		{
			if(_data != 0) CUDA_SAFE_CALL(cudaFreeHost(_data));
			//CUDA_SAFE_CALL(cudaMallocHost((void**) &_data, size * sizeof(T)));
			CUDA_SAFE_CALL(cudaHostAlloc((void**) &_data, size * sizeof(T), cudaHostAllocMapped));
			_capacity = size;
		}
		CUDA_SAFE_CALL(cudaMemcpy(_data, s, size * sizeof(T), cudaMemcpyHostToHost) );
		_size = size;

	}
	void assign(int size)
	{
		if(_capacity < size)
		{
			if(_data != 0) CUDA_SAFE_CALL(cudaFreeHost(_data));
			//CUDA_SAFE_CALL(cudaMallocHost((void**) &_data, size * sizeof(T)));
			CUDA_SAFE_CALL(cudaHostAlloc((void**) &_data, size * sizeof(T), cudaHostAllocMapped));
			_capacity = size;
		}
		//CUDA_SAFE_CALL(cudaMemset(_data, 0, size * sizeof(T));
		_size = size;
	}
private:
	T* _data;
	int _capacity;
	int _size;
};
#endif

