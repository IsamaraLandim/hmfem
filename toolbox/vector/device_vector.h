template<class T> class device_vector
{
public:
	device_vector()
	{
		_data = 0;
		_capacity = 0;
		_size = 0;
	}

	device_vector(int size)
	{
		if(size != 0) CUDA_SAFE_CALL(cudaMalloc((void**) &_data, size * sizeof(T)));
		else _data = 0;
		_capacity = size;
		_size = size;
		//CUDA_SAFE_CALL(cudaMemset(_data, 0, size * sizeof(T));
	}

	device_vector(int size, T val)
	{
		if(size != 0) CUDA_SAFE_CALL(cudaMalloc((void**) &_data, size * sizeof(T)));
		else _data = 0;
		_capacity = size;
		_size = size;
		if(size != 0)
		{
			T* data = new T[size];
			for(int i = 0; i < size; i++) data[i] = val;
			CUDA_SAFE_CALL(cudaMemcpy(_data, data, size * sizeof(T), cudaMemcpyHostToDevice) );
			delete [] data;
		}
	}
	device_vector(const vector<T> &vec)
	{
		int size = vec.size();
		T* data = vec.data();
		if(size != 0) CUDA_SAFE_CALL(cudaMalloc((void**) &_data, size * sizeof(T)));
		else _data = 0;
		_capacity = size;
		_size = size;
		CUDA_SAFE_CALL(cudaMemcpy(_data, data, size * sizeof(T), cudaMemcpyHostToDevice) );
	}
	device_vector(const device_vector<T> &vec)
	{
		int size = vec.size();
		T* data = vec.data();
		if(size != 0) CUDA_SAFE_CALL(cudaMalloc((void**) &_data, size * sizeof(T)));
		else _data = 0;
		_capacity = size;
		_size = size;
		CUDA_SAFE_CALL(cudaMemcpy(_data, data, size * sizeof(T), cudaMemcpyDeviceToDevice) );
	}
	virtual ~device_vector()
	{
		if(_data != 0) CUDA_SAFE_CALL(cudaFree(_data));
	}
	void upload(T* data, int size) const
	{
		//cout << "UP: " << size << endl;
		if(size != _size) cout << "FATAL UP" << endl;
		size = _size;
		CUDA_SAFE_CALL(cudaMemcpy(_data, data, size * sizeof(T), cudaMemcpyHostToDevice) );
	}
	void download(T* data, int size) const
	{
		//cout << "DOWN: " << size << endl;
		if(size != _size) cout << "FATAL DOWN" << endl;
		size = _size;
		CUDA_SAFE_CALL(cudaMemcpy(data, _data, size * sizeof(T), cudaMemcpyDeviceToHost) );
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
			if(_data != 0) CUDA_SAFE_CALL(cudaFree(_data));
			CUDA_SAFE_CALL(cudaMalloc((void**) &_data, size * sizeof(T)));
            _capacity = size;
        }
		CUDA_SAFE_CALL(cudaMemcpy(_data, s, size * sizeof(T), cudaMemcpyHostToDevice) );
        _size = size;
	}
	void assign(int size)
	{
        if(_capacity < size)
        {
			if(_data != 0) CUDA_SAFE_CALL(cudaFree(_data));
			CUDA_SAFE_CALL(cudaMalloc((void**) &_data, size * sizeof(T)));
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
