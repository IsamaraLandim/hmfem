template<class T> class vector
{
public:
	vector()
	{
		_data = 0;
		_capacity = 0;
		_size = 0;
	}
	vector(int size, T val = T())
	{
		if(size != 0) _data = new T[size];
		else _data = 0;
		_capacity = size;
		_size = size;
		for(int i = 0; i < size; i++) _data[i] = val;
	}
	vector(const vector<T> &vec)
	{
		int size = vec.size();
		T* data = vec.data();
		if(size != 0) _data = new T[size];
		else _data = 0;
		_capacity = size;
		_size = size;
		memcpy(_data, data, size * sizeof(T));
	}
	virtual ~vector()
	{
		delete [] _data;
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
			delete [] _data;
			_data = new T[size];
			_capacity = size;
		}
		memcpy(_data, s, size * sizeof(T));
		_size = size;
	}
	void assign(int size, T val = T())
	{
        if(_capacity < size)
		{
			delete [] _data;
			_data = new T[size];
			_capacity = size;
		}
		for(int i = 0; i < size; i++) _data[i] = val;
		_size = size;
	}
    void resize(int size, T val = T())
    {
        if(_capacity < size)
        {
            T* data = _data;
            _data = new T[size];
            memcpy(_data, data, _size * sizeof(T));
            delete [] data;
            _capacity = size;
        }
        for(int i = _size; i < size; i++) _data[i] = val;
        _size = size;
    }

private:
	T* _data;
	int _capacity;
	int _size;
};
