template<class T, class S> class binary_operator
{
public:
	virtual void operator()(const vector<S> &_u, const vector<S> &_v, vector<S> &_w) const = 0;
};

template<class T, class S> class vector_product : public binary_operator<T, S>
{
private:
	communicator<T, S> &_com;
public:
	vector_product(communicator<T, S> &_com)
		: _com(_com)
	{
	}
	void operator()(const vector<S> &_u, const vector<S> &_v, vector<S> &_w) const
	{
		scalar_product(_u, _v, _w);
		_com.collect(_w);
	}
};
