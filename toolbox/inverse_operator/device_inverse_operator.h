template<class T, class S> class device_inverse_operator
{
public:
	virtual void operator()(const device_vector<S> &_f, device_vector<S> &_u) = 0;
};
