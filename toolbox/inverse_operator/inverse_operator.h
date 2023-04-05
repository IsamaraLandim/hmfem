template<class T, class S> class inverse_operator
{
public:
	virtual void operator()(const vector<S> &_f, vector<S> &_u) = 0;
};
