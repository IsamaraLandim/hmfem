/**
* Private procedure, see binary_sort.
*/
template<class T> inline
void _binary_sort(T *_P, T *_Q, T s)
{
	T *P = _P, *Q = _Q;
	T p, q;
	while(P != Q)
	{
		if(((p = P[0]) & s) != 0) break;
		P++;
	}
	while(P != Q)
	{
		if(((q = Q[-1]) & s) == 0) break;
		Q--;
	}
	while(P != Q)
	{
		*P++ = q;
		*--Q = p;
		while(((p = P[0]) & s) == 0) P++;
		while(((q = Q[-1]) & s) != 0) Q--;
	}
	s >>= 1;
	if(s)
	{
		if(_Q - P > 1) _binary_sort(P, _Q, s);
		if(Q - _P > 1) _binary_sort(_P, Q, s);
	}
}

/**
* Private procedure, see binary_sort_copy.
*/
template<class T, class S> inline
void _binary_sort_copy(T *_P, T *_Q, S *_U, S* _V, T s)
{
	T *P = _P, *Q = _Q;
	S *U = _U, *V = _V;
	T p, q;
	while(P != Q)
	{
		if(((p = P[0]) & s) != 0) break;
		P++;
		U++;
	}
	while(P != Q)
	{
		if(((q = Q[-1]) & s) == 0) break;
		Q--;
		V--;
	}
	while(P != Q)
	{
		*P++ = q;
		*--Q = p;
		S u, v;
		u = *U;
		v = *--V;
		*U++ = v;
		*V = u;
		while(((p = P[0]) & s) == 0) P++, U++;
		while(((q = Q[-1]) & s) != 0) Q--, V--;
	}
	s >>= 1;
	if(s)
	{
		if(_Q - P > 1) _binary_sort_copy(P, _Q, U, _V, s);
		if(Q - _P > 1) _binary_sort_copy(_P, Q, _U, V, s);
	}
}

/**
* Private procedure, see binary_sort_copy_copy.
*/
template<class T, class S, class R> inline
void _binary_sort_copy_copy(T *_P, T *_Q, S *_A, S *_B, R *_U, R *_V, T s)
{
	T *P = _P, *Q = _Q;
	S *A = _A, *B = _B;
	R *U = _U, *V = _V;
	T p, q;
	while(P != Q)
	{
		if(((p = P[0]) & s) != 0) break;
		P++;
		A++;
		U++;
	}
	while(P != Q)
	{
		if(((q = Q[-1]) & s) == 0) break;
		Q--;
		B--;
		V--;
	}
	while(P != Q)
	{
		*P++ = q;
		*--Q = p;
		S u, v;
		u = *U;
		v = *--V;
		*U++ = v;
		*V = u;
		S a, b;
		a = *A;
		b = *--B;
		*A++ = b;
		*B = a;
		while(((p = P[0]) & s) == 0) P++, A++, U++;
		while(((q = Q[-1]) & s) != 0) Q--, B--, V--;
	}
	s >>= 1;
	if(s)
	{
		if(_Q - P > 1) _binary_sort_copy_copy(P, _Q, A, _B, U, _V, s);
		if(Q - _P > 1) _binary_sort_copy_copy(_P, Q, _A, B, _U, V, s);
	}
}

/**
* Private procedure, see binary_sort_sort.
*/
template<class T> inline
void _binary_sort_sort(T *_P, T *_Q, T *_A, T *_B, T s, T t)
{
	T *P = _P, *Q = _Q, *A = _A, *B = _B;
	T p, q;
	while(P != Q)
	{
		if(((p = P[0]) & s) != 0) break;
		P++;
		A++;
	}
	while(P != Q)
	{
		if(((q = Q[-1]) & s) == 0) break;
		Q--;
		B--;
	}
	while(P != Q)
	{
		*P++ = q;
		*--Q = p;
		T a, b;
		a = *A;
		b = *--B;
		*A++ = b;
		*B = a;
		while(((p = P[0]) & s) == 0) P++, A++;
		while(((q = Q[-1]) & s) != 0) Q--, B--;
	}
	s >>= 1;
	if(s)
	{
		if(_Q - P > 1) _binary_sort_sort(P, _Q, A, _B, s, t);
		if(Q - _P > 1) _binary_sort_sort(_P, Q, _A, B, s, t);
	}
	else if(t)
	{
		if(_Q - P > 1) _binary_sort_sort(A, _B, P, _Q, t, s);
		if(Q - _P > 1) _binary_sort_sort(_A, B, _P, Q, t, s);
	}
}

/**
* Private procedure, see binary_sort_sort_copy.
*/
template<class T, class S> inline
void _binary_sort_sort_copy(T *_P, T *_Q, T *_A, T *_B, S *_U, S* _V, T s, T t)
{
	T *P = _P, *Q = _Q, *A = _A, *B = _B;
	S *U = _U, *V = _V;
	T p, q;
	while(P != Q)
	{
		if(((p = P[0]) & s) != 0) break;
		P++;
		A++;
		U++;
	}
	while(P != Q)
	{
		if(((q = Q[-1]) & s) == 0) break;
		Q--;
		B--;
		V--;
	}
	while(P != Q)
	{
		*P++ = q;
		*--Q = p;
		T a, b;
		a = *A;
		b = *--B;
		*A++ = b;
		*B = a;
		S u, v;
		u = *U;
		v = *--V;
		*U++ = v;
		*V = u;
		while(((p = P[0]) & s) == 0) P++, A++, U++;
		while(((q = Q[-1]) & s) != 0) Q--, B--, V--;
	}
	s >>= 1;
	if(s)
	{
		if(_Q - P > 1) _binary_sort_sort_copy(P, _Q, A, _B, U, _V, s, t);
		if(Q - _P > 1) _binary_sort_sort_copy(_P, Q, _A, B, _U, V, s, t);
	}
	else if(t)
	{
		if(_Q - P > 1) _binary_sort_sort_copy(A, _B, P, _Q, U, _V, t, s);
		if(Q - _P > 1) _binary_sort_sort_copy(_A, B, _P, Q, _U, V, t, s);
	}
}

/**
* Private procedure, see fractal_sort_sort_copy.
*/
template<class T, class S> inline
void _fractal_sort_sort_copy(T *_P, T *_Q, T *_A, T *_B, S *_U, S* _V, T s, T t)
{
	T *P = _P, *Q = _Q, *A = _A, *B = _B;
	S *U = _U, *V = _V;
	T p, q;
	while(P != Q)
	{
		if(((p = P[0]) & s) != 0) break;
		P++;
		A++;
		U++;
	}
	while(P != Q)
	{
		if(((q = Q[-1]) & s) == 0) break;
		Q--;
		B--;
		V--;
	}
	while(P != Q)
	{
		*P++ = q;
		*--Q = p;
		T a, b;
		a = *A;
		b = *--B;
		*A++ = b;
		*B = a;
		S u, v;
		u = *U;
		v = *--V;
		*U++ = v;
		*V = u;
		while(((p = P[0]) & s) == 0) P++, A++, U++;
		while(((q = Q[-1]) & s) != 0) Q--, B--, V--;
	}
	s >>= 1;
	if(t != 0)
	{
		if(_Q - P > 1) _fractal_sort_sort_copy(A, _B, P, _Q, U, _V, t, s);
		if(Q - _P > 1) _fractal_sort_sort_copy(_A, B, _P, Q, _U, V, t, s);
	}
}

/**
* Private procedure, see binary_sort.
*/
template<class T> inline
T _binary_log(const T *P, const T *Q)
{
	T s = 0;
	while(P != Q)
	{
		s |= *P++;
	}
	T t = ~0;
	while(s & t)
	{
		s &= t;
		t <<= 1;
	}
	if(s < 0)
	{
		cout << "SORT ERROR!" << endl;
		s = 0;
	}
	return s;
}


/**
* The binary_sort procedure sorts a vector of nonnegative integers in place in ascending order.
* \param _V Input: Vector of nonnegative integers to sort. Output: Sorted vector in ascending order.
*/
template<class T> inline
void binary_sort(vector<T> &_V)
{
	if(_V.size() < 2) return;
	_binary_sort(&_V[0], &_V[0]+_V.size(), _binary_log(&_V[0], &_V[0]+_V.size()));
}

/**
* The binary_sort_copy procedure partially sorts pairs in place in ascending order only looking at the first argument.
* \param _V Input: Vector of nonnegative integers to sort, first argument. Output: Partially sorted vector in ascending order, first argument.
* \param _W Input: Vector of arbitrary type to sort, second argument. Output: Partially sorted vector in ascending order, second argument.
*/
template<class T, class S> inline
void binary_sort_copy(vector<T> &_V, vector<S> &_W)
{
	if(_V.size() < 2) return;
	_binary_sort_copy(&_V[0], &_V[0]+_V.size(), &_W[0], &_W[0]+_W.size(), _binary_log(&_V[0], &_V[0]+_V.size()));
}

/**
* The binary_sort_copy_copy procedure partially sorts triples in place in ascending order only looking at the first argument.
* \param _V Input: Vector of nonnegative integers to sort, first argument. Output: Partially sorted vector in ascending order, first argument.
* \param _W Input: Vector of arbitrary type to sort, second argument. Output: Partially sorted vector in ascending order, second argument.
* \param _A Input: Vector of arbitrary type to sort, third argument. Output: Partially sorted vector in ascending order, second argument.
*/
template<class T, class S, class R> inline
void binary_sort_copy_copy(vector<T> &_V, vector<S> &_W, vector<R> &_A)
{
	if(_V.size() < 2) return;
	_binary_sort_copy_copy(&_V[0], &_V[0]+_V.size(), &_W[0], &_W[0]+_W.size(), &_A[0], &_A[0]+_A.size(), _binary_log(&_V[0], &_V[0]+_V.size()));
}

/**
* The binary_sort_sort procedure sorts pairs of nonnegative integers in place in lexicographic order looking at both arguments.
* \param _V Input: Vector of nonnegative integers to sort, first argument. Output: Sorted vector in ascending order, first argument.
* \param _W Input: Vector of nonnegative integers to sort, second argument. Output: Sorted vector in ascending order, second argument.
*/
template<class T> inline
void binary_sort_sort(vector<T> &_V, vector<T> &_W)
{
	if(_V.size() < 2) return;
	_binary_sort_sort(&_V[0], &_V[0]+_V.size(), &_W[0], &_W[0]+_W.size(), 
		_binary_log(&_V[0], &_V[0]+_V.size()), _binary_log(&_W[0], &_W[0]+_W.size()));
}

/**
* The binary_sort_sort_copy procedure partially sorts triples in place in lexicographic order only looking at the first and second argument.
* \param _V Input: Vector of nonnegative integers to sort, first argument. Output: Partially sorted vector in ascending order, first argument.
* \param _W Input: Vector of nonnegative integers to sort, second argument. Output: Partially sorted vector in ascending order, second argument.
* \param _A Input: Vector of arbitrary type to sort, third argument. Output: Partially sorted vector in ascending order, third argument.
*/
template<class T, class S> inline
void binary_sort_sort_copy(vector<T> &_V, vector<T> &_W, vector<S> &_A)
{
	if(_V.size() < 2) return;
	_binary_sort_sort_copy(&_V[0], &_V[0]+_V.size(), &_W[0], &_W[0]+_W.size(), &_A[0], &_A[0]+_A.size(), 
		_binary_log(&_V[0], &_V[0]+_V.size()), _binary_log(&_W[0], &_W[0]+_W.size()));
}

/**
* The fractal_sort_sort_copy procedure partially sorts triples in place in fractal order only looking at the first and second argument.
* \param _V Input: Vector of nonnegative integers to sort, first argument. Output: Partially sorted vector in ascending order, first argument.
* \param _W Input: Vector of nonnegative integers to sort, second argument. Output: Partially sorted vector in ascending order, second argument.
* \param _A Input: Vector of arbitrary type to sort, third argument. Output: Partially sorted vector in ascending order, third argument.
*/
template<class T, class S> inline
void fractal_sort_sort_copy(vector<T> &_V, vector<T> &_W, vector<S> &_A)
{
	if(_V.size() < 2) return;
	T log = max(_binary_log(&_V[0], &_V[0]+_V.size()), _binary_log(&_W[0], &_W[0]+_W.size()));
	_fractal_sort_sort_copy(&_V[0], &_V[0]+_V.size(), &_W[0], &_W[0]+_W.size(), &_A[0], &_A[0]+_A.size(), log, log);
}

/**
* The bucket_sort_count procedure counts how often an index appears in the first vector and stores the counts in the second.
* \param _U Input: Vector of indices.
* \param _W Input: Vector initialized to zero. Output: Sorted Vector of index counts.
*/
template<class T, class S> inline
void bucket_sort_count(const vector<T> &_U, vector<S> &_W)
{
	const T *U = &_U[0], *V = &_U[0] + _U.size();
	S *W = &_W[0];

	while(U != V)
	{
		W[*U++]++;
	}
}

/**
* The bucket_sort_size function calculates the sum over all elements in the vector.
* \param _U Input: Vector of counts.
* \return Output: Sum over all elements in the vector.
*/
template<class T> inline
T bucket_sort_size(const vector<T> &_U)
{
	if(_U.size() == 0) return 0;
	const T *U = &_U[0], *V = &_U[0] + _U.size();

	T t = 0;
	while(U != V)
	{
		t += *U++;
	}
	return t;
}

/**
* The bucket_sort_offset procedure calculates the displacements from a vector of counts.
* \param _U Input: Vector of counts.
* \param _W Output: Vector of displacements.
*/
template<class T> inline
void bucket_sort_offset(const vector<T> &_U, vector<T> &_W)
{
	if(_U.size() == 0) return;
	const T *U = &_U[0], *V = &_U[0] + _U.size();
	T *W = &_W[0];

	T t = 0;
	while(U != V)
	{
		T s = *U++;
		*W++ = t;
		t += s;
	}
}

/**
* The bucket_sort_copy procedure partially sorts pairs in ascending order only looking at the first argument and with partial output.
* \param _U Input: Vector to sort, first argument.
* \param _A Input: Vector to sort, second argument.
* \param _B Output: Sorted vector in ascending order, second argument.
* \param _W Input: Vector of displacements. Output: Vector of final displacements.
* \param _n Input: Size of second argument.
*/
template<class T, class S, class C> inline
void bucket_sort_copy(const vector<T> &_U, const vector<C> &_A,  vector<C> &_B, vector<S> _W, int _n)
{
	const T *U = &_U[0], *V = &_U[0] + _U.size();
	const C *A = &_A[0];
	C *B = &_B[0];
	S *W = &_W[0];

	while(U != V)
	{
		T t = *U++;
		S s = W[t];
		W[t]++;
		B = &_B[_n * s];
		for(int i = 0; i < _n; i++)
		{
			*B++ = *A++;
		}
	}
}

/**
* The unique_resize procedure calculates the set union of an in ascending order sorted vector.
* \param _P Input: Vector sorted in ascending order. Output: Vector sorted in ascending order with only unique elements.
*/
template<class T> inline
void unique_resize(vector<T> &_P)
{
	if(_P.size() < 2) return;
	T *P = &_P[0], *Q = &_P[0] + _P.size();

	if(P != Q)
	{
		T* R = P;
		++P;
		while(P != Q)
		{
			if ((*R != *P))
			{
				*++R = *P;
			}
			++P;
		}
		++R;
		_P.resize(int(R - &_P[0]));
	}
}

/**
* The unique_resize procedure calculates the set union of an in ascending order sorted vector.
* \param _P Input: Vector sorted in ascending order. Output: Vector sorted in ascending order with only unique elements.
* \param _U Input: Vector sorted in ascending order. Output: Vector sorted in ascending order with only unique elements.
*/
template<class T> inline
void unique_resize(vector<T> &_P, vector<T> &_U)
{
	if(_P.size() < 2) return;
	T *P = &_P[0], *Q = &_P[0] + _P.size();
	T *U = &_U[0];

	if(P != Q)
	{
		T* R = P;
		T* W = U;
		++P; ++U;
		while(P != Q)
		{
			if ((*R != *P) || (*W != *U))
			{
				*++R = *P;
				*++W = *U;
			}
			++P; ++U;
		}
		++R;
		++W;
		_P.resize(int(R - &_P[0]));
		_U.resize(int(W - &_U[0]));
	}
}

/**
* The unique_accumulate procedure calculates a partial set union of pairs only looking at the first argument and accumulates the values of the second argument of partially matching pairs.
* \param _P Input: Vector sorted in ascending order, first argument. Output: Vector sorted in ascending order with only unique elements, first argument.
* \param _A Input: Vector to accumulate, second argument. Output: Vector of accumulated values, second argument.
*/
template<class T, class S> inline
void unique_accumulate(vector<T> &_P, vector<S> &_A)
{
	if(_P.size() < 2) return;
	T *P = &_P[0], *Q = &_P[0] + _P.size();
	S *A = &_A[0];

	if(P != Q)
	{
		T* R = P;
		S* C = A;
		++P; ++A;
		while(P != Q)
		{
			if (*R == *P)
			{
				*C += *A;
			}
			else
			{
				*++R = *P;
				*++C = *A;
			}
			++P; ++A;
		}
		++R;
		++C;
		_P.resize(int(R - &_P[0]));
		_A.resize(int(C - &_A[0]));
	}
}

/**
* The unique_accumulate procedure calculates a partial set union of triples only looking at the first and second argument and accumulates the values of the third argument of partially matching triples.
* \param _P Input: Vector sorted in lexicographic order, first argument. Output: Vector sorted in lexicographic order with only unique elements, first argument.
* \param _U Input: Vector sorted in lexicographic order, second argument. Output: Vector sorted in lexicographic order with only unique elements, second argument.
* \param _A Input: Vector to accumulate, third argument. Output: Vector of accumulated values, third argument.
*/
template<class T, class S> inline
void unique_accumulate(vector<T> &_P, vector<T> &_U, vector<S> &_A)
{
	if(_P.size() < 2) return;
	T *P = &_P[0], *Q = &_P[0] + _P.size();
	T *U = &_U[0];
	S *A = &_A[0];

	if(P != Q)
	{
		T* R = P;
		T* W = U;
		S* C = A;
		++P; ++U; ++A;
		while(P != Q)
		{
			if ((*R == *P) && (*W == *U))
			{
				*C += *A;
			}
			else
			{
				*++R = *P;
				*++W = *U;
				*++C = *A;
			}
			++P; ++U; ++A;
		}
		++R;
		++W;
		++C;
		_P.resize(int(R - &_P[0]));
		_U.resize(int(W - &_U[0]));
		_A.resize(int(C - &_A[0]));
	}
}
/**
* The global_intersection procedure calculates the intersection of multiple index sets.
* \param _size Input: Total number of processes.
* \param _rank Input: Rank of the current process.
* \param _P Input: Vector of global node numbers of the current process to intersect.
* \param _U Input: Vector of local node numbers of the current process.
* \param _A Input: Vector of global node numbers for all processes to intersect. Output: Vector of intersecting sets with local node numbers for all processes.
* \param _C Input: Vector of initial counts. Output: Vector of final counts for intersecting sets.
* \param _D Input: Vector of initial displacements. Output: Vector of final displacements for intersecting sets.
*/
template<class T> inline
void global_intersection(const int _size, const int _rank, const vector<T> &_P, const vector<T> &_U, vector<T> &_A, vector<T> &_C, vector<T> &_D)
{
	const T *P = &_P[0], *U = &_U[0];
	T *A = &_A[0], *C = &_C[0], *D = &_D[0];

	int j = 0, k, l;
	int m = (int)_P.size(), n;
	for(int i = 0; i < _size; i++)
	{
		k = 0, l = D[i], n = l + C[i];

		D[i] = j;
		if(_rank != i)
		{
			while((k < m) && (l < n))
			{
				if (P[k] < A[l])
				{
					k++;
				}
				else if (A[l] < P[k])
				{
					l++;
				}
				else
				{
					A[j] = U[k];
					j++;
					k++;
					l++;
				}
			}
		}
		C[i] = j - D[i];
	}
	_A.resize(j);
}
