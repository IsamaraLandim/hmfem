/**
* The scalar_product procedure calculates the scalar product for two vectors.
* \param _x Input: Numerical vector, first argument.
* \param _y Input: Numerical vector, second argument.
* \param _s Output: Scalar product.
*/
template<class S>
void scalar_product(const vector<S> &_x, const vector<S> &_y, S &_s)
{
	S s = 0.0;
	const S* x = &_x[0], *y = &_y[0];
	int m = (int)_x.size();
	for(int i = 0; i < m; i++)
	{
		s += *x++ * *y++;
	}
	_s = s;
}
/**
* The scalar_product procedure calculates the scalar product for two packed vectors.
* \param _x Input: Numerical vector, first argument.
* \param _y Input: Numerical vector, second argument.
* \param _s Output: Vector of scalar products.
*/
template<class S>
void scalar_product(const vector<S> &_x, const vector<S> &_y, vector<S> &_s)
{
	S s0, s1, s2, s3;
	const S* x = &_x[0], *y = &_y[0];
	S* s = &_s[0];
	int m = (int)_x.size(), n = (int)_s.size();
	switch(n)
	{
	case 0:
		break;
	case 1:
		s0 = 0.0;
		for(int i = 0; i < m; i++)
		{
			s0 += *x++ * *y++;
			//s0 += x[i] * y[i];
		}
		s[0] = s0;
		break;
	case 2:
		s0 = 0.0, s1 = 0.0;
		for(int i = 0; i < m; i += 2)
		{
			s0 += x[i + 0] * y[i + 0];
			s1 += x[i + 1] * y[i + 1];
		}
		s[0] = s0, s[1] = s1;
		break;
	case 3:
		s0 = 0.0, s1 = 0.0, s2 = 0.0;
		for(int i = 0; i < m; i += 3)
		{
			s0 += x[i + 0] * y[i + 0];
			s1 += x[i + 1] * y[i + 1];
			s2 += x[i + 2] * y[i + 2];
		}
		s[0] = s0, s[1] = s1, s[2] = s2;
		break;
	case 4:
		s0 = 0.0, s1 = 0.0, s2 = 0.0, s3 = 0.0;
		for(int i = 0; i < m; i += 4)
		{
			s0 += x[i + 0] * y[i + 0];
			s1 += x[i + 1] * y[i + 1];
			s2 += x[i + 2] * y[i + 2];
			s3 += x[i + 3] * y[i + 3];
		}
		s[0] = s0, s[1] = s1, s[2] = s2, s[3] = s3;
		break;
	default:
		memset(&_s[0], 0, sizeof(S) * _s.size());
		for(int i = 0; i < m; i += n)
		{
			for(int j = 0; j < n; j++)
			{
				s[j] += x[i + j] * y[i + j];
			}
		}
	}
}
/**
* The add_scale procedure calculates the sum of two vectors, scaling the second vector. x += y * s
* \param _x Input: Numerical vector, first argument. Output: Scaled vector sum. x += y * s
* \param _y Input: Numerical vector, second argument.
* \param _s Input: Scaling factor.
*/
template<class S>
void add_scale(vector<S> &_x, const vector<S> &_y, const S _s)
{
	S *x = &_x[0];
	const S *y = &_y[0];
	const S s = _s;
	for(int i = 0; i < (int)_x.size(); i++)
	{
		*x++ += *y++ * s;
	}
}
/**
* The add_scale procedure calculates the sum of two packed vectors, scaling the second packed vector. x += y * s
* \param _x Input: Numerical vector, first argument. Output: Scaled vector sum. x += y * s
* \param _y Input: Numerical vector, second argument.
* \param _s Input: Vector of scaling factors.
*/
template<class S>
void add_scale(vector<S> &_x, const vector<S> &_y, const vector<S> &_s)
{
	S s0, s1, s2, s3;
	S *x = &_x[0];
	const S *y = &_y[0], *s = &_s[0];
	int m = (int)_x.size(), n = (int)_s.size();
	switch(n)
	{
	case 0:
		break;
	case 1:
		s0 = s[0];
		for(int i = 0; i < m; i++)
		{
			*x++ += *y++ * s0;
			//x[i] += y[i] * s0;
		}
		break;
	case 2:
		s0 = s[0], s1 = s[1];
		for(int i = 0; i < m; i += 2)
		{
			x[i + 0] += y[i + 0] * s0;
			x[i + 1] += y[i + 1] * s1;
		}
		break;
	case 3:
		s0 = s[0], s1 = s[1], s2 = s[2];
		for(int i = 0; i < m; i += 3)
		{
			x[i + 0] += y[i + 0] * s0;
			x[i + 1] += y[i + 1] * s1;
			x[i + 2] += y[i + 2] * s2;
		}
		break;
	case 4:
		s0 = s[0], s1 = s[1], s2 = s[2], s3 = s[3];
		for(int i = 0; i < m; i += 4)
		{
			x[i + 0] += y[i + 0] * s0;
			x[i + 1] += y[i + 1] * s1;
			x[i + 2] += y[i + 2] * s2;
			x[i + 3] += y[i + 3] * s3;
		}
		break;
	default:
		for(int i = 0; i < m; i += n)
		{
			for(int j = 0; j < n; j++)
			{
				x[i + j] += y[i + j] * s[j];
			}
		}
	}
}
/**
* The sub_scale procedure calculates the difference of two vectors, scaling the second vector. x -= y * s
* \param _x Input: Numerical vector, first argument. Output: Scaled vector sum. x -= y * s
* \param _y Input: Numerical vector, second argument.
* \param _s Input: Scaling factor.
*/
template<class S>
void sub_scale(vector<S> &_x, const vector<S> &_y, const S _s)
{
	S *x = &_x[0];
	const S *y = &_y[0];
	const S s = _s;
	for(int i = 0; i < (int)_x.size(); i++)
	{
		*x++ -= *y++ * s;
	}
}
/**
* The sub_scale procedure calculates the difference of two packed vectors, scaling the second packed vector. x -= y * s
* \param _x Input: Numerical vector, first argument. Output: Scaled vector sum. x -= y * s
* \param _y Input: Numerical vector, second argument.
* \param _s Input: Vector of scaling factors.
*/
template<class S>
void sub_scale(vector<S> &_x, const vector<S> &_y, const vector<S> &_s)
{
	S s0, s1, s2, s3;
	S *x = &_x[0];
	const S *y = &_y[0], *s = &_s[0];
	int m = (int)_x.size(), n = (int)_s.size();
	switch(n)
	{
	case 0:
		break;
	case 1:
		s0 = s[0];
		for(int i = 0; i < m; i++)
		{
			*x++ -= *y++ * s0;
			//x[i] -= y[i] * s0;
		}
		break;
	case 2:
		s0 = s[0], s1 = s[1];
		for(int i = 0; i < m; i += 2)
		{
			x[i + 0] -= y[i + 0] * s0;
			x[i + 1] -= y[i + 1] * s1;
		}
		break;
	case 3:
		s0 = s[0], s1 = s[1], s2 = s[2];
		for(int i = 0; i < m; i += 3)
		{
			x[i + 0] -= y[i + 0] * s0;
			x[i + 1] -= y[i + 1] * s1;
			x[i + 2] -= y[i + 2] * s2;
		}
		break;
	case 4:
		s0 = s[0], s1 = s[1], s2 = s[2], s3 = s[3];
		for(int i = 0; i < m; i += 4)
		{
			x[i + 0] -= y[i + 0] * s0;
			x[i + 1] -= y[i + 1] * s1;
			x[i + 2] -= y[i + 2] * s2;
			x[i + 3] -= y[i + 3] * s3;
		}
		break;
	default:
		for(int i = 0; i < m; i += n)
		{
			for(int j = 0; j < n; j++)
			{
				x[i + j] -= y[i + j] * s[j];
			}
		}
	}
}
/**
* The scale_add procedure calculates the sum of two vectors, scaling the first vector. x = x * s + y
* \param _x Input: Numerical vector, first argument. Output: Scaled vector sum. x = x * s + y
* \param _y Input: Numerical vector, second argument.
* \param _s Input: Scaling factor.
*/
template<class S>
void scale_add(vector<S> &_x, const vector<S> &_y, const S _s)
{
	S *x = &_x[0];
	const S *y = &_y[0];
	const S s = _s;
	for(int i = 0; i < (int)_x.size(); i++)
	{
		*x = *x * s + *y;
		x++; y++;
	}
}
/**
* The scale_add procedure calculates the sum of two packed vectors, scaling the first packed vector. x = x * s + y
* \param _x Input: Numerical vector, first argument. Output: Scaled vector sum. x = x * s + y
* \param _y Input: Numerical vector, second argument.
* \param _s Input: Vector of scaling factors.
*/
template<class S>
void scale_add(vector<S> &_x, const vector<S> &_y, const vector<S> &_s)
{
	S s0, s1, s2, s3;
	S *x = &_x[0];
	const S *y = &_y[0], *s = &_s[0];
	int m = (int)_x.size(), n = (int)_s.size();
	switch(n)
	{
	case 0:
		break;
	case 1:
		s0 = s[0];
		for(int i = 0; i < m; i++)
		{
			*x = *x * s0 + *y;
			x++; y++;
			//x[i + 0] = x[i + 0] * s0 + y[i + 0];
		}
		break;
	case 2:
		s0 = s[0], s1 = s[1];
		for(int i = 0; i < m; i += 2)
		{
			x[i + 0] = x[i + 0] * s0 + y[i + 0];
			x[i + 1] = x[i + 1] * s1 + y[i + 1];
		}
		break;
	case 3:
		s0 = s[0], s1 = s[1], s2 = s[2];
		for(int i = 0; i < m; i += 3)
		{
			x[i + 0] = x[i + 0] * s0 + y[i + 0];
			x[i + 1] = x[i + 1] * s1 + y[i + 1];
			x[i + 2] = x[i + 2] * s2 + y[i + 2];
		}
		break;
	case 4:
		s0 = s[0], s1 = s[1], s2 = s[2], s3 = s[3];
		for(int i = 0; i < m; i += 4)
		{
			x[i + 0] = x[i + 0] * s0 + y[i + 0];
			x[i + 1] = x[i + 1] * s1 + y[i + 1];
			x[i + 2] = x[i + 2] * s2 + y[i + 2];
			x[i + 3] = x[i + 3] * s3 + y[i + 3];
		}
		break;
	default:
		for(int i = 0; i < m; i += n)
		{
			for(int j = 0; j < n; j++)
			{
				x[i + j] = x[i + j] * s[j] + y[i + j];
			}
		}
	}
}

/**
* The element_accumulation procedure accumulates the local finite element matrix and calculates the set of global node numbers for the current process.
* \param _hdr Input: Vector of header information.
* \param _con input: Vector of mesh connectivity data. Output: Vector of global node numbers for the current process.
* \param _row Input: Vector of global row indices. Output: Vector of local row indices.
* \param _col Input: Vector of global column indices. Output: Vector of local column indices.
* \param _ele Input: Vector of element matrix values. Output: Vector of the accumulated matrix values.
*/
template<class T, class S>
void element_accumulation(const vector<T> &_hdr, vector<T> &_con, vector<T> &_row, vector<T> &_col, vector<S> &_ele)
{
	T elemsize = _hdr[3];
	//T gnumnode = _hdr[4];

	_row.resize(_ele.size());
	_col.resize(_ele.size());

	for(int i = 0; i < (int)_con.size(); i += elemsize)
	{
		for(int j = 0; j < elemsize; j++)
		{
			for(int k = 0; k < elemsize; k++)
			{
				_row[i * elemsize + j * elemsize + k] = _con[i + j];
				_col[i * elemsize + j * elemsize + k] = _con[i + k];
			}
		}
	}

	binary_sort_sort_copy(_row, _col, _ele);
	unique_accumulate(_row, _col, _ele);

	binary_sort(_con);
	unique_resize(_con);

	int max = 0;
	for(int i = 0; i < (int)_con.size(); i++) if(max < _con[i]) max = _con[i];
	max++;
	vector<T> glo(max, -1);
	for(int i = 0; i < (int)_con.size(); i++) glo[_con[i]] = i;
	for(int i = 0; i < (int)_row.size(); i++) _row[i] = glo[_row[i]];
	for(int i = 0; i < (int)_row.size(); i++) _col[i] = glo[_col[i]];
}
/**
* The diagonal procedure extracts the diagonal from the local finite element matrix.
* \param _con input: Vector of global node numbers for the current process.
* \param _row Input: Vector of local row indices.
* \param _col Input: Vector of local column indices.
* \param _ele Input: Vector of the accumulated matrix values.
* \param _dia Output: Vector of the diagonal matrix values.
*/
template<class T, class S>
void diagonal(const vector<T> &_con, const vector<T> &_row, const vector<T> &_col, const vector<S> &_ele, vector<S> &_dia)
{
	_dia.resize(_con.size());
	for(int i = 0, j = 0; i < (int)_ele.size(); i++)
	{
		if((_row[i] == j) && (_col[i] == j))
		{
			_dia[j] = _ele[i];
			j++;
		}
	}
}
