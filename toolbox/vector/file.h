/**
* The binary_read procedure reads binary data from a file into a memory block.
* \param _First Input: Iterator pointing to the first element.
* \param _Last Input: Iterator pointing past the last element.
* \param _Filename Input: Name of the input file.
* \param _Seek Input: Offset into the data file
*/
template<class _InIt> inline
void binary_read(_InIt _First, _InIt _Last, const string _Filename, const int _Seek, const bool _Seek_Char = false)
{
	ifstream is(_Filename.c_str(), ios::binary);
	if(_Seek_Char)
	{
		is.seekg(_Seek);
	}
	else
	{
		is.seekg(_Seek * sizeof(*_First));
	}
	is.read((char*)&*_First, streamsize(_Last - _First) * sizeof(*_First));
	is.close();
}

/**
* The binary_write procedure writes binary data to a file from a memory block.
* \param _First Input: Iterator pointing to the first element.
* \param _Last Input: Iterator pointing past the last element.
* \param _Filename Input: Name of the output file.
*/
template<class _InIt> inline
void binary_write(_InIt _First, _InIt _Last, const string _Filename)
{
	ofstream os(_Filename.c_str(), ios::binary);
	os.write((char*)&*_First, streamsize(_Last - _First) * sizeof(*_First));
	os.close();
}

/**
* The read_header procedure reads the header information of the finite element simulation.
* \param _root Input: Directory where the header is located.
* \param _hdr Output: Vector of header information.
*/
template<class T>
void read_header(const string &_root, vector<T> &_hdr)
{
	_hdr.assign(8, 0);
	binary_read(_hdr.begin(), _hdr.end(), _root + "Header.bin", 0);

	T magic = _hdr[0];
	T version = _hdr[1];
	T gnumelem = _hdr[2];
	T elemsize = _hdr[3];
	T gnumnode = _hdr[4];
	T nodesize = _hdr[5];
	T gnumdir = _hdr[6];
	T dirsize = _hdr[7];

	int rank;
	network::rank(rank);
	if(rank == 0)
	{
		cout << "Magic: " << magic << " Version: " << version << endl;
		cout << "Elements: " << gnumelem << " Size: " << elemsize << endl;
		cout << "Nodes: " << gnumnode << " Size: " << nodesize << endl;
		cout << "Dirichlet: " << gnumdir << " Size: " << dirsize << endl;
	}
}

/**
* The read_partition procedure reads the partition data of the finite element simulation.
* \param _root Input: Directory where the header is located.
* \param _hdr Input: Vector of header information.
* \param _par Output: Vector of partition data.
*/
template<class T>
void read_partition(const string &_root, const vector<T> &_hdr, vector<T> &_par)
{
	T gnumelem = _hdr[2];

	int size, rank;
	network::size(size);
	network::rank(rank);
	int poselem = (rank * gnumelem)/size;
	int numelem = (rank * gnumelem + gnumelem)/size - poselem;

	_par.assign(numelem, 0);
	stringstream s;
	s << _root << "Partition" << size << ".bin";
	binary_read(_par.begin(), _par.end(), s.str(), poselem);
}

/**
* The read_connection procedure reads the mesh connectivity data of the finite element simulation.
* \param _root Input: Directory where the header is located.
* \param _hdr Input: Vector of header information.
* \param _par Input: Vector of partition data.
* \param _con Output: Vector of mesh connectivity data.
*/
template<class T>
void read_connection(const string &_root, const vector<T> &_hdr, const vector<T> &_par, vector<T> &_con)
{
	T gnumelem = _hdr[2];
	T elemsize = _hdr[3];

	int size, rank;
	network::size(size);
	network::rank(rank);
	int poselem = (rank * gnumelem)/size;
	int numelem = (rank * gnumelem + gnumelem)/size - poselem;
	vector<T> con(numelem * elemsize);
	binary_read(con.begin(), con.end(), _root + "Connection.bin", poselem * elemsize);

	vector<int> slen(size, 0);
	vector<int> soff(size);
	vector<T> scon(numelem * elemsize);
	bucket_sort_count(_par, slen);
	bucket_sort_offset(slen, soff);
	bucket_sort_copy(_par, con, scon, soff, elemsize);

	vector<int> rlen(size, 0);
	vector<int> roff(size);
	network::alltoall(&slen[0], sizeof(int), &rlen[0], sizeof(int));
	bucket_sort_offset(rlen, roff);
	int rnumelem = bucket_sort_size(rlen);
	//_con.assign(rnumelem * elemsize, 0);
	_con.resize(rnumelem * elemsize);

	for(int i = 0; i < (int)slen.size(); i++) slen[i] *= elemsize * sizeof(T), soff[i] *= elemsize * sizeof(T);
	for(int i = 0; i < (int)slen.size(); i++) rlen[i] *= elemsize * sizeof(T), roff[i] *= elemsize * sizeof(T);
	network::alltoallv(&scon[0], &slen[0], &soff[0], &_con[0], &rlen[0], &roff[0]);
}

/**
* The read_element procedure reads the element matrix data of the finite element simulation.
* \param _root Input: Directory where the header is located.
* \param _hdr Input: Vector of header information.
* \param _par Input: Vector of partition data.
* \param _ele Output: Vector of element matrix data.
*/
template<class T, class S>
void read_element(const string &_root, const vector<T> &_hdr, const vector<T> &_par, vector<S> &_ele)
{
	T gnumelem = _hdr[2];
	T elemsize = _hdr[3];

	int size, rank;
	network::size(size);
	network::rank(rank);
	int poselem = (rank * gnumelem)/size;
	int numelem = (rank * gnumelem + gnumelem)/size - poselem;
	vector<S> ele(numelem * elemsize * elemsize);
	binary_read(ele.begin(), ele.end(), _root + "Element.bin", poselem * elemsize * elemsize);

	vector<int> slen(size, 0);
	vector<int> soff(size);
	vector<S> sele(numelem * elemsize * elemsize);
	bucket_sort_count(_par, slen);
	bucket_sort_offset(slen, soff);
	bucket_sort_copy(_par, ele, sele, soff, elemsize * elemsize);

	vector<int> rlen(size, 0);
	vector<int> roff(size);
	network::alltoall(&slen[0], sizeof(int), &rlen[0], sizeof(int));
	bucket_sort_offset(rlen, roff);
	int rnumelem = bucket_sort_size(rlen);
	//_ele.assign(rnumelem * elemsize * elemsize, 0.0);
	_ele.resize(rnumelem * elemsize * elemsize);

	for(int i = 0; i < (int)slen.size(); i++) slen[i] *= elemsize * elemsize * sizeof(S), soff[i] *= elemsize * elemsize * sizeof(S);
	for(int i = 0; i < (int)slen.size(); i++) rlen[i] *= elemsize * elemsize * sizeof(S), roff[i] *= elemsize * elemsize * sizeof(S);
	network::alltoallv(&sele[0], &slen[0], &soff[0], &_ele[0], &rlen[0], &roff[0]);
}

/**
* The read_element_data procedure reads the element matrix data of the finite element simulation.
* \param _file Input: File where the element data is located.
* \param _hdr Input: Vector of header information.
* \param _par Input: Vector of partition data.
* \param _ele Output: Vector of element matrix data.
*/
template<class T, class S>
void read_element_data(const string &_file, const vector<T> &_hdr, const vector<T> &_par, vector<S> &_ele)
{
	T gnumelem = _hdr[2];
	T elemsize = _hdr[3];

	int size, rank;
	network::size(size);
	network::rank(rank);
	int poselem = (rank * gnumelem)/size;
	int numelem = (rank * gnumelem + gnumelem)/size - poselem;
	vector<S> ele(numelem * elemsize * elemsize);
	binary_read(ele.begin(), ele.end(), _file, poselem * elemsize * elemsize);

	vector<int> slen(size, 0);
	vector<int> soff(size);
	vector<S> sele(numelem * elemsize * elemsize);
	bucket_sort_count(_par, slen);
	bucket_sort_offset(slen, soff);
	bucket_sort_copy(_par, ele, sele, soff, elemsize * elemsize);

	vector<int> rlen(size, 0);
	vector<int> roff(size);
	network::alltoall(&slen[0], sizeof(int), &rlen[0], sizeof(int));
	bucket_sort_offset(rlen, roff);
	int rnumelem = bucket_sort_size(rlen);
	//_ele.assign(rnumelem * elemsize * elemsize, 0.0);
	_ele.resize(rnumelem * elemsize * elemsize);

	for(int i = 0; i < (int)slen.size(); i++) slen[i] *= elemsize * elemsize * sizeof(S), soff[i] *= elemsize * elemsize * sizeof(S);
	for(int i = 0; i < (int)slen.size(); i++) rlen[i] *= elemsize * elemsize * sizeof(S), roff[i] *= elemsize * elemsize * sizeof(S);
	network::alltoallv(&sele[0], &slen[0], &soff[0], &_ele[0], &rlen[0], &roff[0]);
}

/**
* The read_rhs_data procedure reads the element matrix data of the finite element simulation.
* \param _file Input: File where the element data is located.
* \param _hdr Input: Vector of header information.
* \param _par Input: Vector of partition data.
* \param _rhs Output: Vector of right hand side data.
*/
template<class T, class S>
void read_rhs_data(const string &_file, const vector<T> &_hdr, const vector<T> &_par, vector<S> &_rhs)
{
	T gnumelem = _hdr[2];
	T elemsize = _hdr[3];

	int size, rank;
	network::size(size);
	network::rank(rank);
	int poselem = (rank * gnumelem)/size;
	int numelem = (rank * gnumelem + gnumelem)/size - poselem;
	vector<S> rhs(numelem * elemsize);
	binary_read(rhs.begin(), rhs.end(), _file, poselem * elemsize);

	vector<int> slen(size, 0);
	vector<int> soff(size);
	vector<S> srhs(numelem * elemsize);
	bucket_sort_count(_par, slen);
	bucket_sort_offset(slen, soff);
	bucket_sort_copy(_par, rhs, srhs, soff, elemsize);

	vector<int> rlen(size, 0);
	vector<int> roff(size);
	network::alltoall(&slen[0], sizeof(int), &rlen[0], sizeof(int));
	bucket_sort_offset(rlen, roff);
	int rnumelem = bucket_sort_size(rlen);
	//_rhs.assign(rnumelem * elemsize, 0.0);
	_rhs.resize(rnumelem * elemsize);

	for(int i = 0; i < (int)slen.size(); i++) slen[i] *= elemsize * sizeof(S), soff[i] *= elemsize * sizeof(S);
	for(int i = 0; i < (int)slen.size(); i++) rlen[i] *= elemsize * sizeof(S), roff[i] *= elemsize * sizeof(S);
	network::alltoallv(&srhs[0], &slen[0], &soff[0], &_rhs[0], &rlen[0], &roff[0]);
}
