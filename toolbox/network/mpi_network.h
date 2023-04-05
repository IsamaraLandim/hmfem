/**
 * The network class provides global network communication functions using the Message Passing Interface (MPI).
 * \author Manfred Liebmann
 */
class network
{
public:
	/**
	* The init procedure initializes the network infrastructure.
	* \param _argc Input: Passthrough argument argc from main function.
	* \param _argv Input: Passthrough argument argv from main function.
	*/
	static void init(int &_argc, char** &_argv)
	{
		MPI_Init(&_argc, &_argv);
	}
	/**
	* The finalize procedure deinitializes the network infrastructure.
	*/
	static void finalize()
	{
		MPI_Finalize();
	}
	/**
	* The time procedure initializes the network infrastructure.
	* \param _t Output: Time in seconds since an arbitrary time in the past.
	*/
	static void time(double &_t)
	{
		_t = MPI_Wtime();
	}
	/**
	* The rank procedure returns the rank of the current process.
	* \param _rank Output: Rank of the current process.
	*/
	static void rank(int &_rank)
	{
		MPI_Comm_rank(MPI_COMM_WORLD, &_rank);
	}
	/**
	* The size procedure returns the total number of processes.
	* \param _size Output: Total number of processes.
	*/
	static void size(int &_size)
	{
		MPI_Comm_size(MPI_COMM_WORLD, &_size);
	}
	/**
	* The barrier procedure syncronizes all processes.
	*/
	static void barrier()
	{
		MPI_Barrier(MPI_COMM_WORLD);
	}
	/**
	* The alltoall procedure initiates an all-to-all communication using data blocks of the same size.
	* \param _s Input: Address of the send buffer.
	* \param _slen Input: Number of bytes to send to each process.
	* \param _r Input: Address of the receive buffer .
	* \param _rlen input: Number of bytes to receive from any process.
	*/
	static void alltoall(void* _s, int _slen, void* _r, int _rlen)
	{
		//if(_s == 0 || _r == 0) cout << "ALLTOALL NULL PTR! " << _s << " " << _r << endl;
		MPI_Alltoall(_s, _slen, MPI_BYTE, _r, _rlen, MPI_BYTE, MPI_COMM_WORLD);
	}
	/**
	* The alltoallv procedure initiates an all-to-all communication using data blocks of varying size.
	* \param _s Input: Address of the send buffer.
	* \param _slen Input: Array of the number of bytes to send to each process.
	* \param _soff Input: Array of data displacements in the send buffer for each process.
	* \param _r Input: Address of the receive buffer .
	* \param _rlen input: Array of the number of bytes to receive from any process.
	* \param _roff input: Array of data displacements in the receive buffer for any process.
	*/
	static void alltoallv(void* _s, int* _slen, int* _soff, void* _r, int* _rlen, int* _roff)
	{
		//if(_s == 0 || _r == 0) cout << "ALLTOALLV NULL PTR! " << _s << " " << _r << endl;
		MPI_Alltoallv(_s, _slen, _soff, MPI_BYTE, _r, _rlen, _roff, MPI_BYTE, MPI_COMM_WORLD);
	}
	/**
	* The allgather procedure initiates an all-gather communication using data blocks of the same size.
	* \param _s Input: Address of the send buffer.
	* \param _slen Input: Number of bytes to send to each process.
	* \param _r Input: Address of the receive buffer .
	* \param _rlen input: Number of bytes to receive from any process.
	*/
	static void allgather(void* _s, int _slen, void* _r, int _rlen)
	{
		//if(_s == 0 || _r == 0) cout << "ALLGATHER NULL PTR! " << _s << " " << _r << endl;
		MPI_Allgather(_s, _slen, MPI_BYTE, _r, _rlen, MPI_BYTE, MPI_COMM_WORLD);
	}
	/**
	* The allgatherv procedure initiates an all-gather communication using data blocks of varying size.
	* \param _s Input: Address of the send buffer.
	* \param _slen Input: Number of bytes to send to each process.
	* \param _r Input: Address of the receive buffer .
	* \param _rlen input: Array of the number of bytes to receive from any process.
	* \param _roff input: Array of data displacements in the receive buffer for any process.
	*/
	static void allgatherv(void* _s, int _slen, void* _r, int* _rlen, int* _roff)
	{
		//if(_s == 0 || _r == 0) cout << "ALLGATHERV NULL PTR! " << _s << " " << _r << endl;
		MPI_Allgatherv(_s, _slen, MPI_BYTE, _r, _rlen, _roff, MPI_BYTE, MPI_COMM_WORLD);
	}
	static void scatter(void* _s, int _slen, void* _r, int _rlen, int _rank)
	{
		MPI_Scatter(_s, _slen, MPI_BYTE, _r, _rlen, MPI_BYTE, _rank, MPI_COMM_WORLD);
	}
	static void scatterv(void* _s, int* _slen, int* _soff, void* _r, int _rlen, int _rank)
	{
		MPI_Scatterv(_s, _slen, _soff, MPI_BYTE, _r, _rlen, MPI_BYTE, _rank, MPI_COMM_WORLD);
	}

};
