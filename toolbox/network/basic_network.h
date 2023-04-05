/**
 * The network class provides global network communication functions for a single process not using the Message Passing Interface (MPI).
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
	}
	/**
	* The finalize procedure deinitializes the network infrastructure.
	*/
	static void finalize()
	{
	}
	/**
	* The time procedure initializes the network infrastructure.
	* \param _t Output: Time in seconds since an arbitrary time in the past.
	*/
	static void time(double &_t)
	{
		_t = 1.0/double(CLOCKS_PER_SEC)*clock();
	}
	/**
	* The rank procedure returns the rank of the current process.
	* \param _rank Output: Rank of the current process.
	*/
	static void rank(int &_rank)
	{
		_rank = 0;
	}
	/**
	* The size procedure returns the total number of processes.
	* \param _size Output: Total number of processes.
	*/
	static void size(int &_size)
	{
		_size = 1;
	}
	/**
	* The barrier procedure syncronizes all processes.
	*/
	static void barrier()
	{
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
		memcpy(_r, _s, _slen);
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
		memcpy(_r, _s, _slen[0]);
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
		memcpy(_r, _s, _slen);
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
		memcpy(_r, _s, _slen);
	}
    static void scatter(void* _s, int _slen, void* _r, int _rlen, int _rank)
    {
		memcpy(_r, _s, _slen);
    }
    static void scatterv(void* _s, int* _slen, int* _soff, void* _r, int _rlen, int _rank)
    {
		memcpy(_r, _s, _slen[0]);
    }
};
