#include <cstring>
#include <iostream>
using namespace std;
#ifdef NOSSE
#define _mm_prefetch //_mm_prefetch
#define _MM_HINT_NTA 0
#else
#include "xmmintrin.h"
#endif

#ifdef NOSTL
#include "vector/vector.h"
#else
#include <vector>
#endif
#ifdef CUDA
#include <cuda_runtime.h>
#if 0
#define CUDA_SAFE_CALL
#else
#  define CUDA_SAFE_CALL( call) do {                                         \
    cudaError err = call;                                                    \
    if( cudaSuccess != err) {                                                \
        fprintf(stderr, "Cuda error in file '%s' in line %i : %s.\n",        \
                __FILE__, __LINE__, cudaGetErrorString( err) );              \
    } } while (0)
#endif
#include "vector/device_vector.h"
#include "vector/host_vector.h"
#endif

#include "vector/algorithm.h"
#include "vector/sort.h"

#ifdef NOMPI
#include <ctime>
#include "network/basic_network.h"
#else
#define MPICH_IGNORE_CXX_SEEK
#include <mpi.h>
#include "network/mpi_network.h"
#endif

#include "communicator/communicator.h"

#include <fstream>
#include <string>
#include <sstream>
#include "vector/file.h"

#include <cmath>
#include <algorithm>

#include "binary_operator/binary_operator.h"
#include "generic_operator/generic_operator.h"
#include "inverse_operator/inverse_operator.h"
#ifdef SUPERLU
#include "inverse_operator/superlu.h"
#else
#ifdef MUMPS
#include "inverse_operator/mumps.h"
#endif
#endif
#include "inverse_operator/amg.h"
#include "inverse_operator/pcg.h"

#ifdef CUDA
#define L 32
#define M (30*2)
#include "communicator/device_communicator.h"
#include "binary_operator/device_binary_operator.h"
#include "generic_operator/device_generic_operator.h"
#include "inverse_operator/device_inverse_operator.h"
#include "inverse_operator/device_amg.h"
#include "inverse_operator/device_pcg.h"
#endif

