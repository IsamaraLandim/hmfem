#define NOMPI
//#define NOSTL
//#define NOSSE
#include "./toolbox/toolbox.h"
#include "variables.h"

int compute(int* _hdr, int* _cnt, int* _col, double* _ele, double* _b, double* _x)
{
	int size, rank;
	network::size(size);
	network::rank(rank);

	int gnumnode = _hdr[0];
	int gmatsize = _hdr[1];

	vector<int> col(gmatsize);
	vector<double> ele(gmatsize);
	for(int i = 0; i < gmatsize; i++) col[i] = _col[i];
	for(int i = 0; i < gmatsize; i++) ele[i] = _ele[i];

	vector<int> nod(gnumnode);
	for(int i = 0; i < (int)nod.size(); i++) nod[i] = i;

	vector<int> cnt(gnumnode);
	vector<double> b(gnumnode);
	vector<double> x(gnumnode);
	for(int i = 0; i < gnumnode; i++) cnt[i] = _cnt[i];
	for(int i = 0; i < gnumnode; i++) b[i] = _b[i];
	for(int i = 0; i < gnumnode; i++) x[i] = _x[i];

	communicator<int, double> com(nod);
	linear_operator<int, double> lin(cnt, col, ele);
	vector_product<int, double> sca(com);

	double epsilon = 0.0, omega = 1.0;
	int max_level = 16;
	int min_nodes = 1;
	bool gauss_seidel = true;
	amg_solver<int, double> amg(nod, cnt, col, ele, max_level, min_nodes, epsilon, omega, gauss_seidel);

	double error = 1e-08;//1e-8;
	int max_iterations = 10000;
	bool show = false;
	conjugate_gradient<int, double> cg(lin, amg, sca, error, max_iterations, show);
	cg(b, x);
    iter_amg = cg.iterations();

    for(int i = 0; i < gnumnode; i++) _x[i] = x[i];

	return 0;
}
