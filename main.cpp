//
//  main.cpp
//  genInput_SSLP
//
//  Created by Yan Deng on 12/9/14.
//  Copyright (c) 2014 noodle. All rights reserved.
//

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <ctime>
#include <ilcplex/ilocplex.h>
#include <time.h>
#include <algorithm>
#include <vector>
#include <math.h>
using namespace std;
ILOSTLBEGIN
typedef IloArray<IloNumVarArray> IloNumVarArray2;
typedef IloArray<IloIntVarArray> IloIntVarArray2;
typedef IloArray<IloModel> IloModelArray;
typedef IloArray<IloObjective> IloObjArray;
typedef IloArray<IloRangeArray> IloRangeArray2;
const int M = 9999;
class param {
public:
	int num_sce;
	int dimx;
	int dimy;
	int dimb;
	int *c;
	int **b;
	int ***T;
	int ***W;
public:
	param(char*);
	~param();
private:
	void parse1(string, int*, int);
	void parse2(string, int**, int, int);
};
void param::parse1(string line, int* arr, int n){
	string o;
	stringstream s(line);
	s >> o;
	for (int i = 0; i<n; i++)
		s >> arr[i];
}
void param::parse2(string line, int** arr, int m, int n){
	string o;
	stringstream s(line);
	s >> o;
	for (int i = 0; i<m; i++)
		for (int j = 0; j<n; j++)
			s >> arr[i][j];
}
param::param(char* s){
	ifstream in;
	in.open(s, std::ifstream::in);
	if (!in){
		cout << "READING ERROR!" << endl;
		exit(0);
	}
	string line;
	getline(in, line);
	parse1(line, &dimx, 1);
	c = new int[dimx];
	getline(in, line);
	parse1(line, c, dimx);
	getline(in, line);
	parse1(line, &num_sce, 1);
	getline(in, line);
	parse1(line, &dimy, 1);
	getline(in, line);
	parse1(line, &dimb, 1);
	T = new int**[num_sce];
	b = new int*[num_sce];
	for (int w = 0; w<num_sce; w++){
		T[w] = new int*[dimb];
		b[w] = new int[dimb];
		for (int i = 0; i<dimb; i++){
			T[w][i] = new int[dimx];
			getline(in, line);
			parse1(line, T[w][i], dimx);
			getline(in, line);
			parse1(line, &b[w][i], 1);
		}
	}
	in.close();
}
param::~param(){
	delete[] c;
	for (int w = 0; w<num_sce; w++){
		delete[] b[w];
		for (int i = 0; i<dimb; i++)
			delete[] T[w][i];
		delete[] T[w];
	}
	delete[] T;
	delete[] b;
}
class feasoln {
public:
	int obj;
	int dim;
	int *x_sol;
	feasoln(int d){
		obj = 0;
		dim = d;
		x_sol = new int[d];
		for (int i = 0; i<d; i++)
			x_sol[i] = -1;
	};
	feasoln(int* x, int* c, int d){
		x_sol = new int[d];
		dim = d;
		obj = 0;
		for (int i = 0; i<d; i++){
			x_sol[i] = x[i];
			obj += c[i] * x[i];
		}
	};
	feasoln(const feasoln &rhs){
		if (this != &rhs){
			obj = rhs.obj;
			dim = rhs.dim;
			x_sol = new int[dim];
			for (int i = 0; i<dim; i++)
				x_sol[i] = rhs.x_sol[i];
		}
	};
	/*
	~feasoln(){
	delete[] x_sol;
	}
	*/
	int numvio(param* data);
	bool same(feasoln comp);
	bool operator< (const feasoln &sub) const{
		return obj < sub.obj;
	};
};
bool feasoln::same(feasoln comp){
	for (int i = 0; i<dim; i++)
		if (comp.x_sol[i] != x_sol[i])
			return false;
	return true;
}
int feasoln::numvio(param* data){
	int vio = 0;
	for (int n = 0; n<data->num_sce; n++){
		for (int i = 0; i<data->dimb; i++){
			int temp = 0;
			for (int d = 0; d<data->dimx; d++)
				temp += data->T[n][i][d] * x_sol[d];
			if (temp< data->b[n][i]){
				vio++;
				break;
			};
		}
	}
	return vio;
};
class sortedsce {
public:
	int sce;
	double h0;
	double h1;
	double delta;
	int *x0_sol;
	int *x1_sol;
	bool operator< (const sortedsce& sub) const{
		return delta < sub.delta;
	}
	sortedsce(int n, int d){
		sce = n;
		x0_sol = new int[d];
		x1_sol = new int[d];
	};
	/*
	~sortedsce(){
	try{
	delete[] x0_sol;
	} catch (IloException& e) {x0_sol=NULL;}
	try{
	delete[] x1_sol;
	} catch (IloException& e) {x1_sol=NULL;}
	};
	*/
};
class nonanti {
public:
	int dimA;
	int ***A;
	int *r;
private:
	int n_sce;
	int d_x;
public:
	nonanti(int n, int d){
		dimA = n*d;
		n_sce = n;
		d_x = d;
		A = new int**[n];
		r = new int[dimA];
		for (int w = 0; w<n; w++){
			A[w] = new int*[dimA];
			for (int i = 0; i<dimA; i++)
				A[w][i] = new int[d];
		}
	};
	~nonanti(){
		// cout << "delete " << dimA << " " << n_sce << " " << d_x << endl;
		delete[] r;
		for (int w = 0; w<n_sce; w++){
			for (int i = 0; i<dimA; i++)
				delete[] A[w][i];
			delete[] A[w];
		}
		delete[] A;
	}
};
class RadarSBG{
public:
	double* lam;
	double   q;
	int* s;
	double m;
	RadarSBG(){};
	double beta(RadarSBG* rhs, int dim){
		double z = 0;
		double mnk = 0;
		for (int i = 0; i<dim; i++) {
			z += (lam[i] - rhs->lam[i])*rhs->s[i];
			mnk += s[i] * rhs->s[i];
		}
		z += rhs->q - q;
		if (mnk>0)
			return -1;
		if (m - mnk != 0)
			return z / (m - mnk);
		else
			return -1;
	}
};
class cut{
public:
	int sparse;
	string vertex;
	string spv;
public:
	cut(){
		sparse = 0;
		vertex = "";
		spv = "";
	};
	cut(feasoln* soln){
		sparse = soln->dim;
		vertex = "";
		spv = "";
		for (int i = 0; i<sparse; i++){
			spv = spv + " ";
			if (1 - soln->x_sol[i] < soln->x_sol[i])
				vertex = vertex + "1";
			else
				vertex = vertex + "0";
		}
	}
};
bool bysparse(const cut& lhs, const cut& rhs){ return lhs.sparse > rhs.sparse; };
bool byspv(const cut& lhs, const cut& rhs){ return (lhs.spv.compare(rhs.spv)>0); };
cut aggregate(cut* c1, cut* c2){
	cut result;
	int num = 0;
	int place = 0;
	for (int i = 0; i<c1->vertex.length(); i++)
		if (c1->vertex[i] != c2->vertex[i]){
		num++;
		place = i;
		}
	if (num == 1){ // aggregate
		result.vertex = c1->vertex;
		result.vertex[place] = '-';
		result.sparse = c1->sparse - 1;
		result.spv = c1->spv;
		result.spv[place] = '-';
	}
	else{
		result.sparse = -1;
	}
	return result;
}

void buildNAP(int, int, nonanti*);
void subgradient(param*, nonanti*, double);
void solveNAP(param*, nonanti*, double);
void robust(param*, feasoln*);
void genInstance(const int, const int, const int, const int);
ofstream out, detail;
int optimalobj;

const double epsilon = .075;
const int DW = 1000;
const int DX = 20;
const int Db = 20;
const int max_coe = 10;
const int SBG_it = 2;
const int num_rep = 10;

const int maxobjvalue = 100;

int main(){
	// genInstance(DW,DX,Db,max_coe);	

	char inputfile[20] = "input.txt";
	param data(inputfile);
	out.open("out.txt", std::ofstream::out);
	detail.open("detail.txt", std::ofstream::out);

	nonanti nap(data.num_sce, data.dimx);
	buildNAP(data.num_sce, data.dimx, &nap);

	subgradient(&data, &nap, epsilon);
	// solveNAP(&data,&nap,epsilon);

	out.close();
	detail.close();
	return 0;
}
void genInstance(int m, int n, int db, int k){
	char filename[50];
	sprintf(filename, "input.txt");
	ofstream out(filename);
	out << "dimx" << "\t" << n << endl << "c";
	for (int j = 0; j<n; j++)
		out << "\t" << 1;
	out << endl << "numsce\t" << m << endl << "dimy\t" << 0 << endl << "dimb\t" << db << endl;
	int* x;
	x = new int[n];
	for (int j = 0; j<n; j++){
		x[j] = rand() % 2;
	}
	srand(time(NULL));
	for (int i = 0; i<m; i++){
		for (int t = 0; t<db; t++){
			out << "T" << i;
			int b = 0;
			for (int j = 0; j<n; j++){
				int a = rand() % (k + 1);
				out << "\t" << a;
				b += a*x[j];
			}
			b = floor(b*0.8);
			out << endl << "b" << i << "\t" << b << endl;
		}
	}
	out.close();
	delete[] x;
}
void buildNAP(int N, int D, nonanti* result){
	for (int i = 0; i<result->dimA; i++)
		result->r[i] = 0;
	for (int w = 0; w<N; w++)
		for (int i = 0; i<result->dimA; i++)
			for (int d = 0; d<D; d++)
				result->A[w][i][d] = 0;
	for (int d = 0; d<D; d++){
		result->A[0][0 + N*d][d] = 1;
		result->A[0][N - 1 + N*d][d] = -1;
	}
	for (int n = 1; n<N; n++)
		for (int d = 0; d<D; d++){
		result->A[n][n + N*d][d] = 1;
		result->A[n][n - 1 + N*d][d] = -1;
		}
}
void solveNAP(param* data, nonanti* nap, double epsilon){
	IloEnv env;
	IloModel model(env);
	IloIntVarArray2 x(env, data->num_sce);
	IloIntVarArray z(env, data->num_sce, 0, 1);
	for (int n = 0; n<data->num_sce; n++){
		char VarName[10];
		x[n] = IloIntVarArray(env, data->dimx, 0, 1);
		for (int d = 0; d<data->dimx; d++){
			sprintf(VarName, "x(%d)(%d)", n, d);
			x[n][d].setName(VarName);
		}
		sprintf(VarName, "z(%d)", n);
		z[n].setName(VarName);
	}

	IloExpr obj(env);
	for (int n = 0; n<data->num_sce; n++)
		for (int d = 0; d<data->dimx; d++)
			obj += data->c[d] * x[n][d] / data->num_sce;
	model.add(IloMinimize(env, obj));

	// formulate nonanticipativity constraints
	for (int i = 0; i<nap->dimA; i++){
		IloExpr left(env);
		for (int n = 0; n<data->num_sce; n++)
			for (int d = 0; d<data->dimx; d++)
				left += nap->A[n][i][d] * x[n][d];
		model.add(left <= nap->r[i]);
	}

	for (int n = 0; n<data->num_sce; n++)
		for (int i = 0; i<data->dimb; i++){
		IloExpr left(env);
		for (int d = 0; d<data->dimx; d++)
			left += data->T[n][i][d] * x[n][d];
		left += M*z[n];
		model.add(left >= data->b[n][i]);
		}

	if (1 == 1){
		IloExpr left(env);
		for (int n = 0; n<data->num_sce; n++)
			left += z[n];
		model.add(left <= floor(epsilon*data->num_sce));
	}
	IloCplex cplex(model);
	cplex.setOut(env.getNullStream());
	cplex.setWarning(env.getNullStream());
	cplex.setWarning(env.getNullStream());

	double ST = cplex.getCplexTime();
	cplex.solve();
	double ET = cplex.getCplexTime();

	if (cplex.getStatus() == IloAlgorithm::Infeasible){
		cout << "INFEASIBLE!" << endl;
		return;
	}
	else if (cplex.getStatus() == IloAlgorithm::Feasible){
		cout << "TIME OUT!" << endl;
		return;
	}
	else{
		cout << "-------- SOLVE MIP PRIMAL -------" << endl << "OPT OBJ = " << IloRound(cplex.getObjValue()) << ", TIME = " << ET - ST << ", OPT SOLN = ";
		for (int i = 0; i<data->dimx; i++)
			cout << IloRound(cplex.getValue(x[0][i]));
		cout << endl << endl;

		out << "-------- SOLVE MIP PRIMAL -------" << endl << "OPT OBJ = " << IloRound(cplex.getObjValue()) << ", TIME = " << ET - ST << ", OPT SOLN = ";
		for (int i = 0; i<data->dimx; i++)
			out << IloRound(cplex.getValue(x[0][i]));
		out << endl << endl;

		optimalobj = IloRound(cplex.getObjValue());
	}
}
void robust(param* data, feasoln *soln){
	IloEnv env;
	IloModel model(env);
	IloIntVarArray x(env, data->dimx, 0, 1);
	char VarName[10];
	for (int d = 0; d<data->dimx; d++){
		sprintf(VarName, "x(%d)", d);
		x[d].setName(VarName);
	}

	IloExpr obj(env);
	for (int d = 0; d<data->dimx; d++)
		obj += data->c[d] * x[d];
	model.add(IloMinimize(env, obj));
	for (int n = 0; n<data->num_sce; n++)
		for (int i = 0; i<data->dimb; i++){
		IloExpr left(env);
		for (int d = 0; d<data->dimx; d++)
			left += data->T[n][i][d] * x[d];
		model.add(left >= data->b[n][i]);
		}
	IloCplex cplex(model);
	cplex.setOut(env.getNullStream());
	cplex.setWarning(env.getNullStream());
	cplex.setWarning(env.getNullStream());
	cplex.solve();

	if (cplex.getStatus() == IloAlgorithm::Infeasible){
		cout << "ROBUST INFEASIBLE!" << endl;
		exit(2);
	}
	else if (cplex.getStatus() == IloAlgorithm::Feasible){
		cout << "ROBUAT TIME OUT!" << endl;
		exit(2);
	}
	else{
		soln->obj = IloRound(cplex.getObjValue());
		for (int d = 0; d<data->dimx; d++)
			soln->x_sol[d] = IloRound(cplex.getValue(x[d]));
	}
}
void subgradient(param* data, nonanti* nap, double epsilon){

	IloEnv env;

	int num_aggregate = 0;

	// initiate upper bound and lower bound
	feasoln bestfea(data->dimx);
	robust(data, &bestfea);
	int LB = -1;
	int UB = bestfea.obj;

	// initialize lambda
	double* lambda;
	lambda = new double[data->num_sce*data->dimx];

	// initialize h0_n(\lambda), h1_n(\lambda)
	IloModelArray h0(env, data->num_sce);
	IloModelArray h1(env, data->num_sce);
	IloIntVarArray2 x0(env, data->num_sce);
	IloIntVarArray2 x1(env, data->num_sce);
	IloObjArray object0(env, data->num_sce);
	IloObjArray object1(env, data->num_sce);
	IloRangeArray objrange0(env, data->num_sce);
	IloRangeArray objrange1(env, data->num_sce);

	for (int n = 0; n<data->num_sce; n++){
		IloModel model0(env);
		IloModel model1(env);
		h0[n] = model0;
		h1[n] = model1;
		x0[n] = IloIntVarArray(env, data->dimx, 0, 1);
		x1[n] = IloIntVarArray(env, data->dimx, 0, 1);
		for (int d = 0; d<data->dimx; d++){
			char name[15];
			sprintf(name, "xsol%d", d);
			x0[n][d].setName(name);
			sprintf(name, "xsol%d", d);
			x1[n][d].setName(name);
		}
		IloObjective o0 = IloMinimize(env, 0);
		IloObjective o1 = IloMinimize(env, 0);
		object0[n] = o0;
		object1[n] = o1;
		h0[n].add(object0[n]);
		h1[n].add(object1[n]);
		for (int i = 0; i<data->dimb; i++){
			IloExpr left(env);
			for (int d = 0; d<data->dimx; d++)
				left += data->T[n][i][d] * x0[n][d];
			h0[n].add(left >= data->b[n][i]);
		}

		IloExpr obj0(env);
		IloExpr obj1(env);
		for (int d = 0; d<data->dimx; d++){
			obj0 += data->c[d] * x0[n][d];
			obj1 += data->c[d] * x1[n][d];
		}
		IloRange range0(env, LB, UB - 1);
		IloRange range1(env, LB, UB - 1);
		range0.setExpr(obj0);
		range1.setExpr(obj1);
		objrange0[n] = range0;
		objrange1[n] = range1;

		h0[n].add(objrange0[n]);
		h1[n].add(objrange1[n]);
	}

	vector<cut> cuts;
	cut c0(&bestfea);
	cuts.push_back(c0);
	IloRangeArray2 cutconstr0(env, data->num_sce);
	IloRangeArray2 cutconstr1(env, data->num_sce);
	for (int n = 0; n<data->num_sce; n++){
		cutconstr0[n] = IloRangeArray(env);
		cutconstr1[n] = IloRangeArray(env);
	}

	vector<sortedsce> slist;
	vector<feasoln> fealist;
	vector<RadarSBG> radarlist;

	IloCplex cplex(env);
	cplex.setParam(IloCplex::Param::Preprocessing::Presolve, IloFalse);
	cplex.setParam(IloCplex::Threads, 1);
	cplex.setParam(IloCplex::MIPSearch, IloCplex::Traditional);
	
	cplex.setOut(env.getNullStream());
	cplex.setWarning(env.getNullStream());
	cplex.setWarning(env.getNullStream());
	double* gv;
	gv = new double[data->num_sce*data->dimx];
	double gamma = 0;
	for (int d = 0; d<data->dimx; d++)
		gamma += data->c[d];
	gamma = gamma / data->dimx;

	int iteration_totSBG = 0;
	int iteration_BND = 0;
	double totTime = 0;
	double timeAgg = 0;
	
	while (UB - LB >= 1){
		iteration_BND++;
		int iteration_SBG = 0;

		// ----------------- FORMULATE CUTS -----------------
		for (int n = 0; n<data->num_sce; n++){
			cutconstr0[n].clear();
			cutconstr1[n].clear();
			for (int ic = 0; ic<cuts.size(); ic++){
				IloExpr expr0(env);
				IloExpr expr1(env);
				for (int d = 0; d<data->dimx; d++)
					if (cuts.at(ic).vertex.at(d) == '1'){
					expr0 += 1 - x0[n][d];
					expr1 += 1 - x1[n][d];
					}
					else if (cuts.at(ic).vertex.at(d) == '0'){
						expr0 += x0[n][d];
						expr1 += x1[n][d];
					}
					cutconstr0[n].add(expr0 >= 1);
					cutconstr1[n].add(expr1 >= 1);
			}
			h0[n].add(cutconstr0[n]);
			h1[n].add(cutconstr1[n]);
		}
		// ----------------- LOWER BOUNDING -----------------
		fealist.clear();
		radarlist.clear();
		for (int i = 0; i<nap->dimA; i++)
			lambda[i] = gamma;

		double lastV = M;
		double currentV = 0;

		// (0) while(iteration_SBG < iteration_BND + 1){
		// (1) while(iteration_SBG < SBG_it){
		while (!((currentV - lastV >= 0) && (currentV - lastV < 0.05))){

			lastV = currentV;

			slist.clear();

			iteration_SBG++;
			detail << iteration_BND << "\t$\t" << iteration_SBG << "\t$";
			// solve sce subproblem
			for (int n = 0; n<data->num_sce; n++){

				cout << n << ",";
				sortedsce sub(n, data->dimx);

				// reformulate subproblem objective
				IloExpr obj0(env);
				IloExpr obj1(env);
				/*
				double* temp;
				temp = new double[data->dimx];
				for (int d = 0; d<data->dimx; d++){
					obj0 += 1.0 / (double)data->num_sce * x0[n][d] ;
					obj1 += 1.0 / (double)data->num_sce * x1[n][d] ;
					temp[d] = 0.0;
					for (int i = 0; i<nap->dimA; i++)
						temp[d] += lambda[i] * nap->A[n][i][d];
					obj0 += temp[d] * x0[n][d];
					obj1 += temp[d] * x1[n][d];
				}
				delete[] temp;
				*/
				
				for (int d = 0; d < data->dimx; d++){
					obj0 += (1.0 / (double)data->num_sce + lambda[n + data->num_sce*d] - lambda[(n - 1 + data->num_sce) % data->num_sce + data->num_sce*d]) * x0[n][d];
					obj1 += (1.0 / (double)data->num_sce + lambda[n + data->num_sce*d] - lambda[(n - 1 + data->num_sce) % data->num_sce + data->num_sce*d]) * x1[n][d];
				}
				
				object0[n].setExpr(obj0);
				object1[n].setExpr(obj1);

				cplex.extract(h0[n]);
				// char modelName[20];
				// sprintf(modelName,"sub_%d.lp",n);
				// cplex.exportModel(modelName);

				double ST = cplex.getCplexTime();
				cplex.solve();
				double ET = cplex.getCplexTime();
				totTime += ET - ST;
				detail << "\t" << ET - ST << " / ";

				if (cplex.getStatus() == IloAlgorithm::Infeasible){
					sub.h0 = M;
					for (int d = 0; d<data->dimx; d++)
						sub.x0_sol[d] = -1;
				}
				else{
					sub.h0 = cplex.getObjValue();
					for (int d = 0; d<data->dimx; d++){
						try{
							sub.x0_sol[d] = IloRound(cplex.getValue(x0[n][d]));
						}
						catch (IloException& e) { sub.x0_sol[d] = 1; }
					}
				}

				feasoln h0soln(sub.x0_sol, data->c, data->dimx);
				if (h0soln.obj >= LB && h0soln.obj <= (UB - 1)){	// otherwise, the solution will be cut off by the obj cuts. 
					bool flag = true;
					for (int k = 0; k<fealist.size(); k++)
						if (fealist.at(k).same(h0soln)){
						flag = false;
						break;
						}
					if (flag)
						fealist.push_back(h0soln);
				}

				cplex.clear();
				cplex.extract(h1[n]);
				ST = cplex.getCplexTime();
				cplex.solve();
				ET = cplex.getCplexTime();
				totTime += ET - ST;
				detail << ET - ST << "\t$";
				if (cplex.getStatus() == IloAlgorithm::Infeasible){
					sub.h1 = M / 2;
					for (int d = 0; d<data->dimx; d++)
						sub.x1_sol[d] = -1;
				}
				else{
					sub.h1 = cplex.getObjValue();
					for (int d = 0; d<data->dimx; d++){
						try {
							sub.x1_sol[d] = IloRound(cplex.getValue(x1[n][d]));
						}
						catch (IloException& e) { sub.x1_sol[d] = 1; }
					}
				}
				cplex.clear();
				sub.delta = sub.h0 - sub.h1;
				slist.push_back(sub);
			}
			cout << endl;

			std::sort(slist.begin(), slist.end());

			// compute dual objective (lower bound);
			double v = 0;
			for (int i = 0; i<nap->dimA; i++)
				v += -lambda[i] * nap->r[i];
			for (int n = 0; n<data->num_sce; n++)
				v += slist.at(n).h1;
			for (int n = 0; n<(data->num_sce - floor(data->num_sce*epsilon)); n++)
				v += slist.at(n).delta;

			currentV = v;

			// compute gradient		
			double gv_n2 = 0;
			for (int i = 0; i<nap->dimA; i++){
				gv[i] = nap->r[i];
				for (int n = (data->num_sce - floor(data->num_sce*epsilon)); n<data->num_sce; n++)
					for (int d = 0; d<data->dimx; d++)
						gv[i] = gv[i] - (nap->A[slist.at(n).sce][i][d])*(slist.at(n).x1_sol[d]);
				for (int n = 0; n<(data->num_sce - floor(data->num_sce*epsilon)); n++)
					for (int d = 0; d<data->dimx; d++)
						gv[i] = gv[i] - (nap->A[slist.at(n).sce][i][d])*(slist.at(n).x0_sol[d]);
				gv_n2 += gv[i] * gv[i];
			}

			// update the best v
			if (v>LB){
				LB = ceil(v);
				if (LB >= UB)
					break;
				std::sort(fealist.begin(), fealist.end());
				if (1 == 1){
					int i;
					for (i = 0; i<fealist.size(); i++)
						if (fealist.at(i).obj >= LB)
							break;
					fealist.erase(fealist.begin(), fealist.begin() + i);
				}
			}

			// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ Radar for deciding step length 		
			/*
			RadarSBG radar;
			radar.lam = new double[nap->dimA];
			radar.s = new int[nap->dimA];
			for (int i=0;i<nap->dimA;i++){
			radar.lam[i] = lambda[i];
			radar.s[i] = -gv[i];
			}
			radar.q = v;
			radar.m = gv_n2;
			double alpha = M;
			for (int p=0;p<radarlist.size();p++){
			double z = radar.beta(&radarlist.at(p),nap->dimA);
			if ((z>0) && (z<alpha))
			alpha = z;
			}
			if (alpha == M)
			alpha = ((UB+LB)/2 - v)/(gv_n2);

			radarlist.push_back(radar);
			*/
			// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

			// update lambda
			// double alpha = 2/(9*data.num_sce*data.dimx);			                 //  (0) constant step size	
			// double alpha = gamma/iteration_SBG;						             //  (1) unnormalized divergent square-summable.  
			// double alpha = gamma/(iteration_SBG*sqrt(gv_n2));		                 // *(2) normalizd divergent square-summable.		(converge well for input-problematic.txt)
			double alpha = ((UB + LB) / 2 - v) / gv_n2;
			// double alpha = ((UB+LB)/2 + gamma/iteration_SBG  - v)/gv_n2;	         // *(3) Polyak (converge not as well as (2))
			// double alpha = (optimalobj - v)/gv_n2;                                   // *(3') Polyak, use the real optimal obj as the esti.
			// double alpha = (UB - v)/(gv_n2);	 //  (4) Polyak' (not converge) 

			/* // Use if want to enforce lambda >= 0
			for (int i=0;i<nap.dimA;i++)
			if (alpha*gv[i]>lambda[i])
			alpha = (lambda[i]/gv[i])*0.95;
			*/

			for (int i = 0; i<nap->dimA; i++)
				lambda[i] = lambda[i] - alpha*gv[i];
			double change = sqrt(gv_n2);

			std::cout << iteration_BND << "\t$\t" << iteration_SBG << "\t$\t" << "UB: " << UB << "; LB: " << LB << "; curV: " << v << "; |grad|: " << change << endl;
			out << iteration_BND << "\t$\t" << iteration_SBG << "\t$\t" << "UB: " << UB << "; LB: " << LB << "; curV: " << v << "; |grad|: " << change << endl;

			detail << endl;
		}

		for (int n = 0; n<data->num_sce; n++){
			h0[n].remove(cutconstr0[n]);
			h1[n].remove(cutconstr1[n]);
		}

		iteration_totSBG += iteration_SBG;

		bool stop = false;
		if (LB <= UB - 1){

			// --------------- UPPER BOUNDING ------------------
			std::sort(fealist.begin(), fealist.end());
			for (int i = 0; i<fealist.size(); i++){
				if (fealist.at(i).numvio(data) <= floor(data->num_sce*epsilon)){
					if (fealist.at(i).obj <= UB - 1){
						bestfea.obj = fealist.at(i).obj;
						for (int j = 0; j<data->dimx; j++)
							bestfea.x_sol[j] = fealist.at(i).x_sol[j];
						UB = bestfea.obj;
						break;
					}
				}
			}

			if (LB <= UB - 1){
				int num_cut = 0;
				// ---------------- CUTTING -------------------
				for (int j = 0; j<fealist.size(); j++){					// add LL cuts
					if (fealist.at(j).obj <= UB - 1){
						num_cut++;
						cut newcut(&fealist.at(j));
						cuts.push_back(newcut);
					}
					else
						break;
				}

				double ST = cplex.getCplexTime();
				sort(cuts.begin(), cuts.end(), bysparse);
				int head = 0;
				int tail;
				while (head<cuts.size() - 1){
					int sp = cuts.at(head).sparse;
					tail = head + 1;
					while (cuts.at(tail).sparse == sp){
						tail++;
						if (tail >= cuts.size())
							break;
					}
					sort(cuts.begin() + head, cuts.begin() + tail, byspv);

					bool* deletecut = new bool[tail - head];
					for (int del = 0; del<tail - head; del++)
						deletecut[del] = false;

					vector<int> index;
					index.push_back(head);
					for (int pt = head + 1; pt<tail; pt++){
						if (cuts.at(pt).spv.compare(cuts.at(pt - 1).spv) != 0)
							index.push_back(pt);
					}
					index.push_back(tail);

					for (int sec = 0; sec<index.size() - 1; sec++){
						for (int start = index.at(sec); start<index.at(sec + 1) - 1; start++)
							if (!deletecut[start - head])
								for (int end = start + 1; end<index.at(sec + 1); end++)
									if (!deletecut[end - head]){
							cut temp = aggregate(&cuts.at(start), &cuts.at(end));
							if (temp.sparse >= 0){
								num_aggregate++;
								cuts.push_back(temp);
								deletecut[end - head] = true;
								deletecut[start - head] = true;
								break;
							}
									}
					}

					for (int del = tail - head - 1; del >= 0; del--)
						if (deletecut[del]){
						tail--;
						cuts.erase(cuts.begin() + head + del, cuts.begin() + head + del + 1);
						}

					sort(cuts.begin() + tail, cuts.end(), bysparse);
					head = tail;
				}
				double ET = cplex.getCplexTime();
				timeAgg += ET-ST;
				for (int n = 0; n<data->num_sce; n++){
					objrange0[n].setBounds(LB, UB - 1);
					objrange1[n].setBounds(LB, UB - 1);
				}

				out << " - " << endl;
				detail << " - \t # added LL cuts:" << num_cut << endl;

			}
			else
				stop = true;
		}
		else
			stop = true;

		if (stop){
			cout << endl << "-------- SOLVE NOO SUBGRADIENT -------" << endl;
			cout << "Iterations: " << iteration_BND << " (BND), " << iteration_totSBG << " (SBG)." << " # Aggregation: " << num_aggregate << endl;
			cout << "OPT OBJ = " << bestfea.obj << ", TIME = " << totTime + timeAgg << ", AggTime = " << timeAgg << ", OPT SOLN = ";
			for (int s = 0; s<data->dimx; s++)
				cout << bestfea.x_sol[s];
			cout << endl << endl << endl;

			out << endl << "-------- SOLVE NOO SUBGRADIENT -------" << endl;
			out << "Iterations: " << iteration_BND << " (BND), " << iteration_totSBG << " (SBG)." << " # Aggregation: " << num_aggregate << endl;
			out << "OPT OBJ = " << bestfea.obj << ", TIME = " << totTime + timeAgg << ", AggTime = " << timeAgg << ", OPT SOLN = ";
			for (int s = 0; s<data->dimx; s++)
				out << bestfea.x_sol[s];
			out << endl;
			break;
		}
	}

	// relea
	delete[] lambda;
	delete[] gv;
}
