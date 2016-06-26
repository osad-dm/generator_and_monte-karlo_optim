// Optim_monte-karlo.cpp : Defines the entry point for the console application.
//������ ���������� ������ ������ ������� ��������� ���� �� ������ ������������ ������, ��������� �� �� ������� ��������
//� ����� ��� ��������������� ���� ������� ���������� ������� ���������� ���������� ����� (����., ����������� ��������������),
//��� ������� ����� �������� � ���������� ������������.

#include "stdafx.h"

using namespace std;

#define func function<complex<double>()>
#define M_PI 3.14159265

//������������ ����� ������������� ����������
enum operators_types
{
	beamsplitter,
	directCoupler,
	waveplate
};

//�������� ��� ���������� �� ��������� �����
struct graph
{
	/*
	p(orts) - ����� ������ �����-������
	cp (coupPlates) - ����� ���������������� ���������
	dc (directionCoupler) - ����� ������������ ��������������
	w(aveplates) - ����� �������� ��������� (��������������)
	*/
	graph(char ports = 0, char beamsplitters = 0, char directCouplers = 0, char waveplates = 0)
	{
		p = ports;
		bs = beamsplitters;
		dc = directCouplers;
		w = waveplates;
		q = bs + dc + w;
		size = p + 2 * (bs + dc + w);
		edges.resize(size);
		busy.resize(size, false);
		comb.resize(q);
		efficiency = 0;
		//������������� comb
		{
			for (char i = 0; i < bs; i++)			comb[i] = beamsplitter;
			for (char i = bs; i < bs + dc; i++)		comb[i] = directCoupler;
			for (char i = bs + dc; i < q; i++)			comb[i] = waveplate;
		}
		var.resize(bs + dc + 2 * w, 0);
	};

	//������������ ����
	vector<char> edges;
	//������ � ���� ���� � ���, � ����� �� ������ ��� ������� �����
	vector<bool> busy;
	//������ � ���� ���������� ����� ������������� ����������
	vector<operators_types> comb;
	
	char size;//������ ������������� �����
	char p;//p(orts) - ����� ������ �����-������
	char q;//����� ������������� ����������
	char bs;//(beamsplitters) - ����� ���������������� ���������
	char dc;//(direction couplers) - ����� ������������ ��������������
	char w;//(waveplates) - ����� �������� ��������� (��������������)
	double efficiency;//������������� ������ �����

	vector<double> var;//���������� ��������� �����
	
	//���������� ��������� �� ���������� �� ������� var[], ��������������� �������������� ��������� comb->op[oper_num]
	//���� ��� ������� �������������� ��������� ����� ��������� ���������� - �� ����� ���������� ������ �� ������ �� ���
	double* var_num(char oper_num)
	{
		//�������� ����� ������������ ��������� �� ������� var[]
		char var_num = 0;
		for (char i = 0; i < oper_num; i++) 
			switch (comb[i])
			{
			case beamsplitter:
			case directCoupler: var_num++; break;
			case waveplate: var_num += 2; break;
			};
		return &(var[var_num]);
	};
	
	//���������� ������ ���� �������������� ��������� ����������� ���������� � ������� var_num
	operators_types oper_type(char var_num)
	{
		char oper_num = 0;
		while (var_num > 0)
			switch (comb[oper_num])
			{
			case beamsplitter:
			case directCoupler:
				var_num--; break;
			case waveplate:
				var_num -= 2;
			}
		return comb[oper_num];
	};
	void operator= (graph &other)
	{
		graph(other.p, other.bs, other.dc, other.w);
		edges = other.edges;
		var = other.var;
		busy = other.busy;
		comb = other.comb;
	};
};

//���������� ������� ��� ��������� ����
void paths(char start, graph &g, vector<vector<set<vector<char>>>> &matrix, vector<char> way)
{
	way.push_back(start);
	way.push_back(g.edges[start]);
	if (g.edges[start] < g.q * 2) //���� ��� ���� ������� � ������������� ��������
	{
		paths((g.edges[start] / 2) * 2, g, matrix, way);
		paths((g.edges[start] / 2) * 2 + 1, g, matrix, way);
	}
	else//���� ��� ���� ������� � ���� ������
	{
		//������� ������������ ���������� � ������
		matrix[way.front() - g.q * 2][way.back() - g.q * 2].insert(way);
	}
};

//���������� �������, ��������� ��� ��������� �� ������ oper_num �� ����� g, � ������ ������ �����-������ �� in � out
func get_func(graph &g, char oper_num, char in, char out)
{
	switch (g.comb[oper_num])
	{
	case beamsplitter:
	{
		double *tmp = g.var_num(oper_num);
		if (in == 0 && out == 0) return [tmp](void)->complex<double> {return sqrt(*tmp);}; else
		if (in == 0 && out == 1) return [tmp](void)->complex<double> {return sqrt(complex<double>::complex(1,0) - *tmp);}; else
		if (in == 1 && out == 0) return [tmp](void)->complex<double> {return sqrt(complex<double>::complex(1,0) - *tmp);}; else
		if (in == 1 && out == 1) return [tmp](void)->complex<double> {return -sqrt(*tmp);};
		break;
	}
	case directCoupler:
	{
		double *tmp = g.var_num(oper_num);

		if (in == 0 && out == 0) return [tmp](void)->complex<double> {return sqrt(*tmp);}; else
		if (in == 0 && out == 1) return [tmp](void)->complex<double> {return sqrt((complex<double>)1 - *tmp)*exp(complex<double>::complex(0, M_PI / 2));}; else
		if (in == 1 && out == 0) return [tmp](void)->complex<double> {return sqrt((complex<double>)1 - *tmp)*exp(complex<double>::complex(0, M_PI / 2));}; else
		if (in == 1 && out == 1) return [tmp](void)->complex<double> {return sqrt(*tmp);};
		break;
	}
	case waveplate:
	{
		double *phi = g.var_num(oper_num);
		double *alpha = phi + 1;

		if (in == 0 && out == 0) return [phi, alpha](void)->complex<double> 
			{return exp(complex<double>::complex(0, *phi))*pow(cos(*alpha), 2) + pow(sin(*alpha),2);}; else
		if (in == 0 && out == 1) return [phi, alpha](void)->complex<double> 
			{return (exp(complex<double>::complex(0,*phi)) - (complex<double>)1)*cos(*alpha)*sin(*alpha);}; else
		if (in == 1 && out == 0) return [phi, alpha](void)->complex<double> 
			{return (exp(complex<double>::complex(0, *phi)) - (complex<double>)1)*cos(*alpha)*sin(*alpha);}; else
		if (in == 1 && out == 1) return [phi, alpha](void)->complex<double> 
			{return exp(complex<double>::complex(0, *phi))*pow(sin(*alpha), 2) + pow(cos(*alpha), 2);};
		break;
	}
	default: return []()->complex<double> {return 0;};
	}
}

//���������� ������� �������� �� ������� ����������. ���� g ����� ��� ����������� ���� �������������� �������� � ����� ��� ��� ���� ���������
vector<vector<func>> make_matrix_amplitude(vector<vector<set<vector<char>>>> &matrix, graph &g)
{
	vector<vector<func>> ans(g.p, vector<func>(g.p));;//������� ����� ����������� �����
	
	for (char i = 0; i < g.p; i++)
		for (char j = 0; j < g.p; j++)
		{
			if (!matrix[i][j].empty())
			{
				while (!matrix[i][j].empty())
				{
					vector<char> one_path = *(matrix[i][j].begin());//��������� ����
					matrix[i][j].erase(matrix[i][j].begin());
					func summand = nullptr;//��������� ���������
					for (char k = 1; k < (char)one_path.size() - 1; k += 2)
					{
						char oper_num = one_path[k] / 2;//���������� ����� �������������� ��������� �� g.comb
						char in = one_path[k] % 2;//����� ����� �����
						char out = one_path[k + 1] % 2;//����� ����� ������
						if (summand == nullptr)
						{
							//� ���� ����� �� ����� ����� �������������� ����� � ��� ����� �����-������
							summand = get_func(g, oper_num, in, out);
						}
						else
						{
							func old_func = summand;
							func new_func = get_func(g, oper_num, in, out);
							summand = [old_func, new_func]()->complex<double> {return old_func() * new_func();};
						}
					}//end for(k)

					//� ������ ����� �� ������������ �������, ������� ������������ ����� ��������� ���������
					if (ans[i][j] == nullptr) ans[i][j] = summand;
					else
					{
						func old_func = ans[i][j];
						ans[i][j] = [old_func, summand]()->complex<double> {return old_func() + summand();};
					}
				}
			}
			else ans[i][j] = []()->complex<double> {return 0;};
		}
	return ans;
};

//���������� ������ �����-�����
//���������� ������������ �������� ������� ������� aim
//g - ������������ ����
//L - ���������� ������� �� �������
//restrictions - ������ ������� ���������
//eps - �������� �������������� �������� ���������
//satisfy - ����� �����, ������� ������ ������������� ���� �������� ��������� - ������� ������ �� �����
double monte_karlo(graph &g, func aiming_function, vector<func> &restrictions, double eps, int satisfy)
{
	//���������� ������� N ��������� ����� � ������������ g.var[]
	//� ���� ��� � ��������� eps ������������� �������� ���������
	//���������� ������ �����-����� - ������ ��������� ������� ����� � ��������� ��� ��� �������
	
	graph g_best(g);//��������� ���� - � ������ ������� ��� �� ���� ����� ������ ���������� � �������� �������������
	bool nonzero = false;//���� ����, ��� ����� �� ������������� ������ �� �������
	srand((unsigned)time(0));

	unsigned int satisfy_restrictions = 0;//����� �����, ��������������� ���� ��������

	//������ ����� �������� ��������� �����
	do {
		//����������� ��������� ����� � ������� ���������� ��������
		for (char i = 0; i < g.var.size(); i++)
			switch (g.oper_type(i))
			{
			case beamsplitter:
			case directCoupler:
				g.var[i] = (double)rand() / RAND_MAX; break;
			case waveplate:
				g.var[i] = (double)rand() / RAND_MAX * 2 * M_PI; break;
			}

		//�������� �� ������� ���������
		for (size_t i = 0; i < restrictions.size(); i++)
			if (abs(restrictions[i]()) > eps)
			{
				nonzero = true;
				break;
			};

		if (nonzero)//���� ���� ���� ������� ��������� �� �����������
		{
			nonzero = false;
			continue;//���� ������ ����� �� ������������� ������ �� ������� ���������, �� ������� ���������
		}

		//� ������ ����� � ��� ���� �����, ��������������� ���� �������� ���������
		satisfy_restrictions++;

		//������ �������� ��, � ����� �������������� �������� ����� � ������������ ����������� ���������� ����������
		double efficiency = abs(aiming_function());
		if (g_best.efficiency < efficiency)
		{
			g_best = g;
			g_best.efficiency = efficiency;
		}
	} while (satisfy_restrictions < satisfy);//����� ���� �������� ����� �� ����� ��������������� �����
	//����� ���������� ������ �����-�����

	g = g_best;
	return g_best.efficiency;
}

//������ �������� ����� �������� �� ����� ������ ����, ��� ��� ���������
//���������� ��������� ���� �� ���� ��������� ��� ������ ���������
//me - ����, � �������� �� ������ ��������
void choose_best_graph(char me, graph g, graph &g_best)
{
	//���������� ��������� ������������ ���� g
	if (me < g.q * 2)//���� �� ������ �������� �� �������������� ���������
	{
		for (char i = (me / 2) * 2 + 2; i < g.edges.size(); i++)//����, ���� ������� ���� ������ �������������� ���������
		{
			if (!g.busy[i])
			{
				g.busy[i] = true;
				g.edges[me] = i;
				choose_best_graph(me + 1, g, g_best);
				g.busy[i] = false;
			}
		}
	}
	else//� ���� ����� �� �������� �� ����� �����
	{
		if (me < g.edges.size())
			//���� ����� ����� �������� � ����� ���� �����
			for (int i = 0; i < g.edges.size(); i++)
			{
				if (!g.busy[i])
				{
					g.busy[i] = true;
					g.edges[me] = i;
					choose_best_graph(me + 1, g, g_best);
					g.busy[i] = false;
				}
			}
		else
		{
			//� ���� ����� �� ��� �������� �� ���� ��������� �������� �����.
			//������ � g ��������� ���������� ������, �������������� ����� ������������ ����. �� �� ��� ���������� �� �������������� � ��������������

			//����������� �������, ������� ������ ��� ���� �� ����� ����� ������ ������ ����� � ������(�) ��������,
			//� �� ���� �������� ������� ������������ �����
			{
				vector<vector<set<vector<char>>>> matrix(g.p, vector<set<vector<char>>>(g.p));//������� ����������

				//�������� ������� ����������
				//�������� �������� ������ ���������� �� ���� ������ �����
				//��������� ����� � matrix
				for (int i = g.q * 2; i < g.size; i++) paths(i, g, matrix, vector<char>::vector());

				//������ � ��� ���� ������� ���������� ��� ���������������� �����. ������ ��������� ������������� �� ���� ���� �������.
				if (matrix[0][2].empty() && matrix[0][3].empty() && //������ ������
					!matrix[1][2].empty() && !matrix[1][3].empty() &&
					!matrix[0][0].empty() && matrix[0][1].empty() && //������ ������ ������
					matrix[1][0].empty() && !matrix[1][1].empty() &&
					matrix[2][0].empty() && !matrix[2][1].empty() && //������ �������� ������
					!matrix[2][2].empty() && !matrix[2][3].empty() &&
					matrix[3][0].empty() && !matrix[3][1].empty() &&
					!matrix[3][2].empty() && !matrix[3][3].empty())
				{
					//� ���� ����� � ��� ���� ���� � �������� ����������, ������� ������������� ���� ��������
					//������ �������� ��� ���������� ����� ������������� ����������

					set<vector<operators_types>> memory_oper_comb;//��� ���������� ������������� ����������

					do
					{
						if (memory_oper_comb.find(g.comb) == memory_oper_comb.end())
						{
							//� ������ ����� � ��� ���������� ����������� ���� g, � �������� ���������� matrix, � ���������� ����������� ������������� ����������
							memory_oper_comb.insert(g.comb);//�������� ��� ����������, ����� � ������� �� �����������

							//������ ������� ������� �������� �� �������
							vector<vector<func>> MAmpl = make_matrix_amplitude(matrix, g);

							//������ �� ���� ������� �������� ���������� �������� ���������� �������, ������� �������� 4x4
							func L[4][4];
							{
								L[0][0] = [&MAmpl]() {return MAmpl[1][1]() * MAmpl[3][3]() + MAmpl[1][3]() * MAmpl[3][1]();};
								L[0][1] = [&MAmpl]() {return MAmpl[1][1]() * MAmpl[2][3]() + MAmpl[1][3]() * MAmpl[2][1]();};
								L[0][2] = [&MAmpl]() {return MAmpl[0][1]() * MAmpl[3][3]() + MAmpl[0][3]() * MAmpl[3][1]();};
								L[0][3] = [&MAmpl]() {return MAmpl[0][1]() * MAmpl[2][3]() + MAmpl[0][3]() * MAmpl[2][1]();};

								L[1][0] = [&MAmpl]() {return MAmpl[1][1]() * MAmpl[3][2]() + MAmpl[1][2]() * MAmpl[3][1]();};
								L[1][1] = [&MAmpl]() {return MAmpl[1][1]() * MAmpl[2][2]() + MAmpl[1][3]() * MAmpl[2][1]();};
								L[1][2] = [&MAmpl]() {return MAmpl[0][1]() * MAmpl[3][2]() + MAmpl[0][2]() * MAmpl[3][1]();};
								L[1][3] = [&MAmpl]() {return MAmpl[0][1]() * MAmpl[2][2]() + MAmpl[0][2]() * MAmpl[2][1]();};

								L[2][0] = [&MAmpl]() {return MAmpl[1][0]() * MAmpl[3][3]() + MAmpl[1][3]() * MAmpl[3][0]();};
								L[2][1] = [&MAmpl]() {return MAmpl[1][0]() * MAmpl[2][3]() + MAmpl[1][3]() * MAmpl[2][0]();};
								L[2][2] = [&MAmpl]() {return MAmpl[0][0]() * MAmpl[3][3]() + MAmpl[0][3]() * MAmpl[3][0]();};
								L[2][3] = [&MAmpl]() {return MAmpl[0][0]() * MAmpl[2][3]() + MAmpl[0][3]() * MAmpl[2][0]();};

								L[3][0] = [&MAmpl]() {return MAmpl[1][0]() * MAmpl[3][2]() + MAmpl[1][2]() * MAmpl[3][0]();};
								L[3][1] = [&MAmpl]() {return MAmpl[1][0]() * MAmpl[2][2]() + MAmpl[1][2]() * MAmpl[2][0]();};
								L[3][2] = [&MAmpl]() {return MAmpl[0][0]() * MAmpl[3][2]() + MAmpl[1][2]() * MAmpl[3][0]();};
								L[3][3] = [&MAmpl]() {return MAmpl[0][0]() * MAmpl[2][2]() + MAmpl[0][2]() * MAmpl[2][0]();};
							}

							//������ ���������� ��������� ������ �������:
							vector<func> restrictions;//����������� ��������� " = 0 " - ��� ��� ������ ���� ����� ����
							{
								//�����������
								{
									//������� ���� ��� ��������.
									for (int k1 = 0; k1 < g.p - 1; k1++)//����� ������� ������� ����
										for (int k2 = k1 + 1; k2 < g.p; k2++)//����� ������� ������� ����
										{
											restrictions.push_back(
												[&MAmpl, k1, k2, g]()->complex<double> {
												complex<double> val;
												for (char i = 0; i < g.p; i++)
													val += MAmpl[i][k1]() * MAmpl[i][k2]();
												return val;
											});//end vector::push_back()
										}
								}//����� �����������

								 //������� �������� � ���������� �������
								{
									restrictions.push_back(L[0][0]);
									restrictions.push_back(L[0][2]);
									restrictions.push_back(L[0][3]);

									restrictions.push_back(L[1][1]);
									restrictions.push_back(L[1][2]);
									restrictions.push_back(L[1][3]);

									restrictions.push_back(L[2][0]);
									restrictions.push_back(L[2][1]);
									restrictions.push_back(L[2][3]);

									restrictions.push_back(L[3][0]);
									restrictions.push_back(L[3][1]);
									restrictions.push_back(L[3][2]);
								}//����� ������� �������� ���������� �������

								 //��������� ����������� ��������� ���������� �������
								{
									restrictions.push_back([&L]()->complex<double> {return L[0][1]() - L[1][0]();});
									restrictions.push_back([&L]()->complex<double> {return L[1][0]() - L[2][2]();});
									restrictions.push_back([&L]()->complex<double> {return L[2][2]() - L[3][3]();});
								}//����� ��������� ����������� ��������� ����� �����
							}

							//� ������ ����� � ��� ���� ������� ��� ���� ������� ��������� � ��������� ��������� ���� - �������� ������ ����� ���������� �������
							if (g_best.efficiency < monte_karlo(g, L[3][3], restrictions, 1e-1, 5))
							{
								g_best = g;
							}
						}
					} while (next_permutation(g.comb.begin(), g.comb.end()));

					//� ������ ����� � ��� ���� ���� g_best � ������������ ����������� �����������, ��� ������� �� �������� � ������������ best
				}//����� ��������� �����, ���������������� �������� ������� ����������
			}
		}
	}//����� ��������� ���������������� ������������� �����
}

//������� ���������� ��������� ������, ��������� ��� ���������� ����, ���� ������� ������ count-������ �����
void make_templates_graphs(char count, graph g, vector<graph> &answer, char deep = 0)
{
	if (deep < count)
	{
		for (char to = 0; to < g.edges.size(); to++)//����� ����, ���� ������� ���� �����
		{
			if (!g.busy[to])
			{
				g.edges[g.q * 2 + deep] = to;
				g.busy[to] = true;
				make_templates_graphs(count, g, answer, deep + 1);
				g.busy[to] = false;
			}
		}
	}
	else answer.push_back(g);
}

int main(void)
{
	//����� ����, ��������� ��� ������� ����������� ���� ����
	if (false)
	{
		graph g(6, 5, 0, 0);
		{
			char tmp_graph[] = { 10, 15, 4, 6, 8, 11, 9, 14, 12, 13, 0, 5, 2, 3, 7, 1 };
			for (char i = 0; i < g.size; i++)
			{
				g.edges[i] = tmp_graph[i];
				g.busy[i] = true;
			}
		}

		vector<vector<set<vector<char>>>> matrix_traj(g.p, vector<set<vector<char>>>(g.p));		
		for (char i = 0; i < g.p; i++) paths(i + g.q * 2, g, matrix_traj, vector<char>::vector());
		
		g.var = {1./3, 1./2, 1./3, 1./3, 1./2};

		vector<vector<func>> MAmpl = make_matrix_amplitude(matrix_traj, g);

		for (char i = 0; i < g.p; i++)
		{
			cout << endl;
			cout << setw(3) << real(MAmpl[i][0]());
			for (char j = 1; j < g.p; j++) cout << '\t' << setw(3) << real(MAmpl[i][j]());
		}

		func L[4][4];
		{
			L[0][0] = [&MAmpl]() {return MAmpl[1][1]() * MAmpl[3][3]() + MAmpl[1][3]() * MAmpl[3][1]();};
			L[0][1] = [&MAmpl]() {return MAmpl[1][1]() * MAmpl[2][3]() + MAmpl[1][3]() * MAmpl[2][1]();};
			L[0][2] = [&MAmpl]() {return MAmpl[0][1]() * MAmpl[3][3]() + MAmpl[0][3]() * MAmpl[3][1]();};
			L[0][3] = [&MAmpl]() {return MAmpl[0][1]() * MAmpl[2][3]() + MAmpl[0][3]() * MAmpl[2][1]();};

			L[1][0] = [&MAmpl]() {return MAmpl[1][1]() * MAmpl[3][2]() + MAmpl[1][2]() * MAmpl[3][1]();};
			L[1][1] = [&MAmpl]() {return MAmpl[1][1]() * MAmpl[2][2]() + MAmpl[1][3]() * MAmpl[2][1]();};
			L[1][2] = [&MAmpl]() {return MAmpl[0][1]() * MAmpl[3][2]() + MAmpl[0][2]() * MAmpl[3][1]();};
			L[1][3] = [&MAmpl]() {return MAmpl[0][1]() * MAmpl[2][2]() + MAmpl[0][2]() * MAmpl[2][1]();};

			L[2][0] = [&MAmpl]() {return MAmpl[1][0]() * MAmpl[3][3]() + MAmpl[1][3]() * MAmpl[3][0]();};
			L[2][1] = [&MAmpl]() {return MAmpl[1][0]() * MAmpl[2][3]() + MAmpl[1][3]() * MAmpl[2][0]();};
			L[2][2] = [&MAmpl]() {return MAmpl[0][0]() * MAmpl[3][3]() + MAmpl[0][3]() * MAmpl[3][0]();};
			L[2][3] = [&MAmpl]() {return MAmpl[0][0]() * MAmpl[2][3]() + MAmpl[0][3]() * MAmpl[2][0]();};

			L[3][0] = [&MAmpl]() {return MAmpl[1][0]() * MAmpl[3][2]() + MAmpl[1][2]() * MAmpl[3][0]();};
			L[3][1] = [&MAmpl]() {return MAmpl[1][0]() * MAmpl[2][2]() + MAmpl[1][2]() * MAmpl[2][0]();};
			L[3][2] = [&MAmpl]() {return MAmpl[0][0]() * MAmpl[3][2]() + MAmpl[1][2]() * MAmpl[3][0]();};
			L[3][3] = [&MAmpl]() {return MAmpl[0][0]() * MAmpl[2][2]() + MAmpl[0][2]() * MAmpl[2][0]();};
		}

		cout << endl;
		for (char i = 0; i < 4; i++)
		{
			cout << endl << setw(3) << real(L[i][0]());
			for (char j = 1; j < 4; j++) cout << '\t' << setw(3) << real(L[i][j]());
		}

		cout << endl << endl;
		graph g_best(g);
		choose_best_graph(g.size, g, g_best);
		cout << endl << "Variables:" << endl;
		for (int i = 0; i < g_best.var.size(); i++) cout << "var[" << i << "]=" << g_best.var[i] << endl;

		cout << "\nLogic matrix:";
		for (char i = 0; i < 4; i++)
		{
			cout << endl << setw(3) << abs(L[i][0]());
			for (char j = 1; j < 4; j++) cout << '\t' << setw(3) << abs(L[i][j]());
		}
		cout << endl;
	}

	//�������� ������ ��������
	{
		for (char p = 6; p <= 6; p++)
			for (char bs = 5; bs <= 5; bs++)
				for (char dc = 0; dc <= 0; dc++)
					for (char w = 0; w <= 0; w++)
					{
						//�������� ������ ��������� - �������� ��� �������� ���� ����� �������� ������ ��� ����� �����
						//��� 16!/13! = 3360 ��� p = 6 � q = 5
						{
							graph g(p, bs, dc, w);
							vector<graph> templates_graphs;
							make_templates_graphs(3, g, templates_graphs);
							
							//� ������ ����� � ��� ���� ���� ��������� ������ - ���������� ��������� �� �� ��������� ������� � ���-�� ������� � ������� �� ��� ����
							int n_threads;//����� ��������� �������
							vector<graph> best_graphs(n_threads, g);//��������� ��������� ������ �� ���� �������
#pragma omp parallel for
							for (int i = 0; i < templates_graphs.size(); i++)//�������� �� ��������� �����
							{
								int my_num_thread;//����� �������� ������
								graph loc_template = templates_graphs[i];
								graph loc_best_graph = best_graphs[my_num_thread];
								choose_best_graph(0, loc_template, loc_best_graph);
								best_graphs[my_num_thread] = loc_best_graph;
#pragma omp master
								if (i % n_threads == 0 && i != 0)//������ ���, ����� �� ��������� �� ��� ������ ����� ���������
								{
									//�� ���������� ������� �������� � ��������� ����
									ofstream f("current_progress", ios::trunc);
									f << "Persentage: " << (double)i / templates_graphs.size() * 100 << "%" << endl;

									//������� ����������� ����� �� �������� �� �������������
									vector<graph> sorted_best_graphs(best_graphs);
									sort(sorted_best_graphs.begin(), sorted_best_graphs.end(),
										[](graph &left, graph &right)->bool {return left.efficiency > right.efficiency;});
									for (int i = 0; i < best_graphs.size(); i++)
									{
										f << "Graph " << i << endl;
										f << "\tEfficiency: " << best_graphs[i].efficiency << endl;
										f << "\tg = {" << best_graphs[i].edges[0];
										for (char j = 1; j < best_graphs[i].edges.size(); j++) f << ", " << best_graphs[i].edges[j];
										f << "}" << endl;
										for (char v = 0; v < best_graphs[i].var.size(); v++)
											f << "var[" << v << "]=" << best_graphs[i].var[v] << endl;
									}
									f.close();
								}
#pragma omp master end
							}
#pragma omp end
						}
					}
	}//����� ���������� ��������� ������� ���������

	return 0;
}