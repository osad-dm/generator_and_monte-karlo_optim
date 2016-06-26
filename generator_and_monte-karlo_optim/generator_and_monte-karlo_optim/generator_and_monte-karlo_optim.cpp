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
	//������ � ���� ���� � ���, ����� �� ���� ���� ����������������
	vector<bool> busy;
	//������ � ���� ���������� ����� ������������� ����������
	vector<operators_types> comb;
	
	char size;//������ ������������� �����
	char p;//p(orts) - ����� ������ �����-������
	char q;//����� ������������� ����������
	char bs;//(beamsplitters) - ����� ���������������� ���������
	char dc;//(direction couplers) - ����� ������������ ��������������
	char w;//(waveplates) - ����� �������� ��������� (��������������)

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

//������ �������� ����� �������� �� ����� ������ ����, ��� ��� ���������
//���������� ��������� ���� �� ���� ��������� ��� ������ ���������
//me - ����, � �������� �� ������ ��������
graph choose_best_graph(char me, graph g)
{
	double best = 0;//����������� ������ ��������� ����� � ������������ �����������
	graph g_best = g;//����, ����������� ������ �������� ���������

	//���������� ��������� ���������� ������������ ����
	if (me < g.q * 2)//���� �� ������ �������� �� �������������� ���������
	{
		//�������� �� ���� ������������� ���������
		for (int i = (me / 2) * 2 + 2; i < g.size; i++)
		{
			if (!g.busy[i])
			{
				g.busy[i] = true;
				g.edges[me] = i;
				choose_best_graph(me + 1, g);
				g.busy[i] = false;
			}
		}
	}
	else//� ���� ����� �� �������� �� ����� �����
	{
		if (me < g.size)
			//���� ����� ����� �������� � ����� ���� �����
			for (int i = 0; i < g.size; i++)
			{
				if (!g.busy[i])
				{
					g.busy[i] = true;
					g.edges[me] = i;
					choose_best_graph(me + 1, g);
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

							//������ ���������� ������ �����-����� - ���������� ������� N ��������� ����� � ������������ g.var[]
							//� ���� ��� � ��������� eps ������������� �������� ���������
							//���������� ������ �����-����� - ������ ��������� ������� ����� � ��������� ��� ��� �������
							{
								double eps = 1e-1;//�������� ��� ������� ���������
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
									double efficiency = abs(L[3][3]());
									if (best < efficiency)
									{
										best = efficiency;
										g_best = g;
									}
								} while (satisfy_restrictions < 5);//����� ���� �������� ����� �� ����� ��������������� �����
							}//����� ���������� ������ �����-�����
						}
					} while (next_permutation(g.comb.begin(), g.comb.end()));

					//� ������ ����� � ��� ���� ���� g_best � ������������ ����������� �����������, ��� ������� �� �������� � ������������ best
				}//����� ��������� �����, ���������������� �������� ������� ����������
			}
		}
	}//����� ��������� ���������������� ������������� �����

	return g_best;
}

int main(void)
{
	//����� ����, ��������� ��� ������� ����������� ���� ����
	if (true)
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
		g = choose_best_graph(g.size, g);
		cout << endl << "Variables:" << endl;
		for (int i = 0; i < g.var.size(); i++) cout << "var[" << i << "]=" << g.var[i] << endl;

		cout << "\nLogic matrix:";
		for (char i = 0; i < 4; i++)
		{
			cout << endl << setw(3) << abs(L[i][0]());
			for (char j = 1; j < 4; j++) cout << '\t' << setw(3) << abs(L[i][j]());
		}
		cout << endl;
	}

	return 0;
}