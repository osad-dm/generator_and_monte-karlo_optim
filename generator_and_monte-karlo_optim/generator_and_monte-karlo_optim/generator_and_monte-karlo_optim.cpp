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

//������ � ���� ���������� ����� ������������� ����������
struct operators_comb
{
	//q - ����� ������������� ����������
	operators_comb(char q)
	{
		size = q;
		op = new operators_types[size];
	}
	~operators_comb()
	{
		delete op;
	}
	operators_types* op;
	char size;//����� ������������� ����������

	bool operator== (operators_comb other)
	{
		bool equal = true;
		for (char i = 0; i < size; i++)
			if (op[i] != other.op[i])
			{
				equal = false; break;
			}
		return equal;
	}
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
	graph(char ports, char beamsplitters, char directCouplers, char waveplates)
	{
		p = ports;
		bs = beamsplitters;
		dc = directCouplers;
		w = waveplates;
		q = bs + dc + w;
		size = p + 2 * (bs + dc + w);
		edges = new char[size];
		busy = new bool[size];
		for (char i = 0; i < size; i++) busy[i] = false;
		comb = new operators_comb(q);
		var_size = bs + dc + 2 * w;
		var = new double[var_size];
	};
	/*graph(graph other)
	{
		p = other.p;
		bs = other.bs;
		dc = other.dc;
		w = other.w;
		q = other.q;
		size = other.size;
		edges = new char[size];
		for (char i = 0; i < size; i++) edges[i] = other.edges[i];
		busy = new bool[size];
		for (char i = 0; i < size; i++) busy[i] = other.busy[i];
		comb = new operators_comb(q);
		for (char i = 0; i < comb->size; i++) comb->op[i] = other.comb->op[i];
		var_size = bs + dc + 2 * w;
		var = new double[var_size];
		for (char i = 0; i < var_size; i++) var[i] = other.var[i];
	};*/
	~graph()
	{
		delete edges;
		delete busy;
		delete comb;
		delete var;
	};

	//������������ ����
	char *edges;
	//������ � ���� ���� � ���, ����� �� ���� ���� ����������������
	bool *busy;
	//������ � ���� ���������� ����� ������������� ����������
	operators_comb *comb;
	
	char size;//������ ������������� ����� - ����� ������ ��� �� ���������� ��������� ���������
	char p;//p(orts) - ����� ������ �����-������
	char q;//����� ������������� ����������
	char bs;//(beamsplitters) - ����� ���������������� ���������
	char dc;//(direction couplers) - ����� ������������ ��������������
	char w;//(waveplates) - ����� �������� ��������� (��������������)
	char var_size;//����� ����������

	complex<double> *var;//���������� ��������� �����
	
	//���������� ����� ���������� �� ������� var[], ��������������� �������������� ��������� comb->op[oper_num]
	//���� ��� ������� �������������� ��������� ����� ��������� ���������� - �� ����� ���������� ������ �� ������� �� ���
	complex<double>& var_num(char oper_num)
	{
		//�������� ����� ������������ ��������� �� ������� var[]
		char var_num = 0;
		for (char i = 0; i < oper_num; i++) 
			switch (comb->op[i])
			{
			case beamsplitter:
			case directCoupler: var_num++; break;
			case waveplate: var_num += 2; break;
			};
		return var[var_num];
	};
	
	//���������� ������ ���� �������������� ��������� ����������� ���������� � ������� var_num
	operators_types oper_type(char var_num)
	{
		char oper_num = 0;
		while (var_num > 0)
			switch (comb->op[oper_num])
			{
			case beamsplitter:
			case directCoupler:
				var_num--; break;
			case waveplate:
				var_num -= 2;
			}
		return comb->op[oper_num];
	};
};

//���������� ������� ��� ��������� ����
void paths(char start, graph g, set<vector<char>> **matrix, vector<char> way)
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
	switch (g.comb->op[oper_num])
	{
	case beamsplitter:
	{
		complex<double> &tmp = g.var_num(oper_num);
		if (in == 0 && out == 0) return [&tmp](void)->complex<double> {return sqrt(tmp);}; else
		if (in == 0 && out == 1) return [&tmp](void)->complex<double> {return sqrt((complex<double>)1 - tmp);}; else
		if (in == 1 && out == 0) return [&tmp](void)->complex<double> {return sqrt((complex<double>)1 - tmp);}; else
		if (in == 1 && out == 1) return [&tmp](void)->complex<double> {return -sqrt(tmp);};
		break;
	}
	case directCoupler:
	{
		complex<double> &tmp = g.var_num(oper_num);

		if (in == 0 && out == 0) return [&tmp](void)->complex<double> {return sqrt(tmp);}; else
		if (in == 0 && out == 1) return [&tmp](void)->complex<double> {return sqrt((complex<double>)1 - tmp)*exp(complex<double>::complex(0, M_PI / 2));}; else
		if (in == 1 && out == 0) return [&tmp](void)->complex<double> {return sqrt((complex<double>)1 - tmp)*exp(complex<double>::complex(0, M_PI / 2));}; else
		if (in == 1 && out == 1) return [&tmp](void)->complex<double> {return sqrt(tmp);};
		break;
	}
	}
}

//���������� ������� �������� �� ������� ����������. ���� g ����� ��� ����������� ���� �������������� �������� � ����� ��� ��� ���� ���������
func** make_matrix_amplitude(set<vector<char>> **matrix, graph &g)
{
	func **ans = new func*[g.p];
	for (char i = 0; i < g.p; i++) ans[i] = new func[g.p];

	for (char i = 0; i < g.p; i++)
		for (char j = 0; j < g.p; j++)
		{
			ans[i][j] = nullptr;
			if (!matrix[i][j].empty())
			{
				vector<char> one_path = *(matrix[i][j].begin());//��������� ����
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
					} else
					{
						func old_func = summand;
						summand = [old_func, &g, oper_num, in, out]()->complex<double> {return old_func() * get_func(g, oper_num, in, out)();};
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
	return ans;
};

//������ �������� ����� �������� �� ����� ������ ����, ��� ��� ���������
//���������� ��������� ���� �� ���� ��������� ��� ������ ���������
//me - ����, � �������� �� ������ ��������
pair<graph,double> choose_best_graph(char me, graph g)
{
	double best = 0;//����������� ������ ��������� ����� � ������������ �����������
	graph g_best(g);//����, ����������� ������ �������� ���������

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
				set<vector<char>> **matrix = new set<vector<char>>*[g.p];//������� ����������
				for (char i = 0; i < g.p; i++) matrix[i] = new set<vector<char>>[g.p];

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

					set<operators_comb> memory_oper_comb;//��� ���������� ������������� ����������
					for (char i = 0; i < g.bs; i++)				g.comb->op[i] = beamsplitter;
					for (char i = g.bs; i < g.bs + g.dc; i++)	g.comb->op[i] = directCoupler;
					for (char i = g.dc; i < g.q; i++)			g.comb->op[i] = waveplate;

					do
					{
						if (memory_oper_comb.find(*g.comb) == memory_oper_comb.end())
						{
							//� ������ ����� � ��� ���������� ����������� ���� g, � �������� ���������� matrix, � ���������� ����������� ������������� ����������
							memory_oper_comb.insert(*g.comb);//�������� ��� ����������, ����� � ������� �� �����������

							//������ ������� ������� �������� �� �������
							func **MAmpl = make_matrix_amplitude(matrix, g);

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
									for (int k1 = 0; k1 < 5; k1++)//����� ������� ������� ����
										for (int k2 = k1 + 1; k2 < 6; k2++)//����� ������� ������� ����
										{
											restrictions.push_back(
												[&MAmpl, k1, k2]()->complex<double> {
												complex<double> val;
												for (int i = 0; i < 6; i++)
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
								srand(time(0));

								unsigned int total = 0;//����� ����� ��������������� �����
								unsigned int satisfy_restrictions = 0;//����� �����, ��������������� ���� ��������

								
								for (char i = 0; i < g.var_size; i++) g.var[i] = 0;

								//������ ����� �������� ��������� �����
								do {
									//����������� ��������� ����� � ������� ���������� ��������
									for (char i = 0; i < g.var_size; i++)
										switch (g.oper_type(i))
										{
										case beamsplitter:
										case directCoupler:
											g.var[i] = (double)rand() / RAND_MAX;
										case waveplate:
											g.var[i] = (double)rand() / RAND_MAX * 2 * M_PI;
										}

									total++;

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

								for (char i = 0; i < g.p; i++) delete MAmpl;
								delete MAmpl;
							}//����� ���������� ������ �����-�����
						}
					} while (next_permutation(g.comb->op, g.comb->op + g.q));

					//� ������ ����� � ��� ���� ���� g_best � ������������ ����������� �����������, ��� ������� �� �������� � ������������ best
					//���������� ��������� ��� ��� ���� �� ������� ���
				}//����� ��������� �����, ���������������� ������� ����������

				for (char i = 0; i < g.p; i++) delete matrix[i];
				delete matrix;
			}
		}//����� ��������� ���������������� ������������� �����

		return pair<graph, double>::pair(g_best, best);
}

int main(void)
{
	

	return 0;
}