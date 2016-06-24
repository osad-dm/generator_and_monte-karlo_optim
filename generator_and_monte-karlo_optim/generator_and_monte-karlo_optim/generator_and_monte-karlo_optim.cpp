// Optim_monte-karlo.cpp : Defines the entry point for the console application.
//������ ���������� ������ ������ ������� ��������� ���� �� ������ ������������ ������, ��������� �� �� ������� ��������
//� ����� ��� ��������������� ���� ������� ���������� ������� ���������� ���������� ����� (����., ����������� ��������������),
//��� ������� ����� �������� � ���������� ������������.

#include "stdafx.h"

using namespace std;

#define func function<complex<double>()>
#define usint unsigned short int

usint cp;//����� ���������������� ���������
usint wp;//����� �������� ���������
usint p;//����� ������ �����-������
usint graph_size;//������ �����

//������ ��� �������� edges
struct edges_status
{
	usint p;
	usint c;
	usint w;
	vector<usint> arr;
	vector<bool> busy;
};

//��������� ���� - �������� ����� ������������������ ����
struct graph
{
	graph(short int )
};

//���������� ������� ��� ��������� ����
//TODO: ���������� �� ������������� �����������
void paths(int start, edges_status stat, vector<vector<set<vector<int>>>> &matrix, vector<int> way)
{
	way.push_back(start);
	way.push_back(stat.arr[start]);
	if (stat.arr[start] < (stat.c + stat.w) * 2) //���� ��� ���� ������� � ������������� ��������
	{
		paths((stat.arr[start] / 2) * 2, stat, matrix, way);
		paths((stat.arr[start] / 2) * 2 + 1, stat, matrix, way);
	}
	else//���� ��� ���� ������� � ���� ������
	{
		//������� ������������ ���������� � ������
		matrix[way.front() - (stat.c + stat.w) * 2][way.back() - (stat.c + stat.w) * 2].insert(way);
	}
};

//��������� ��� �������� ���� ��������������� ������
struct vertice
{
	vertice *left, *right;
	func f;
	vertice() { left = right = nullptr; }
	~vertice() { delete left; delete right; };
};

//���������, ������� ������ ������������ ���� � ������� �������� �� �������
struct graph_and_Mampl
{
	short int *graph = nullptr;
	func **MAmpl = nullptr;
	//p - ����� ������ �����-������
	//c - ����� ���������������� ���������
	//w - ����� �������� ���������
	graph_and_Mampl(void)
	{
		graph = new short int[graph_size];
		MAmpl = new func*[p];
		for (usint i = 0; i < p; i++) MAmpl[i] = new func[p];
	};
	~graph_and_Mampl()
	{
		for (usint i = 0; i < p; i++) delete MAmpl[i];
		delete MAmpl;
		delete graph;
	}
};

int main(void)
{
	
	unsigned int graph[16] = { 10, 15, 4, 6, 8, 11, 9, 14, 12, 13, 0, 5, 2, 3, 7, 1 };//������������ ����, ����������� �����

	//���� ���, ������� ������� ���������� ������� ���������� ������������� �����
	{
		vector <vector<set<vector<int>>>> matrix;//������� ���������� - ������ � ���� ������ ����� �����, ����� ������� �������� �����

		//���� ����� ���� �������� ������������� � ��������� ������� ����������
		{
			edges_status stat;
			int p = stat.p = 6;
			int c = stat.c = 5;
			int w = stat.w = 0;

			for (int i = 0; i < 16; i++) stat.arr[i] = graph[i];
			stat.busy.resize((c + w) * 2 + p, false);

			matrix.resize(stat.p);
			for (int i = 0; i < stat.p; i++) matrix[i].resize(stat.p);

			for (int i = (stat.c + stat.w) * 2; i < (stat.c + stat.w) * 2 + stat.p; i++)
			{
				vector<int> way;
				paths(i, stat, matrix, way);
			}
		}

		double eta[5] = { 1. / 3, 1. / 2, 1. / 3, 1. / 3, 1. / 2 };
		func couplers[5][2][2];

		for (int num = 0; num < 5; num++)
		{
			couplers[num][0][0] = [&eta, num]()->complex<double> {return sqrt(eta[num]);};
			couplers[num][0][1] = [&eta, num]()->complex<double> {return sqrt(1 - eta[num]);};
			couplers[num][1][0] = [&eta, num]()->complex<double> {return sqrt(1 - eta[num]);};
			couplers[num][1][1] = [&eta, num]()->complex<double> {return -sqrt(eta[num]);};
		}

		vertice *M[6][6];
		{
			for (int i = 0; i < 6; i++)
				for (int j = 0; j < 6; j++)
					M[i][j] = nullptr;
		}

		//�������� ��������
		{
			for (int i = 0; i < 6; i++)
				for (int j = 0; j < 6; j++)
					do {//���� ���� ���������� ��� ���� � ������ ������� matrix[i][j]
						//���� � ��������� ����
						if (!matrix[i][j].empty())
						{
							vertice *path = nullptr;//��������� ���� - ������ ���� ������������ ���� ���������, �/� ������� �������� �����
							vector<int> tmp = *matrix[i][j].begin();//����������� ��������� ����
							for (int k = 1; k < (int)tmp.size() - 1; k += 2)//����� �������������� ��������� � ���� tmp
							{
								pair<int, int> coup = { tmp[k], tmp[k + 1] };//��������� ������������� ��������
								if (path == nullptr)//���� ��� ������ ������������� ��������
								{//,�� ������ ����� ���� ������
									path = new vertice;
									path->f = couplers[coup.first / 2][coup.first % 2][coup.second % 2];//���� ���� ������������ ���� � ����� ����� ��������������� ������
								}
								else
								{//�����:
									vertice *newcenter = new vertice;
									newcenter->left = path;//������� ������ ������ �����
														   //� ������ ������ ����� ������������� ��������
									newcenter->right = new vertice;
									newcenter->right->f = couplers[coup.first / 2][coup.first % 2][coup.second % 2];
									//� ����� ����������� �������
									func fl = newcenter->left->f, fr = newcenter->right->f;
									newcenter->f = [fl, fr]()->complex<double> { return fl() * fr(); };

									path = newcenter;
								}
							}
							matrix[i][j].erase(matrix[i][j].begin());//������� ����, ������� �� ������ ��� ����������

							if (M[i][j] == nullptr)//���� ��� ������ �� ����� � ������ ������
							{
								M[i][j] = path;
							}
							else
							{
								//���� ��� �� ������ �� �����, �� ����� ���� ����� ������� � ���, ������� ��� ������� �  M[i][j]
								vertice *newM = new vertice;

								newM->left = M[i][j];
								newM->right = path;
								func fl = newM->left->f;
								func fr = newM->right->f;
								newM->f = [fl, fr]()->complex<double> { return fl() + fr(); };

								M[i][j] = newM;
							}
						}
						else//����� ������� ������� matrix[i][j] �����/���� ������
						{
							if (M[i][j] == nullptr)
							{
								M[i][j] = new vertice;
								M[i][j]->f = []() {return 0;};
							}
							break;
						}
					} while (true);//����� �������� ���� ����� ����� ������
		}

		/*
		� ���� ����� �� �������� �� ���� ������� ������� ���������� matrix[i][j]
		� �������� ��� ���������� �� �������� ������������� ���������. ��� ��������� ����������
		����� �������� ��������, ������������ ��������� �� ������� M
		*/

		//��� ����� ������� �� ����� ���������� ������� �������� M
		if (true)
		{
			cout << endl;
			for (int i = 0; i < 6; i++)
			{
				cout << setprecision(3) << M[i][0]->f().real();
				for (int j = 1; j < 6; j++)
					cout << '\t' << setprecision(3) << M[i][j]->f().real();
				cout << endl;
			}
		}

		//���������� ���������� �������
		func L[4][4];//���������� �������
		func A[6][6];//������� �������� - ���� ������� M, �� ������� ����� ������ �������. ��� ��� ��������, ����� �� ������ ������ ��� M[i][j]->f
		{
			for (int i = 0; i < 6; i++)
				for (int j = 0; j < 6; j++)
					A[i][j] = M[i][j]->f;

			L[0][0] = [&A]() {return A[1][1]() * A[3][3]() + A[1][3]() * A[3][1]();};
			L[0][1] = [&A]() {return A[1][1]() * A[2][3]() + A[1][3]() * A[2][1]();};
			L[0][2] = [&A]() {return A[0][1]() * A[3][3]() + A[0][3]() * A[3][1]();};
			L[0][3] = [&A]() {return A[0][1]() * A[2][3]() + A[0][3]() * A[2][1]();};

			L[1][0] = [&A]() {return A[1][1]() * A[3][2]() + A[1][2]() * A[3][1]();};
			L[1][1] = [&A]() {return A[1][1]() * A[2][2]() + A[1][3]() * A[2][1]();};
			L[1][2] = [&A]() {return A[0][1]() * A[3][2]() + A[0][2]() * A[3][1]();};
			L[1][3] = [&A]() {return A[0][1]() * A[2][2]() + A[0][2]() * A[2][1]();};

			L[2][0] = [&A]() {return A[1][0]() * A[3][3]() + A[1][3]() * A[3][0]();};
			L[2][1] = [&A]() {return A[1][0]() * A[2][3]() + A[1][3]() * A[2][0]();};
			L[2][2] = [&A]() {return A[0][0]() * A[3][3]() + A[0][3]() * A[3][0]();};
			L[2][3] = [&A]() {return A[0][0]() * A[2][3]() + A[0][3]() * A[2][0]();};

			L[3][0] = [&A]() {return A[1][0]() * A[3][2]() + A[1][2]() * A[3][0]();};
			L[3][1] = [&A]() {return A[1][0]() * A[2][2]() + A[1][2]() * A[2][0]();};
			L[3][2] = [&A]() {return A[0][0]() * A[3][2]() + A[1][2]() * A[3][0]();};
			L[3][3] = [&A]() {return A[0][0]() * A[2][2]() + A[0][2]() * A[2][0]();};
		}

		//����� ���������� ������� L �� �����
		if (true)
		{
			cout << endl;
			for (int i = 0; i < 4; i++)
			{
				cout.setf(ios::fixed);
				cout << setprecision(3) << L[i][0]().real();
				for (int j = 1; j < 4; j++)
					cout << '\t' << setprecision(3) << L[i][j]().real();
				cout << endl;
			}
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
							[&A, k1, k2]()->complex<double> {
							complex<double> val;
							for (int i = 0; i < 6; i++)
								val += A[i][k1]() * A[i][k2]();
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

		//��������, ������������� �� ����������� �������� �������������� ���� ���� ����� ��������
		if (true)
		{
			bool satisfy = true;//���� ����, ��� ����������� ��� ������� ���������
			double eps = 1e-2;
			for (size_t i = 0; i < restrictions.size(); i++)
				if (abs(restrictions[i]()) > eps) {
					satisfy = false; break;
				}
			if (satisfy) cout << "Ideal parameters satisfying all restrictions" << endl;
			else cout << "Ideal parameters NOT satisfying all restrictions" << endl;

			//������� �������� ���� ������� ��������� ����
			if (!satisfy)
			{
				for (size_t i = 0; i < restrictions.size(); i++)
					cout << i << ": " << restrictions[i]() << endl;
			}
		}

		//���������� ������ �����-����� - ������ ��������� ������� ����� � ��������� ��� ��� �������
		if (true)
		{
			cout << endl;
			double eps = 1e-1;//�������� ��� ������� ���������
			bool nonzero = false;//���� ����, ��� ����� �� ������������� ������ �� �������
			srand(time(0));

			unsigned int total = 0;//����� ����� ��������������� �����
			unsigned int satisfy_restrictions = 0;//����� �����, ��������������� ���� ��������
			clock_t start = clock();//��� ��������� �������, ������������ ��� ���������� ����������� ���������

			vector<unsigned int> statistic_total, statistic_satis;
			vector<clock_t> statistic_clock;

			do {//���� �������� ��������� �����
				//����������� ��������� ����� � ������� ���������� ��������
				for (int i = 0; i < 5; i++)
					eta[i] = (double)rand() / RAND_MAX;
				total++;

				//if (total % 10000 == 0) cout << "Total = " << total << "; satisfy(%) = " << setprecision(2) << (double)satisfy_restrictions/total * 100 << endl;

				for (size_t i = 0; i < restrictions.size(); i++)
					if (abs(restrictions[i]()) > eps)
					{
						nonzero = true;
						break;
					};
				if (nonzero)
				{
					nonzero = false;
					continue;//���� ������ ����� �� ������������� ������ �� ������� ���������, �� ������� ���������
				}

				//� ������ ����� � ��� ���� �����, ��������������� ���� �������� ���������
				satisfy_restrictions++;

				//������ ��������� ����� �����, �������� � ���������� �������� ������� �������
				if (abs(L[3][3]()) > 0.3)
				{
					statistic_total.push_back(total);
					statistic_satis.push_back(satisfy_restrictions);
					statistic_clock.push_back(clock() - start);
					cout << statistic_total.size() <<
						": Tot = " << total <<
						";\tsatis_restr = " << satisfy_restrictions <<
						";\ttime = " << setprecision(3) << (double)(statistic_clock.back() * 1000) / CLOCKS_PER_SEC << "ms" << endl;
					total = satisfy_restrictions = 0;
					start = clock();
				}

			} while (statistic_total.size() < 1000);//����� ���� �������� �����

													//������� ���������� � ���� ��� ���������
			if (true)
			{
				ofstream file("data_total", ios::trunc);
				sort(statistic_total.begin(), statistic_total.end());
				for (size_t i = 0; i < statistic_total.size(); i++)
					file << statistic_total[i] << '\t' << (double)i / statistic_total.size() << endl;
				file.close();

				sort(statistic_satis.begin(), statistic_satis.end());
				file.open("data_satis", ios::trunc);
				for (size_t i = 0; i < statistic_satis.size(); i++)
					file << statistic_satis[i] << '\t' << (double)i / statistic_satis.size() << endl;
				file.close();

				sort(statistic_clock.begin(), statistic_clock.end());
				file.open("data_clock", ios::trunc);
				for (size_t i = 0; i < statistic_clock.size(); i++)
					file << (double)(statistic_clock[i] * 1000) / CLOCKS_PER_SEC << '\t' << (double)i / statistic_clock.size() << endl;
				file.close();
			}
		}

		for (int i = 0; i < 6; i++)
			for (int j = 0; j < 6; j++)
				delete M[i][j];
	}

	return 0;
}