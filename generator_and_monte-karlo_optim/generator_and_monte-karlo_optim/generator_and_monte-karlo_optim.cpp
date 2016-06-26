// Optim_monte-karlo.cpp : Defines the entry point for the console application.
//Данное приложение делает полный перебор квантовых схем на основе направленных графов, фильтрует их по матрице амплитуд
//и далее для отфильтрованных схем находит глобальный оптимум внутренних параметров схемы (напр., соотношения светоделителей),
//при которых схема работает с наибольшей вероятностью.

#include "stdafx.h"

using namespace std;

#define func function<complex<double>()>
#define M_PI 3.14159265

//Перечисление типов однокубитовых операторов
enum operators_types
{
	beamsplitter,
	directCoupler,
	waveplate
};

//Содержит всю информацию об отдельном графе
struct graph
{
	/*
	p(orts) - число портов ввода-вывода
	cp (coupPlates) - число светоделительных пластинок
	dc (directionCoupler) - число направленных светоделителей
	w(aveplates) - число волновых пластинок (фазовращателей)
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
		//Инициализация comb
		{
			for (char i = 0; i < bs; i++)			comb[i] = beamsplitter;
			for (char i = bs; i < bs + dc; i++)		comb[i] = directCoupler;
			for (char i = bs + dc; i < q; i++)			comb[i] = waveplate;
		}
		var.resize(bs + dc + 2 * w, 0);
	};

	//Направленный граф
	vector<char> edges;
	//Хранит в себе инфу о том, в какие из вершин уже смотрит ребро
	vector<bool> busy;
	//Хранит в себе комбинацию типов однокубитовых операторов
	vector<operators_types> comb;
	
	char size;//Размер направленного графа
	char p;//p(orts) - число портов ввода-вывода
	char q;//Число однокубитовых операторов
	char bs;//(beamsplitters) - число светоделительных пластинок
	char dc;//(direction couplers) - число направленных светоделителей
	char w;//(waveplates) - число волновых пластинок (фазовращателей)
	double efficiency;//Эффективность работы графа

	vector<double> var;//Внутренние параметры графа
	
	//Возвращает указатель на переменную из массива var[], соответствующей однокубитовому оператору comb->op[oper_num]
	//Если для данного однокубитового оператора нужно несколько переменных - то будет возвращёна ссылка на первую из них
	double* var_num(char oper_num)
	{
		//Вычислим номер необходимого оператора из массива var[]
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
	
	//Определяет какому типу однокубитового оператора принадлежит переменная с номером var_num
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

//Рекурсивно находит все возможные пути
void paths(char start, graph &g, vector<vector<set<vector<char>>>> &matrix, vector<char> way)
{
	way.push_back(start);
	way.push_back(g.edges[start]);
	if (g.edges[start] < g.q * 2) //Если наш порт смотрит в однокубитовый оператор
	{
		paths((g.edges[start] / 2) * 2, g, matrix, way);
		paths((g.edges[start] / 2) * 2 + 1, g, matrix, way);
	}
	else//Если наш порт смотрит в порт вывода
	{
		//Запишем получившуюся траекторию в массив
		matrix[way.front() - g.q * 2][way.back() - g.q * 2].insert(way);
	}
};

//Возвращает функцию, определяя тип оператора по номеру oper_num из графа g, а номера портов ввода-вывода по in и out
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

//Возвращает матрицу амплитуд по матрице траекторий. Граф g нужен для определения типа однокубитового элемента и числа тех или иных элементов
vector<vector<func>> make_matrix_amplitude(vector<vector<set<vector<char>>>> &matrix, graph &g)
{
	vector<vector<func>> ans(g.p, vector<func>(g.p));;//Функции будем формировать здесь
	
	for (char i = 0; i < g.p; i++)
		for (char j = 0; j < g.p; j++)
		{
			if (!matrix[i][j].empty())
			{
				while (!matrix[i][j].empty())
				{
					vector<char> one_path = *(matrix[i][j].begin());//Отдельный путь
					matrix[i][j].erase(matrix[i][j].begin());
					func summand = nullptr;//Отдельное слагаемое
					for (char k = 1; k < (char)one_path.size() - 1; k += 2)
					{
						char oper_num = one_path[k] / 2;//Порядковый номер однокубитового оператора из g.comb
						char in = one_path[k] % 2;//Номер порта ввода
						char out = one_path[k + 1] % 2;//Номер порта вывода
						if (summand == nullptr)
						{
							//В этой точке мы знаем номер однокубитового порта и его порты ввода-вывода
							summand = get_func(g, oper_num, in, out);
						}
						else
						{
							func old_func = summand;
							func new_func = get_func(g, oper_num, in, out);
							summand = [old_func, new_func]()->complex<double> {return old_func() * new_func();};
						}
					}//end for(k)

					//В данной точке мы сформировали функцию, которая представляет собой отдельное слагаемое
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

//Реализация метода Монте-Карло
//Возвращает максимальное значение целевой функции aim
//g - направленный граф
//L - логическая матрица из функций
//restrictions - массив условий равенства
//eps - точность удовлетворения условиям равенства
//satisfy - число точек, которые должны удовлетворять всем условиям равенства - условие выхода из цикла
double monte_karlo(graph &g, func aiming_function, vector<func> &restrictions, double eps, int satisfy)
{
	//Необходимо бросить N случайных точек в пространстве g.var[]
	//И если они с точностью eps удовлетворяют условиям равенства
	//Реализация метода Монте-Карло - кидаем случайным образом точки и посмотрим что нам подойдёт
	
	graph g_best(g);//Наилучший граф - в данной функции нам от него нужны только переменные и конечная эффективность
	bool nonzero = false;//Флаг того, что точка не удовлетворяет одному из условий
	srand((unsigned)time(0));

	unsigned int satisfy_restrictions = 0;//Число точек, удовлетворяющих всем условиям

	//Начало цикла перебора случайных точек
	do {
		//Сгенерируем случайную точку в области допустимых значений
		for (char i = 0; i < g.var.size(); i++)
			switch (g.oper_type(i))
			{
			case beamsplitter:
			case directCoupler:
				g.var[i] = (double)rand() / RAND_MAX; break;
			case waveplate:
				g.var[i] = (double)rand() / RAND_MAX * 2 * M_PI; break;
			}

		//Проверка на условия равенства
		for (size_t i = 0; i < restrictions.size(); i++)
			if (abs(restrictions[i]()) > eps)
			{
				nonzero = true;
				break;
			};

		if (nonzero)//Если хоть одно условие равенства не выполнилось
		{
			nonzero = false;
			continue;//Если данная точка не удовлетворяет одному из условий равенства, то пробуем следующую
		}

		//В данной точке у нас есть точка, удовлетворяющая всем условиям равенства
		satisfy_restrictions++;

		//Теперь вычислим то, с какой эффективностью работает схема с получившейся комбинацией внутренних параметров
		double efficiency = abs(aiming_function());
		if (g_best.efficiency < efficiency)
		{
			g_best = g;
			g_best.efficiency = efficiency;
		}
	} while (satisfy_restrictions < satisfy);//конец цикл перебора точек по числу сгенерированных точек
	//Конец реализации метода Монте-Карло

	g = g_best;
	return g_best.efficiency;
}

//Данная рекурсия может получить на входе пустой граф, или его заготовку
//Возвращает наилучший граф из всех возможных для данной заготовки
//me - узел, с которого мы сейчас стартуем
void choose_best_graph(char me, graph g, graph &g_best)
{
	//Необходимо достроить направленный граф g
	if (me < g.q * 2)//Если мы сейчас стартуем из однокубитового оператора
	{
		for (char i = (me / 2) * 2 + 2; i < g.edges.size(); i++)//Узел, куда смотрит порт выхода однокубитового оператора
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
	else//В этой точке мы стартуем из порта ввода
	{
		if (me < g.edges.size())
			//Порт ввода может смотреть в любой узел графа
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
			//В этой ветке мы уже прошлись по всем исходящим вершинам графа.
			//Теперь в g полностью запоненный массив, представляющий собой направленный граф. Но он ещё абстрактен от светоделителей и фазовращателей

			//Накостыляем условия, которые скажут нам есть ли связь между нужным портом ввода и нужным(и) выводами,
			//и по этим условиям отсеять неподходящие графы
			{
				vector<vector<set<vector<char>>>> matrix(g.p, vector<set<vector<char>>>(g.p));//Матрица траекторий

				//Заполним матрицу траекторий
				//Стартанём алгоритм поиска траекторий из всех портов ввода
				//Результат будет в matrix
				for (int i = g.q * 2; i < g.size; i++) paths(i, g, matrix, vector<char>::vector());

				//Теперь у нас есть матрица траекторий для сгенерированного графа. Теперь определим удовлетворяет ли этот граф условию.
				if (matrix[0][2].empty() && matrix[0][3].empty() && //Первая версия
					!matrix[1][2].empty() && !matrix[1][3].empty() &&
					!matrix[0][0].empty() && matrix[0][1].empty() && //Начало второй версии
					matrix[1][0].empty() && !matrix[1][1].empty() &&
					matrix[2][0].empty() && !matrix[2][1].empty() && //Начало четвёртой версии
					!matrix[2][2].empty() && !matrix[2][3].empty() &&
					matrix[3][0].empty() && !matrix[3][1].empty() &&
					!matrix[3][2].empty() && !matrix[3][3].empty())
				{
					//В этой точке у нас есть граф с матрицей траекторий, который удовлетворяет всем условиям
					//Теперь переберём все комбинации типов однокубитовых операторов

					set<vector<operators_types>> memory_oper_comb;//Все комбинации однокубитовых операторов

					do
					{
						if (memory_oper_comb.find(g.comb) == memory_oper_comb.end())
						{
							//В данной точке у нас уникальный заполненный граф g, с матрицей траекторий matrix, с уникальной комбинацией однокубитовых операторов
							memory_oper_comb.insert(g.comb);//Сохраним эту комбинацию, чтобы в будущем не повторяться

							//Теперь получим матрицу амплитуд из функций
							vector<vector<func>> MAmpl = make_matrix_amplitude(matrix, g);

							//Теперь из этой матрицы амплитуд необходимо получить логическую матрицу, которая размером 4x4
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

							//Теперь необходимо составить список условий:
							vector<func> restrictions;//Ограничения равенства " = 0 " - все они должны быть равны нулю
							{
								//унитарность
								{
									//Перебор всех пар столбцов.
									for (int k1 = 0; k1 < g.p - 1; k1++)//Номер первого столбца пары
										for (int k2 = k1 + 1; k2 < g.p; k2++)//Номер второго столбца пары
										{
											restrictions.push_back(
												[&MAmpl, k1, k2, g]()->complex<double> {
												complex<double> val;
												for (char i = 0; i < g.p; i++)
													val += MAmpl[i][k1]() * MAmpl[i][k2]();
												return val;
											});//end vector::push_back()
										}
								}//конец унитарность

								 //нулевые элементы в логической матрице
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
								}//конец нулевые элементы логической матрицы

								 //Равенство действующих элементов логической матрицы
								{
									restrictions.push_back([&L]()->complex<double> {return L[0][1]() - L[1][0]();});
									restrictions.push_back([&L]()->complex<double> {return L[1][0]() - L[2][2]();});
									restrictions.push_back([&L]()->complex<double> {return L[2][2]() - L[3][3]();});
								}//конец равенство действующих элементов между собой
							}

							//В данной точке у нас есть функции для всех условий равенства и полностью собранный граф - осталось только найти глобальный оптимум
							if (g_best.efficiency < monte_karlo(g, L[3][3], restrictions, 1e-1, 5))
							{
								g_best = g;
							}
						}
					} while (next_permutation(g.comb.begin(), g.comb.end()));

					//В данной точке у нас есть граф g_best с оптимальными внутренними параметрами, при который он работает с вероятностью best
				}//Конец обработки графа, удовлетворяющего условиям матрицы траекторий
			}
		}
	}//Конец обработки сгенерированного направленного графа
}

//Функция возвращает заготовки графов, перебирая все комбинации того, куда смотрят первые count-портов ввода
void make_templates_graphs(char count, graph g, vector<graph> &answer, char deep = 0)
{
	if (deep < count)
	{
		for (char to = 0; to < g.edges.size(); to++)//Номер узла, куда смотрит порт ввода
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
	//Кусок кода, созданный для отладки написанного выше кода
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

	//Основной боевой алгоритм
	{
		for (char p = 6; p <= 6; p++)
			for (char bs = 5; bs <= 5; bs++)
				for (char dc = 0; dc <= 0; dc++)
					for (char w = 0; w <= 0; w++)
					{
						//Создадим массив заготовок - переберём все варианты куда могут смотреть первые три порта ввода
						//это 16!/13! = 3360 для p = 6 и q = 5
						{
							graph g(p, bs, dc, w);
							vector<graph> templates_graphs;
							make_templates_graphs(3, g, templates_graphs);
							
							//В данной точке у нас есть куча заготовок графов - необходимо раскидать их по отдельным потокам и как-то собрать с каждого из них инфу
							int n_threads;//Число отдельных потоков
							vector<graph> best_graphs(n_threads, g);//Коллекция наилучших графов со всех потоков
#pragma omp parallel for
							for (int i = 0; i < templates_graphs.size(); i++)//Итератор на заготовку графа
							{
								int my_num_thread;//Номер текущего потока
								graph loc_template = templates_graphs[i];
								graph loc_best_graph = best_graphs[my_num_thread];
								choose_best_graph(0, loc_template, loc_best_graph);
								best_graphs[my_num_thread] = loc_best_graph;
#pragma omp master
								if (i % n_threads == 0 && i != 0)//Каждый раз, когда мы раскидаем на все потоки новые заготовки
								{
									//Мы записываем текущий прогресс в текстовый файл
									ofstream f("current_progress", ios::trunc);
									f << "Persentage: " << (double)i / templates_graphs.size() * 100 << "%" << endl;

									//Сначала отсортируем графы по убыванию их эффективности
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
	}//Конец реализации основного боевого алгоритма

	return 0;
}