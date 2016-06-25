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

//Хранит в себе комбинацию типов однокубитовых операторов
struct operators_comb
{
	//q - число однокубитовых операторов
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
	char size;//Число однокубитовых операторов

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

//Содержит всю информацию об отдельном графе
struct graph
{
	/*
	p(orts) - число портов ввода-вывода
	cp (coupPlates) - число светоделительных пластинок
	dc (directionCoupler) - число направленных светоделителей
	w(aveplates) - число волновых пластинок (фазовращателей)
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

	//Направленный граф
	char *edges;
	//Хранит в себе инфу о том, какие из рёбер были инициализированы
	bool *busy;
	//Хранит в себе комбинацию типов однокубитовых операторов
	operators_comb *comb;
	
	char size;//Размер направленного графа - чтобы каждый раз не складывать отдельные слагаемые
	char p;//p(orts) - число портов ввода-вывода
	char q;//Число однокубитовых операторов
	char bs;//(beamsplitters) - число светоделительных пластинок
	char dc;//(direction couplers) - число направленных светоделителей
	char w;//(waveplates) - число волновых пластинок (фазовращателей)
	char var_size;//Число переменных

	complex<double> *var;//Внутренние параметры графа
	
	//Возвращает номер переменной из массива var[], соответствующей однокубитовому оператору comb->op[oper_num]
	//Если для данного однокубитового оператора нужно несколько переменных - то будет возвращёна ссылка на первого из них
	complex<double>& var_num(char oper_num)
	{
		//Вычислим номер необходимого оператора из массива var[]
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
	
	//Определяет какому типу однокубитового оператора принадлежит переменная с номером var_num
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

//Рекурсивно находит все возможные пути
void paths(char start, graph g, set<vector<char>> **matrix, vector<char> way)
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

//Возвращает матрицу амплитуд по матрице траекторий. Граф g нужен для определения типа однокубитового элемента и числа тех или иных элементов
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
				vector<char> one_path = *(matrix[i][j].begin());//Отдельный путь
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
					} else
					{
						func old_func = summand;
						summand = [old_func, &g, oper_num, in, out]()->complex<double> {return old_func() * get_func(g, oper_num, in, out)();};
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
	return ans;
};

//Данная рекурсия может получить на входе пустой граф, или его заготовку
//Возвращает наилучший граф из всех возможных для данной заготовки
//me - узел, с которого мы сейчас стартуем
pair<graph,double> choose_best_graph(char me, graph g)
{
	double best = 0;//Вероятность работы наилучшей схемы с оптимальными параметрами
	graph g_best(g);//Граф, вероятность работы которого наивысшая

	//Необходимо достроить полученный направленный граф
	if (me < g.q * 2)//Если мы сейчас стартуем из однокубитового оператора
	{
		//Пройдёмся по всем однокубитовым элементам
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
	else//В этой точке мы стартуем из порта ввода
		if (me < g.size)
			//Порт ввода может смотреть в любой узел графа
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
			//В этой ветке мы уже прошлись по всем исходящим вершинам графа.
			//Теперь в g полностью запоненный массив, представляющий собой направленный граф. Но он ещё абстрактен от светоделителей и фазовращателей

			//Накостыляем условия, которые скажут нам есть ли связь между нужным портом ввода и нужным(и) выводами,
			//и по этим условиям отсеять неподходящие графы
			{
				set<vector<char>> **matrix = new set<vector<char>>*[g.p];//Матрица траекторий
				for (char i = 0; i < g.p; i++) matrix[i] = new set<vector<char>>[g.p];

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

					set<operators_comb> memory_oper_comb;//Все комбинации одонкубитовых операторов
					for (char i = 0; i < g.bs; i++)				g.comb->op[i] = beamsplitter;
					for (char i = g.bs; i < g.bs + g.dc; i++)	g.comb->op[i] = directCoupler;
					for (char i = g.dc; i < g.q; i++)			g.comb->op[i] = waveplate;

					do
					{
						if (memory_oper_comb.find(*g.comb) == memory_oper_comb.end())
						{
							//В данной точке у нас уникальный заполненный граф g, с матрицей траекторий matrix, с уникальной комбинацией однокубитовых операторов
							memory_oper_comb.insert(*g.comb);//Сохраним эту комбинацию, чтобы в будущем не повторяться

							//Теперь получим матрицу амплитуд из функций
							func **MAmpl = make_matrix_amplitude(matrix, g);

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
									for (int k1 = 0; k1 < 5; k1++)//Номер первого столбца пары
										for (int k2 = k1 + 1; k2 < 6; k2++)//Номер второго столбца пары
										{
											restrictions.push_back(
												[&MAmpl, k1, k2]()->complex<double> {
												complex<double> val;
												for (int i = 0; i < 6; i++)
													val += MAmpl[i][k1]() * MAmpl[i][k2]();
												return val;
											});//end vector::push_back()
										}
								}//конец униатрность

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

							//Теперь реализация метода Монте-Карло - необходимо бросить N случайных точек в пространстве g.var[]
							//И если они с точностью eps удовлетворяют условиям равенства
							//Реализация метода Монте-Карло - кидаем случайным образом точки и посмотрим что нам подойдёт
							{
								double eps = 1e-1;//Точность для условий равенства
								bool nonzero = false;//Флаг того, что точка не удовлетворяет одному из условий
								srand(time(0));

								unsigned int total = 0;//Общее число сгенерированных точек
								unsigned int satisfy_restrictions = 0;//Число точек, удовлетворяющих всем условиям

								
								for (char i = 0; i < g.var_size; i++) g.var[i] = 0;

								//Начало цикла перебора случайных точек
								do {
									//Сгенерируем случайную точку в области допустимых значений
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
									double efficiency = abs(L[3][3]());
									if (best < efficiency)
									{
										best = efficiency;
										g_best = g;
									}
								} while (satisfy_restrictions < 5);//конец цикл перебора точек по числу сгенерированных точек

								for (char i = 0; i < g.p; i++) delete MAmpl;
								delete MAmpl;
							}//Конец реализации метода Монте-Карло
						}
					} while (next_permutation(g.comb->op, g.comb->op + g.q));

					//В данной точке у нас есть граф g_best с оптимальными внутренними параметрами, при который он работает с вероятностью best
					//Необходимо выгрузить всю эту инфу во внешний мир
				}//Конец обработки графа, удовлетворяющего матрице траекторий

				for (char i = 0; i < g.p; i++) delete matrix[i];
				delete matrix;
			}
		}//Конец обработки сгенерированного направленного графа

		return pair<graph, double>::pair(g_best, best);
}

int main(void)
{
	

	return 0;
}