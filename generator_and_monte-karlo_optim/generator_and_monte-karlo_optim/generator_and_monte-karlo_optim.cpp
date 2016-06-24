// Optim_monte-karlo.cpp : Defines the entry point for the console application.
//Данное приложение делает полный перебор квантовых схем на основе направленных графов, фильтрует их по матрице амплитуд
//и далее для отфильтрованных схем находит глобальный оптимум внутренних параметров схемы (напр., соотношения светоделителей),
//при которых схема работает с наибольшей вероятностью.

#include "stdafx.h"

using namespace std;

#define func function<complex<double>()>
#define usint unsigned short int

usint cp;//Число светоделительных пластинок
usint wp;//Число волновых пластинок
usint p;//Число портов ввода-вывода
usint graph_size;//Размер графа

//Память для рекурсии edges
struct edges_status
{
	usint p;
	usint c;
	usint w;
	vector<usint> arr;
	vector<bool> busy;
};

//Отдельный граф - содержит флаги инициализированных рёбер
struct graph
{
	graph(short int )
};

//Рекурсивно находит все возможные пути
//TODO: Избавиться от использования контейнеров
void paths(int start, edges_status stat, vector<vector<set<vector<int>>>> &matrix, vector<int> way)
{
	way.push_back(start);
	way.push_back(stat.arr[start]);
	if (stat.arr[start] < (stat.c + stat.w) * 2) //Если наш порт смотрит в однокубитовый оператор
	{
		paths((stat.arr[start] / 2) * 2, stat, matrix, way);
		paths((stat.arr[start] / 2) * 2 + 1, stat, matrix, way);
	}
	else//Если наш порт смотрит в порт вывода
	{
		//Запишем получившуюся траекторию в массив
		matrix[way.front() - (stat.c + stat.w) * 2][way.back() - (stat.c + stat.w) * 2].insert(way);
	}
};

//Структура для хранения узла синтаксического дерева
struct vertice
{
	vertice *left, *right;
	func f;
	vertice() { left = right = nullptr; }
	~vertice() { delete left; delete right; };
};

//Структура, которая хранит направленный граф и матрицу амплитуд из функций
struct graph_and_Mampl
{
	short int *graph = nullptr;
	func **MAmpl = nullptr;
	//p - число портов ввода-вывода
	//c - число светоделительных пластинок
	//w - число волновых пластинок
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
	
	unsigned int graph[16] = { 10, 15, 4, 6, 8, 11, 9, 14, 12, 13, 0, 5, 2, 3, 7, 1 };//Направленный граф, описывающий схему

	//Ниже код, который находит глобальный оптимум отдельного направленного графа
	{
		vector <vector<set<vector<int>>>> matrix;//Матрица траекторий - хранит в себе номера узлов графа, через которые проходит фотон

		//Этот кусок кода проводит инициализацию и формирует матрицу траекторий
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

		//Основной алгоритм
		{
			for (int i = 0; i < 6; i++)
				for (int j = 0; j < 6; j++)
					do {//Этот цикл перебирает все пути в ячейке матрицы matrix[i][j]
						//Вход в отдельный путь
						if (!matrix[i][j].empty())
						{
							vertice *path = nullptr;//Отдельный путь - должно быть перемножение всех элементов, ч/з которые проходит фотон
							vector<int> tmp = *matrix[i][j].begin();//Вытаскиваем отдельный путь
							for (int k = 1; k < (int)tmp.size() - 1; k += 2)//Номер однокубитового оператора в пути tmp
							{
								pair<int, int> coup = { tmp[k], tmp[k + 1] };//Отдельный однокубитовый оператор
								if (path == nullptr)//Если это первый однокубитовый оператор
								{//,то создаём новый узел дерева
									path = new vertice;
									path->f = couplers[coup.first / 2][coup.first % 2][coup.second % 2];//Этот узел впоследствии уйдёт в самый конец синтаксического дерева
								}
								else
								{//Иначе:
									vertice *newcenter = new vertice;
									newcenter->left = path;//убираем старое дерево влево
														   //а справа ставим новый однокубитовый оператор
									newcenter->right = new vertice;
									newcenter->right->f = couplers[coup.first / 2][coup.first % 2][coup.second % 2];
									//И задаём центральную функцию
									func fl = newcenter->left->f, fr = newcenter->right->f;
									newcenter->f = [fl, fr]()->complex<double> { return fl() * fr(); };

									path = newcenter;
								}
							}
							matrix[i][j].erase(matrix[i][j].begin());//Удаляем путь, который мы только что обработали

							if (M[i][j] == nullptr)//Если это первый из путей в данной ячейке
							{
								M[i][j] = path;
							}
							else
							{
								//Если это не первый из путей, то новый путь нужно сложить с тем, который уже сохранён в  M[i][j]
								vertice *newM = new vertice;

								newM->left = M[i][j];
								newM->right = path;
								func fl = newM->left->f;
								func fr = newM->right->f;
								newM->f = [fl, fr]()->complex<double> { return fl() + fr(); };

								M[i][j] = newM;
							}
						}
						else//Когда элемент матрицы matrix[i][j] стала/была пустой
						{
							if (M[i][j] == nullptr)
							{
								M[i][j] = new vertice;
								M[i][j]->f = []() {return 0;};
							}
							break;
						}
					} while (true);//Конец перебора всех путей одной ячейки
		}

		/*
		В этой точке мы прошлись по всем ячейкам матрицы траекторий matrix[i][j]
		и заменили эти траектории на значения однокубитовых элементов. При изменении переменных
		будут меняться значения, возвращаемые функциями из матрицы M
		*/

		//Для теста выведем на экран содержимое матрицы амплитуд M
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

		//Сформируем логическую матрицу
		func L[4][4];//Логическая матрица
		func A[6][6];//Матрица амплитуд - есть матрица M, из которой взяты только функции. Ввёл для удобства, чтобы не писать каждый раз M[i][j]->f
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

		//Вывод логической матрицы L на экран
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
							[&A, k1, k2]()->complex<double> {
							complex<double> val;
							for (int i = 0; i < 6; i++)
								val += A[i][k1]() * A[i][k2]();
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

		//Проверим, удовлетворяют ли оптимальные значения светоделителей всем этим нашим условиям
		if (true)
		{
			bool satisfy = true;//Флаг того, что выполняются все условия равенства
			double eps = 1e-2;
			for (size_t i = 0; i < restrictions.size(); i++)
				if (abs(restrictions[i]()) > eps) {
					satisfy = false; break;
				}
			if (satisfy) cout << "Ideal parameters satisfying all restrictions" << endl;
			else cout << "Ideal parameters NOT satisfying all restrictions" << endl;

			//Выведем значнеия всех условий равенства нулю
			if (!satisfy)
			{
				for (size_t i = 0; i < restrictions.size(); i++)
					cout << i << ": " << restrictions[i]() << endl;
			}
		}

		//Реализация метода Монте-Карло - кидаем случайным образом точки и посмотрим что нам подойдёт
		if (true)
		{
			cout << endl;
			double eps = 1e-1;//Точность для условий равенства
			bool nonzero = false;//Флаг того, что точка не удовлетворяет одному из условий
			srand(time(0));

			unsigned int total = 0;//Общее число сгенерированных точек
			unsigned int satisfy_restrictions = 0;//Число точек, удовлетворяющих всем условиям
			clock_t start = clock();//Для измерения времени, необходимого для нахождения глобального максимума

			vector<unsigned int> statistic_total, statistic_satis;
			vector<clock_t> statistic_clock;

			do {//цикл перебора случайных точек
				//Сгенерируем случайную точку в области допустимых значений
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
					continue;//Если данная точка не удовлетворяет одному из условий равенства, то пробуем следующую
				}

				//В данной точке у нас есть точка, удовлетворяющая всем условиям равенства
				satisfy_restrictions++;

				//Теперь посчитаем число точек, попавших в допустимый интервал целевой функции
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

			} while (statistic_total.size() < 1000);//конец цикл перебора точек

													//Выведем статистику в файл для отрисовки
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