#include "mpi.h"
#include <iostream>
#include <fstream>
#include <windows.h>
#include <ctime>
#include <cmath>
#include <vector>
#include <iomanip>
using namespace std;

template<typename T> 
T* generate_matrix (int n) {
	T* matrix = new T[n * n];
	return matrix;
}

template<typename T> 
void init_matrix (T* matrix, int n) {
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			if (i == j) {
				matrix[i * n + i] = 100.;
			} else {
				matrix[i * n + j] = (float)(2 * i + j) / 100000.f;
			}
		}
	}
}

template<typename T> 
T* multiply_matrix_vector (T* matrix, T* vector, int n) {
	T* result = new T[n];
	for (int i = 0; i < n; ++i) {
		result[i] = 0;
		for (int j = 0; j < n; ++j) {
			result[i] += matrix[i * n + j] * vector[j];
		}
	}
	return result;
}

template<typename T> 
void print_matrix (T* matrix, int n1, int n2) {
	for (int i = 0; i < n1; i++) {
		for (int j = 0; j < n2; ++j) {

			cout << setw(5) << matrix[i * n2 + j] << " ";
		}
		cout << endl;
	}
}

template<typename T> 
void print_vector (T* vector, int n) {
	for (int i = 0; i < n; ++i) {
		cout << vector[i] << " ";
	}
	cout << endl;
}

template<typename T> 
T* hstack_matrix_vector(T* matrix, T* vector, int n) {
	T* extended_matrix = new T[n * (n + 1)];
	for (int i = 0; i < n; ++i) {
		extended_matrix[i * (n + 1) + n] = vector[i];
		for (int j = 0; j < n; ++j) {
			extended_matrix[i * (n + 1) + j] = matrix[i * n + j];
		}
	}
	return extended_matrix;
}

template<typename T> 
T* transpose_matrix(T* matrix, int n1, int n2) {
	T* transposed_matrix = new T[n1 * n2];
	for (int i = 0; i < n1; ++i) {
		for (int j = 0; j < n2; ++j) {
			transposed_matrix[j * n1 + i] = matrix[i * n2 + j];
		}
	}
	return transposed_matrix;
}

template<typename T> 
inline void gauss_forward_step(T* a, int n) {
	T l_k;
	for (int k = 0; k < n - 1; ++k) {
		for (int i = k + 1; i < n; ++i) {
			l_k = a[i * (n + 1) + k] / a[k * (n + 1) + k];
			for (int j = k + 1; j < n + 1; ++j) {
				a[i * (n + 1) + j] -= l_k * a[k * (n + 1) + j] ;
			}
		}
	}
}

template<typename T> 
inline T* gauss_backward_step(T* a, int n) {
	T* x = new T[n];
	x[n - 1] = a[(n - 1) * (n + 1)  + n] / a[(n - 1) * (n + 1) + n - 1];
	for (int i = n - 2; i >= 0; --i) {
		x[i] = a[i * (n + 1) + n];
		for (int j = n - 1; j > i; --j) {
			x[i] -= a[i * (n + 1) + j] * x[j];
		}
		x[i] /= a[i * (n + 1) + i];
	}
	return x;
}

template<typename T> 
inline T* gauss_method_sequential(T* a, int n) {
	gauss_forward_step(a, n);
	T* result = gauss_backward_step(a, n);
	return result;
}

template<typename T> 
T check_solution(T* result, T* exact_solution, int n) {
	T difference = 0;
	for (int i = 0; i < n; ++i) {
		difference += (result[i] - exact_solution[i]) * (result[i] - exact_solution[i]);
	}
	difference = sqrt(difference);
	return difference;
}

vector<int> read_sequential_test(string file_name) throw (exception){
	ifstream test(file_name);
	int size;
	vector<int> tests;
	if (test.is_open()) {
		while (test >> size) {
			tests.push_back(size);
		}
	} else {
		throw exception("File not found.");
	}
	test.close();
	return tests;
}

template <typename T>
void write_sequential_result_to_file (string file_name, int n, T difference, float time) throw (exception){
	ofstream seq_result(file_name, std::ios::app);

	seq_result << "Size of matrix is " << n << endl;
	seq_result << "Difference between solutions by L_2 norm is " << difference << endl;
	seq_result << "Sequential algorith took " << time << " seconds." << endl;
	seq_result << "__________________________________________" << endl;
		
	seq_result.close();
}

template <typename T>
void write_sequential_statistic_to_console (int n, T difference, float time) throw (exception){
	cout << "Size of matrix is " << n << endl;
	cout << "Difference between solutions by L_2 norm is " << difference << endl;
	cout << "Sequential algorith took " << time << " seconds." << endl;
	cout << "__________________________________________" << endl;
}

template <typename T>
void write_parallel_result_to_file (string file_name, int n, int q3, T difference, float time) throw (exception){
	ofstream seq_result(file_name, std::ios::app);

	seq_result << "Size of matrix is " << n << endl;
	seq_result << "Number of processes is " << q3 << endl;
	seq_result << "Difference between solutions by L_2 norm is " << difference << endl;
	seq_result << "Parallel algorith took " << time << " milliseconds." << endl;
	seq_result << "__________________________________________" << endl;

	seq_result.close();
}

int main(int argc, char *argv[]) {

	/*
*******************sequential algorithm************************
	*/
	
	//vector<int> tests;
	//try{
	//	 tests = read_sequential_test("sequence_test.txt");
	//} catch (exception const & e){
	//	cout << e.what() << endl;
	//}

	//for(int n : tests) {
	//	float* a_matrix = generate_matrix<float>(n);
	//	float* exact_solution = new float[n];
	//	fill(exact_solution, exact_solution + n, 1.);

	//	init_matrix(a_matrix, n);

	//	float* b_vector = multiply_matrix_vector(a_matrix, exact_solution, n);

	//	float* a_extended = hstack_matrix_vector(a_matrix, b_vector, n);

	//	//cout<<*(a_extended+0)<< " " << *(a_extended+1) << " " << *(a_extended+3) << endl;

	//	float* a_transposed = transpose_matrix(a_extended, n, n + 1);

	//	clock_t start_time = clock();
	//	float* result = gauss_method_sequential(a_extended, n);
	//	clock_t end_time = clock();

	//	float difference = check_solution(result, exact_solution, n);
	//
	//	write_sequential_statistic_to_console(n, difference, (float)(end_time - start_time) / CLOCKS_PER_SEC);

	//	write_sequential_result_to_file("sequence_result.txt", n, difference, (float)(end_time - start_time) / CLOCKS_PER_SEC);
	//}
	//std::system("pause");
	
//***********Paraller algorithm*******************

	// number of processes
	int q3;
	// process number
	int p;

	if (int rc = MPI_Init(&argc, &argv)) { 
		cout << "Error while starting MPI." << endl;
		MPI_Abort(MPI_COMM_WORLD, rc);
	}

	int n, r1, r2;

	if (argc == 4) {
		n = atoi(argv[1]);
		r1 = atoi(argv[2]);
		r2 = atoi(argv[3]);
	} else {
		n = 6;
		r1 = 2;
		r2 = 2;
	}

	MPI_Comm_size(MPI_COMM_WORLD, &q3);
	MPI_Comm_rank(MPI_COMM_WORLD, &p);

	int q1 = ceil((double) n / r1);
	int q2 = ceil((double) n / r2);
	int r3 = ceil((double) (n + 1) / q3);

	float* exact_solution;

	float a_extended [6 * 7] = {100, 1, 2, 3, 4, 5, 115, 6, 100, 7, 8, 9,1, 131, 2,3,100,4,5,6,120,7,8,9,100,1,2,127, 3,4,5,6,100,7, 125, 8,9,1,2,3,100,123};
		//= new float[];
	double start_wtime;
	// initializing matrix
	if(p == 0) {
		exact_solution = new float[n];
		fill(exact_solution, exact_solution + n, 1.);
		
		float* a_matrix = generate_matrix<float>(n);
		
		init_matrix(a_matrix, n);

		float* b_vector = multiply_matrix_vector(a_matrix, exact_solution, n);

		//a_extended = hstack_matrix_vector(a_matrix, b_vector, n);

		//delete[] b_vector;
		//delete[] a_matrix;
		start_wtime = MPI_Wtime();
	}

	// sending part a_p
	// matrix n*r3
	float* a_p = new float[n * r3];
	if(p == 0){
		for(int i = q3 - 1; i >= 0; --i) {
			for(int j = 0; j < n; ++j) {
				for(int k = i * r3; k < min((i + 1) * r3, n + 1); ++k) {
					int k_p = k - i * r3;
					a_p[j * r3 + k_p] = a_extended[j * (n + 1) + k];
				}
			}
			if( i != 0 ) {
				MPI_Send(a_p, n * r3, MPI_FLOAT, i, 0, MPI_COMM_WORLD);
			}
		}
	}
	if (p != 0) {
		MPI_Recv(a_p, n * r3, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}

	if(p==0){
		print_matrix(a_p, n, r3);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	if(p == 1) {
		print_matrix(a_p, n, r3);
	}
	MPI_Barrier(MPI_COMM_WORLD);

	for(int k_gl = 0; k_gl < q1; ++k_gl) {
		cout << "k_gl = " << k_gl << endl;
		for(int i_gl = 0; i_gl < q2; ++i_gl) {
			cout << "i_gl = " << i_gl << endl;
			float* l_p = new float[r2 * r1]; // может выделить под размер треугольной матрицы
			//float l_p;
			for(int k = k_gl * r1; k <= min((k_gl + 1) * r1 - 1, n - 2); ++k) { 
				int k_p = k - k_gl * r1;
				for(int i = max(i_gl * r2, k + 1); i <= min((i_gl + 1) * r2 - 1, n - 1); ++i) {
					int i_p = i - i_gl * r2;
					//int pr_num = floor((double) k / r3);
					if( p == 0 ){
						///????????????????????????????????????????????????
						l_p[i_p * r1 + k_p] = a_p[i * (n + 1) + k%r3] / a_p[k * (n + 1) + k%r3];
						//l_p[i_p * r1 + k_p] = a_extended[i * (n + 1) + k] / a_extended[k * (n + 1) + k];
						//l_p = a_p[i * (n + 1) + k_p] / a_p[k * (n + 1) + k_p];
						cout << "l(" << i << ", " << k << ") = " <<  l_p[i_p * r1 + k_p] << endl;
						//cout<<p<< endl;
						if(q3 > 1){
							MPI_Send(l_p, r2 * r1, MPI_FLOAT, 1, 0, MPI_COMM_WORLD);
							//MPI_Send(&l_p, 1, MPI_FLOAT, 1, 0, MPI_COMM_WORLD);
						}
					}
					if(p > 0){
						MPI_Recv(l_p, r2 * r1, MPI_FLOAT, p - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
						//MPI_Recv(&l_p, 1, MPI_FLOAT, p - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
						if( p != q3 - 1 ) {
							MPI_Send(l_p, r2 * r1, MPI_FLOAT, p + 1, 0, MPI_COMM_WORLD);
							//MPI_Send(&l_p, 1, MPI_FLOAT, p + 1, 0, MPI_COMM_WORLD);
						}
					}
					for(int j = max(p * r3, k + 1); j <= min((p + 1) * r3 - 1, n); ++j) {
						int j_p = j - p * r3;
						a_p[i * r3 + j_p] -= l_p[i_p * r1 + k_p] * a_p[k * r3 + j_p];
						//a_p[i * r3 + j_p] -= l_p * a_p[k * r3 + j_p];
					}
				}
				if(p==0){
		print_matrix(a_p, n, r3);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	if(p == 1) {
		print_matrix(a_p, n, r3);
	}
	MPI_Barrier(MPI_COMM_WORLD);
			}
		}
		cout<<"asdasd" <<endl;
	}
	//MPI_Barrier(MPI_COMM_WORLD);
	if(p != 0) {
		for(int i = 1; i < q3; ++i) {
			MPI_Send(a_p, n * r3, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);
		}
	}

	if(p == 0) {
		cout << "A" << endl;
		for(int j = 0; j < n; ++j) {
			for(int k = 0; k < r3; ++k) {
				a_extended[j * (n + 1) + k] = a_p[j * r3 + k];
			}
		}
		for(int i = 1; i < q3; ++i) {
			MPI_Recv(a_p, n * r3, MPI_FLOAT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			for(int j = 0; j < n; ++j) {
				for(int k = i * r3; k < min((i + 1) * r3, n + 1); ++k) {
					int k_p = k - i * r3;
					a_extended[j * (n + 1) + k] = a_p[j * r3 + k_p];
				}
			}
		}
		float* solution = gauss_backward_step(a_extended, n);
		double time_spent = MPI_Wtime() - start_wtime;
		float check = check_solution(solution, exact_solution, n);
		print_vector(solution, n);
		write_parallel_result_to_file("parallel_result.txt", n, q3, check, time_spent);
		std::system("pause");
	}
	delete[] exact_solution;
	delete[] a_extended;
	delete[] a_p;
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
}

/*if(p == 0) {
				for(int k = k_gl * r1; k <= min((k_gl + 1) * r1 - 1, n - 2); ++k) {
					int k_p = k - k_gl * r1;
					for(int i = max(i_gl * r2, k + 1); i <= min((i_gl + 1) * r2 - 1, n - 1); ++i) {
						int i_p = i - i_gl * r2;
						l_p[i_p * r1 + k_p] = a_extended[i * (n + 1) + k] / a_extended[k * (n + 1) + k];
					}
				}
				if(q3 > 1){
					MPI_Send(l_p, r2 * r1, MPI_FLOAT, 1, 0, MPI_COMM_WORLD);
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
			if(p != 0){
				MPI_Recv(l_p, r2 * r1, MPI_FLOAT, p - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				if( p != q3 - 1 ) {
					MPI_Send(l_p, r2 * r1, MPI_FLOAT, p + 1, 0, MPI_COMM_WORLD);
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);

			for(int k = k_gl * r1; k <= min((k_gl + 1) * r1 - 1, n - 2); ++k) {
				int k_p = k - k_gl * r1;
				for(int i = max(i_gl * r2, k + 1); i <= min((i_gl + 1) * r2 - 1, n - 1); ++i) {
					int i_p = i - i_gl * r2;
					for(int j = max(p * r3, k + 1); j <= min((p + 1) * r3 - 1, n); ++j) {
						int j_p = j - p * r3;
						a_p[i * r3 + j_p] -= l_p[i_p * r1 + k_p] * a_p[k * r3 + j_p];
					}
				}
			}*/