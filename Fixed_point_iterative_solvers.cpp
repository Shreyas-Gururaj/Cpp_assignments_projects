#include <iostream>
#include <vector>
#include <math.h>
# define M_PI           3.14159265358979323846  /* pi */
#include <chrono>
#include <algorithm>
#include <fstream>

//Struct contains all variables required for storing FD_Linear system matrix in CRS_Format.
struct FD_Linear_System{

    /// Stores number of inner points.
    int total_inner_points;
    std::vector<double> values;
    std::vector<int> col_indices;
    std::vector<int> row_indices;
};

//Function to compute the Euclidean norm of the vector.
double vector_euclidean_norm(const std::vector<double> &x, const std::vector<double> &y)
{
    double result_norm_prod = 0.0;

    for(unsigned int elem_index = 0; elem_index < x.size(); elem_index++)
        {
            result_norm_prod += x[elem_index] * y[elem_index];
        }
    result_norm_prod = sqrt(result_norm_prod);
    return result_norm_prod;
}

///// Function to get FD_Linear matrix_A in CRS format.
void get_crs_matrix(FD_Linear_System &obj_FD)
{
    const int nr_non_zero_entries(5*obj_FD.total_inner_points*obj_FD.total_inner_points - 4*obj_FD.total_inner_points);
    const int dim = obj_FD.total_inner_points * obj_FD.total_inner_points;
    
    ///resize vectors
    obj_FD.values.resize(nr_non_zero_entries);
    obj_FD.col_indices.resize(nr_non_zero_entries);
    obj_FD.row_indices.resize(dim+1);    


    if(obj_FD.total_inner_points == 1)
    {
        obj_FD.values[0] = 4.0;
        obj_FD.col_indices[0] = 0;
        obj_FD.row_indices[1] = 1;
    }
    else
    {
        int counter = 0;
        obj_FD.row_indices[0] = 0;

        /// for-loop over rows of Matrix
        for(int current_row = 0; current_row < dim; current_row++)
        {   
            if(current_row >=  obj_FD.total_inner_points)
            {
                obj_FD.values[counter] = -1.0;
                obj_FD.col_indices[counter] = current_row - obj_FD.total_inner_points;
                counter++;
            }

            if(current_row % obj_FD.total_inner_points != 0)
            {
                obj_FD.values[counter] = -1.0;
                obj_FD.col_indices[counter] = current_row -1;
                counter++;
            }

            obj_FD.values[counter] = 4.0;
            obj_FD.col_indices[counter] = current_row;
            counter++;

            if(current_row % obj_FD.total_inner_points != obj_FD.total_inner_points-1)
            {
                obj_FD.values[counter] = -1.0;
                obj_FD.col_indices[counter] = current_row +1;
                counter++;
            }

            if(current_row < dim-obj_FD.total_inner_points)
            {
                obj_FD.values[counter] = -1.0;
                obj_FD.col_indices[counter] = current_row + obj_FD.total_inner_points;
                counter++;
            }

            obj_FD.row_indices[current_row+1] = counter;
        }
    }
}

// RHS vector
void get_b(std::vector<double> &b, const int dimension)
{
    const double m_h = 1/static_cast<double>(dimension);
    b.resize(dimension*dimension, 0.0);
    double x_1(0);
    double x_2(0);
    double coefficient = 1.25 * m_h*m_h * M_PI*M_PI;
    double x = 0.0;

    for(int i=0;i<dimension;i++)
    {
        for(int j=0; j<dimension;j++)
        {
            x_1 = m_h * static_cast<double>(j+1);                        
            x_2 = m_h * static_cast<double>(i+1);
            b[j+i*dimension] = coefficient * sin(M_PI*x_1) * cos(0.5*M_PI*x_2); 
        }
    }

    //Add Boundary Condition
    for(int i=0;i<dimension;i++)
    {
        x = m_h * static_cast<double>(i+1);
        b[i] += sin(M_PI*x);
    }
}

// True solution
void true_solution(std::vector<double> &u, const int dimension)
{
    const double m_h = 1/static_cast<double>(dimension);
    u.resize(dimension*dimension, 0.0);
    double x_1(0);
    double x_2(0);

    for(int i=0;i<dimension;i++)
    {
        for(int j=0; j<dimension;j++)
        {
            x_1 = m_h * static_cast<double>(j+1);                        
            x_2 = m_h * static_cast<double>(i+1);
            u[j+i*dimension] = sin(M_PI*x_1) * cos(0.5*M_PI*x_2); 
        }
    }
}

// Function to compute residue vector
void residue(const FD_Linear_System &obj_CRS_Matrix, std::vector<double> &x_temp, const std::vector <double> &b_RHS_vector,
             std::vector <double> &residue_vector)
{
    residue_vector = b_RHS_vector;
    for(unsigned int i = 0; i < obj_CRS_Matrix.row_indices.size() - 1; i++)
    {
        for(int k = obj_CRS_Matrix.row_indices[i]; k < obj_CRS_Matrix.row_indices[i + 1]; k++)
            {
                {
                    residue_vector[i] += -(obj_CRS_Matrix.values[k] * x_temp[obj_CRS_Matrix.col_indices[k]]);
                }
            }
    }
}

// Function to compute one step of jacobi iteration
void jacobi_step(const FD_Linear_System &obj_CRS_Matrix, std::vector<double> &x_jacobi,
                 const std::vector <double> &b_RHS_vector)
{
    std::vector<double> x_jacobi_temp = x_jacobi;
    const double omega  = 1.0;
    for(unsigned int i = 0; i < obj_CRS_Matrix.row_indices.size() - 1; i++)
    {
        x_jacobi[i] = 0.25 * b_RHS_vector[i]  * omega;
        for(int k = obj_CRS_Matrix.row_indices[i]; k < obj_CRS_Matrix.row_indices[i + 1]; k++)
            {
                if(obj_CRS_Matrix.col_indices[k] !=i)
                { 
                    x_jacobi[i] += -(obj_CRS_Matrix.values[k] * x_jacobi_temp[obj_CRS_Matrix.col_indices[k]]) * 0.25 * omega;
                }
            }
        x_jacobi[i] += (1 - omega) * x_jacobi_temp[i];
    }
}

// Function to compute one step of Gauss-Seidel iteration(Forward)
void Gauss_Seidel_step(const FD_Linear_System &obj_CRS_Matrix, std::vector<double> &x_GS,
                 const std::vector <double> &b_RHS_vector)
{
    std::vector<double> x_GS_temp = x_GS;
    const double omega  = 1.0;
    for(unsigned int i = 0; i < obj_CRS_Matrix.row_indices.size() - 1; i++)
    {
        x_GS[i] = 0.25 * b_RHS_vector[i] * omega;
        for(int k = obj_CRS_Matrix.row_indices[i]; k < obj_CRS_Matrix.row_indices[i + 1]; k++)
            {
                if(obj_CRS_Matrix.col_indices[k] < i)
                {
                    x_GS[i] += -(obj_CRS_Matrix.values[k] * x_GS[obj_CRS_Matrix.col_indices[k]]) * 0.25 * omega;
                }

                else if(obj_CRS_Matrix.col_indices[k] > i)
                {
                    x_GS[i] += -(obj_CRS_Matrix.values[k] * x_GS_temp[obj_CRS_Matrix.col_indices[k]]) * 0.25 * omega;
                }

                else
                {
                    continue;
                }
            }
        x_GS[i] += (1 - omega) * x_GS_temp[i];
    }
}

// Function to compute one step of Gauss-Seidel iteration(Backward)
void Gauss_Seidel_step_backward(const FD_Linear_System &obj_CRS_Matrix, std::vector<double> &x_GS_backward,
                 const std::vector <double> &b_RHS_vector)
{
    std::vector<double> x_GS_backward_temp = x_GS_backward;
    const double omega  = 1.0;
    for(unsigned int i = obj_CRS_Matrix.row_indices.size() - 1; i--;)
    {
        x_GS_backward[i] = 0.25 * b_RHS_vector[i] * omega;
        for(int k = obj_CRS_Matrix.row_indices[i]; k < obj_CRS_Matrix.row_indices[i + 1]; k++)
            {
                if(obj_CRS_Matrix.col_indices[k] < i)
                {
                    x_GS_backward[i] += -(obj_CRS_Matrix.values[k] * x_GS_backward_temp[obj_CRS_Matrix.col_indices[k]]) * 0.25 * omega;
                }

                else if(obj_CRS_Matrix.col_indices[k] > i)
                {
                    x_GS_backward[i] += -(obj_CRS_Matrix.values[k] * x_GS_backward[obj_CRS_Matrix.col_indices[k]]) * 0.25 * omega;
                }

                else
                {
                    continue;
                }
            }
        x_GS_backward[i] += (1 - omega) * x_GS_backward_temp[i];
    }
}

// Function to compute one step of symmetric Gauss-Seidel iteration
void Gauss_Seidel_step_symmetric(const FD_Linear_System &obj_CRS_Matrix, std::vector<double> &x_GS_symmetric,
                                 const std::vector <double> &b_RHS_vector)
{
    Gauss_Seidel_step(obj_CRS_Matrix, x_GS_symmetric, b_RHS_vector);
    Gauss_Seidel_step_backward(obj_CRS_Matrix, x_GS_symmetric, b_RHS_vector);
}

void result(const int inner_points, unsigned int maximun_iterations)
{
    FD_Linear_System obj_FD_LS_A;
    obj_FD_LS_A.total_inner_points = inner_points;    // Very important, else will go into an infinte loop
    get_crs_matrix(obj_FD_LS_A);

    std::vector <double> b_RHS_vector;
    get_b(b_RHS_vector, inner_points);

    double relative_accuracy = 1e-4;
    unsigned int k = 0;     //For the while loop(jacobi).
    unsigned int j = 0;     //For the while loop(Gauss-Seidel).
    unsigned int l = 0;     //For the while loop(Symmetric Gauss-Seidel).

    std::vector<double> x_jacobi(inner_points*inner_points, 0.0);
    std::vector<double> x_GS(inner_points*inner_points, 0.0);
    std::vector<double> x_GS_symmetric(inner_points*inner_points, 0.0);
    std::vector<double> residue_vector_jacobi;
    std::vector<double> residue_vector_GS;
    std::vector<double> residue_vector_GS_symmetric;
    residue(obj_FD_LS_A, x_jacobi, b_RHS_vector, residue_vector_jacobi);
    residue(obj_FD_LS_A, x_GS, b_RHS_vector, residue_vector_GS);
    double initial_residue = vector_euclidean_norm(residue_vector_jacobi, residue_vector_jacobi);
    std::cout << " the initial residue is ::   " << initial_residue << std::endl;

//Checking Jacobi method 
auto start_time = std::chrono::high_resolution_clock::now();
    do
    {
        jacobi_step(obj_FD_LS_A, x_jacobi, b_RHS_vector);
        residue(obj_FD_LS_A, x_jacobi, b_RHS_vector, residue_vector_jacobi);
        k = k + 1;

    } while ((vector_euclidean_norm(residue_vector_jacobi, residue_vector_jacobi) / initial_residue > relative_accuracy) && (k <= maximun_iterations));

std::cout << " the number of iterations for jacobi is ::   " << k << std::endl;
auto end_time = std::chrono::high_resolution_clock::now();
auto run_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
std::cout << "Time in milliseconds for Jacobi is :  " << run_time.count() << std::endl;

std::vector<double> true_u;
true_solution(true_u, inner_points);

auto biggest_numerical_jacobi = std::max_element(std::begin(x_jacobi), std::end(x_jacobi));
auto biggest_analytical = std::max_element(std::begin(true_u), std::end(true_u));

std::cout << "Discrete maximum norm b/n anlytical and numerical solution using Jacobi is ::   " << abs(*biggest_numerical_jacobi - *biggest_analytical) << std::endl;

//Checking Gauss-Seidel method
start_time = std::chrono::high_resolution_clock::now();
do
    {
        Gauss_Seidel_step(obj_FD_LS_A, x_GS, b_RHS_vector);
        residue(obj_FD_LS_A, x_GS, b_RHS_vector, residue_vector_GS);
        j = j + 1;

    } while ((vector_euclidean_norm(residue_vector_GS, residue_vector_GS) / initial_residue > relative_accuracy) && (j <= maximun_iterations));

std::cout << " the number of iterations for Gauss-Seidel is ::   " << j << std::endl;
end_time = std::chrono::high_resolution_clock::now();
run_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
std::cout << "Time in milliseconds for Gauss-Seidel is :  " << run_time.count() << std::endl;

auto biggest_numerical_GS = std::max_element(std::begin(x_GS), std::end(x_GS));
std::cout << "Discrete maximum norm b/n anlytical and numerical solution using Gauss Seidel is ::   " << abs(*biggest_numerical_GS - *biggest_analytical) << std::endl;

//Checking symmetric Gauss-Seidel method
start_time = std::chrono::high_resolution_clock::now();
do
    {
        Gauss_Seidel_step_symmetric(obj_FD_LS_A, x_GS_symmetric, b_RHS_vector);
        residue(obj_FD_LS_A, x_GS_symmetric, b_RHS_vector, residue_vector_GS_symmetric);
        l = l + 1;

    } while ((vector_euclidean_norm(residue_vector_GS_symmetric, residue_vector_GS_symmetric) / (initial_residue) > relative_accuracy) && (l <= maximun_iterations));

std::cout << " the number of iterations for symmetric Gauss-Seidel is ::   " << l << std::endl;
end_time = std::chrono::high_resolution_clock::now();
run_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
std::cout << "Time in milliseconds for symmetric Gauss-Seidel is :  " << run_time.count() << std::endl;

auto biggest_numerical_symmetric_GS = std::max_element(std::begin(x_GS_symmetric), std::end(x_GS_symmetric));
std::cout << "Discrete maximum norm b/n anlytical and numerical solution using Gauss Seidel is ::   " << abs(*biggest_numerical_symmetric_GS - *biggest_analytical) << std::endl;

}

int main(int argc, char *argv[])
{
    int inner_points = 20;    //default inner_points
    unsigned int maximun_iterations = 100000;  //for while loop exit condition.

    if(argc == 2)
	{
        inner_points = std::stoi(argv[1]);
	}

    std::cout << "The number of inner points are ::   " << inner_points << std::endl;

    // Calling results of all the methods
    result(inner_points, maximun_iterations);

    return 0;
}