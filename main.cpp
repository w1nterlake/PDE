#include <iostream>
#include <cmath>
#include <ctime>
#include <fstream>

void task(int N, int M, double *u){
    int i = 0, j = 0;
    double del_t = 1. / (N - 1);
    double del_x = 1. / (M - 1);
    std::ofstream fout;
    fout.open("text.txt");
    for(i = 0; i < M; i++){
        u[i] = 0;
        fout << 0 << " " << i * del_x << " " << u[i] << "\n";
    }
    u[M] = sin(3 * del_t) * sin(3 * del_t);
    fout << del_t << " " << 0 << " " << u[M] << "\n";
    //std::cout << u[M] << "\n";
    for(i = 1; i < M - 1; i++){
        u[M + i] = 0;
        // u[M + i] = del_t * del_t * (u[i + 1] - 2 * u[i]  + u[i - 1]) / (2 * del_x * del_x);
        fout << del_t << " " << i * del_x << " " << u[M + i] << "\n";
        //std::cout << u[M + i] << "\n";
    }


    u[2 * M - 1] = u[2 * M - 2] + del_t * del_x / (1 + del_t * del_t);
    //std::cout << u[2 * M - 1] << "\n";
    fout << del_t << " " << 1 << " " << u[2 * M - 1] << "\n";
    for(i = 2; i < N; i++){
        u[i * M] =  sin(3 * i * del_t) * sin(3 * i * del_t);
        //std::cout  << u[i * M];
        fout << i * del_t << " " << 0 << " " << u[i * M] << "\n";
        for(j = 1; j < M - 1; j++){
          
            u[i * M + j] = 2 * u[(i - 1) * M + j] - u[(i - 2) * M + j] + del_t * del_t * (u[(i - 1) * M + j + 1] - 2 * u[(i - 1) * M + j] + u[(i - 1) * M + j - 1]) / (del_x * del_x) + del_t * del_t  * (u[(i - 1) * M + j + 1] - u[(i - 1) * M + j - 1]) /del_x + del_t * del_t * u[(i - 1) * M + j];
           // std::cout << u[i * M + j] << "\n";
           fout << i * del_t << " " << j * del_x << " " << u[i * M + j] << "\n"; 
        }
        u[i * M + j] = u[i * M + j - 1] + i * del_t * del_t / (1 + i * i * del_t * del_t); 
       // std::cout << u[i * M + j] << "\n";
       fout << i * del_t << " " << 1 << " " << u[i * M + j] << "\n";
    }
    fout.close();
}

int main(){
    int N = 101, M = 301;
    double *u1, *u2, *u3;
    u1 = (double*) malloc(201 * 101 * sizeof(double));
    u2 = (double*)malloc(201 * 401 * sizeof(double));
    u3 = (double*) malloc(401 * 801 *sizeof(double));
    task(201, 101, u1);
    task(401, 201, u2);
    task(801, 401, u3);
    double p = (u1[101 * (201 - 1) + 25] - u2[201 * (401 - 1) + 50]) ;
    double q = (u2[201 * (401 - 1) + 50] - u3[401 * (801 - 1) + 100]); 
    std::cout << "(u(1/4, 1, 1/ 100) - u(1/4, 1, 1/200)) / (u(1/4, 1, 1/200) - u(1/4, 1, 1/800)) = " << p/ q<< "\n";
    p = (u1[101 * (201 - 1) + 50] - u2[201 * (401 - 1) + 100]) ;
    std::cout << "(u(1/2, 1, 1/ 100) - u(1/2, 1, 1/200)) / (u(1/2, 1, 1/200) - u(1/2, 1, 1/800)) = " << p / (u2[201 * (401 - 1) + 100] - u3[401 * (801 - 1) + 200]) << "\n";
    std::cout << "(u(3/4, 1, 1/ 100) - u(3/4, 1, 1/200)) / (u(3/4, 1, 1/200) - u(3/4, 1, 1/800)) = " << (u1[101 * (201 - 1) + 75] - u2[201 * (401 - 1) + 150]) / (u2[201 * (401 - 1) + 150] - u3[401 * (801 - 1) + 300]) << "\n";
    std::cout << "(u(1, 1, 1/ 100) - u(1, 1, 1/200)) / (u(1, 1, 1/200) - u(1, 1, 1/800)) = " << (u1[101 * (201 - 1) + 100] - u2[201 * (401 - 1) + 200]) / (u2[201 * (401 - 1) + 200] - u3[401 * (801 - 1) + 400]) << "\n";

    
    // for(int i = 0; i < N; i++){
    //     for(int j = 0; j < M; j++){
    //         std::cout << u[i * M + j] << " ";
    //     }
    //     std::cout << "\n";
    // }

    free(u1);
    free(u2);
    free(u3);
    return 0;
}