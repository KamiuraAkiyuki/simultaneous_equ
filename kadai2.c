#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 7
#define TRUE 1
#define FALSE 0

// -------------------------------------------------------------------------
//     プロトタイプ宣言
// -------------------------------------------------------------------------

void print_array(double a[N][N]);
void print_vector(double v[N]);
int is_lower_triangler_matrix(double mat[N][N]);
int is_upper_triangler_matrix(double mat[N][N]);
void mat_mlt(double in1[N][N], double in2[N][N], double out[N][N]);
void forward_substitution(double coefficient_matrix[N][N], double right_hand_vector[N], double solution_vector[N]);
void backward_substitution(double coefficient_matrix[N][N], double right_hand_vector[N], double solution_vector[N]);
void lu_decomposition(double a_mat[N][N], double l_mat[N][N], double u_mat[N][N]);
void solute_SLE_with_n_variables(double coefficient_matrix[N][N], double right_hand_vector[N], double solution_vector[N]);
void evaluate_reverse_matrix(double a_matrix[N][N], double reverse_matrix[N][N]);


// -------------------------------------------------------------------------
//     関数の定義
// -------------------------------------------------------------------------

// NxN配列をコマンドラインに表示する
void print_array(double a[N][N]){
    int i, j;
    for (i = 0; i < N; i++){
        for (j = 0; j < N; j++){
            if (j == 0) printf("| ");
            printf("%f ", a[i][j]);
            if (j == N-1) printf("|");
        }
        puts("");
    }
    puts("");
}

// N次元ベクトルをコマンドラインに表示する
void print_vector(double v[N]){
    int i;
    printf("| ");
    for (i = 0; i < N; i++){
        printf("%f ", v[i]);
    }
    puts("|");
}

// NxN配列が下三角行列か判断する
int is_lower_triangler_matrix( double mat[N][N] ){
    int i, j;
    for (i = 0; i < N; i++){
        for(j = 0; j < N; j++){
            if (i == j){ // 対角成分
                if (mat[i][j] != 1.0) return FALSE; 
            } else if (i < j){ // 上の部分
                if (mat[i][j] != 0.0) return FALSE;
            }
        }
    }
    return TRUE;
}

// NxN配列が上三角行列か判断する
int is_upper_triangler_matrix( double mat[N][N] ){
    int i, j;
    for (i = 0; i < N; i++){
        for(j = 0; j < N; j++){
            if (i > j){ // 下の部分
                if (mat[i][j] != 0.0) return FALSE;
            }
        }
    }
    return TRUE;
}

// 行列の積
void mat_mlt( double in1[N][N], double in2[N][N], double out[N][N] ) {
    int i, j, k;
    for (i = 0; i < N; i++){
        for (j = 0; j < N; j++){
            out[i][j] = 0;
        }
    }

    for (i = 0; i < N; i++){
        for (j = 0; j < N; j++){
            for (k = 0; k < N; k++){
                out[i][k] += in1[i][j] * in2[j][k];
            }
        }
    }
}

// 前進代入
void forward_substitution(double coefficient_matrix[N][N], double right_hand_vector[N], double solution_vector[N]){

    int i, j;
    for (i = 0; i < N; i++){
        double sum = 0;
        for (j = 0; j < i; j++){
            sum += coefficient_matrix[i][j]*solution_vector[j];
        }
        solution_vector[i] = right_hand_vector[i] - sum;
    }
}

// 後退代入
void backward_substitution(double coefficient_matrix[N][N], double right_hand_vector[N], double solution_vector[N]){

    int i, j;
    for (i = N-1; i >= 0; i--){
        double sum = 0;
        for (j = N-1; j > i; j--){
            sum += coefficient_matrix[i][j]*solution_vector[j];
        }
        solution_vector[i] = (right_hand_vector[i] - sum)/ coefficient_matrix[i][i];
    }
}  

// LU分解　
void lu_decomposition(double a_mat[N][N], double l_mat[N][N], double u_mat[N][N]){
    int i, j, k;
    double sum;
    for (i = 0; i < N; i++){
        for (j = 0; j < N; j++){
            if (i <= j){ // 下三角成分
                sum = 0;
                if (i != 0){
                    for (k = 0; k <= j-1; k++){
                        sum += l_mat[i][k]*u_mat[k][j];
                    } 
                }
                u_mat[i][j] = a_mat[i][j] - sum;
                if (i == j) l_mat[i][j] = 1.0;
            } else if(i > j){ // 上三角成分
                sum = 0;
                if (j != 0){
                    for (k = 0; k <= j-1; k++){
                        sum += l_mat[i][k]*u_mat[k][j];
                    }
                }
                l_mat[i][j] = (a_mat[i][j] - sum)/u_mat[j][j];
            }
        }
    }
}

// simultaneous linear equation 
void solute_SLE_with_n_variables(double coefficient_matrix[N][N], double right_hand_vector[N], double solution_vector[N]){
    double l_matrix[N][N] = {0};
    double u_matrix[N][N] = {0};
    double mid_solution_vector[N];

    lu_decomposition(coefficient_matrix, l_matrix, u_matrix);
    forward_substitution(l_matrix, right_hand_vector, mid_solution_vector);
    backward_substitution(u_matrix, mid_solution_vector, solution_vector);
}

// 逆行列を求める
void evaluate_reverse_matrix(double a_matrix[N][N], double reverse_matrix[N][N]){
    
    double b_matrix[N][N];;

    int i, j;
    for (i = 0; i < N; i++){
        for (j = 0; j < N; j++){
            if (i == j){
                b_matrix[i][j] = 1.0;
            } else {
                b_matrix[i][j] = 0.0;
            }
        }
    }

    for (i = 0; i < N; i++){
        double x[N] = {0};
        solute_SLE_with_n_variables(a_matrix, b_matrix[i], x);
        print_vector(x);
        for (j = 0; j < N; j++){
            reverse_matrix[j][i] = x[j];
        }
    }
}

// -------------------------------------------------------------------------
//    main
// -------------------------------------------------------------------------

int main(){
    int kadai;

    // 課題番号指定
    while(1){
        printf("課題番号: ");
        scanf("%d", &kadai);
        if (kadai < 1 || kadai > 8 || kadai == 7){
            printf("1~8!!\n");
        } else {
            break;
        }
    }

    if (kadai == 1){

        double coefficient_mat[N][N] = {{1, 0, 0}, {3, 1, 0}, {-2, 2, 1}};
        double right_hand_vec[N] = {2, 3, -1};
        double solution_vec[N];

        printf("--------------- L行列 ---------------\n");
        print_array(coefficient_mat);
        printf("------------- 右辺ベクトル ------------\n");
        print_vector(right_hand_vec);

        forward_substitution(coefficient_mat, right_hand_vec, solution_vec);

        printf("------------- 解ベクトル --------------\n");
        print_vector(solution_vec);

    } 
    else if (kadai == 2)
    {

        double coefficient_mat[N][N] = {{2, 1, -1}, {0, 3, 2}, {0, 0, -3}};
        double right_hand_vec[N] = {2, -3, 9};
        double solution_vec[N];

        printf("--------------- U行列 ----------------\n");
        print_array(coefficient_mat);
        printf("------------- 右辺ベクトル -------------\n");
        print_vector(right_hand_vec);

        backward_substitution(coefficient_mat, right_hand_vec, solution_vec);

        printf("-------------- 解ベクトル -------------\n");
        print_vector(solution_vec);

    }
    else if (kadai == 3)
    {

        double coefficient_mat[N][N];
        double l_matrix[N][N] = {{1, 0, 0}, {3, 1, 0}, {-2, 2, 1}};
        double u_matrix[N][N] = {{2, 1, -1}, {0, 3, 2}, {0, 0, -3}};

        mat_mlt(l_matrix, u_matrix, coefficient_mat);

        printf("---------- 係数行列 ----------\n");
        print_array(coefficient_mat);

    }
    else if (kadai == 4)
    {

        double coefficient_mat[N][N] = {{2, 1, -1}, {6, 6, -1}, {-4, 4, 3}};
        double l_matrix[N][N] = {0};
        double u_matrix[N][N] = {0};

        lu_decomposition(coefficient_mat, l_matrix, u_matrix);

        printf("---------- L行列 ----------\n");
        print_array(l_matrix);
        printf("---------- U行列 ----------\n");
        print_array(u_matrix);

    }
    else if (kadai == 5)
    {

        double coefficient_mat[N][N] = {{2, 1, -1}, {6, 6, -1}, {-4, 4, 3}};
        double right_hand_vec[N] = {2, 3, -1};
        double solution_vec[N];

        printf("---------- 係数行列 ----------\n");
        print_array(coefficient_mat);
        printf("---------- 右辺ベクトル ----------\n");
        print_vector(right_hand_vec);

        solute_SLE_with_n_variables(coefficient_mat, right_hand_vec, solution_vec);

        printf("---------- 解ベクトル ----------\n");
        print_vector(solution_vec);

    }
    else if (kadai == 6)
    {

        double h_coefficient_mat[N][N];
        double right_hand_vec[N];
        double solution_vec[N];

        int i, j;
        for (i = 0; i < N; i++){
            for (j = 0; j < N; j++){
                h_coefficient_mat[i][j] = pow(0.5, abs(i-j));
            }
        }
        printf("--------------- 行列H ---------------\n");
        print_array(h_coefficient_mat);

        for (i = 0; i < N; i++){
            right_hand_vec[i] = 3 - pow(2, i-N+1) - pow(2, -i);
        }
        printf("------------- 右辺ベクトル -------------\n");
        print_vector(right_hand_vec);

        solute_SLE_with_n_variables(h_coefficient_mat, right_hand_vec, solution_vec);

        printf("---------- 解ベクトル ----------\n");
        print_vector(solution_vec);

    }
    else if (kadai == 8)
    {

        double h_mat[N][N];
        double h_reverse_mat[N][N] = {0};
        double solution_mat[N][N];

        int i, j;
        for (i = 0; i < N; i++){
            for (j = 0; j < N; j++){
                h_mat[i][j] = pow(0.5, abs(i-j));
            }
        }

        evaluate_reverse_matrix(h_mat, h_reverse_mat);

        mat_mlt(h_mat, h_reverse_mat, solution_mat);
        printf("--------------- 行列H ---------------\n");
        print_array(h_mat);
        printf("------------ 行列Hの逆行列 ------------\n");
        print_array(h_reverse_mat);
        printf("---------- 解行列 ----------\n");
        print_array(solution_mat);
    } 

    return 0;
}