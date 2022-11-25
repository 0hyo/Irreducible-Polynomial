#include <stdio.h>
#include <math.h>
#include <time.h>
#include "flint.h"
#include "fq_nmod_poly.h"
#include "fq_nmod_mat.h"
#include "fq_nmod.h"
#include "nmod_poly.h"

#include "fmpz.h"

/**
 * gcc -Wall -O2 *.c -L/opt/homebrew/Cellar/mpfr/4.1.0/lib -L/opt/homebrew/Cellar/gmp/6.2.1_1/lib -I/opt/homebrew/Cellar/gmp/6.2.1_1/include -I/opt/homebrew/Cellar/mpfr/4.1.0/include -I/usr/local/include/flint/ -lflint -lmpfr -lgmp -lpthread
 */

/**
 * @brief 
 * 
 * @param f  판별 대상 다항식
 * @param ctx 
 * @param irr 
 * @return int 1 or 0 
 */
int ben_or(fq_nmod_poly_t f, fq_nmod_ctx_t ctx)// 이진 유한체에 대한 기약다항식 판별만 가능
{
    // printf("\n@@@@@@@@@@@@@@@@  ben or  start @@@@@@@@@@@@@@@@\n");
    fq_nmod_poly_t gcd, x_q, x_q_i, x_q_i_1, poly_x; 
    fq_nmod_t zero, one;   
    fq_nmod_poly_init(x_q_i, ctx);// a = 0, x_q_i = x^q^n - x;
    fq_nmod_poly_init(gcd, ctx);fq_nmod_poly_init(x_q_i_1, ctx);
    fq_nmod_poly_init(poly_x, ctx);  fq_nmod_poly_init(x_q,ctx);
    long m = ctx->j[ctx->len - 1];   // m : 기약다항식의 차수
    long q = pow(2,m); // 2 ^ m = q
    // F_q의 원소 0, 1
    fq_nmod_init(zero, ctx); 
    fq_nmod_init(one, ctx); 
    fq_nmod_zero(zero, ctx);
    fq_nmod_one(one, ctx);

    int n = fq_nmod_poly_degree(f, ctx);
    fq_nmod_poly_set_coeff(x_q, q, one, ctx);   //x^q
    fq_nmod_poly_set_coeff(x_q_i_1, q, one, ctx);//x^q
    fq_nmod_poly_set_coeff(poly_x, 1, one, ctx); //x

    fmpz_t i, d, fmpz_1, fmpz_q, fmpz_q_i;
    fmpz_init(i); fmpz_init(d); fmpz_init(fmpz_1); fmpz_init(fmpz_q); fmpz_init(fmpz_q_i);
    fmpz_one(i); fmpz_one(fmpz_1);
    fmpz_set_si(fmpz_q, q);
    fmpz_init_set_ui(d, floor(n/2));


    while(fmpz_cmp(i,d) <= 0)
    {
        // printf("i = "); fmpz_print(i); printf("\n");
        // printf("d = "); fmpz_print(d); printf("\n");

        fmpz_pow_fmpz(fmpz_q_i,fmpz_q,i); //q^i
        // printf("q^i = "); fmpz_print(fmpz_q_i); printf("\n");

        fq_nmod_poly_powmod_fmpz_binexp(x_q_i, poly_x, fmpz_q_i , f, ctx); // x^(q^i) 
        fq_nmod_poly_add(x_q_i_1, x_q_i, poly_x, ctx);
        // printf("x^q^i = "); fq_nmod_poly_print_pretty(x_q_i,"X",ctx); printf("\n");
        // printf("x^q^i-x = "); fq_nmod_poly_print_pretty(x_q_i_1,"X",ctx); printf("\n");
        // printf("gcd start = \n");

        // fq_nmod_poly_gcd(gcd, f, x_q_i_1, ctx);
        fq_nmod_poly_gcd_euclidean(gcd, f, x_q_i_1, ctx);

        // printf("gcd end = \n"); 
        // printf("1x_q_n : "); fq_nmod_poly_print_pretty(x_q_i,"X",ctx); printf("\n");
        // printf("f : "); fq_nmod_poly_print_pretty(f,"X",ctx); printf("\n");
        // printf("gcd : "); fq_nmod_poly_print_pretty(gcd,"X",ctx); printf("\n");
        // // printf("\nis one = %d\n",fq_nmod_poly_is_one(gcd, ctx));        
        if(!fq_nmod_poly_is_one(gcd, ctx)){
            
            fq_nmod_poly_clear(gcd, ctx);
            fq_nmod_poly_clear(x_q, ctx);
            fq_nmod_poly_clear(x_q_i, ctx);
            fq_nmod_poly_clear(x_q_i_1, ctx);
            fq_nmod_poly_clear(poly_x, ctx);

            fq_nmod_clear(one, ctx);
            fq_nmod_clear(zero, ctx);
  


            fmpz_clear(i);
            fmpz_clear(d);
            fmpz_clear(fmpz_1);
            fmpz_clear(fmpz_q);
            fmpz_clear(fmpz_q_i);
            return 0;
        }
        // printf("2x_q_n : "); fq_nmod_poly_print_pretty(x_q_i,"X",ctx); printf("\n");
        // printf("3x_q_n : "); fq_nmod_poly_print_pretty(x_q_i,"X",ctx); printf("\n");
        fmpz_add_ui(i, i, 1);
    }
    // printf("m = %ld\n",m);
    // printf("q = %ld\n",q);

    fq_nmod_poly_clear(gcd,ctx);
    fq_nmod_poly_clear(x_q,ctx);
    fq_nmod_poly_clear(x_q_i,ctx);
    fq_nmod_poly_clear(x_q_i_1,ctx);
    fq_nmod_poly_clear(poly_x,ctx);

    fq_nmod_clear(one, ctx);
    fq_nmod_clear(zero, ctx);

    fmpz_clear(i);
    fmpz_clear(d);
    fmpz_clear(fmpz_1);
    fmpz_clear(fmpz_q);
    fmpz_clear(fmpz_q_i);
    return 1;
}

void func_Test()
{
    fq_nmod_poly_t x, y, z;
    fq_nmod_ctx_t ctx;
    fq_nmod_t f, g, irr;
    fmpz_t p;
    fmpz_init(p);
    fmpz_set_si(p, 2);
    nmod_poly_init(f,  2);
    nmod_poly_init(g,  2);
    nmod_poly_init(irr,  2);
    nmod_poly_set_str(f, "2 2 1 1");
    nmod_poly_set_str(g, "2 2 0 1");
    nmod_poly_set_str(irr, "13 2 1 0 0 1 0 0 0 0 0 0 0 0 1");
    fq_nmod_ctx_init_modulus(ctx, irr, "z");
    printf("irr = ");    nmod_poly_print_pretty(irr,"z"); printf("\n");
    fq_nmod_poly_init(x, ctx);    
    fq_nmod_poly_init(y, ctx);    
    fq_nmod_poly_init(z, ctx);    
    fq_nmod_poly_set_coeff(x, 1, f, ctx);
    fq_nmod_poly_set_coeff(x, 0, g, ctx);
    fq_nmod_poly_set_coeff(y, 1, g, ctx);


/**
    nmod_poly_init(f,  2);
    nmod_poly_init(g,  5);
    fq_nmod_poly_init(x,ctx);
    fq_nmod_poly_init(y,ctx);
    nmod_poly_set_coeff_ui (f , 3 , 1) ;
    nmod_poly_set_coeff_ui (f , 0 , 1) ;
    nmod_poly_set_str(g, "2 2 1 1");
    nmod_poly_mul(f, f, g);
 */


    // fq_nmod_poly_set_coeff(y, 0, g, ctx);
    printf("ctx \n");    fq_nmod_ctx_print(ctx);
    // fq_nmod_poly_print_pretty(x,"x",ctx); printf("\n");
    // void fq_nmod_poly_pow(fq_nmod_poly_t rop, const fq_nmod_poly_t op, ulong e, const fq_nmod_ctx_t ctx)
    for(int i = 0; i < 13; i++)
    {
        // fq_nmod_poly_mul_classical(x,x,y,ctx); printf("x = ");    fq_nmod_poly_print_pretty(x,"x",ctx); printf("\n");
        // fq_nmod_poly_gcd(y,x,y,ctx); printf("y = "); fq_nmod_poly_print_pretty(y,"x",ctx); printf("\n");
    }

    fq_nmod_poly_clear(x, ctx);
    fq_nmod_poly_clear(y, ctx);    
    fq_nmod_poly_clear(z, ctx);    
    // flint_randclear(state);

    nmod_poly_clear(f);
    nmod_poly_clear(g);

    fq_nmod_ctx_clear(ctx);
    fmpz_clear(p);
}

void ben_or_test(int Case, int times)
{
    printf("\n@@@@@@@@@@@@@@@@  ben or test start @@@@@@@@@@@@@@@@\n");
    fq_nmod_poly_t x, rand_poly;   // F_q[X]
    fq_nmod_ctx_t ctx;             // ctx
    nmod_poly_t irr;
    fmpz_t p;
    flint_rand_s* rand_state;
    rand_state = flint_rand_alloc();
    flint_randinit(rand_state);
    flint_randseed(rand_state, rand(), rand());

    fmpz_init(p);
    fmpz_set_ui(p, 2);
    nmod_poly_init(irr,  2);
    if(Case == 1)
        fq_nmod_ctx_init(ctx, p, 12, "z");   
    if(Case == 2 ||Case == 3 ||Case == 4 ||Case == 5 )
        fq_nmod_ctx_init(ctx, p, 13, "z");   
    fq_nmod_init(irr,  ctx);


    slong deg = 0;
    if(Case == 1)
    {
        nmod_poly_set_str(irr, "13 2 1 0 0 1 0 0 0 0 0 0 0 0 1 ");
        fq_nmod_ctx_init_modulus(ctx, irr, "z");
        // printf("\ncase 1 ctx : ");    fq_nmod_ctx_print(ctx);printf("\n");
        deg = 65;
    }
    else if(Case == 2 ||Case == 3 ||Case == 4 ||Case == 5 )
    {
        nmod_poly_set_str(irr, "14 2 1 1 0 1 1 0 0 0 0 0 0 0 0 1 ");
        fq_nmod_ctx_init_modulus(ctx, irr, "z");
        // printf("\ncase %d ctx : ",Case);    fq_nmod_ctx_print(ctx);printf("\n");
        if(Case == 2)   //deg는 항의 갯수로 실제 차수 + 1을 저장해야한다.
            deg = 129;
        else if(Case == 3)
            deg = 120;
        else if(Case == 4)
            deg = 129;
        else if(Case == 5)
            deg = 129;
    }

    // printf("\nirr = ");    nmod_poly_print_pretty(irr, "x"); printf("\n");
    // printf("\nctx : ");    fq_nmod_ctx_print(ctx);printf("\n");

    fq_nmod_poly_init(x, ctx);     
    fq_nmod_poly_init(rand_poly, ctx);
    
    fq_nmod_poly_randtest(rand_poly, rand_state, deg, ctx);
    // printf("\nrand poly = ");   fq_nmod_poly_print_pretty(rand_poly, "X", ctx);
    // printf("deg of poly = %ld",fq_nmod_poly_degree(rand_poly, ctx));

    /*********      시간측정        ******/    
    clock_t start, end;
    start = clock();
    int k;  int i = 0;
    for (k = 0; k < times; k++)
    {
        // printf("k = %d\n",k);
        fq_nmod_poly_randtest(rand_poly,rand_state,deg,ctx);
        // ben-or start
        while(!ben_or(rand_poly, ctx))
        {
            // printf("함수 동작 횟수 = %d\n", i++);
            i++;
            fq_nmod_poly_randtest(rand_poly,rand_state,deg,ctx);
        }
        // /** 기약다항식 출력 */
        // printf("irr by ben or = ");    fq_nmod_poly_print_pretty(rand_poly, "X", ctx); printf("\n");
        
        // printf("%d번 만에 찾은 기약다항식\n", i);
        /** 기약다항식 확인(내장함수로 오류 판별) */
        // if(!fq_nmod_poly_is_irreducible(rand_poly, ctx))
        // {
        //     printf("************** error test case **************"); 
        //     fq_nmod_poly_print_pretty(rand_poly, "X", ctx); printf("\n");
        // }
    }
    end = clock();
    printf("\n평균 %d회 테스트 후 기약다항식 생성", i/times);
    printf("\n%d회 측정결과 1회 평균 시간 : %lu초\n",k, ((end-start)/(CLOCKS_PER_SEC))/k);

    // printf("Is irr %d\n", fq_nmod_poly_is_irreducible(rand_poly,ctx));
    // printf("Irr deg = %ld\n", fq_nmod_poly_degree(rand_poly, ctx));

    /** 변수 clear */
    flint_rand_free(rand_state);
    fq_nmod_poly_clear(rand_poly, ctx);
    fq_nmod_poly_clear(x,ctx);
    fq_nmod_ctx_clear(ctx);
    nmod_poly_clear(irr);
}

int rabin(fq_nmod_poly_t f, fq_nmod_ctx_t ctx, long* poly_deg_factors, long poly_deg_factors_num)
{
    // printf("\n@@@@@@@@@@@@@@@@  rabin  start @@@@@@@@@@@@@@@@\n");
    fq_nmod_poly_t gcd, x_q, x_q_i, x_q_i_1, poly_x; 
    fq_nmod_t zero, one;   
    fq_nmod_poly_init(x_q_i, ctx);// a = 0, x_q_i = x^q^n - x;
    fq_nmod_poly_init(gcd, ctx);fq_nmod_poly_init(x_q_i_1, ctx);
    fq_nmod_poly_init(poly_x, ctx);  fq_nmod_poly_init(x_q,ctx);
    long m = ctx->j[ctx->len - 1];   // m : 기약다항식의 차수
    long q = pow(2,m); // 2 ^ m = q
    // F_q의 원소 0, 1
    fq_nmod_init(zero, ctx); 
    fq_nmod_init(one, ctx); 
    fq_nmod_zero(zero, ctx);
    fq_nmod_one(one, ctx);

    int n = fq_nmod_poly_degree(f, ctx);        //f의 차수
    fq_nmod_poly_set_coeff(x_q, q, one, ctx);   //x^q
    fq_nmod_poly_set_coeff(poly_x, 1, one, ctx); //x

    fmpz_t  d, fmpz_1, fmpz_q, fmpz_q_i, fmpz_n_over_p;
    fmpz_init(d); fmpz_init(fmpz_1); fmpz_init(fmpz_q); fmpz_init(fmpz_q_i); fmpz_init(fmpz_n_over_p);
    fmpz_one(fmpz_1);
    fmpz_set_si(fmpz_q, q);
    fmpz_init_set_ui(d, n); // d == deg of f ( fq_nmod_poly )

    for (int i = 0; i < poly_deg_factors_num; i++)
    {
        // printf("\nDeg  = %d\n",n);
        // printf("\n n / p_[%d] = %ld",i, (n/poly_deg_factors[i]));
        fmpz_set_ui(fmpz_n_over_p, (n/poly_deg_factors[i]) );
        fmpz_pow_fmpz(fmpz_q_i,fmpz_q, fmpz_n_over_p); //q^ (n/p)

        // printf("\nfmpz_ n / p = "); fmpz_print(fmpz_n_over_p); 
        // printf("\nfmpz_q = "); fmpz_print(fmpz_q); 


        fq_nmod_poly_powmod_fmpz_binexp(x_q_i, poly_x, fmpz_q_i , f, ctx); // x^(q^(n/p)) 
        fq_nmod_poly_add(x_q_i_1, x_q_i, poly_x, ctx);
        // printf("\n x_q_i_1 = "); fq_nmod_poly_print_pretty(x_q_i_1,"X", ctx); printf("\n");
        // printf("\n f = "); fq_nmod_poly_print_pretty(f,"X", ctx); printf("\n");
        fq_nmod_poly_gcd_euclidean(gcd, f, x_q_i_1, ctx);
        // printf("\n gcd = "); fq_nmod_poly_print_pretty(gcd,"X", ctx); printf("\n");
        if(!fq_nmod_poly_is_one(gcd, ctx)){
            // printf("\nis one\n");
            fq_nmod_poly_clear(gcd, ctx);
            fq_nmod_poly_clear(x_q, ctx);
            fq_nmod_poly_clear(x_q_i, ctx);
            fq_nmod_poly_clear(x_q_i_1, ctx);
            fq_nmod_poly_clear(poly_x, ctx);

            fq_nmod_clear(one, ctx);
            fq_nmod_clear(zero, ctx);
  
            fmpz_clear(d);
            fmpz_clear(fmpz_1);
            fmpz_clear(fmpz_q);
            fmpz_clear(fmpz_q_i);
            fmpz_clear(fmpz_n_over_p);
            return 0;
        }
    }
    fmpq_pow_fmpz(fmpz_q_i, fmpz_q, d);    
    fq_nmod_poly_powmod_fmpz_binexp(x_q_i, poly_x, fmpz_q_i , f, ctx); // x^(q^d) 
    fq_nmod_poly_add(x_q_i_1, x_q_i, poly_x, ctx);
    // printf("\ng = "); fq_nmod_poly_print_pretty(x_q_i_1,"X",ctx);
    if(fq_nmod_poly_is_zero(x_q_i_1, ctx))
    {
        fq_nmod_poly_clear(gcd,ctx);
        fq_nmod_poly_clear(x_q,ctx);
        fq_nmod_poly_clear(x_q_i,ctx);
        fq_nmod_poly_clear(x_q_i_1,ctx);
        fq_nmod_poly_clear(poly_x,ctx);

        fq_nmod_clear(one, ctx);
        fq_nmod_clear(zero, ctx);

        fmpz_clear(d);
        fmpz_clear(fmpz_1);
        fmpz_clear(fmpz_q);
        fmpz_clear(fmpz_q_i);
        fmpz_clear(fmpz_n_over_p);

        return 1;
    }


    fq_nmod_poly_clear(gcd,ctx);
    fq_nmod_poly_clear(x_q,ctx);
    fq_nmod_poly_clear(x_q_i,ctx);
    fq_nmod_poly_clear(x_q_i_1,ctx);
    fq_nmod_poly_clear(poly_x,ctx);

    fq_nmod_clear(one, ctx);
    fq_nmod_clear(zero, ctx);

    fmpz_clear(d);
    fmpz_clear(fmpz_1);
    fmpz_clear(fmpz_q);
    fmpz_clear(fmpz_q_i);
    fmpz_clear(fmpz_n_over_p);
    return 0;
}

void rabin_test(int Case, int times)
{
    printf("\n@@@@@@@@@@@@@@@@  rabin test start @@@@@@@@@@@@@@@@\n");
    fq_nmod_poly_t  rand_poly;   // F_q[X]
    fq_nmod_ctx_t ctx;             // ctx
    fq_nmod_t h;  //F_2[X]  
    nmod_poly_t irr;
    fmpz_t p;
    flint_rand_s* rand_state;
    rand_state = flint_rand_alloc();
    flint_randinit(rand_state);
    flint_randseed(rand_state, rand(), rand());

    fmpz_init(p);
    fmpz_set_ui(p, 2);

    fq_nmod_ctx_init(ctx, p, 13, "z");
    fq_nmod_init(h,  ctx);
    fq_nmod_init(irr,  ctx);
    
    nmod_poly_init(irr,  2);
    if(Case == 1)
        fq_nmod_ctx_init(ctx, p, 12, "z");   
    if(Case == 2 ||Case == 3 ||Case == 4 ||Case == 5 )
        fq_nmod_ctx_init(ctx, p, 13, "z");   
    fq_nmod_init(irr,  ctx);

    slong deg = 0;
    if(Case == 1)
    {
        nmod_poly_set_str(irr, "13 2 1 0 0 1 0 0 0 0 0 0 0 0 1 ");
        fq_nmod_ctx_init_modulus(ctx, irr, "z");
        // printf("\ncase 1 ctx : ");    fq_nmod_ctx_print(ctx);printf("\n");
        deg = 65;
    }
    else if(Case == 2 ||Case == 3 ||Case == 4 ||Case == 5 )
    {
        nmod_poly_set_str(irr, "14 2 1 1 0 1 1 0 0 0 0 0 0 0 0 1 ");
        fq_nmod_ctx_init_modulus(ctx, irr, "z");
        // printf("\ncase %d ctx : ",Case);    fq_nmod_ctx_print(ctx);printf("\n");
        if(Case == 2)   //deg는 항의 갯수로 실제 차수 + 1을 저장해야한다.
            deg = 97;
        if(Case == 3)
            deg = 129;
        if(Case == 4)
            deg = 120;
        if(Case == 5)
            deg = 129;
    }    
    

    /** rand poly의 차수에 따른 약수와 그 갯수 */
    long n[2], num_factor;
    if(((deg - 1) == 128) || ((deg - 1) == 64) ){
            n[0] = 2;
            num_factor = 1;
        }
    else if ((deg - 1) == 96){
            n[0] = 2;
            n[1] = 3;
            num_factor = 2;
        }
    else if ((deg - 1) == 119){
            n[0] = 7;
            n[1] = 19;
            num_factor = 2;
        }


    fq_nmod_poly_init(rand_poly, ctx);    
    /*********      시간측정        ******/    
    clock_t start, end;
    start = clock();
    int cnt_i;  int i = 0;
    for(cnt_i = 0; cnt_i< times; cnt_i++)
    {
        // printf("k = %d\n",cnt_i);
        fq_nmod_poly_randtest(rand_poly, rand_state, deg, ctx);
        while(!rabin(rand_poly, ctx, n, num_factor))
        {
            // printf("함수 동작 횟수 = %d\n", i++);
            i++;
            fq_nmod_poly_randtest(rand_poly, rand_state, deg, ctx);
            // if(fq_nmod_poly_is_irreducible(rand_poly, ctx)){// 오류 확인.
            //     printf("\nis irreducible %d\n",rabin(rand_poly, ctx, n, num_factor)); 
            //     break;
            // }
        }
        // fq_nmod_poly_print_pretty(irr,"X", ctx);
        // /** 기약다항식 출력 */
        // printf("irr by rabin = ");    fq_nmod_poly_print_pretty(rand_poly, "X", ctx); printf("\n");
        // /** 기약다항식 확인 */
        // printf("%d 번만에 찾은 기약다항식\n",i);
        // printf("is irreducible %d\n",fq_nmod_poly_is_irreducible(rand_poly, ctx)); 
    }
    end = clock();
    printf("\n평균 %d회 테스트 후 기약다항식 생성", i/times);
    printf("\n%d회 측정결과 1회 평균 시간 : %lu 초\n",cnt_i, ((end-start)/(CLOCKS_PER_SEC))/cnt_i);

    // printf("\ndeg of rand poly %ld", fq_nmod_poly_degree(rand_poly, ctx));
    // printf("\nrand poly = ");    fq_nmod_poly_print_pretty(rand_poly, "X", ctx);

    /** 변수 clear */
    flint_rand_free(rand_state);
    fq_nmod_poly_clear(rand_poly, ctx);
    fq_nmod_ctx_clear(ctx);
    nmod_poly_clear(h);
    nmod_poly_clear(irr);
}

int ClassicMC_Irr_Gen(fq_nmod_poly_t ret_Irr, fq_nmod_poly_t CM_Parameter, fq_nmod_ctx_t ctx)
{
    fq_nmod_mat_t Matrix;

    fq_nmod_poly_t Beta;
    nmod_poly_t coeff_1;
    slong deg;

    /**  random   */
    flint_rand_s* rand_state;
    rand_state = flint_rand_alloc();
    flint_randinit(rand_state);
    flint_randseed(rand_state, rand(), rand());
 
    deg = fq_nmod_poly_degree(CM_Parameter, ctx);
    // printf("\ndeg = %lu \n",deg);
    
    fq_nmod_mat_init(Matrix, deg, deg + 1, ctx);
    fq_nmod_poly_init(Beta, ctx);
    nmod_poly_init(coeff_1, 2);

    nmod_poly_set_coeff_ui(coeff_1, 0, 1);
    fq_nmod_mat_entry_set(Matrix, 0, 0, coeff_1, ctx);    
    // printf("\n MAT\n"); fq_nmod_mat_print_pretty(Matrix, ctx);

    // printf("CM Parameter = ");    fq_nmod_poly_print_pretty(CM_Parameter,"X", ctx);
    fq_nmod_poly_randtest(Beta, rand_state, deg, ctx);

    // fmpz_t temp;    fmpz_init(temp);
    // for(int i = 0; i <= deg; i++)
    // {
    //     for(int j = 0; j <= fq_nmod_ctx_degree(ctx); j++)
    //     {
    //         fmpz_randbits(temp, rand_state, fq_nmod_ctx_degree(ctx));
    //     }
    // }

    // printf("Beta = ");    fq_nmod_poly_print_pretty(Beta, "x", ctx);    printf("\n");

    fq_nmod_t temp;
    fq_nmod_init(temp, ctx);
    fq_nmod_poly_t temp_poly1, temp_poly2;
    fq_nmod_poly_init(temp_poly1, ctx);    fq_nmod_poly_init(temp_poly2, ctx);
    fq_nmod_poly_set(temp_poly1, Beta, ctx);
    for(int j = 1; j < deg + 1; j++)
    {
        // printf("\ntemp_poly1 = "); fq_nmod_poly_print_pretty(temp_poly1, "x", ctx);
        for (int i = 0; i < deg; i++)
        {
            fq_nmod_poly_get_coeff(temp, temp_poly1, i, ctx);
            fq_nmod_mat_entry_set(Matrix, i, j , temp, ctx);    
        }
        // printf("\n\tbeta = "); fq_nmod_poly_print_pretty(Beta, "x", ctx); printf("\n");
        // fq_nmod_poly_mulmod(temp_poly1, temp_poly1, Beta, CM_Parameter , ctx);
        // printf("print(  (() "); fq_nmod_poly_print_pretty(temp_poly1,"x", ctx); printf(")*("); fq_nmod_poly_print_pretty(Beta, "x", ctx); printf(")")


        fq_nmod_poly_mul(temp_poly2, temp_poly1, Beta, ctx);
        // printf("print(  ( "); fq_nmod_poly_print_pretty(Beta, "x", ctx); printf(") * (");   fq_nmod_poly_print_pretty(temp_poly1,"x",ctx);printf(") == ("); fq_nmod_poly_print_pretty(temp_poly2, "x", ctx); printf("))\n");
        fq_nmod_poly_rem(temp_poly1, temp_poly2, CM_Parameter, ctx);
        // printf("print(  ( "); fq_nmod_poly_print_pretty(temp_poly2, "x", ctx); printf(")%%(");fq_nmod_poly_print_pretty(CM_Parameter, "x", ctx);printf(") == (");fq_nmod_poly_print_pretty(temp_poly1,"x", ctx); printf(") )\n");
    }
    // printf("print(  ("); fq_nmod_poly_print_pretty(Beta, "x", ctx); printf(")*(");fq_nmod_poly_print_pretty(Beta,"x",ctx);printf(") == "); fq_nmod_poly_print_pretty(temp_poly2,"x",ctx); printf(")");



    // printf("\n MAT\n"); fq_nmod_mat_print_pretty(Matrix, ctx);
    fq_nmod_mat_rref(Matrix, ctx);
    // printf("element = ");    fq_nmod_print_pretty(fq_nmod_mat_entry(Matrix,deg-1,deg-1),ctx); printf("\n");

#if 1
    /*        Full Rank가 아닌 경우 다시           **/
    if(fq_nmod_is_zero(fq_nmod_mat_entry(Matrix,deg-1,deg-1),ctx))
    {
        return 0;
        ClassicMC_Irr_Gen(ret_Irr, CM_Parameter, ctx);
    }
    // printf("\n MAT\n"); fq_nmod_mat_print_pretty(Matrix, ctx);
#endif

    for(int i = 0; i < deg; i++)
    {
        fq_nmod_poly_set_coeff(ret_Irr, i, fq_nmod_mat_entry(Matrix, i, deg), ctx);
    }
    fq_nmod_poly_set_coeff(ret_Irr, deg, coeff_1, ctx);
    // printf("\n irr = "); fq_nmod_poly_print_pretty(ret_Irr, "X", ctx);
    // printf("\n Is irr = %d",fq_nmod_poly_is_irreducible(ret_Irr, ctx)); 

    flint_randclear(rand_state);
    nmod_poly_clear(coeff_1);
    fq_nmod_poly_clear(temp_poly1, ctx);
    fq_nmod_poly_clear(temp_poly2, ctx);
    fq_nmod_poly_clear(Beta, ctx);
    fq_nmod_mat_clear(Matrix, ctx);
    fq_nmod_clear(temp, ctx);
    return 1;
}

void CM_Irr_Gen_Test(int Case, int times)
{
    printf("-------CM Irr Test Start--------\n");
    fq_nmod_poly_t CM_Parameter, ret_irr;
    nmod_poly_t Coeff_1, Coeff_z;
    fq_nmod_ctx_t ctx;
    nmod_poly_t irr;
    
    fmpz_t p;
    fmpz_init(p);
    fmpz_set_ui(p, 2);
    
    nmod_poly_init(Coeff_1, 2);    nmod_poly_init(Coeff_z, 2);

    /**  random   */
    flint_rand_s* rand_state;
    rand_state = flint_rand_alloc();
    flint_randinit(rand_state);
    flint_randseed(rand_state, rand(), rand());

    if(Case == 1)
        fq_nmod_ctx_init(ctx, p, 12, "z");   
    if(Case == 2 ||Case == 3 ||Case == 4 ||Case == 5 )
        fq_nmod_ctx_init(ctx, p, 13, "z");   
    fq_nmod_init(irr,  ctx);


    // fq_nmod_ctx_init(ctx, p, 14, "z");
    // fq_nmod_poly_init(irr, ctx);
    // printf("ctx = = "); fq_nmod_ctx_print(ctx);
    if(Case == 1)
    {
        nmod_poly_set_str(irr, "13 2 1 0 0 1 0 0 0 0 0 0 0 0 1 ");
        fq_nmod_ctx_init_modulus(ctx, irr, "z");
        // printf("\ncase 1 ctx : ");    fq_nmod_ctx_print(ctx);printf("\n");
    }
    else if(Case == 2 ||Case == 3 ||Case == 4 ||Case == 5 )
    {
        nmod_poly_set_str(irr, "14 2 1 1 0 1 1 0 0 0 0 0 0 0 0 1 ");
        fq_nmod_ctx_init_modulus(ctx, irr, "z");
        // printf("\ncase %d ctx : ",Case);    fq_nmod_ctx_print(ctx);printf("\n");
    }    
    nmod_poly_set_str(Coeff_1, "1 2 1");
    nmod_poly_set_str(Coeff_z,"2 2 0 1");
    // printf("\nnmod 1 = "); nmod_poly_print_pretty(Coeff_1,"z");
    // printf("\nnmod z = "); nmod_poly_print_pretty(Coeff_z,"z"); printf("\n");
    /**    Classic Mceliece에 있는 식으로 바꿀 것!   */
    fq_nmod_poly_init(CM_Parameter, ctx);
    if(Case == 1)
    {
        fq_nmod_poly_set_coeff(CM_Parameter, 0, Coeff_z, ctx);
        fq_nmod_poly_set_coeff(CM_Parameter, 1, Coeff_1, ctx);
        fq_nmod_poly_set_coeff(CM_Parameter, 3, Coeff_1, ctx);
        fq_nmod_poly_set_coeff(CM_Parameter, 64, Coeff_1, ctx);
    }
    else if(Case == 2)
    {
        fq_nmod_poly_set_coeff(CM_Parameter, 0, Coeff_1, ctx);
        fq_nmod_poly_set_coeff(CM_Parameter, 6, Coeff_1, ctx);
        fq_nmod_poly_set_coeff(CM_Parameter, 9, Coeff_1, ctx);
        fq_nmod_poly_set_coeff(CM_Parameter, 10, Coeff_1, ctx);
        fq_nmod_poly_set_coeff(CM_Parameter, 96, Coeff_1, ctx);
    }
    else if(Case == 3)
    {
        fq_nmod_poly_set_coeff(CM_Parameter, 0, Coeff_1, ctx);
        fq_nmod_poly_set_coeff(CM_Parameter, 1, Coeff_1, ctx);
        fq_nmod_poly_set_coeff(CM_Parameter, 2, Coeff_1, ctx);
        fq_nmod_poly_set_coeff(CM_Parameter, 7, Coeff_1, ctx);
        fq_nmod_poly_set_coeff(CM_Parameter, 128, Coeff_1, ctx);
    }
    else if(Case == 4)
    {
        fq_nmod_poly_set_coeff(CM_Parameter, 0, Coeff_1, ctx);
        fq_nmod_poly_set_coeff(CM_Parameter, 8, Coeff_1, ctx);
        fq_nmod_poly_set_coeff(CM_Parameter, 119, Coeff_1, ctx);
    }
    else if(Case == 5)
    {
        fq_nmod_poly_set_coeff(CM_Parameter, 0, Coeff_1, ctx);
        fq_nmod_poly_set_coeff(CM_Parameter, 1, Coeff_1, ctx);
        fq_nmod_poly_set_coeff(CM_Parameter, 2, Coeff_1, ctx);
        fq_nmod_poly_set_coeff(CM_Parameter, 7, Coeff_1, ctx);
        fq_nmod_poly_set_coeff(CM_Parameter, 128, Coeff_1, ctx);
    }
    // fq_nmod_poly_randtest_irreducible(CM_Parameter, rand_state, 13, ctx);
    // printf("CM Parameter = ");    fq_nmod_poly_print_pretty(CM_Parameter,"x", ctx);
    fq_nmod_poly_init(ret_irr, ctx);
    /*********      시간측정        ******/    
    clock_t start, end;
    start = clock();
    int cnt_i;
    for(cnt_i = 0; cnt_i < times; cnt_i++)
    {
        ClassicMC_Irr_Gen(ret_irr, CM_Parameter, ctx);

        // printf("\n Is irr = %d",fq_nmod_poly_is_irreducible(ret_irr, ctx)); 
        // if(!fq_nmod_poly_is_irreducible(ret_irr, ctx))
        //     fq_nmod_poly_print_pretty(ret_irr,"X",ctx);
    }
    end = clock();
    printf("\n%d회 측정결과 1회 평균 시간 : %lu 미리초\n",cnt_i, ((end-start)*1000/(CLOCKS_PER_SEC))/cnt_i);
    fq_nmod_poly_clear(CM_Parameter, ctx);
    fq_nmod_poly_clear(ret_irr, ctx);
    // nmod_poly_clear(Coeff_1);
    // nmod_poly_clear(Coeff_z);
    nmod_poly_clear(irr);
    flint_randclear(rand_state);
    fmpz_clear(p);
    fq_nmod_ctx_clear(ctx);
}


void fq_nmod_poly_mul_timetest()
{
    fq_nmod_poly_t poly1, poly2, mod_poly;
    fmpz_t p;
    fq_nmod_ctx_t ctx;
    fmpz_init(p);
    fmpz_set_ui(p, 2);
    fq_nmod_ctx_init(ctx, p, 13, "z");
    fq_nmod_poly_init(poly1, ctx);
    fq_nmod_poly_init(poly2, ctx);
    fq_nmod_poly_init(mod_poly, ctx);

    /**  random   */
    flint_rand_s* rand_state;
    rand_state = flint_rand_alloc();
    flint_randinit(rand_state);
    flint_randseed(rand_state, rand(), rand());
    // printf("ctx : ");   fq_nmod_ctx_print(ctx);
    fq_nmod_poly_randtest(poly1, rand_state, 65,ctx);
    fq_nmod_poly_randtest(poly2, rand_state, 65,ctx);
    fq_nmod_poly_randtest_irreducible(mod_poly, rand_state, 65, ctx);

    fq_nmod_poly_t temp;
        /*********      시간측정        ******/    
    clock_t start, end;
    start = clock();
    int cnt_i = 0;
    fmpz_set_si(p, 2<<30);
    fmpz_mul(p,p,p);    fmpz_mul(p,p,p);    fmpz_mul(p,p,p);    fmpz_mul(p,p,p);
    fmpz_print(p    );
    for(cnt_i = 0; cnt_i < 3; cnt_i++)
    {
        // fq_nmod_poly_randtest_irreducible(poly1,rand_state, 65, ctx);
        // ben_or(poly1, ctx);
        // long asd[1] = {2};
        // rabin(poly1, ctx, asd, 1);

        // fq_nmod_poly_randtest(poly2, rand_state, 65,ctx);
        // fq_nmod_poly_powmod_fmpz_binexp(poly1, poly2, p, mod_poly, ctx);
        // fq_nmod_poly_gcd_euclidean(poly1, poly2, mod_poly, ctx);
        // fq_nmod_poly_mulmod(poly1,poly1,poly2,mod_poly,ctx);
        fq_nmod_poly_init(temp, ctx);
        fq_nmod_poly_clear(temp, ctx);
    }
    end = clock();
    printf("\n%d회 측정결과 1회 평균 시간 : %lu 미리초\n",cnt_i, ((end-start)*1000/(CLOCKS_PER_SEC))/cnt_i);

    // fq_nmod_poly_print_pretty(poly1,"X",ctx);   printf("\n\n");
    // fq_nmod_poly_print_pretty(poly2,"X",ctx);   printf("\n\n");
    // fq_nmod_poly_print_pretty(mod_poly,"X",ctx);printf("\n\n");

    flint_rand_free(rand_state);
    fmpz_clear(p);
    fq_nmod_poly_clear(poly1, ctx);
    fq_nmod_poly_clear(poly2, ctx);
    fq_nmod_poly_clear(mod_poly, ctx);
    fq_nmod_ctx_clear(ctx);
}

void cm_irr_proba_Test(int deg_t, int times)
{
    fq_nmod_poly_t rand_poly, irr_poly;
    fq_nmod_ctx_t ctx;
    fmpz_t p;
    fmpz_init(p);
    fmpz_set_ui(p, 2);
    fq_nmod_ctx_init(ctx, p, 13, "z");
    fq_nmod_poly_init(rand_poly, ctx);
    fq_nmod_poly_init(irr_poly, ctx);

    // nmod_poly_t irr;
    // nmod_poly_init(irr, ctx);
    // nmod_poly_set_str(irr, "13 2 1 0 0 1 0 0 0 0 0 0 0 0 1 ");
    // fq_nmod_ctx_init_modulus(ctx, irr, "z");

    flint_rand_s* rand_state;
    rand_state = flint_rand_alloc();
    flint_randinit(rand_state);
    flint_randseed(rand_state, rand(), rand());
    fq_nmod_poly_randtest_irreducible(irr_poly, rand_state, deg_t, ctx);

    int true_check, false_check ;   true_check = 0; false_check = 0;

    for(int i = 0; i < times; i++)
    {
        if(ClassicMC_Irr_Gen(rand_poly, irr_poly, ctx))
            true_check++;
        else
            false_check++;
    }
    printf("성공 = %d\n",true_check);
    printf("실패 = %d\n",false_check);
    // printf("확률 = %d\n",(true_check/times));

    flint_randclear(rand_state);
    fq_nmod_poly_clear(rand_poly, ctx);
    fq_nmod_poly_clear(irr_poly, ctx);
    fmpz_clear(p);
    fq_nmod_ctx_clear(ctx);
}

int main()
{
    srand(time(NULL));

    // printf("-       Case 1 \n");
    // ben_or_test(1, 5);
    // rabin_test(1, 5);
    // CM_Irr_Gen_Test(1, 50);

    // printf("-       Case 2 \n");
    // ben_or_test(2, 50);
    // rabin_test(2, 50);
    // CM_Irr_Gen_Test(2, 50);

    // printf("-       Case 3 \n");
    // ben_or_test(3, 50);
    // rabin_test(3, 50);
    // CM_Irr_Gen_Test(3, 50);

    // printf("-       Case 4 \n");
    // ben_or_test(4, 50);
    // rabin_test(4, 50);
    // CM_Irr_Gen_Test(4, 50);

    // fq_nmod_poly_mul_timetest();


    for(int t = 2; t <=16; t++)
    {
        printf("%d차에 대한 결과\n", t);
        cm_irr_proba_Test(t, 100000);        
    }


    return 0;
}