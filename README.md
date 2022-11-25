## Irreducible-Polynomial

#Classic McEliece에서 사용하는 최소다항식 생성 알고리듬과, Ben-or, Rabin 기약다항식 생성 알고리듬
/** 컴파일 명령어
* gcc -Wall -O2 *.c -L/opt/homebrew/Cellar/mpfr/4.1.0/lib -L/opt/homebrew/Cellar/gmp/6.2.1_1/lib -I/opt/homebrew/Cellar/gmp/6.2.1_1/include -I/opt/homebrew/Cellar/mpfr/4.1.0/include -I/usr/local/include/flint/ -lflint -lmpfr -lgmp -lpthread
*/


# 구현 환경

| 하드웨어  | MacBook Air. Apple M2. 8GB RAM  | 
|:—————————:|:—————————:|
| **컴파일러** | **gcc 13.1.6(-O2)** |
| **정수 연산 라이브러리** | **FLINT 2.9.0** |
