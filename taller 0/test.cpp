#include <stdio.h>
int main(){
    float i=10.2;
    printf("value of i -> %.1f \n", i);
    float *P=&i;
    float *H=P;
    *H=23.0;
    printf("value of i -> %.1f pointer %.1f \n", i, *P);

    return 0;
}