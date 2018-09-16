/* exp example */
#include <stdio.h>      /* printf */
#include <math.h>       /* exp */
int potencia(int base, int exponente){
    if(exponente=0){
        return 1;
    }
    printf("entrando con exponen %i\n", exponente);
    for(int i=1; i<exponente; i++){
        base=base*base;
        printf("%i", base);
    }
    return base;
}

int main ()
{
    int bast=2, expo=14;
  printf("%i", potencia(bast,expo));
  return 0;
}