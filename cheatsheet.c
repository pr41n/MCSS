/*  --------VARIABLES-------- (32 bits)
TIPO                Bits  Range
char                8     -128 a 127
unsigned char       8     0 a 255
signed char         8     -128 a 127
short               16     -32768 a 32767
int                 16     -32768 a 32767
unsigned int        16     0 a 65535
signed int          16     -32768 a 32767
short int           16     -32768 a 32767
unsigned short int  16     0 a 65535
signed short int    16     -32768 a 32767
long int            32     -2147483648 a 2147483647
signed long int     32     -2147483648 a 2147483647
unsigned long int   32     0 a 4294967295
long                32     -2147483648 a 2147483647
unsigned long       32     0 a 4294967295
float               32     3.4E-38 a 3.4E+38
double              64     1.7E-308 a 1.7E+308
long double         64     3.4E-4932 a 1.1E+4932
*/

/*  --------VARIABLES-------- (64 bits)
char                  1 bytes    8 bits                 from -128 to 127
unsigned char         1 bytes    8 bits                 from 0 to 255
signed char           1 bytes    8 bits                 from -128 to 127
short                 2 bytes   16 bits                 from -32768 to 32767
short int             2 bytes   16 bits                 from -32768 to 32767
unsigned short int    2 bytes   16 bits                 from 0 to 65535
signed short int      2 bytes   16 bits                 from -32768 to 32767
int                   4 bytes   32 bits                 from -2147483648 (-2.15e+09) to 2147483647 (2.15e+09)
unsigned int          4 bytes   32 bits                 from 0 to 4294967295 (4.29e+09)
signed int            4 bytes   32 bits                 from -2147483648 (-2.15e+09) to 2147483647 (2.15e+09)
long                  8 bytes   64 bits                 from -9223372036854775808 (-9.22e+18) to 9223372036854775807 (9.22e+18)
long int              8 bytes   64 bits                 from -9223372036854775808 (-9.22e+18) to 9223372036854775807 (9.22e+18)
unsigned long int     8 bytes   64 bits                 from 0 to 18446744073709551615 (1.84e+19)
signed long int       8 bytes   64 bits                 from -9223372036854775808 (-9.22e+18) to 9223372036854775807 (9.22e+18)
float                 4 bytes   32 bits                 from 1.175494351e-38 (1.18e-38) to 3.402823466e+38 (3.40e+38)
double                8 bytes   64 bits                 from 2.2250738585072014e-308 (2.23e-308) to 1.7976931348623158e+308 (1.79e+308)
long double          16 bytes   128 bits                from --------------------------------------------
*/

/* --------BASIC FORMAT SPECIFIERS--------
Format Specifier	   Type
%c	                 Character
%d	                 Signed integer
%e or %E	         Scientific notation of floats
%f	                 Float values
%g or %G	         Similar as %e or %E
%hi	                 Signed integer (short)
%hu	                 Unsigned Integer (short)
%i	                 Unsigned integer
%l or %ld or %li	 Long
%lf	                 Double
%Lf	                 Long double
%lu	                 Unsigned int or unsigned long
%lli or %lld	     Long long
%llu	             Unsigned long long
%o	                 Octal representation
%p	                 Pointer
%s	                 String
%u	                 Unsigned int
%x or %X	         Hexadecimal representation
%n	                 Prints nothing
%%	                 Prints % character
%-20.5s              left align and print string up to 5 character
*/

/* --------FUNCTION RETURN TYPES--------
integral datatype                                       (_Bool[C99]/char/short/int/long/long long and signed/unsigned variants)
floating-point datatype                                 (float/double/long double and _Complex variants[C99])
structure and union values                              (values of type struct ... or union ...)
enumeration values                                      (values of type enum ...)
pointers to any of the above and to any pointer
function pointers
void
*/


#include <stdio.h>

#define LOWER 0             // doesnt use memory for storage (precompilation change)
#define NAME "nombre"
const double PI = 3.14;     // uses memory for storage


int main(){
/* basic variables */
    int integer1, integer2;    // 2 bytes or 1 word
    float decimal = 1.9;    // 4 bytes, 1 word
    double double_precission_floating;  // 8 bytes, 2 words
    char a;     // 1 byte

    double_precission_floating = 1.2;
    a = 'k';

/* arrays */
    float arr1[5] = {3.2, 6.55};     // missing values set to 0
    arr1[2] = arr1[1];
    int arr2[2][3] = {
        {3, 2, 6},
        {4, 5, 20}
        };

/* memory size of variables */
    printf("int: %ld \n", sizeof(int));     // ld for long decimalleft align and print string up to 5 character
    printf("float: %ld \n", sizeof(float));
    printf("double: %ld \n", sizeof(double));
    printf("char: %ld \n", sizeof(char));

    printf("short: %ld \n", sizeof(short));
    printf("long: %ld \n", sizeof(long));
    printf("long long: %ld \n", sizeof(long long));

/* input & output */
    char b = getchar();
    char str[100];
    gets(str);          // EXTREMELY UNSECURE

    printf("CHAR: %-c \nSTRING: %-5s", b, str);
    putchar(b);
    puts(str);

    int x; float y; char text[20];

    scanf("%d %f %s", &x, &y, text);        // pointers needed
    printf("%d %f %s", x, y, text);

    scanf("%2d %d %*f %5s", &x, &y, text);        // * skips the input field
    printf("%d %f %s", x, y, text);

    printf("%3d%9.2f \n", 1234, 123.51546);
    printf("%c\t%s\bTEXT \n", 'a', "text");
    printf("%o\t0x%x \n", 16, 0x6a);

/* arithmetic operators */
    printf("%d %d %f %f", 5/9, 10%3, 5.0/9.0, (float) 5/9);
    x = 2; x+=7-2; x-=3*2; x*=6; x/=3; x%=3;        // x=2

    int x1, x2, x3; x1=3; x2=4;
    printf("%d %d %d \n", x1, x2, x3);
    x3 = x1--; printf("%d %d %d \n", x1, x2, x3);
    x3 = ++x2; printf("%d %d %d \n", x1, x2, x3);

/* logical operators */
    if ( !('a'=='b' && 1>3 || 2!=3) ) printf("true1");
    if ( !( 'a'=='b' && (1>3 || 2!=3) ) ) printf("true2");
    if (3) printf("true3");
    if (1^0) printf("true4 \n");   // XOR operator

/* conditionals & loops */
    x = 4; printf("x=%d \n", x);

    y = (x < 3) ? 5.1 : 7.2;        // y = 5.1 if x<3 else 7.2

    if (x>3) printf("x>3 \n");
    else if (x==3) printf("x==3 \n");
    else printf("x<3 \n");

    switch (x+2) {
        case 1:
            printf("case1 \n");
            break;
        case 2:
        case 4:
            printf("case2");
        default:
            printf("default \n");
    }

    x = 5;
    while (x >= 0)
        printf("%d", x--);
    printf("\n");

    do {
        printf("%d", x);
        x++;
    } while (x < 7);
    printf("\n");

    int i = 0;
    for (; i <= 3; i++) {           // initvalue, condition and/or increment skipable
        printf("Cuenta1\n");
    }
    int j, k;
    for (j=0, k=3; j<k; j++, k--) {     // initvalue, condition and/or increment can be expanded
        printf("Cuenta2\n");
    }

    for (i=0; i<1000; i++) {            // goto example
        for (j=0; j<1000; j++) {
            for (k=0; k<1000; k++) {
                if (i==0 && j==0 && k==2) goto found;
            }
        }
    }
    found:
        printf("Found!");

/* pointers */


}

