print(getwd())
sourceCpp("fibonacci.cpp")

# sourceCpp(code='
#   #include <Rcpp.h>
#  
#   // [[Rcpp::export]]
#   int fibonacci(const int x) {
#     if (x == 0) return(0);
#     if (x == 1) return(1);
#     return (fibonacci(x - 1)) + fibonacci(x - 2);
#   }'
# )

myFibo = fibonacci(7)

print(myFibo)