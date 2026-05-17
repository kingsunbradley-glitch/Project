#include <iostream>

int main()
{
    int a {}, b {}; // local variable, no linkage
	std::cin >> a >> b;
	if (a > b)
	 {
			   std::cout << "a("<<a<<" )is bigger than b(" << b << ")" << std::endl;

    }
	else 
	{
			   std::cout << "b("<<b<<" )is bigger than a(" << a << ")" << std::endl;
			   std::cout << "swapping..." << std::endl;
			
			   a = (b+a)/2; // a now holds the sum of a and b
			   std::cout << "a: " << a << ", b: " << b << std::endl;
			   b = (a-b);
			   std::cout << "a: " << a << ", b: " << b << std::endl;
			   a = a+b; // a now holds the original value of b
			   std::cout << "a: " << a << ", b: " << b << std::endl;
			   b = 2*a+b; // b now holds the original value of a
			   std::cout << "a: " << a << ", b: " << b << std::endl;

	}
    return 0;
}