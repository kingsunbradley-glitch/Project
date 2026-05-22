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
			
			  a = a + b;  // a 变为和
	b = a - b;  // b 变为原 a
	a = a - b;  // a 变为原 b

	}
    return 0;
}