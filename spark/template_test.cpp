#include <string>
#include <iostream>

using namespace std;

template <class T>
class MyClass
{
  public:
    std::wstring msg = L"hey";
    MyClass(){};
};

template <class S = float, class T = shared_ptr<MyClass<S>>>
class MyClass2
{
  public:
    MyClass2(T m, int i = 0);
};
template <class S, class T>
MyClass2<S, T>::MyClass2(T m, int i) { std::wcout << m->msg << std::endl; }

int main(int argc, char **argv)
{
    shared_ptr<MyClass<float>> mc(new MyClass<float>);
    wcout << mc->msg << endl;
    MyClass2<> mc2(mc, 1);
    return 0;
}