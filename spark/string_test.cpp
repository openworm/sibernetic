#include <string>
#include <iostream>
#include <cmath>
#include <bitset>
using namespace std;
union C{
    uint64_t i;
    double d;
};
int main(int argc, char ** argw){
    C a;
    a.d = 1223.23213;
    double e = 1223.23213;
    bitset<64> b(a.i);
    uint64_t lfractionalPart = (*(reinterpret_cast<uint64_t *>(&e)) & ~0xF7000000);
    //b>>=12;//(b.size() >> 2);
    //b<<=12;
    //cout << int(a.d) << endl <<b.to_string() << " " << uint(b.to_ulong()) << " " << b.size()<< endl;
    cout << lfractionalPart << endl;
    return 0;
}