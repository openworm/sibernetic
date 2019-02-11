#include <iostream>

using namespace std;

int main(int argc, char ** argv){
	int i;
	float f;
	do{
		
		cout << "In f" << endl;
		cin >> f;
		if(cin.fail()){
			cout <<"Try again"<< endl; 
			cin.clear();
			cin.ignore(numeric_limits<streamsize>::max(),'\n');
		}
	}while(cin.fail());
	cout << f << endl;
	return 0;
}
