#include <iostream>
#include <sstream>
#include <vector>
#include <list>
#include <algorithm>

using namespace std;

int main()
{

	list<int> * v = new list<int>({0 ,1, 2, 8 ,9 ,10,16,17 ,18 ,24,25 ,26});
	list<int> * w = new list<int>({0,1,2,3,8,9,10,11,16,17,18,19,24,25,26,27});
	list<int> x  = {1,5,3};

	// sort(v.begin(), v.end());
	// sort(w.begin(), w.end());
	// sort(x.begin(), x.end());
	v->sort();
	w->sort();
	x.sort();


	list<int> diff;
	list<int> inter;
	list<int> uni;
	
	set_difference(
		v->begin(), v->end(),
		w->begin(), w->end(),
		back_inserter(diff)
		)	;

	set_intersection(
		v->begin(), v->end(),
		w->begin(), w->end(),
		back_inserter(inter)

		);

	set_union(
		v->begin(), v->end(),
		x.begin(), x.end(),
		back_inserter(uni)

		);
	cout << "diff: ";
	for(auto d: diff)
		cout << d << " ";

	cout << endl;

	cout << "inter: ";
	for(auto d: inter)
		cout << d << " ";
	cout << endl;

	cout << "uni: ";
	for(auto d: uni)
		cout << d << " ";
	cout << endl;





}