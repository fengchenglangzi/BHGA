/*#include"stdlib.h"
#include"string.h"

template <class AType>
class array
{
public:
	array(int size);
	~array()
	{delete []a;}
	AType &operator [](int i);
private:
	int length;
	AType *a;
};

template<class AType>
array<AType>::array(int size)
{
	register int i;
	length=size;
	a=new AType[size];
	if(!a)
	{
		cout<<"can't allocate array.\n";
		exit(1);
		
	}
	for(i=0;i<size;i++)
		a[i]=0;
}

template<class AType>
array<AType>::array(int size)
{
	if(i<0||i>length-1)
	{
		cout<<"\nIndx value of";
		cout<<i<<"is out of bounds.\n";
		exit(2);

		
	}
	return a[i];

}




void main2()
{
	//printf("hello.\n");

	system("pause");
}*/