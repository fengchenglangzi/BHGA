#include "BasinHopping.h"

class BH_Manage
{
public:
	BH_Manage(void);
	~BH_Manage(void);
	void checkin();
	void process();
	bool EndingCondition(BasinHopping &);
	bool EndingCondition1(BasinHopping &,int FeNum);
	void set_objbest(double);
	double get_objbest();
	void disand(int N,int FeNum,int F);
private:
	
	int popNum;
	int FeNum;
	int iterations;
	int repeat;
	int current;
	int atomNum;
	double objbest;

};