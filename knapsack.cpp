 #include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/timeb.h>
#include <string.h>
#include <string> 
#include "knapsack.h" 

#define TIMEB struct timeb
#define FTIME ftime
#define UDT_TIME long
#define MAX_SIZE_TO_PRINT 10

UDT_TIME gRefTime = 0;


UDT_TIME GetMilliSecondTime(TIMEB timeBuf);
void SetTime(void);
UDT_TIME GetTime(void);

int main(int argc, char* argv[])
{
	UDT_TIME time, BFTime;
    float speedup;
	int itemCnt;
	KnapsackInstance* inst; //a Knapsack instance object
	KnapsackDPSolver DPSolver; //dynamic programming solver
	KnapsackBFSolver BFSolver; //brute-force solver
	KnapsackBTSolver BTSolver; //backtracking solver
	KnapsackBBSolver BBSolver1(UB1); //branch-and-bound solver with UB1
	KnapsackBBSolver BBSolver2(UB2); //branch-and-bound solver with UB2
	KnapsackBBSolver BBSolver3(UB3); //branch-and-bound solver with UB3
	KnapsackSolution *DPSoln, *BFSoln, *BTSoln, *BBSoln1, *BBSoln2, *BBSoln3;

    if(argc != 2) {
        printf("Invalid Number of command-line arguments\n");
	    exit(1);
    }
    itemCnt = atoi(argv[1]);
    if(itemCnt < 1) {
        printf("Invalid number of items\n");
        exit(1);
    }

	inst = new KnapsackInstance(itemCnt);
	DPSoln = new KnapsackSolution(inst);
	BFSoln = new KnapsackSolution(inst);
	BTSoln = new KnapsackSolution(inst);
	BBSoln1 = new KnapsackSolution(inst);
	BBSoln2 = new KnapsackSolution(inst);
	BBSoln3 = new KnapsackSolution(inst);
	
	inst->Generate();
	inst->Print();

	SetTime();
	DPSolver.Solve(inst,DPSoln);
	time = GetTime();
	printf("\n\nSolved using dynamic programming (DP) in %ld ms. Optimal value = %d",time, DPSoln->GetValue());
	if(itemCnt <= MAX_SIZE_TO_PRINT)
		DPSoln->Print("Dynamic Programming Solution");

	SetTime();
	BFSolver.Solve(inst,BFSoln);
	BFTime = time = GetTime();
	printf("\n\nSolved using brute-force enumeration (BF) in %ld ms. Optimal value = %d",time, BFSoln->GetValue());
	if(itemCnt <= MAX_SIZE_TO_PRINT)
		BFSoln->Print("Brute-Force Solution");
	if(*DPSoln == *BFSoln)
		printf("\nSUCCESS: DP and BF solutions match");
	else
		printf("\nERROR: DP and BF solutions mismatch");

	SetTime();
	BTSolver.Solve(inst,BTSoln);
	time = GetTime();
	printf("\n\nSolved using backtracking (BT) in %ld ms. Optimal value = %d",time, BTSoln->GetValue());
	if(itemCnt <= MAX_SIZE_TO_PRINT)
		BTSoln->Print("Backtracking Solution");
	if(*BFSoln == *BTSoln)
		printf("\nSUCCESS: BF and BT solutions match");
	else
		printf("\nERROR: BF and BT solutions mismatch");
    speedup = time==0? 0 : 100.0*(BFTime-time)/(float)BFTime;
	printf("\nSpeedup of BT relative to BF is %.2f%c",speedup,'%');

	SetTime();
	BBSolver1.Solve(inst,BBSoln1);
	time = GetTime();
	printf("\n\nSolved using branch-and-bound (BB) with UB1 in %ld ms. Optimal value = %d",time, BBSoln1->GetValue());
	if(itemCnt <= MAX_SIZE_TO_PRINT)
		BBSoln1->Print("BB-UB1 Solution");
	if(*BFSoln == *BBSoln1)
		printf("\nSUCCESS: BF and BB-UB1 solutions match");
	else
		printf("\nERROR: BF and BB-UB1 solutions mismatch");
    speedup = time==0? 0 : 100.0*(BFTime-time)/(float)BFTime;
	printf("\nSpeedup of BB-UB1 relative to BF is %.2f%c",speedup,'%');

	SetTime();
	BBSolver2.Solve(inst,BBSoln2);
	time = GetTime();
	printf("\n\nSolved using branch-and-bound (BB) with UB2 in %ld ms. Optimal value = %d",time, BBSoln2->GetValue());
	if(itemCnt <= MAX_SIZE_TO_PRINT)
		BBSoln2->Print("BB-UB2 Solution");
	if(*BFSoln == *BBSoln2)
		printf("\nSUCCESS: BF and BB-UB2 solutions match");
	else
		printf("\nERROR: BF and BB-UB2 solutions mismatch");
    speedup = time==0? 0 : 100.0*(BFTime-time)/(float)BFTime;
	printf("\nSpeedup of BB-UB2 relative to BF is %.2f%c",speedup,'%');

	SetTime();
	BBSolver3.Solve(inst,BBSoln3);
	time = GetTime();
	printf("\n\nSolved using branch-and-bound (BB) with UB3 in %ld ms. Optimal value = %d",time, BBSoln3->GetValue());
	if(itemCnt <= MAX_SIZE_TO_PRINT)
		BBSoln3->Print("BB-UB3 Solution");
	if(*BFSoln == *BBSoln3)
		printf("\nSUCCESS: BF and BB-UB3 solutions match");
	else
		printf("\nERROR: BF and BB-UB3 solutions mismatch");
    speedup = time==0? 0 : 100.0*(BFTime-time)/(float)BFTime;
	printf("\nSpeedup of BB-UB3 relative to BF is %.2f%c",speedup,'%');
	
	delete inst;
	delete DPSoln;
	delete BFSoln;
	delete BTSoln;
	delete BBSoln1;
	delete BBSoln2;
	delete BBSoln3;

	printf("\n\nProgram Completed Successfully\n");

	return 0;
}
/********************************************************************/

KnapsackInstance::KnapsackInstance(int itemCnt_)
{
	itemCnt = itemCnt_;

	weights = new int[itemCnt+1];
	values = new int[itemCnt+1];
	cap = 0;
}
/********************************************************************/
	
KnapsackInstance::~KnapsackInstance()
{
	delete [] weights;
	delete [] values;
}
/********************************************************************/

void KnapsackInstance::Generate()
{
	int i, wghtSum;  

	weights[0] = 0;
	values[0] = 0;
wghtSum = 0;

	for(i=1; i<= itemCnt; i++)
    {
        weights[i] = rand()%100 + 1;
        values[i] = weights[i] + 10;
        wghtSum += weights[i]; 
	}
cap = wghtSum/2;

}
/********************************************************************/

int KnapsackInstance::GetItemCnt()
{
	return itemCnt;
}
/********************************************************************/

int KnapsackInstance::GetItemWeight(int itemNum)
{
	return weights[itemNum];
}
/********************************************************************/

int KnapsackInstance::GetItemValue(int itemNum)
{
	return values[itemNum];
}
/********************************************************************/

int KnapsackInstance::GetCapacity()
{
	return cap;
}
/********************************************************************/

void KnapsackInstance::Sort(){
	
	int ratio[GetItemCnt() + 1];
	int i;
	for(i = 1;i <= GetItemCnt(); i++){
		ratio[i] = values[i] / weights[i];
	}
	
	int j,key, temp1,temp2;
	for(j = 2;j <= GetItemCnt() ; j++){
		key = ratio[j];
		temp1 = weights[j];
		temp2 = values[j];
		i=j-1;
		while(i >= 1 && ratio[i] < key){
			ratio[i+1] = ratio[i];
			weights[i+1] = weights[i];
			values[i+1] = values[i];
			i=i-1;
		}
		ratio[i+1] = key;
		weights[i+1] = temp1;
		values[i+1] = temp2;
	}
	
	
	
}
/*******************************************************************/
void KnapsackInstance::Print()
{
	int i;

	printf("Number of items = %d, Capacity = %d\n",itemCnt, cap);
	printf("Weights: ");
	for(i=1; i<=itemCnt; i++)
	{
		printf("%d ",weights[i]);
	}
	printf("\nValues: ");
	for(i=1; i<=itemCnt; i++)
	{
		printf("%d ",values[i]);
	}
	printf("\n");
}
/*****************************************************************************/

KnapsackSolution::KnapsackSolution(KnapsackInstance* inst_)
{
	int i, itemCnt = inst_->GetItemCnt();

	inst = inst_;
	isTaken = new bool[itemCnt+1];
	value = INVALID_VALUE;

	for(i=1; i<=itemCnt; i++)
	{
		isTaken[i] = false;
	}
}
/********************************************************************/

KnapsackSolution::~KnapsackSolution()
{
	delete [] isTaken;
}
/********************************************************************/

bool KnapsackSolution::operator == (KnapsackSolution& otherSoln)
{
	return value == otherSoln.value;
}
/********************************************************************/

void KnapsackSolution::TakeItem(int itemNum)
{
	isTaken[itemNum] = true;
	CurWeight += inst->GetItemWeight(itemNum);
	
	//  Branch and Bound
	
	TakenValue += inst->GetItemValue(itemNum);
	RemainingCap -= inst->GetItemWeight(itemNum);
	
}
/********************************************************************/

void KnapsackSolution::UndoTake(int itemNum){
	isTaken[itemNum] = false;
	CurWeight -= inst->GetItemWeight(itemNum);
	RemainingCap += inst->GetItemWeight(itemNum);
	TakenValue -= inst->GetItemValue(itemNum);
}

/*******************************************************************/

	
void KnapsackSolution::DontTakeItem(int itemNum)
{
	isTaken[itemNum] = false;
	
	//  Branch and Bound
	
	UntakenValue += inst->GetItemValue(itemNum);	
	
}
/********************************************************************/

void KnapsackSolution::UndoDontTake(int itemNum)
{
	isTaken[itemNum] = false;
	
	//  Branch and Bound
	
	UntakenValue -= inst->GetItemValue(itemNum);	
	
}
/********************************************************************/


/********************************************************************/


int KnapsackSolution::ComputeValue()
{
	int i, itemCnt = inst->GetItemCnt(), weight = 0;

	value = 0;
	for(i=1; i<=itemCnt; i++)
	{
		if(isTaken[i] == true)
		{
			weight += inst->GetItemWeight(i);
			if(weight > inst->GetCapacity())
			{
				value = INVALID_VALUE;
				break;
			}
			value += inst->GetItemValue(i);
		}
	}
	return value;
}
/********************************************************************/

int KnapsackSolution::GetValue()
{
	return value;
}
/********************************************************************/

void KnapsackSolution::Copy(KnapsackSolution* otherSoln)
{
	int i, itemCnt = inst->GetItemCnt();

	for(i=1; i<=itemCnt; i++)
	{
		isTaken[i] = otherSoln->isTaken[i];
	}
	value = otherSoln->value;
}
/********************************************************************/

void KnapsackSolution::Print(std::string title)
{
	int i, itemCnt = inst->GetItemCnt();

	printf("\n%s: ",title.c_str());
	for(i=1; i<=itemCnt; i++)
	{
		if(isTaken[i] == true)
			printf("%d ",i);
	}
	printf("\nValue = %d\n",value);
	
}
/*****************************************************************************/

KnapsackBFSolver::KnapsackBFSolver()
{
	crntSoln = NULL;
}
/********************************************************************/

KnapsackBFSolver::~KnapsackBFSolver()
{
	
}
/********************************************************************/

void KnapsackBFSolver::Solve(KnapsackInstance* inst_, KnapsackSolution* soln_)
{
	inst = inst_;	
	bestSoln = soln_;
	crntSoln = new KnapsackSolution(inst);
	FindSolns(1);
}
/********************************************************************/

void KnapsackBFSolver::FindSolns(int itemNum)
{
	int itemCnt = inst->GetItemCnt();

	if(itemNum == itemCnt + 1)
	{
		CheckCrntSoln();
		return;
	}
	
	crntSoln->DontTakeItem(itemNum);
	FindSolns(itemNum+1);
	crntSoln->TakeItem(itemNum);
	FindSolns(itemNum+1);
}
/********************************************************************/

void KnapsackBFSolver::CheckCrntSoln()
{
	int crntVal = crntSoln->ComputeValue();

#ifdef KNAPSACK_DEBUG
	printf("\nChecking solution ");
	crntSoln->Print(" ");
#endif

	if(crntVal == INVALID_VALUE)
		return;

	if(bestSoln->GetValue() == INVALID_VALUE) //The first solution is initially the best solution
		bestSoln->Copy(crntSoln);
	else
	{
		if(crntVal > bestSoln->GetValue())
			bestSoln->Copy(crntSoln);
	}
}
/********************************************************************/

// Write code below to implement the DP solver, the backtracking (BT) solver
// and the Branch-and-Bound (BB) solver. 
// Note that the BT and BB solvers inherit from the brute-force solver.
// You may add any private data members or functions that you need.
// Your solve() function takes an object of class KnapsackInstance as input
// and produces an object of class KnapsackSolution as output. 
// See how the given KnapsackBFSolver::Solve() writes its result into the 
// KnapsackSolution object and make the solvers that you write do the same.

KnapsackDPSolver::KnapsackDPSolver()
{

}

KnapsackDPSolver::~KnapsackDPSolver()
{

}

void KnapsackDPSolver::Solve(KnapsackInstance* inst_, KnapsackSolution* soln_)
{
	inst = inst_;
	soln = soln_;
	
	int itemCnt = inst->GetItemCnt();
	int cap = inst->GetCapacity();
	
	
	int matrix[itemCnt+1][cap+1] ;
	int i,j;
	
	for(j = 0;j <= cap;j++)
		matrix[0][j] = 0;
	
	for(i = 1; i <= itemCnt;i++){
		matrix[i][0] = 0;
		for(j = 1;j <= cap;j++){
			if(j < inst->GetItemWeight(i)){
				matrix[i][j] = matrix[i-1][j];				
			}
			else if(matrix[i-1][j] > (inst->GetItemValue(i) + matrix[i-1][j-inst->GetItemWeight(i)])){
				matrix[i][j] = matrix[i-1][j];
				
			}
			else{
					matrix[i][j] = (inst->GetItemValue(i) + matrix[i-1][j-inst->GetItemWeight(i)]);
					
			}
		}	
	}
	
	int n = inst->GetItemCnt();
	int c = inst->GetCapacity();
	
	for(n = itemCnt; n >= 1; n--){
		if(matrix[n][c] != matrix[n-1][c]){
			soln->TakeItem(n);
			c = c - inst->GetItemWeight(n);
		}
	} 
	soln->ComputeValue();
	
	
	
			
}
/*****************************************************************************/

KnapsackBTSolver::KnapsackBTSolver()
{
	crntSoln = NULL;
}
/********************************************************************/

KnapsackBTSolver::~KnapsackBTSolver()
{

}
/********************************************************************/

void KnapsackBTSolver::Solve(KnapsackInstance* inst_, KnapsackSolution* soln_)
{
	inst = inst_;	
	bestSoln = soln_;
	crntSoln = new KnapsackSolution(inst);
	crntSoln->CurWeight = 0;
	FindSolns(1);
	
}


void KnapsackBTSolver::FindSolns(int itemNum)
{
	int itemCnt = inst->GetItemCnt();

	if(crntSoln->CurWeight > inst->GetCapacity()){
		return;
	}

	if(itemNum == itemCnt + 1)
	{
		CheckCrntSoln();
		return;
	}
	crntSoln->DontTakeItem(itemNum);
	FindSolns(itemNum+1);
	crntSoln->TakeItem(itemNum);
	FindSolns(itemNum+1);
	crntSoln->UndoTake(itemNum);
}

/********************************************************************************************************/



/*****************************************************************************/

KnapsackBBSolver::KnapsackBBSolver(enum UPPER_BOUND ub_)
{
	ub = ub_;
	

}
/********************************************************************/

KnapsackBBSolver::~KnapsackBBSolver()
{

}
/********************************************************************/

void KnapsackBBSolver::Solve(KnapsackInstance* inst_, KnapsackSolution* soln_)
{
	inst = inst_;	
	bestSoln = soln_;
	crntSoln = new KnapsackSolution(inst);
	
	if (ub == UB3){
		inst->Sort();
	}
	
	

	
	crntSoln->CurWeight = 0;
	
	crntSoln->TotalValue = 0 ;
	crntSoln->TakenValue = 0 ;
	crntSoln->UntakenValue = 0;
	crntSoln->UndecidedValue = 0; 
	crntSoln->RemainingCap = inst->GetCapacity();
	
	
	int itemCnt = inst->GetItemCnt();
	int i;


	// Calculate Total sum of values
	
	for(i = 1; i <= itemCnt; i++){
		crntSoln->TotalValue += inst->GetItemValue(i);
	}

	FindSolns(1);
	
}
/*****************************************************************************/

void KnapsackBBSolver::FindSolns(int itemNum)
{
	int itemCnt = inst->GetItemCnt();
	
	if(crntSoln->CurWeight > inst->GetCapacity()){
		return;
	}
	else if(ub == UB1){
		if((crntSoln->TotalValue - crntSoln->UntakenValue) <= bestSoln->GetValue())
			return;
	}
	
	else if(ub == UB2){
	
		crntSoln->UndecidedValue = 0;
		for(int i = itemNum; i <= itemCnt; i++ ){
	
			if(inst->GetItemWeight(i) <= crntSoln->RemainingCap){
			
				crntSoln->UndecidedValue += inst->GetItemValue(i);
			}
			
	}
		
	
		
		if((crntSoln->TakenValue + crntSoln->UndecidedValue) <= bestSoln->GetValue())
			return;
	}
	
	else if(ub == UB3){
		 crntSoln->FractionalRemaining = crntSoln->RemainingCap ;
		
		
		 crntSoln->UndecidedValue = 0;
		 
		 
	
		 for(int i = itemNum; i <= itemCnt; i++){
		 	if(inst->GetItemWeight(i) <= (crntSoln->FractionalRemaining)){

		 		crntSoln->UndecidedValue += inst->GetItemValue(i);
		 		crntSoln->FractionalRemaining -= inst->GetItemWeight(i);
			 }
			 
			 else {

			 	crntSoln->UndecidedValue += ((float) inst->GetItemValue(i) * ((crntSoln->FractionalRemaining) / (float) inst->GetItemWeight(i)));
			 	crntSoln->FractionalRemaining =0;
				break;
			 }
		 }

		 if((crntSoln->TakenValue + crntSoln->UndecidedValue) <= bestSoln->GetValue()){


		 	return;
		 }
			
	}
	

	if(itemNum == itemCnt + 1)
	{
		CheckCrntSoln();
		return;
	}
	crntSoln->DontTakeItem(itemNum);
	FindSolns(itemNum+1);
	crntSoln->UndoDontTake(itemNum);
	crntSoln->TakeItem(itemNum);
	FindSolns(itemNum+1);
	crntSoln->UndoTake(itemNum);
}

/***************************************************************************/

UDT_TIME GetCurrentTime(void)
{
	UDT_TIME crntTime=0;

	TIMEB timeBuf;
	FTIME(&timeBuf);
    crntTime = GetMilliSecondTime(timeBuf);

	return crntTime;
}
/********************************************************************/

void SetTime(void)
{
	gRefTime = GetCurrentTime();
}
/********************************************************************/

UDT_TIME GetTime(void)
{
	UDT_TIME crntTime = GetCurrentTime();

	return (crntTime - gRefTime);
}
/********************************************************************/

UDT_TIME GetMilliSecondTime(TIMEB timeBuf)
{
	UDT_TIME mliScndTime;

	mliScndTime = timeBuf.time;
	mliScndTime *= 1000;
	mliScndTime += timeBuf.millitm;
	return mliScndTime;
}
/*****************************************************************************/
