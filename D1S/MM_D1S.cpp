#include <ilcplex/ilocplex.h>
#include <stdlib.h>   
#include <vector>
#include <time.h>
//#include <random>
#include <cmath>
using namespace std;

ILOSTLBEGIN

typedef IloArray<IloNumVarArray> VarArray2;
typedef IloArray<VarArray2>      VarArray3;
typedef IloArray<IloRangeArray>  EquationArray;
typedef IloArray<EquationArray>  EquationArray2;
typedef IloArray<IloNumArray2>   IloNumArray3;
typedef IloArray<IloExpr>        ObjArray;  
typedef IloArray<IloCplex>       CplexArray;            
typedef IloArray<IloModel>       ModelArray;

void merge(IloNumArray a, IloNumArray b, IloInt low, IloInt high)
{
    int h,i,j,k;
	int pivot = floor((high + low) / 2);
    h=low;
    i=low;
    j=pivot+1;
	
    while((h<=pivot)&&(j<=high)){
        if(a[h]<=a[j]){
            b[i]=a[h];
            h++;
        }
        else{
            b[i]=a[j];
            j++;
        }
        i++;
    }
    if(h>pivot){
        for(k=j; k<=high; k++){
            b[i]=a[k];
            i++;
        }
    }
    else{
        for(k=h; k<=pivot; k++){
            b[i]=a[k];
            i++;
        }
    }
    for(k=low; k<=high; k++){
		a[k]=b[k];
	}
}
void merge_sort(IloNumArray a, IloNumArray b, IloInt low, IloInt high)
{
    int pivot;
    if(low<high)
    {
        pivot=floor((low+high)/2);
        merge_sort(a,b,low,pivot);
        merge_sort(a,b, pivot+1,high);
        merge(a, b, low, high);
    }
}

void insertion_sort(IloInt size, IloNumArray item){
	IloInt i, j;
	IloBool ordered;
	IloNum temp;
	ordered=IloFalse;
	for (i=size; i>0; i--){
		if (!ordered){
			ordered=IloTrue;
			for (j=0; j<i-1; j++){
				if (item[j]> item[j+1]){
					temp=item[j];
					item[j]=item[j+1];
					item[j+1]=temp;
					ordered=IloFalse;
				}
			}
		}
	}
	
					
}


int main(int argc, char **argv) {

 IloEnv env1;

 try {

	  //read differenet values of lambda and alpha
	string line;
	ifstream riskpar ("config.txt");
	vector<double> vect;
	int kk;
	if (riskpar.is_open()){
		while ( getline (riskpar,line) ){
			kk++;
			stringstream ss(line);

			double i;

			while (ss >> i){
				vect.push_back(i);
				if (ss.peek() == ' ')
					ss.ignore();
			}
			
		}
	}
	riskpar.close();
	const int No_instance=kk;
	/*for(int i=0; i<vect.size(); i++){
		cout<<vect[i]<<endl;
	}*/
		
	double config[8][2];
	kk=0;
	for(int i=0; i<vect.size(); i++){
		if ((i+1) % 2!=0){
			//read lambda
			config[kk][0]=vect[i];		
		}
		else{
			//read alpha
			config[kk][1]=vect[i];
			kk++;
		}
	}
	
		 IloNum toler = 0.000001, disRate = 0.02;
		 IloInt i, j, k, t, w, w2, s;
		 const IloInt No_type = 10, No_node = 62, No_stage=4; 
	
		 // The data below read from input file
		 IloNum discRate = 0.04; //discount Rate
		 IloIntArray type(env1), No_year(env1);                      //The type of corresponding node
		 IloIntArray startNode(env1), endNode(env1);
		 IloInt No_link; 	//Number of year in each stage 
		 IloNumArray  cost(env1);   //cost on each link
		 IloNum   returnRate;
		 IloNumArray   loss(env1), capacity(env1), Storage0(env1),  storageUB(env1);
	///////////////// DATA FILE READING ////////////////////////////////
		   
		  const char* filename1  = "Data\\input1_800.dat";	      
		
		  ifstream file1(filename1);
		  if (!file1) {
				cerr << "ERROR: can not open file '" << filename1 << "' for reading" << endl;
				cerr << "usage:   " << argv[0] << " <file>" << endl;
				throw(-1);
			}
		   
          file1 >> No_year;
		  file1 >> type >> startNode >> endNode >> cost ;
	      file1 >> capacity >>  returnRate; 
		  file1 >>  loss >> Storage0>>storageUB;

		  No_link = startNode.getSize();
		  IloInt No_yearMax = 15 ; 
	 
 		  IloBool consistentData = ( No_year.getSize() == No_stage && endNode.getSize() == No_link 
			  && cost.getSize() == No_link && type.getSize() == No_node );
           
		  if (!consistentData) { 
 			  cerr << "ERROR: Inconsistent data1!" << endl;
			  throw(-1);
		  }
	   
		  file1.close();
		  cout<<"Data1 reading down"<<endl;
  ///////////////////Identify type group////////////////////////
		  //identify role of each node in the water network, and add node number to their corresponding role vector
		  //it acts to know which constraints are needed for each node

		  vector<int> userID;
		  vector<int> potUserID;
		  vector<int> rechargeID;
		  vector<int> balanceID;
		  vector<int> capacityID;

		  for(i = 0; i<No_node; i++){
			  if(type[i] == 4){
				  rechargeID.push_back(i);
			  }
			  else if(type[i] == 0 || type[i] == 2 || type[i] == 5 || type[i] == 6 || type[i] == 7){
				  balanceID.push_back(i); 
			  }
			  else if(type[i] == 1){
				  capacityID.push_back(i); 
			  }
			  else if(type[i] == 8 ){
				  userID.push_back(i);
				  potUserID.push_back(i);
			  }
		  }
          for(i = 0; i<No_node; i++){
			  if(type[i] == 9 ){
				  userID.push_back(i);
			  }
			  else if(type[i] == 4){
				  capacityID.push_back(i); 
			  }
		  }
		  const IloInt No_user = int(userID.size());
		  const IloInt No_potUser = int(potUserID.size());
		  const IloInt No_balance = int(balanceID.size());
          const IloInt No_capacity = int(capacityID.size());
          const IloInt No_recharge = int(rechargeID.size());
		  cout<<No_user<<"  "<<No_potUser<<endl;
		 
///////////////////read demeand, capacity, and intial storage//////////////////
          
		  IloNumArray3   demand(env1,No_stage);
		  IloNumArray2   population(env1,No_stage), TucsonPop(env1,No_stage);
		  IloNumArray    demUnit(env1,No_user);
		  for(s = 0; s< No_stage; s++){
			  demand[s] = IloNumArray2(env1,No_yearMax);
			  population[s] = IloNumArray(env1,No_yearMax);
			  TucsonPop[s] = IloNumArray(env1,No_yearMax);
			  for(t=0; t<No_yearMax; t++)
				  demand[s][t] = IloNumArray(env1,No_user);
		  }
		  const char* filename2  =  "Data\\input2.dat";
		  ifstream file2(filename2);
		  if (!file2) {
				cerr << "ERROR: can not open file '" << filename2 << "' for reading" << endl;
				cerr << "usage:   " << argv[0] << " <file>" << endl;
				throw(-1);
			}
           
		  file2 >> population[0]>>population[1]>>population[2]>>population[3]; 
		  file2 >> demUnit; 
		  file2 >> TucsonPop[0]>> TucsonPop[1]>> TucsonPop[2]>> TucsonPop[3];
		  
		  file2.close();
		  cout<<"Data part 2 read done"<<endl;
		 
		  const char* filename3  =  "Data\\DemPortion.txt";
		  ifstream file3(filename3);
		  if (!file3) {
				cerr << "ERROR: can not open file '" << filename3 << "' for reading" << endl;
				cerr << "usage:   " << argv[0] << " <file>" << endl;
				throw(-1);
			}
		  IloNumArray2 DembyNode(env1,No_year[0]+No_year[1]+No_year[2]+No_year[3]);
		  for(i= 0; i<No_year[0]+No_year[1]+No_year[2]+No_year[3];i++){
			  DembyNode[i]=IloNumArray(env1,No_potUser);
			  /*file3>>DembyNode[i];*/
		  }
		  for(i= 0; i<No_year[0]+No_year[1]+No_year[2]+No_year[3];i++){ 
			  file3>>DembyNode[i];
		  }
       
		  file3.close(); 
		  cout<<"Data part 3 read done"<<endl; 
          
		  for(t=0; t<No_year[0]; t++){
			  for(i=0; i<No_potUser;i++){
			     demand[0][t][i] = double(135*0.00112*population[0][t]*DembyNode[t][i]*0.8);
			  }
		  } 
		  for(t=0; t<No_year[1]; t++){
			  for(i=0; i<No_potUser;i++){
			     demand[1][t][i] = double(135*0.00112*population[1][t]*DembyNode[t + No_year[0]][i]*0.8);
			  }
		  }
		  for(t=0; t<No_year[2]; t++){
			  for(i=0; i<No_potUser;i++){
			     demand[2][t][i] = double(135*0.00112*population[2][t]*DembyNode[t+ No_year[0]+ No_year[1]][i]*0.8);
			  }
		  }
		  for(t=0; t<No_year[3]; t++){
			  for(i=0; i<No_potUser;i++){
			     demand[3][t][i] = double(135*0.00112*population[3][t]*DembyNode[t+ No_year[0]+ No_year[1]+ No_year[2]][i]*0.8);
			  }
		  }
		  for(s=0;s<No_stage;s++){
			  for(t=0; t<No_yearMax; t++){
				  for(i=No_potUser; i<No_user;i++){
					 demand[s][t][i] = demand[s][t][i-No_potUser]*0.25;
				  }
			  }
		  }
		   
		  IloNumArray2  CAPamt(env1,No_stage);
		  for(s = 0; s< No_stage; s++){
			  CAPamt[s] = IloNumArray(env1,No_yearMax);
			  for(t=0; t<No_yearMax; t++){
				  CAPamt[s][t] = 144000*population[s][t]/TucsonPop[s][t] ;
			  }
		  }
          population.end();TucsonPop.end(); demUnit.end();DembyNode.end();
		  
        
		  ///////////////Define scenario tree nodes 
		  IloInt No_scen =15; //Number of scenarios in each stage
		  IloNum prob = (double)1/No_scen;
		  
		  IloInt No_scenNode = 1 + No_scen+ No_scen*No_scen + No_scen*No_scen*No_scen; //Total number of nodes in the scenario tree 
		  IloInt No_scenElm = No_yearMax+2;
		  
		  
		 IloNumArray2 scenFac(env1);
		const char* filename4  =  "Data\\scenFac_15.txt";
		ifstream file4(filename4);
		file4>>scenFac;
		file4.close();

		  //determine stage of each node in the tree
		  IloIntArray stageMap(env1, No_scenNode);
		  stageMap[0] = 0;
          for(w=1;w<No_scenNode; w++){
			  if((w>=1) && (w<No_scen+1)){
				  stageMap[w] = 1;}
			  else if((w>=No_scen+1) && (w<1+No_scen+ No_scen*No_scen)){
				  stageMap[w] = 2;
			  }
			  else 
				  stageMap[w] = 3;
		  }

		  IloIntArray2 ancestor(env1, No_scenNode);
		  for(w=0; w<No_scenNode; w++){
			  ancestor[w] = IloIntArray(env1, No_scenNode);
			  for(w2=0; w2<No_scenNode; w2++){
				  ancestor[w][w2] = 0;
			  }
		  }
		  for(w=1; w<No_scen+1;w++)
		     ancestor[0][w] = 1;
		  for(i=1; i<No_scen+1;i++){
			  for(j= 1+No_scen+No_scen*(i-1); j<1+No_scen+No_scen*i;j++){
				  //i is ancestor for j, i is in stage 2
				  ancestor[i][j] = 1;
				  for(k = 1+No_scen+No_scen*(j-1);k<1+No_scen+No_scen*j;k++){
					  //j is ancestor for k, j is in stage 3
					  ancestor[j][k] = 1;
				  }
			  }
		  }
         
		  //denotes decendant of each node in the scenario tree by an array, whose elements are the child node number
		  IloIntArray2 decendant(env1, No_scenNode);
		  for(w=0;w<No_scenNode;w++){
			  decendant[w]= IloIntArray(env1,No_scen);
			  for(j=0; j<No_scen;j++){
				  decendant[w][j] = 0;
			  }
		  }
	      
		  //decendant of root
		  for(i=0;i<No_scen;i++){
			
			  decendant[0][i] = i+1;
		  }
		  
		  //th next two blocks could be combined
		  //decendant of nodes in stage 2
		  for(i=1;i<No_scen+1;i++){
			  for(j=0;j<No_scen;j++){
				  decendant[i][j] = 1+ No_scen+ No_scen*(i-1)+j;
			  }
		  }
		  //decendant of nodes in stage 3
		  for(i=No_scen+1; i<1+No_scen+No_scen*No_scen;i++){
			  for(j=0;j<No_scen;j++){
				  decendant[i][j] = 1+No_scen*i+j;
			  }
		  }

 
//		  IloInt No_scenElm = No_yearMax+2;
		  
		//  IloNumArray2  scenFac(env,No_scenNode); 
		//  for(w=0; w<No_scenNode; w++) {
		//	  scenFac[w] = IloNumArray(env,No_scenElm);
		//	  for(t = 0;t<No_scenElm;t++){
		//		  scenFac[w][t] = 1;
		//	  }
		//  }

		//  //generate \Xi_1 for every node in stage 2, end of for must be No_scen however it has been taken care later
		//  for(w=1; w<No_scenNode; w++) {
		//	  double randomNo = (double)rand()/(double)RAND_MAX;
		//	  if(randomNo <= 0.1){ 
		//		 scenFac[w][0] = 0.9;
		//	  }
		//	  else {
		//		 scenFac[w][0] = 1; 
		//	  } 
	 //
  //            //it is supposed to be between 0.9 and 1.1
		//	  for(t = 1;t<No_scenElm;t++){
		//		 double randomNo2 = 0.9 + 0.2*(double)rand()/(double)RAND_MAX;
		//		 scenFac[w][t] = randomNo2; 
		//	  }	 
		//  }
		//  
		//  //Change scenFac for second, third stage 
		//  for(w=No_scen+1;w<1+No_scen+No_scen*No_scen;w++){
		//	  double randomNo = (double)rand()/(double)RAND_MAX;
		//	  if(randomNo <= 0.25){ 
		//		 scenFac[w][0] = 0.9;
		//	  }
		//	  else {
		//		 scenFac[w][0] = 1; 
		//	  }
		//  }

		//  for(w=1+No_scen+No_scen*No_scen; w<1+No_scen+No_scen*No_scen+No_scen*No_scen*No_scen;w++){
		//	  double randomNo = (double)rand()/(double)RAND_MAX;
		//	  if(randomNo <= 0.35){ 
		//		 scenFac[w][0] = 0.9;
		//	  }
		//	  else {
		//		 scenFac[w][0] = 1; 
		//	  }
		//  }

		///// ******
  // 		const char* Scenfilename  = "scenFac_512.txt";
		//ofstream scenTree(Scenfilename);
		//scenTree<<'[';
		//for(w=0; w<No_scenNode; ++w) {
		//	scenTree<<'[';
		//	for(t = 0;t<No_scenElm;++t){
		//		 if(t+1< No_scenElm)
		//		     scenTree<<scenFac[w][t] <<',' ;
		//		 else
		//			 scenTree<<scenFac[w][t]; 
		//	 }

		//	if (w+1 < No_scenNode){
		//		scenTree<<']'<< ','<< endl;
		//	}
		//	else{
		//		scenTree<<']' << ']' << endl;
		//	}
		//}
		//scenTree.close();


		

///******
 for (kk=0; kk<No_instance; kk++){
		IloEnv env;
		 IloNum lambda=config[kk][0];
		 IloNum alpha=config[kk][1];
		 IloInt CVaR_No;
		 for (j=1; j<No_scen+1; j++){
			  if(j*prob>= alpha && (j-1)*prob< alpha)
				CVaR_No=j-1;
		  }		  
 

		 
 /////////////////// DECISION VARIABLES WITH NAMES  /////////////////////////////
		  
		IloInt nCut = 0;
		IloInt MaxCut =50;

		//Upper letters -- Variables
		
		char varName[100];

		//define x variables
		VarArray3 Q(env, No_link);
		for(i =0; i<No_link; i++){
			Q[i] = VarArray2(env, No_stage);
			for(s = 0; s<No_stage; s++){
				Q[i][s] = IloNumVarArray(env, No_yearMax, 0, IloInfinity, ILOFLOAT); 
				for(t=0;t<No_yearMax;t++){
					sprintf_s(varName, "Q_%d_%d_%d",(int) i, (int) s, (int) t);
					Q[i][s][t].setName(varName);
				}
			}
		}
		//define y variables
		VarArray3 Storage(env, No_recharge);
		for(i =0; i<No_recharge; i++){
			Storage[i] = VarArray2(env, No_stage);
			for(s = 0; s<No_stage; s++){
				Storage[i][s]= IloNumVarArray(env, No_yearMax, 0, IloInfinity, ILOFLOAT);
				for(t=0;t<No_yearMax;t++){
					sprintf_s(varName, "S_%d_%d_%d",(int) i, (int) s, (int)t);
					Storage[i][s][t].setName(varName);
				}
			}
		}

		VarArray2	Theta1(env,No_stage);
		VarArray2	Theta2(env,No_stage);
		for(s=0; s<No_stage; s++){
			Theta1[s] = IloNumVarArray(env, No_scen, 0, IloInfinity, ILOFLOAT);
			Theta2[s] = IloNumVarArray(env, No_scen, 0, IloInfinity, ILOFLOAT);
			for(j=0;j<No_scen;j++){
			   sprintf_s(varName, "Theta1_%d_%d",(int) s, (int) j);
			   Theta1[s][j].setName(varName);
			}
			for(j=0;j<No_scen;j++){
			   sprintf_s(varName, "Theta2_%d_%d",(int) s, (int) j);
			   Theta2[s][j].setName(varName);
			}
		}


		VarArray2 VaR(env, MaxCut);
        VarArray2 V(env, MaxCut);
		for(nCut=0; nCut<MaxCut;nCut++){
		   VaR[nCut] = IloNumVarArray(env,No_stage, 0, IloInfinity, ILOFLOAT);
		   V[nCut] = IloNumVarArray(env, No_scenNode, 0, IloInfinity, ILOFLOAT);
		   for(w=0;w<No_scenNode;w++){
			   sprintf_s(varName, "V_%d_%d",(int) nCut, (int) w);
			   V[nCut][w].setName(varName);
		   }
		} 

		//Lower case letters -- passed values for the variables 
		IloNumArray2 storage_hat(env, No_recharge);
		for(i =0; i<No_recharge; i++){
			storage_hat[i] = IloNumArray(env, No_scenNode); 
		} 

        //Dual values 
		IloNumArray3 pi_balance(env, No_yearMax);  //dual variable for nodal flow balance constraints
		IloNumArray3 pi_capacity(env, No_yearMax);  //dual variable for pump/surface water capacity constraints
		IloNumArray3 pi_demand(env, No_yearMax);  //dual variable for demand satisfication constraints
		IloNumArray3 pi_return(env, No_yearMax);  //dual variable for potable return constraints
		IloNumArray3 pi_safe(env, No_yearMax);  //dual variable for safe yield
		IloNumArray2 pi_storage(env, No_scenNode);  //dual variable for storage balance constraints
		IloNumArray2 pi_flowBound(env, No_scenNode);  //dual variable for RF outflow bound constraints 

		for(t=0;t<No_yearMax;t++){
			pi_balance[t] = IloNumArray2(env, No_scenNode);
			pi_capacity[t] = IloNumArray2(env, No_scenNode);
			pi_demand[t] = IloNumArray2(env, No_scenNode);
			pi_return[t] = IloNumArray2(env, No_scenNode);
			pi_safe[t] = IloNumArray2(env, No_scenNode);
			for(w=0;w<No_scenNode;w++){
				pi_balance[t][w] = IloNumArray(env, No_balance);
				pi_capacity[t][w] = IloNumArray(env, No_capacity);
				pi_demand[t][w] = IloNumArray(env, No_user);
				pi_return[t][w] = IloNumArray(env, No_potUser);
				pi_safe[t][w] = IloNumArray(env, No_recharge);
				pi_storage[w] = IloNumArray(env, No_recharge);
				pi_flowBound[w] = IloNumArray(env, No_recharge); 

				for(i=0; i<No_balance; i++) 
 	   				pi_balance[t][w][i] = 0; 
				
				for(i=0; i<No_capacity; i++) 
 	   				pi_capacity[t][w][i] = 0;
				
				for(i=0; i<No_user; i++) 
 	   				pi_demand[t][w][i] = 0;
				
				for(i=0; i<No_potUser; i++) 
 	   				pi_return[t][w][i] = 0;

				for(i=0; i<No_recharge; i++) 
 	   				pi_safe[t][w][i] = 0;
		   		
				for(i=0; i<No_recharge; i++) 
 	   				pi_storage[w][i] = 0;
		        
				for(i=0; i<No_recharge; i++) 
 	   				pi_flowBound[w][i] = 0;
			}
		}

		IloNumArray3  pi_Cut1(env, No_scenNode);
		IloNumArray3  pi_Cut3(env, No_scenNode);
		
        for(w=0;w<No_scenNode;w++){
			pi_Cut1[w] = IloNumArray2(env,No_scen);
			pi_Cut3[w] = IloNumArray2(env,No_scen);
			for(j=0;j<No_scen;j++){
				pi_Cut1[w][j] = IloNumArray(env,MaxCut);
				pi_Cut3[w][j] = IloNumArray(env,MaxCut);
				for(nCut=0;nCut<MaxCut;nCut++){
					pi_Cut1[w][j][nCut] = 0;
					pi_Cut3[w][j][nCut] = 0;
				}
			}
		
		} 
        
////////////BUILD MASTER AND SUBPROBLEM MODELS  //////////////////////////

		cout<<"building master and subproblems"<<endl;
        // arrays for model and lp
        CplexArray CPX(env, No_stage);
        ModelArray MODEL(env, No_stage);

		//Objective function array
		ObjArray      OBJ(env, No_stage);

		//Constraints array
		 EquationArray2 FlowBalance(env, No_stage);
		 EquationArray2 CapBound(env,  No_stage);
         EquationArray2 MeetDemand(env, No_stage);
		 EquationArray2 ReturnFlow(env,  No_stage); 
		 EquationArray2 SafeYield(env, No_stage); 
         EquationArray2 StorageBalance(env, No_stage); 
		 EquationArray2 RFoutflowBound(env,  No_stage);


		 //Optimality Cut
		 EquationArray2   Cut1(env, No_stage);
		 EquationArray2   Cut2(env, No_stage);
		 EquationArray2  Cut3(env, No_stage);

		 for(s = 0; s<No_stage; s++){
			 MODEL[s] = IloModel(env);
			 CPX[s]   = IloCplex(env);
			 OBJ[s]   = IloExpr(env);
			 FlowBalance[s] = EquationArray(env, No_yearMax);
			 CapBound[s]    = EquationArray(env, No_yearMax);
			 MeetDemand[s]  = EquationArray(env, No_yearMax);
			 ReturnFlow[s]  = EquationArray(env, No_yearMax);
			 SafeYield[s]  =  EquationArray(env, No_yearMax);
			 StorageBalance[s] = EquationArray(env, No_yearMax);
			 RFoutflowBound[s] = EquationArray(env, No_yearMax); 
 		     //Obj definition at each stage s 
			 for(i =0; i<No_link; i++){
				for(t=0;t<No_year[s];t++){
					if(s==1){
					OBJ[s] += pow(1+ disRate,-double(No_year[0]+ t))*cost[i]*Q[i][s][t] ;	
					}
					else if(s==2){
					OBJ[s] += pow(1+ disRate,-double(No_year[0]+ No_year[1]+ t))*cost[i]*Q[i][s][t] ;	
					}
					else if(s==3){
					OBJ[s] += pow(1+ disRate,-double(No_year[0]+ No_year[1] + No_year[2]+t))*cost[i]*Q[i][s][t] ;	
					}
					else{
					OBJ[s] += cost[i]*Q[i][s][t] ;	
					}
				}
			 }
			 
			 if( s<No_stage-1){  
				 for(j=0; j<No_scen; j++){
					 OBJ[s] += prob*lambda*Theta2[s+1][j]+prob*(1-lambda)*Theta1[s+1][j];
				 }
			 }
			 

			 
			 MODEL[s].add(IloMinimize(env, OBJ[s])); 
			 CPX[s].extract(MODEL[s]); 
			 OBJ[s].end();

			 char name1[100];
			 for(t = 0;t<No_yearMax;t++){
				 FlowBalance[s][t]     = IloRangeArray(env);
			     CapBound[s][t]        = IloRangeArray(env);
				 MeetDemand[s][t]      = IloRangeArray(env); 
				 ReturnFlow[s][t]      = IloRangeArray(env);
				 SafeYield[s][t]       = IloRangeArray(env);
				 StorageBalance[s][t]  = IloRangeArray(env); 
				 RFoutflowBound[s][t] = IloRangeArray(env); 
                 //define balnace
				 for(i=0;i<No_balance;i++){
					sprintf_s(name1, "flowbalance_%d_%d_%d",(int) t,(int) s, (int)i);
					IloExpr LHS(env), RHS(env);
					for(j=0;j<No_link;j++){
						if(endNode[j] == balanceID[i])
							LHS += Q[j][s][t]*(1-loss[j]);      //inflow
					} 
					for(j=0;j<No_link;j++){
						if(startNode[j] == balanceID[i])       
							RHS += Q[j][s][t] ;                 //outflow
					}
					if(t<No_year[s]){
						FlowBalance[s][t].add(LHS - RHS == 0); 
					    FlowBalance[s][t][i].setName(name1);
					}
					LHS.end(); RHS.end();
				 }
				

                 // CAP allocation + recharge facility capacity
				 //which constraints? Surface Water Supply Allotment and Treatment Plant Capacity Bounds?
                 for(i=0;i<No_capacity; i++){  
					 sprintf_s(name1, "CapBound_%d_%d_%d",(int) t,(int) s, (int)i);
					 IloExpr LHS(env);
					 IloNum RHS ; 
					 for(j=0;j<No_link;j++){
						if(startNode[j] == capacityID[i])
							LHS += Q[j][s][t];      //outflow
				     } 
					 RHS = capacity[i];
					 if(t<No_year[s]){
					    CapBound[s][t].add(LHS <= RHS);
					    CapBound[s][t][i].setName(name1);
					 }

					 LHS.end();
				 }
			

                 for(i=0;i<No_recharge; i++){ // recharge facility storage safe yield
					sprintf_s(name1, "SafeYield_%d_%d_%d",(int) i, (int) s, (int)t);
					if(t<No_year[s]){
						SafeYield[s][t].add(Storage[i][s][t] <= storageUB[i]);
						///>=???
						SafeYield[s][t][i].setName(name1);
					}
				 }
                
				 //potable and non-potable users 
				 for(i=0;i<No_user; i++){  
					sprintf_s(name1, "MeetDemand_%d_%d_%d",(int) s, (int) t, (int)i);
					IloExpr LHS(env);
					IloNum RHS; 
					RHS = demand[s][t][i];// right hand side is the demand; inflow = demand
					for(j=0;j<No_link;j++){
						if(endNode[j] == userID[i])
							LHS += Q[j][s][t]*(1-loss[j]); 
					}
					if(t<No_year[s]){
						MeetDemand[s][t].add(LHS == RHS);
					    MeetDemand[s][t][i].setName(name1);
					}
					LHS.end(); 
				 }
                 //  Return flow from potable users
                 for(i=0;i<No_potUser; i++){  
					sprintf_s(name1, "ReturnFlow_%d_%d_%d",(int) s,(int) t, (int)i);
					IloExpr LHS(env);
					IloNum RHS;	 
					RHS = returnRate*demand[s][t][i];                // 98% of the potable uses return to the facility
					for(j=0;j<No_link;j++){
						//why sum???
						if(startNode[j] == potUserID[i])
							LHS += Q[j][s][t];
					}  
					if(t<No_year[s]){
						ReturnFlow[s][t].add(LHS == RHS);
					    ReturnFlow[s][t][i].setName(name1);
					}
					LHS.end();
				 }
				
				 //Storage balance
                 for(i=0;i<No_recharge; i++){                  
					IloExpr LHS(env), RHS(env);
					sprintf_s(name1, "StorageBalance_%d_%d_%d",(int) s,(int) t, (int)i);
					for(j=0;j<No_link;j++){
						if(endNode[j] == rechargeID[i])
							LHS += Q[j][s][t]*(1-loss[j]);      //inflow
					} 
					for(j=0;j<No_link;j++){
						if(startNode[j] == rechargeID[i])       
							RHS += Q[j][s][t] ;                 //outflow
					} 
					if( t == 0){
						//inflow -Storage[i][t] - outflow  =  -Storage0[i]
						LHS += - Storage[i][s][t] ;
						if(s==0){ 
							StorageBalance[s][t].add(LHS - RHS == - Storage0[i]);
						}
						else{
							//I believe it has been taken care of later since it needs ancestor value (storage)??? later in the code it has been taken care of
							StorageBalance[s][t].add(LHS - RHS == 0);
						}
						StorageBalance[s][t][i].setName(name1);
					}
					else{
						//why it is not taking care of ancestor???!! (4.20)? I beleive the report must change and here is OK! like 4.24 later
						//-Storage[i][t] + Storage[i][t-1] + inflow -outflow=  0
						LHS += - Storage[i][s][t] + Storage[i][s][t-1];
						if(t<No_year[s]){
							StorageBalance[s][t].add(LHS - RHS == 0);
							StorageBalance[s][t][i].setName(name1);}
					} 
					LHS.end(); RHS.end();
				 }
				
                 //Recharge facility outflow bound
                 for(i=0;i<No_recharge; i++){ 
					IloExpr LHS(env);
					for(j=0;j<No_link;j++){
						if(startNode[j] == rechargeID[i])
							LHS += Q[j][s][t];
					}
					if(s>0){
						if(t==0)
							RFoutflowBound[s][t].add(LHS  <= 0);
						else{
							if(t<No_year[s])
								//why is it not having ancestor? later in the code it is clear that there is no need to ancestor, means that
								//equation (4.24) in the report must change
								RFoutflowBound[s][t].add(LHS - Storage[i][s][t-1] <= 0);
						}
					}
					else{
						if(t==0)
							RFoutflowBound[s][t].add(LHS - Storage0[i] <= 0);
						else{
							if(t<No_year[s])
								RFoutflowBound[s][t].add(LHS - Storage[i][s][t-1] <= 0);
						}
					}
					LHS.end(); 
				 }
				 if(t<No_year[s]){
					 MODEL[s].add(FlowBalance[s][t]);
					 MODEL[s].add(CapBound[s][t]);
					 MODEL[s].add(MeetDemand[s][t]);
					 MODEL[s].add(ReturnFlow[s][t]);
					 MODEL[s].add(SafeYield[s][t]);
					 MODEL[s].add(StorageBalance[s][t]);
					 MODEL[s].add(RFoutflowBound[s][t]);  
				 }
				}//end year 
				 

             Cut1[s] =EquationArray(env,No_scen);
			 Cut2[s] = EquationArray(env,No_scen);
			 Cut3[s] = EquationArray(env,No_scen);
			 if(s<No_stage-1){
				for(w = 0; w<No_scen;w++){
					Cut1[s][w] = IloRangeArray(env);
					Cut2[s][w] = IloRangeArray(env);
					Cut3[s][w] = IloRangeArray(env);				 
					for(nCut = 0; nCut<MaxCut; nCut++){
                         IloExpr OptCutExpr1(env) ;
                         for(i=0;i<No_recharge; i++){ 
							 OptCutExpr1 += Storage[i][s][No_year[s]-1]; 
						 }
						 sprintf_s(name1, "Cut1_%d_%d_%d",(int) s, (int)w,(int) nCut);
						 Cut1[s][w].add(Theta1[s+1][w] + OptCutExpr1 >= -IloInfinity );
						 Cut1[s][w][nCut].setName(name1);
						 OptCutExpr1.end();

                         IloExpr OptCutExpr2(env) ;
						 sprintf_s(name1, "Cut2_%d_%d_%d",(int) s, (int)w,(int) nCut);
						 OptCutExpr2 = VaR[nCut][s];
						 for(w2 = 0; w2<No_scenNode;w2++){
							 if(stageMap[w2]==s+1){ 
								 OptCutExpr2 += V[nCut][w2];
							 }
						 }
						 Cut2[s][w].add(Theta2[s+1][w] - OptCutExpr2 >= 0 );
						 Cut2[s][w][nCut].setName(name1);
						 OptCutExpr2.end(); 

						 IloExpr OptCutExpr3(env);
						 OptCutExpr3 = VaR[nCut][s];
						 for(w2 = 0; w2<No_scenNode;w2++){
							 if(stageMap[w2]==s+1){ 
								 OptCutExpr3 += V[nCut][w2];
							 }
						 }
						 for(i=0;i<No_recharge; i++){ 
							 OptCutExpr3 += Storage[i][s][No_year[s]-1]; 
						 }
						 sprintf_s(name1, "Cut3_%d_%d_%d",(int) s, (int)w,(int) nCut);
						 Cut3[s][w].add(OptCutExpr3 >= -IloInfinity);
						 Cut3[s][w][nCut].setName(name1);
						 OptCutExpr3.end();
					 }				
				}
			}
			 
			CPX[s].extract(MODEL[s]);  

		 }//end stage
		 for(s=0;s<No_stage;s++){
			 CPX[s].extract(MODEL[s]); 
		 }
		

		 
		//////////BEGIN ITERATIONS/////////////////////////////////
        cout<<"begin iteration..."<<endl;
		cout.precision(10);
        //Initialize LB, UB
		IloNum rel_Gap = IloInfinity;
		IloNum LB = 0;
		IloNum UB = IloInfinity;
        //Scenario cost
		IloNumArray subObj_hat(env, No_scenNode);
		//IloNumArray subObj_cx(env, No_scenNode);
        IloNum z_hat=0;
        		char resName[100];
		int l, a;
		if (lambda==0.3){
			l=3;
		}
		else if (lambda==0.5){
			l=5;
		}
		else if (lambda==0.7){
			l=7;
		}
		else if (lambda==0.8){
			l=8;
		}
		else {
			l=9;
		}

		if (alpha==0.8){
			a=8;
		}
		else if (alpha==0.9){
			a=9;
		}
		else {
			a=95;
		}
		sprintf_s(resName, "%d_%d.txt",(int) l, (int) a);

        const char* Resfilename  = resName;
		ofstream Result(Resfilename);
        Result << "Iter  \t    Z_hat    \t     LB      \t UB    \t rel_Gap"<<endl;


	


		IloNumArray3 h(env,No_scenNode);
        IloNumArray3 CutRhs1(env,No_scenNode); 
		IloNumArray4 CutLhs1(env,No_scenNode);
		IloNumArray3 CutRhs3(env,No_scenNode);
		IloNumArray4 CutLhs3(env,No_scenNode);
		for(w2=0; w2<No_scenNode; w2++){
			h[w2] = IloNumArray2(env,No_scen); 

			CutRhs1[w2] = IloNumArray2(env,No_scen); 
			CutLhs1[w2] = IloNumArray3(env,No_scen);
			CutRhs3[w2] = IloNumArray2(env,No_scen);
			CutLhs3[w2] = IloNumArray3(env,No_scen);  
			for(j = 0;j<No_scen;j++){

				h[w2][j] = IloNumArray(env,MaxCut);
				CutRhs1[w2][j] = IloNumArray(env,MaxCut);
				CutLhs1[w2][j] = IloNumArray2(env,MaxCut);
				CutRhs3[w2][j] = IloNumArray(env,MaxCut);
				CutLhs3[w2][j] = IloNumArray2(env,MaxCut);
				for(nCut=0;nCut<MaxCut;nCut++){
					h[w2][j][nCut] = 0;

					CutRhs1[w2][j][nCut] = 0;
					CutRhs3[w2][j][nCut] = 0;
					CutLhs1[w2][j][nCut] = IloNumArray(env,No_recharge);
					CutLhs3[w2][j][nCut] = IloNumArray(env,No_recharge);
					for(i=0;i<No_recharge;i++){
						CutLhs1[w2][j][nCut][i] = 0;
					    CutLhs3[w2][j][nCut][i] = 0;
					}
				}
			}
		}

		///////////////////////////////////////////////////////////
        const clock_t begin_time = clock();
        for(nCut = 0; nCut<MaxCut; nCut++){ //main loop
			cout<<"nCut ="<<nCut<<endl;
			for(i=0;i<No_recharge;i++){
				for(w=0;w<No_scenNode;w++){
					storage_hat[i][w] = 0;
				}
			}
			for(w=0;w<No_scenNode;w++){
				subObj_hat[w] = 0;
				//subObj_cx[w] = 0; 
			}

            //***Forward pass***///
            //Step 1: solve master problem -- First stage relaxed problem 
			 

            CPX[0].solve();  //solve first stage relaxed problem
		    if( !CPX[0].solve () ){
			   env.error() << "Master problem infeaasible" << endl;
			   throw(-1);
		    } 

			//record storage_hat[i][0]
			//to pass to next stage
			for(i=0;i<No_recharge;i++){
			    storage_hat[i][0] = CPX[0].getValue(Storage[i][0][No_year[0] - 1]);
		    }
			//record obj_master = CPX[0].getObjValue();
			LB = CPX[0].getObjValue();
			

			subObj_hat[0] = CPX[0].getObjValue() ;
			for(j=0;j<No_scen;j++){
				subObj_hat[0] -=  prob*lambda*CPX[0].getValue(Theta2[1][j]);
				subObj_hat[0] -= prob*(1-lambda)*CPX[0].getValue(Theta1[1][j]); 
			}
			
			//subObj_hat[0] = subObj_cx[0];
			
			//Step 2: For each node at each stage, update the sceanrio demands
            for(s=1; s<No_stage;s++){ 
				//update subproblem
	            for(w=1; w<No_scenNode; w++){
					if(stageMap[w]==s){ 
						for(t=0;t<No_year[s];t++){
							//randomness for RHS of supply, for capacityID=0
							//scenFac[w][0] is reserved for \Xi_1
						   CapBound[s][t][0].setBounds(-IloInfinity, CAPamt[s][t]*scenFac[w][0]);
		                   for(i=0;i<No_user;i++)
							   //scenFac[w][t+1] corresponds to \Xi_2 (t)
							   MeetDemand[s][t][i].setBounds(scenFac[w][t+1]*demand[s][t][i], scenFac[w][t+1]*demand[s][t][i]);
						   for(i=0;i<No_potUser; i++)
							   ReturnFlow[s][t][i].setBounds(returnRate*scenFac[w][t+1]*demand[s][t][i],returnRate*scenFac[w][t+1]*demand[s][t][i]);
						}
						
						 
						for(w2=0; w2<No_scenNode; w2++){
							if(ancestor[w2][w]==1){
								//w2 is ancestor of w
							   for(i=0;i<No_recharge; i++){
								   StorageBalance[s][0][i].setBounds(-storage_hat[i][w2],-storage_hat[i][w2]);
								   RFoutflowBound[s][0][i].setBounds(-IloInfinity, storage_hat[i][w2]);
								}
							   //note if ancestor of w is found, there is no need to continue loop, since every node has only one ancestor
							}
						}
						//update cuts
						if(s<No_stage-1 && nCut>0){ 
							for(k=0;k<nCut;k++){ 
								for(j=0;j<No_scen;j++){
									Cut1[s][j][k].setBounds(CutRhs1[w][j][k], IloInfinity);
									Cut2[s][j][k].setBounds(0,IloInfinity);
									Cut3[s][j][k].setBounds(CutRhs3[w][j][k], IloInfinity);
									for(i=0;i<No_recharge; i++){ 
									   Cut3[s][j][k].setLinearCoef(Storage[i][s][No_year[s]-1],-CutLhs3[w][j][k][i]);
									   Cut1[s][j][k].setLinearCoef(Storage[i][s][No_year[s]-1], -CutLhs1[w][j][k][i]);
									}
									for(w2=0;w2<No_scenNode;w2++){
										if(w2 == decendant[w][j]){
											Cut2[s][j][k].setLinearCoef(V[k][w2], -1/(1-alpha));
											Cut3[s][j][k].setLinearCoef(V[k][w2],1);
										}
										else{
											Cut2[s][j][k].setLinearCoef(V[k][w2],0);
											Cut3[s][j][k].setLinearCoef(V[k][w2],0);
										}
									}
								}
							}
						}
						
			          
						CPX[s].solve();
		 
						for(i=0;i<No_recharge;i++){
							storage_hat[i][w] = CPX[s].getValue(Storage[i][s][No_year[s] - 1]);
						} 
						//store dual variables at last stage	
						if(s==No_stage-1){
							//they are the same since there is no Theta
							subObj_hat[w] = CPX[s].getObjValue();
							//subObj_cx[w] = subObj_hat[w];
							
							for(t=0;t<No_year[s];t++){ 
								for(i=0; i<No_capacity; i++) {
 	   								pi_capacity[t][w][i] = CPX[s].getDual(CapBound[s][t][i]); 
								}
								
								for(i=0; i<No_user; i++) {
 	   								pi_demand[t][w][i] = CPX[s].getDual(MeetDemand[s][t][i]); 
								}
								
								for(i=0; i<No_potUser; i++) {
 	   								pi_return[t][w][i] = CPX[s].getDual(ReturnFlow[s][t][i]); 
								}

								for(i=0; i<No_recharge; i++){ 
 	   								pi_safe[t][w][i] = CPX[s].getDual(SafeYield[s][t][i]);
								}
								
							}
							for(i=0; i<No_recharge; i++){ 
								//we need it to obtain coefficent of the Storage[i][s][No_year[s] - 1] in the cut
 	   							pi_storage[w][i] = CPX[s].getDual(StorageBalance[s][0][i]); 
							}
						        
							for(i=0; i<No_recharge; i++){ 
								//we need it to obtain coefficent of the Storage[i][s][No_year[s] - 1] in the cut
 	   							pi_flowBound[w][i] = CPX[s].getDual(RFoutflowBound[s][0][i]); 
							}
							//note these cuts are not appearing in last stage
							
							for(j=0;j<No_scen;j++){
							   pi_Cut3[w][j][nCut] =  0; 
							   pi_Cut1[w][j][nCut] = 0;
							}
						}
						else{   
							//for other stages than last stage
							subObj_hat[w] = CPX[s].getObjValue() ;
							for(j=0;j<No_scen;j++){
								subObj_hat[w] -=  prob*lambda*CPX[s].getValue(Theta2[s+1][j]);
								subObj_hat[w] -= prob*(1-lambda)*CPX[s].getValue(Theta1[s+1][j]); 
							}
			                //subObj_hat[w] = subObj_cx[w]; 
						}//end if
					}
				}//end scenario
			}//end stage 
			//end forward pass

			cout<<"calculating gap"<<endl;
			///*****Calculate UB*****////

			//*****************************************
				///Need to form CVaR problem for UB
			//*****************************************
			for(s = No_stage-1; s>0; s--){
                for(w=0; w<No_scenNode; w++){
					IloNumArray z_w(env, No_scen);
					IloNumArray temp(env, No_scen);
					if(stageMap[w] == s-1){
						for(j = 0; j<No_scen; j++){
							z_w[j] = subObj_hat[decendant[w][j]]; 
							subObj_hat[w] +=  prob*subObj_hat[decendant[w][j]]*(1-lambda); 
							
						} 
						
						merge_sort(z_w, temp, 0, No_scen-1);
						//insertion_sort(No_scen, z_w);
						subObj_hat[w] += lambda* z_w[CVaR_No];
						for(i = CVaR_No+1; i<No_scen; i++){
							subObj_hat[w]+= prob*lambda*(z_w[i] - z_w[CVaR_No])/(1-alpha);
						}
					}
					
				}  
			}

			z_hat = subObj_hat[0];


			            

			//*****Calculate gap****////
			if(z_hat < UB)
				UB = z_hat;
			rel_Gap = (UB - LB)/LB;
            Result<< nCut+1 << '\t' << z_hat << '\t' << LB <<'\t' << UB <<'\t' << rel_Gap << endl;
			
			// Convergence criterion
			if( rel_Gap <= toler ||nCut ==MaxCut-1 || float( clock () - begin_time ) /  CLOCKS_PER_SEC >7200){ 
			   cout<<"DONE"<<endl;
			   Result << float( clock () - begin_time ) /  CLOCKS_PER_SEC <<endl;
			   std::cout << "time =" << float( clock () - begin_time ) /  CLOCKS_PER_SEC<<endl; 
			   break;
			}


			////////////////////////// 

            for(s = No_stage-1; s>0; s--){
				
				for(j=0;j<No_scen;j++){
					MODEL[s-1].add(Cut1[s-1][j][nCut]);
					MODEL[s-1].add(Cut2[s-1][j][nCut]);
					MODEL[s-1].add(Cut3[s-1][j][nCut]);
				}
                for(w=0; w<No_scenNode; w++){
					if(stageMap[w]==s-1){ 
						for(k=0;k<nCut;k++){  
							for(j=0;j<No_scen;j++){
								Cut1[s-1][j][k].setBounds(CutRhs1[w][j][k], IloInfinity);
								Cut2[s-1][j][k].setBounds(0,IloInfinity);
								Cut3[s-1][j][k].setBounds(CutRhs3[w][j][k], IloInfinity);
								for(i=0;i<No_recharge; i++){ 
								   Cut1[s-1][j][k].setLinearCoef(Storage[i][s-1][No_year[s-1]-1], -CutLhs1[w][j][k][i]);
								   Cut3[s-1][j][k].setLinearCoef(Storage[i][s-1][No_year[s-1]-1],-CutLhs3[w][j][k][i]);
								}
								for(w2=0;w2<No_scenNode;w2++){
									if(w2 == decendant[w][j]){
										Cut2[s-1][j][k].setLinearCoef(V[k][w2], -1/(1-alpha));
										Cut3[s-1][j][k].setLinearCoef(V[k][w2],1);
									}
									else{
										Cut2[s-1][j][k].setLinearCoef(V[k][w2],0);
										Cut3[s-1][j][k].setLinearCoef(V[k][w2],0);
									}
								}
							}
						}		
				
						
							
						
						
                       //you have done that before
                        for(t=0;t<No_year[s-1];t++){
							CapBound[s-1][t][0].setBounds(-IloInfinity, CAPamt[s-1][t]*scenFac[w][0]);
		                   for(i=0;i<No_user;i++)
							   MeetDemand[s-1][t][i].setBounds(scenFac[w][t+1]*demand[s-1][t][i], scenFac[w][t+1]*demand[s-1][t][i]);
						   for(i=0;i<No_potUser; i++)
							   ReturnFlow[s-1][t][i].setBounds(returnRate*scenFac[w][t+1]*demand[s-1][t][i],returnRate*scenFac[w][t+1]*demand[s-1][t][i]);
						}
						for(w2=0; w2<No_scenNode; w2++){
							if(s-1>0){
								if(ancestor[w2][w]==1){
								   for(i=0;i<No_recharge; i++){
									   //update rhs using storage of previous stage
									   StorageBalance[s-1][0][i].setBounds(-storage_hat[i][w2],-storage_hat[i][w2]);
									   RFoutflowBound[s-1][0][i].setBounds(-IloInfinity, storage_hat[i][w2]);
									}
								}
							}
						}
						//*******************
						for(j=0; j<No_scen; j++){
						    for(t =0; t<No_year[s];t++){
				               for(i=0; i<No_user; i++){
					              h[w][j][nCut] += pi_demand[t][decendant[w][j]][i]*scenFac[decendant[w][j]][t+1]*demand[s][t][i];
							   }
				               for(i=0; i<No_potUser; i++){
					              h[w][j][nCut] += pi_return[t][decendant[w][j]][i]*returnRate*scenFac[decendant[w][j]][t+1]*demand[s][t][i];
							   }
						       for(i=1;i<No_capacity; i++){ 
					              h[w][j][nCut] += pi_capacity[t][decendant[w][j]][i]*capacity[i];
							   }
							   h[w][j][nCut] += pi_capacity[t][decendant[w][j]][0]*CAPamt[s][t]*scenFac[decendant[w][j]][0];
							   for(i=0;i<No_recharge; i++){ 
					              h[w][j][nCut] += pi_safe[t][decendant[w][j]][i]*storageUB[i];
							   }
							}
					   }
						if (s< No_stage-1){
							for(j=0; j<No_scen; j++){ 
								for(k=0;k<nCut+1;k++){
									for(w2 = 0; w2<No_scen; w2++){
									
										h[w][j][nCut] += pi_Cut1[decendant[w][j]][w2][k]*CutRhs1[decendant[w][j]][w2][k];
									
										h[w][j][nCut] += pi_Cut3[decendant[w][j]][w2][k]*CutRhs3[decendant[w][j]][w2][k];
								  
									} 
								} 
							}
						}
						 for(j=0; j<No_scen; j++){
							
							CutRhs3[w][j][nCut]=h[w][j][nCut];
							CutRhs1[w][j][nCut]=CutRhs3[w][j][nCut];
							Cut3[s-1][j][nCut].setBounds(CutRhs3[w][j][nCut], IloInfinity);
							
						    Cut1[s-1][j][nCut].setBounds(CutRhs1[w][j][nCut], IloInfinity);
						}

						for(i=0;i<No_recharge; i++){
							for(j=0; j<No_scen; j++){
								
								CutLhs3[w][j][nCut][i] = -pi_storage[decendant[w][j]][i]+ pi_flowBound[decendant[w][j]][i];
								Cut3[s-1][j][nCut].setLinearCoef(Storage[i][s-1][No_year[s-1]-1], -CutLhs3[w][j][nCut][i]);
								CutLhs1[w][j][nCut][i]=CutLhs3[w][j][nCut][i];
								Cut1[s-1][j][nCut].setLinearCoef(Storage[i][s-1][No_year[s-1]-1], -CutLhs1[w][j][nCut][i]); 
								 
							}
						}  

						for(j=0; j<No_scen; j++){  
							Cut2[s-1][j][nCut].setBounds(0,IloInfinity);
							for(w2=0;w2<No_scenNode;w2++){
								if(w2 == decendant[w][j]){
									Cut2[s-1][j][nCut].setLinearCoef(V[nCut][w2],-1/(1-alpha));
									Cut3[s-1][j][nCut].setLinearCoef(V[nCut][w2],1);
								}
								else{
									Cut2[s-1][j][nCut].setLinearCoef(V[nCut][w2],0);
									Cut3[s-1][j][nCut].setLinearCoef(V[nCut][w2],0);
								}
							}
						}
						///*******************

	
					  
						
						if(s<=No_stage-1 && s>1){
							CPX[s-1].solve();
							for(t=0;t<No_year[s-1];t++){ 
								for(i=0; i<No_capacity; i++) {
   									pi_capacity[t][w][i] = CPX[s-1].getDual(CapBound[s-1][t][i]); 
								}
								
								for(i=0; i<No_user; i++) {
   									pi_demand[t][w][i] = CPX[s-1].getDual(MeetDemand[s-1][t][i]); 
								}
								
								for(i=0; i<No_potUser; i++) {
   									pi_return[t][w][i] = CPX[s-1].getDual(ReturnFlow[s-1][t][i]); 
								}

								for(i=0; i<No_recharge; i++) 
 	   								pi_safe[t][w][i] = CPX[s-1].getDual(SafeYield[s-1][t][i]);
							}
							for(i=0; i<No_recharge; i++) {
   								pi_storage[w][i] = CPX[s-1].getDual(StorageBalance[s-1][0][i]); 
							}
						        
							for(i=0; i<No_recharge; i++) {
   								pi_flowBound[w][i] = CPX[s-1].getDual(RFoutflowBound[s-1][0][i]); 
							}

							

						    //Dual value for cut3
							for(k=0;k<nCut+1;k++){
								for(j=0;j<No_scen;j++){
									pi_Cut1[w][j][k] =  CPX[s-1].getDual(Cut1[s-1][j][k]); 
									pi_Cut3[w][j][k] =  CPX[s-1].getDual(Cut3[s-1][j][k]); 
								}
							}
						} 
					}
				}
			}//end backward  
			
		}//end main loop
env.end();		
}//end loop on instances
	}

    catch(IloException &e) {
	env1.out() << "ERROR: " << e << endl;
	}
	catch(...){
	env1.out() << "Unknown exception" << endl;
	}
	env1.end();
	
	getchar();
	return 0;

}