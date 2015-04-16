///////////////////////////////////////////////////////////////////
//    RANDOM BAMBOO    |     v0.5.1    |      April 16, 2015     //
//---------------------------------------------------------------//
//              (C) 2014 Han Zhang, Yifan Yang                   //
//              GNU General Public License  V3                   //
//---------------------------------------------------------------//
//   For documentation, citation and bug-report instructions:    //
//           http://www.hanzhang.name/softwares/rb               //
//---------------------------------------------------------------//
//                                                               //
//   Update: 0.5.x                                               //
//                                                               //
//     October 9, 2014                                           //
//     (1) remove --neighbor and --searchback introduced         //
//         in v0.4.x. As expected, these two options don't help  //
//         to improve the power due to multiple-comparison issue //
//                                                               //
//     April 16, 2015                                            //
//     (1) add --swt to allow individual spcific sample weights  //
//         in bootstrap                                          //
///////////////////////////////////////////////////////////////////


#include "config.h"
#include "constant.h"
#include "bamboo.h"


int main(int argc, char **argv){
	
	cout << endl;
	cout << "+------------------------+-------------------+---------------------+" << endl;
	cout << "|     Random Bamboo      |     " << setw(8) << RANDOM_BAMBOO_VERSION << "      |      04/16/2015     |" << endl;
	cout << "+------------------------+-------------------+---------------------+" << endl;
	cout << "|                   (C) 2015 Han Zhang, Yifan Yang                 |" << endl;
	cout << "|                   GNU General Public License, V3                 |" << endl;
	cout << "+------------------------------------------------------------------+" << endl;
	cout << "|     For documentation, citation and bug-report instructions:     |" << endl;
	cout << "|               http://www.hanzhang.name/softwares/rb              |" << endl;
	cout << "+------------------------------------------------------------------+" << endl;
	cout << endl;
	
	char *file = NULL;//f
	char *out = NULL;//o
	char *cont = NULL;//c
	char *cate = NULL;//a
	char *pred = NULL;//p
	char *bam = NULL;//b
	char *trainid = NULL;//y
	char *testid = NULL;//z
	char *snpid = NULL;//S
	char *varid = NULL;//j
	char *samplewt = NULL;//W
	
	int ntree = 1;//t
	int mtry = 0;//m
	int seed = 1;//s
	int max_nleaf = 1000000;//l
	int min_leaf_size = 1;//r
	int imp_measure = 1;//i
	int nthread = sysconf( _SC_NPROCESSORS_ONLN );//h
	//if class_weight < 0, reweight to balance the data
	double class_weight = -1.0;//w
	double cutoff = 1.0;//u
	
	bool flip = false;//g
	bool output_prox = false;//x
	bool output_imp = true;//n
	bool balance = false;//B //if it is turned on and class_weight is not set (< .0), then class_weight by balancing the data (according to case/control ratio in data)
	bool trace = false;//r
	bool output_bamboo = true;//N
	
	bool file_spec = false;
	bool out_spec = false;
	bool pred_spec = false;
	bool bam_spec = false;
	bool ntree_spec = false;
	
	int c;
	while(true){
		static struct option long_options[] = {
			{"file", 1, NULL, 'f'},
			{"out", 1, NULL, 'o'},
			{"cont", 1, NULL, 'c'},
			{"cate", 1, NULL, 'a'},
			{"pred", 1, NULL, 'p'},
			{"bam", 1, NULL, 'b'},
			{"trainid", 1, NULL, 'y'},
			{"testid", 1, NULL, 'z'},
			{"snpid", 1, NULL, 'S'},
			{"varid", 1, NULL, 'j'},
			{"ntree", 1, NULL, 't'},
			{"mtry", 1, NULL, 'm'},
			{"seed", 1, NULL, 's'},
			{"maxnleaf", 1, NULL, 'l'},
			{"minleafsize", 1, NULL, 'e'},
			{"imp", 1, NULL, 'i'},
			{"nthread", 1, NULL, 'd'},
			{"classwt", 1, NULL, 'w'},
			{"swt", 1, NULL, 'W'}, 
			{"cutoff", 1, NULL, 'u'},
			{"flip", 0, NULL, 'g'},
			{"prox", 0, NULL, 'x'},
			{"noimp", 0, NULL, 'n'},
			{"balance", 0, NULL, 'B'},
			{"trace", 0, NULL, 'r'},
			{"nobam", 0, NULL, 'N'},
			{"help", 0, NULL, 'h'},
			{"version", 0, NULL, 'v'},
			{0, 0, 0, 0}
		};
		
		//jkq
		
		int option_index = 0;
		c = getopt_long(argc, argv, "f:o:c:a:p:b:y:z:S:j:t:m:s:l:e:i:d:w:W:u:gxnBrNhv", long_options, &option_index);
		
		if(c == -1){
			break;
		}
		
		switch(c){
			case 'f':
				file = optarg;
				file_spec = true;
				break;
			case 'o':
				out = optarg;
				out_spec = true;
				break;
			case 'c':
				cont = optarg;
				break;
			case 'a':
				cate = optarg;
				break;
			case 'p':
				pred = optarg;
				pred_spec = true;
				break;
			case 'b':
				bam = optarg;
				bam_spec = true;
				break;
			case 'y':
				trainid = optarg;
				break;
			case 'z':
				testid = optarg;
				break;
			case 'S':
				snpid = optarg;
				break;
			case 'j':
				varid = optarg;
				break;
			case 't':
				ntree = atoi(optarg);
				ntree_spec = true;
				break;
			case 'm':
				mtry = atoi(optarg);
				if(mtry <= 0){
					cout << "Error: The option --mtry must be a positive integer" << endl;
					return 0;
				}
				break;
			case 's':
				seed = atoi(optarg);
				break;
			case 'l':
				max_nleaf = atoi(optarg);
				if(max_nleaf <= 0){
					cout << "Error: The option --maxnleaf must be a positive integer" << endl;
					return 0;
				}
				break;
			case 'e':
				min_leaf_size = atoi(optarg);
				if(min_leaf_size <= 0){
					cout << "Error: The option --minleafsize must be a positive integer" << endl;
					return 0;
				}
				break;
			case 'i':
				imp_measure = atoi(optarg);
				if(imp_measure != IMP_GINI && imp_measure != IMP_BREIMAN_CUTLER && imp_measure != IMP_LIAW_WIENER){
					cout << "Error: The option --imp only support " << IMP_GINI << " (default), " << IMP_BREIMAN_CUTLER << " and " << IMP_LIAW_WIENER << " right now. Program terminates" << endl;
					return 0;
				}
				if(imp_measure == IMP_NULL){
					imp_measure = 1;
				}
				break;
			case 'd':
				nthread = atoi(optarg);
				if(nthread <= 0){
					cout << "Error: The option --nthread must be a positive integer" << endl;
					return 0;
				}
				
				if(nthread > sysconf( _SC_NPROCESSORS_ONLN )){
					nthread = sysconf( _SC_NPROCESSORS_ONLN );
					cout << "Warning: --nthread is reset to be " << sysconf( _SC_NPROCESSORS_ONLN ) << " (# CPUs)" << endl;
				}
				break;
			case 'w':
				class_weight = atof(optarg);
				if(class_weight <= .0){
					cout << "Error: The option --classwt must be a positive number" << endl;
					return 0;
				}
				break;
			case 'W':
				samplewt = optarg;
				break;
			case 'u':
				cutoff = atof(optarg);
				if(cutoff <= .0){
					cout << "Error: The option --cutoff must be a positive number" << endl;
					return 0;
				}
				break;
			case 'g':
				flip = true;
				break;
			case 'x':
				output_prox = true;
				break;
			case 'n':
				output_imp = false;
				break;
			case 'B':
				balance = true;
				break;
			case 'r':
				trace = true;
				break;
			case 'N':
				output_bamboo = false;
				break;
			case 'h':
				//print help info here
				return 0;
			case 'v':
				return 0;
			default:
				cout << "Error: Unknow options. Program terminates" << endl;
				exit(1);
		}
	}
	
	
	if(!file_spec){//no plink file is specified
		if(!pred_spec && !bam_spec){
			cout << "Error: The option --file must be specified" << endl;
			return 0;
		}else if(pred_spec && !bam_spec){
			cout << "Error: The option --bam must be specified for prediction purpose" << endl;
			return 0;
		}else if(!pred_spec && bam_spec){//print the information of the specified model
			BAMBOO bb(bam, nthread);
			bb.PrintModelInfo();
			return 0;
		}else{
			;//good
		}
	}else{//train model from specified data
		bam = NULL;
		bam_spec = false;
	}
	
	if(!out_spec){//if the name of output file is not specified
		if(file_spec && pred_spec){//fit model and predict new data
			cout << "Error: When options --file and --pred are specified simultaneouslythe option --out must be spcified by the user" << endl;
			exit(1);
		}else if(file_spec && !pred_spec){//fit model but don't predict new data
			out = file;
		}else if(!file_spec && pred_spec){//predict new data by given model
			out = pred;
		}else{//impossible
			cout << "Error: debug" << endl;
			exit(1);
		}
	}
	
	if(balance && class_weight > .0){
		cout << "Error: Please don't use the option --balance with option --classwt simultaneously" << endl;
		return 0;
	}
	
	//if balance = true, the program will not check the value of class_weight and will determine it automatically
	//so don't worry about the following. If balance = false and class_weight is not set by the user (class_weight < .0), 
	//then set it as 1.0 (means no weight). So it is important to initialize class_weight as -1.0 at the very beginning
	if(class_weight <= .0){
		class_weight = 1.0;
	}
	
	
	if(!file_spec && pred_spec){//predicting testing data (--pred) using models specified in --bam
		if(!ntree_spec){
			ntree = -1;
		}
		BAMBOO bb (out, pred, bam, testid, nthread, ntree);
		bb.PredictTestingSample();
		return 0;
	}else{//training model, and then predicting testing data if --pred is on
		BAMBOO bb (file, out, cont, cate, pred, trainid, testid, snpid, samplewt, ntree, mtry, max_nleaf, 
		min_leaf_size, imp_measure, seed, nthread, class_weight, cutoff, flip, output_prox, output_imp, 
		output_bamboo, balance, trace);
		bb.GrowForest();
		bb.PredictTestingSample();
		return 0;
	}
	
	return 0;
	
}
