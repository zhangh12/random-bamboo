///////////////////////////////////////////////////////////////////
//    RANDOM BAMBOO    |     v0.2.3    |   September 9, 2014     //
//---------------------------------------------------------------//
//      (C) 2014 Han Zhang, GNU General Public License  V3       //
//---------------------------------------------------------------//
//   For documentation, citation and bug-report instructions:    //
//           http://www.hanzhang.name/softwares/rb               //
//---------------------------------------------------------------//
//                                                               //
//   Update: v0.1.x                                              //
//                                                               //
//     June 2, 2014                                              //
//     (1) Bugs fixed in computing OOB error                     //
//     (2) Bugs fixed in computing Gini importance               //
//         The sum of decrease in impurity from splitting on     //
//         a variable should be weighted by the sample size      //
//         in its parent node, i.e.                              //
//         DEC = N_P ( tau(P) - p_L tau(L) - p_R tau(R) )        //
//         The bug is because I missed N_P in the formula        //
//         which was not mentioned in many literatures           //
//                                                               //
//     June 3, 2014                                              //
//     (1) Revise and simplify the data structure to handle      //
//         continuous covariates                                 //
//                                                               //
//     June 6, 2014                                              //
//     (1) Fix bugs when continuous covariates are involved      //
//     (2) handle categorical covariates                         //
//                                                               //
//     June 18, 2014                                             //
//     (1) Use reference in loops in SplitNode() to speedup      //
//         the program. x2 accelerate                            //
//                                                               //
//     June 20, 2014                                             //
//     (1) For the nodes with small sample size, a better way    //
//         to compute the criteria is NOT to use the Boolean     //
//         operations to construct the contingency tables.       //
//         x2 accelerate                                         //
//                                                               //
//     June 23, 2014                                             //
//     (1) unify the data structure used for genotype and        //
//         categorical variables                                 //
//                                                               //
//     July 7, 2014                                              //
//     (1) Re-write the whole program. The result is slightly    //
//         different from v0.1.5. The reason is that I           //
//         re-define the left and right children of a tree       //
//         and thus different random series (selected SNPs       //
//         and/or covariates) are applied to the child leaves    //
//                                                               //
//     July 8, 2014                                              //
//     (1) Add CompPermutationImportance()                       //
//     (2) Fix bugs in prediction function PutDownSampleToLeaf   //
//                                                               //
//     July 14, 2014                                             //
//     (1) Re-wirte the prediction function. Predict a sample    //
//         by using the tree structure rather than screening     //
//         all the leaves                                        //
//     (2) Add option --noflip. By default, flip is on and the   //
//         genotype is recoded if MAF > 0.5. Note that the       //
//         tree built on a flipped genotypes can be different    //
//         from the one built on the original genotypes          //
//         (compared with last version)                          //
//     (3) Enable multi-Processing by adding option --nthread    //
//     (4) Fix bug in calling openMP                             //
//     (5) Fix bug in computing OOB error                        //
//     (6) Output various permuted variable importances          //
//                                                               //
//     July 17, 2014                                             //
//     (1) Add option --classwt to handle imbalanced data.       //
//         The modifications occured in two places. One is in    //
//         function SplitNode when computing the Gini impurity.  //
//         The other one is in determining the class prediction  //
//         in a terminal node. Some of the members of struct     //
//         NODE was modified accordingly                         //
//     (2) Add option --cutoff. It only works in the final vote, //
//         i.e., aggregatting the weighted votes from each       //
//         individual tree                                       //
//     (3) Compute confusion matrix using oob samples            //
//     (4) Modify the definition of gini impurity.               //
//     (5) Modify the way to define the class prediction in a    //
//         leaf                                                  //
//     (6) Modify the way to define the final class prediction   //
//         of the whole forest                                   //
//     (7) (1, 4-6) are based on rpart's manual when class       //
//         weights are taken into account                        //
//                                                               //
//     July 18, 2014                                             //
//     (1) Fix bugs. In the former version, I only updated the   //
//         gini impurity when splitting with SNPs. Now I made it //
//         for categorical and continuous covariates             //
//     (2) Reorganized the file structure                        //
//     (3) Turn off automatical reweight and remove options      //
//         --noreweight, since the default weight (naively       //
//         balance two classes) doesn't give result close to     //
//         the optimal one. The option --classwt is still        //
//         available. Instead, add an option --balance to adjust //
//         the data to be 1:1 balanced                           //
//     (4) Output confusion matrix to local file *.cof           //
//     (5) Delete path2node in struct NODE                       //
//                                                               //
//   Update: v0.2.x                                              //
//                                                               //
//     July 18, 2014                                             //
//     (1) Simplify the struct NODE to facilitate the saving     //
//         procedure of the forest. Potential memory issues      //
//         fixed by redesigning the copy constructor             //
//                                                               //
//     July 19, 2014                                             //
//     (1) Add feature to save the forest by default             //
//     (2) Add option --nobamboo to switch off the feature of    //
//         saving forest                                         //
//                                                               //
//     July 21, 2014                                             //
//     (1) Fix bug when computing permutation importance, i.e.,  //
//         if the SD is zero, then the scaled importance is      //
//         marked as special number                              //
//     (2) Serialize the forest to local file for future         //
//         prediction                                            //
//                                                               //
//     July 26, 2014                                             //
//     (1) Load test data                                        //
//     (2) Add prediction feature                                //
//     (3) Mute option --flip                                    //
//                                                               //
//     September 9, 2014                                         //
//     (1) Load specified model and testing data                 //
//     (2) Modify prediction method                              //
//     (3) Interface extended. Saved model can be used directly  //
//         in prediction                                         //
//     (4) Add progress bar                                      //
//                                                               //
///////////////////////////////////////////////////////////////////

#include "config.h"
#include "constant.h"
#include "bamboo.h"

int main(int argc, char **argv){
	
	cout << endl;
	cout << "+------------------------+-------------------+---------------------+" << endl;
	cout << "|     Random Bamboo      |     " << setw(8) << RANDOM_BAMBOO_VERSION << "      |      09/09/2014     |" << endl;
	cout << "+------------------------+-------------------+---------------------+" << endl;
	cout << "|        (C) 2014 Han Zhang, GNU General Public License, V3        |" << endl;
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
	char *model = NULL;//b
	
	int ntree = 1;//t
	int mtry = 0;//m
	int seed = 0;//s
	int max_nleaf = 1000000;//l
	int min_leaf_size = 1;//r
	int imp_measure = 1;//i
	int nthread = 1;//h
	//if class_weight < 0, reweight to balance the data
	double class_weight = -1.0;//w
	double cutoff = 1.0;//u
	
	bool flip = true;//g
	bool output_prox = false;//x
	bool output_imp = true;//n
	bool balance = false;//B //if it is turned on and class_weight is not set (< .0), then class_weight by balancing the data (according to case/control ratio in data)
	bool trace = false;//r
	bool output_bamboo = true;//N
	
	bool file_spec = false;
	bool out_spec = false;
	bool pred_spec = false;
	bool model_spec = false;
	
	int c;
	while(true){
		static struct option long_options[] = {
			{"file", 1, NULL, 'f'},
			{"out", 1, NULL, 'o'},
			{"cont", 1, NULL, 'c'},
			{"cate", 1, NULL, 'a'},
			{"pred", 1, NULL, 'p'},
			{"bam", 1, NULL, 'b'},
			{"ntree", 1, NULL, 't'},
			{"mtry", 1, NULL, 'm'},
			{"seed", 1, NULL, 's'},
			{"maxnleaf", 1, NULL, 'l'},
			{"minleafsize", 1, NULL, 'e'},
			{"imp", 1, NULL, 'i'},
			{"nthread", 1, NULL, 'd'},
			{"classwt", 1, NULL, 'w'},
			{"cutoff", 1, NULL, 'u'},
			{"noflip", 0, NULL, 'g'},
			{"prox", 0, NULL, 'x'},
			{"noimp", 0, NULL, 'n'},
			{"balance", 0, NULL, 'B'},
			{"trace", 0, NULL, 'r'},
			{"nobamboo", 0, NULL, 'N'},
			{0, 0, 0, 0}
		};
		int option_index = 0;
		c = getopt_long(argc, argv, "f:o:c:a:p:b:t:m:s:l:e:i:d:w:u:gxnBrN", long_options, &option_index);
		
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
				model = optarg;
				model_spec = true;
				break;
			case 't':
				ntree = atoi(optarg);
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
			case 'u':
				cutoff = atof(optarg);
				if(cutoff <= .0){
					cout << "Error: The option --cutoff must be a positive number" << endl;
					return 0;
				}
				break;
			case 'g':
				flip = false;
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
			default:
				cout << "Error: Unknow options. Program terminates" << endl;
				exit(1);
		}
	}
	
	if(!file_spec){//no plink file is specified
		if(!pred_spec || !model_spec){
			cout << "Error: The option --file must be specified. Program terminates" << endl;
			return 0;
		}else{
			if(pred_spec && !model_spec){
				cout << "Error: The bamboo must be specified by option --bam for prediction purpose" << endl;
				return 0;
			}else if(!pred_spec && model_spec){
				cout << "Error: Please specify test dataset by option --pred" << endl;
				return 0;
			}else{
				;//good
			}
		}
	}else{//train model from specified data
		model = NULL;
		model_spec = false;
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
	
	if(!file_spec && pred_spec){
		BAMBOO bb(out, pred, model);
		bb.GrowForest();
		bb.PredictFromBamboo();
	}else{
		BAMBOO bb(file, out, cont, cate, pred, ntree, mtry, max_nleaf, min_leaf_size, imp_measure, seed, 
		nthread, class_weight, cutoff, flip, output_prox, output_imp, output_bamboo, balance, trace);
		bb.GrowForest();
		bb.PredictFromBamboo();
	}
	
	
	return 0;
	
}
