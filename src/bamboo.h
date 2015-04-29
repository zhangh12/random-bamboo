
#ifndef BAMBOO_H
#define BAMBOO_H

#include <omp.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/string.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/filesystem.hpp>
#include <iomanip>
#include <math.h>
#include <cmath>
#include <algorithm>
#include <assert.h>
#include <float.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>


using namespace std;
using namespace boost::filesystem;

typedef unsigned char uint8;
typedef unsigned short uint16;
typedef unsigned long long uint64;
typedef vector<uint64> bitvec;
typedef vector<bitvec > bitmat;

inline uint64_t rdtsc() {
    uint32_t lo, hi;
    __asm__ __volatile__ (
      "xorl %%eax, %%eax\n"
      "cpuid\n"
      "rdtsc\n"
      : "=a" (lo), "=d" (hi)
      :
      : "%ebx", "%ecx");
    return (uint64_t)hi << 32 | lo;
}

namespace boost {
namespace serialization {
class access;
}
}

class NODE{
	
	public: 
	vector<bitvec > sample_id64; // indicate that which samples are belong to the node
	vector<int> sample_id;
	bool large_sample_size;//if false, don't use the boolean operation for this node and its children
	
	vector<int> cont_lower;
	vector<int> cont_upper;
	int node_id; // id of the node
	int leaf_id; //id of the node, if it is a leaf_id
	int tree_id;
	
	//string split_snp;
	int split_snp_id; // the index of the snp indicator finally chosen to split the node.
	int split_snp_left_val; // subjects with SNP <= split_snp_left_val will fall into left. Range: 0, 1
	
	//string split_cont;
	int split_cont_id;
	int split_cont_lower;
	int split_cont_upper;
	double split_cont_thr;
	
	//string split_cate;
	int split_cate_id;
	vector<int> split_cate_left_code;
	vector<int> split_cate_right_code;
	
	int ncase_all;
	int nctrl_all;
	
	int ncase; // number of case in the node
	int nctrl; // number of control in the node
	int N; // total sample size in the node
	int split_by;
	bool processed; // indicate that whether the node is needed to be process
	bool terminal; // indicate that whether the node is set as terminal
	double split_stat; // the largest reduced gini if splitted
	int node_class;
	double node_prob[2]; //to save memory, consider to discard this mumber as it can be computed from node_risk_case and node_risk_ctrl
	double node_risk_case;//Pr(y & c), where y is case/control, c is a leaf
	double node_risk_ctrl;
	
	double class_weight;
	
	//if the node can be further split
	//the sample sizes in its children
	int ncase_left;
	int nctrl_left;
	int N_left;
	int ncase_right;
	int nctrl_right;
	int N_right;
	
	//some pointers
	NODE *parent;
	NODE *child1;//left child
	NODE *child2;//right child
	
	//if saved as model, the vector index
	int iparent;
	int ichild1;
	int ichild2;
	
	NODE(){
		
		parent = NULL;
		child1 = NULL;
		child2 = NULL;
	}

	NODE(const int input_ncase_all, const int input_nctrl_all, const vector<bitvec > &input_sample_id64, 
	const vector<int> &input_sample_id, const vector<int> &input_cont_lower, 
	const vector<int> &input_cont_upper, const int input_tree_id, const int input_node_id, 
	const bool input_processed, const bool input_terminal, const double input_class_weight, 
	const int input_ncase, const int input_nctrl, const int input_sample_threshold, 
	NODE *const par, NODE *const chd1, NODE *const chd2){
		
		ncase_all = input_ncase_all;
		nctrl_all = input_nctrl_all;
		
		sample_id64 = input_sample_id64;
		sample_id = input_sample_id;
		cont_lower = input_cont_lower;
		cont_upper = input_cont_upper;
		tree_id = input_tree_id;
		node_id = input_node_id;
		split_snp_id = -1;
		split_snp_left_val = -1;
		split_cont_id = -1;
		split_cont_lower = -1;
		split_cont_upper = -1;
		split_cate_id = -1;
		split_cate_left_code.clear();
		split_cate_right_code.clear();
		
		class_weight = input_class_weight;
		ncase = input_ncase;
		nctrl = input_nctrl;
		N = ncase + nctrl;
		large_sample_size = (N >= input_sample_threshold) ? true : false;
		
		processed = input_processed;
		terminal = input_terminal;
		leaf_id = -1;//for all non-leaf, set as -1
		split_stat = -1.0;
		
		//the sums of node_risk_case (or node_risk_ctrl) over all trees in the forest are used to produce the final determinant prediction of the whole forest
		//large sum means larger risk, and hence the prediction will be the other class
		//e.g. if sum of node_risk_case > sum of node_risk_ctrl, then predict is as control and versa vise.
		//see function CompOOBErrorSingleProc()
		node_risk_case = nctrl * 1.0 / nctrl_all;//Pr(y=0, A): risk of predicting as 1, proportional to P(y=0 | A), A is a leaf
		node_risk_ctrl = ncase * class_weight / ncase_all;//Pr(y=1, A): risk of predicting as 0, proportional to P(y=1 | A)
		
		//to get a probability prediction from the forest, we can average node_prob over all trees, 
		//or, calculate the sum of node_risk_case (or node_risk_ctrl) over all trees first, 
		//and then use sum of node_risk_ctrl / (sum of node_risk_ctrl + sum of node_risk_case) 
		//as the estimated posterior probability of a case
		//we will compare these two in later version
		node_class = (node_risk_case > node_risk_ctrl) ? 0 : 1;
		node_prob[0] = node_risk_case / (node_risk_case + node_risk_ctrl);//P(y=0 | A)
		node_prob[1] = 1.0 - node_prob[0];//P(y=1 | A)
		
		if(!ncase || !nctrl){
			terminal = true;
		}
		
		split_by = CODE_NULL;
		
		ncase_left = 0;
		nctrl_left = 0;
		N_left = 0;
		ncase_right = 0;
		nctrl_right = 0;
		N_right = 0;

		parent = par;
		child1 = chd1;
		child2 = chd2;
		
	}
	
};

class MINI_NODE{
	
	public:
	int node_id;
	//int leaf_id;
	int tree_id;
	
	int split_snp_id;
	int split_snp_left_val;
	
	int split_cont_id;
	double split_cont_thr;
	
	int split_cate_id;
	vector<int> split_cate_left_code;
	vector<int> split_cate_right_code;
	
	int split_by;
	double split_stat;
	bool terminal;
	
	//int node_class;
	double node_risk_case;
	double node_risk_ctrl;
	
	//double node_prob[2];//discard this member as it can be computed by node_risk_case and node_risk_ctrl
	
	int N;
	
	int iparent;
	int ichild1;
	int ichild2;
	
	MINI_NODE(){};
	
	MINI_NODE(const NODE& node){
		
		node_id = node.node_id;
		//leaf_id = node.leaf_id;
		tree_id = node.tree_id;
		
		split_snp_id = node.split_snp_id;
		split_snp_left_val = node.split_snp_left_val;
		
		split_cont_id = node.split_cont_id;
		split_cont_thr = node.split_cont_thr;
		
		split_cate_id = node.split_cate_id;
		split_cate_left_code = node.split_cate_left_code;
		split_cate_right_code = node.split_cate_right_code;
		
		split_by = node.split_by;
		split_stat = node.split_stat;
		terminal = node.terminal;
		//node_class = node.node_class;
		node_risk_case = node.node_risk_case;
		node_risk_ctrl = node.node_risk_ctrl;
		//node_prob[0] = node.node_prob[0];
		//node_prob[1] = node.node_prob[1];
		
		N = node.N;
		
		iparent = node.iparent;
		ichild1 = node.ichild1;
		ichild2 = node.ichild2;
		
	}
	
	MINI_NODE(const MINI_NODE& mini_node){
		
		node_id = mini_node.node_id;
		//leaf_id = mini_node.leaf_id;
		tree_id = mini_node.tree_id;
		
		split_snp_id = mini_node.split_snp_id;
		split_snp_left_val = mini_node.split_snp_left_val;
		
		split_cont_id = mini_node.split_cont_id;
		split_cont_thr = mini_node.split_cont_thr;
		
		split_cate_id = mini_node.split_cate_id;
		split_cate_left_code = mini_node.split_cate_left_code;
		split_cate_right_code = mini_node.split_cate_right_code;
		
		split_by = mini_node.split_by;
		split_stat = mini_node.split_stat;
		terminal = mini_node.terminal;
		//node_class = mini_node.node_class;
		node_risk_case = mini_node.node_risk_case;
		node_risk_ctrl = mini_node.node_risk_ctrl;
		//node_prob[0] = mini_node.node_prob[0];
		//node_prob[1] = mini_node.node_prob[1];
		
		N = mini_node.N;
		
		iparent = mini_node.iparent;
		ichild1 = mini_node.ichild1;
		ichild2 = mini_node.ichild2;
		
	}
	
	
	MINI_NODE& operator=(const MINI_NODE &mini_node){
		
		node_id = mini_node.node_id;
		//leaf_id = mini_node.leaf_id;
		tree_id = mini_node.tree_id;
		
		//split_snp = node.split_snp;
		split_snp_id = mini_node.split_snp_id;
		split_snp_left_val = mini_node.split_snp_left_val;
		
		//split_cont = node.split_cont;
		split_cont_id = mini_node.split_cont_id;
		split_cont_thr = mini_node.split_cont_thr;
		
		//split_cate = node.split_cate;
		split_cate_id = mini_node.split_cate_id;
		split_cate_left_code = mini_node.split_cate_left_code;
		split_cate_right_code = mini_node.split_cate_right_code;
		
		split_by = mini_node.split_by;
		split_stat = mini_node.split_stat;
		terminal = mini_node.terminal;
		//node_class = mini_node.node_class;
		node_risk_case = mini_node.node_risk_case;
		node_risk_ctrl = mini_node.node_risk_ctrl;
		//node_prob[0] = mini_node.node_prob[0];
		//node_prob[1] = mini_node.node_prob[1];
		
		N = mini_node.N;
		
		iparent = mini_node.iparent;
		ichild1 = mini_node.ichild1;
		ichild2 = mini_node.ichild2;
		
		return *this;
		
	}
	
	bool operator==(const MINI_NODE &mini_node) const;
	
	friend class BAMBOO;
	
	private:
	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive & ar, const unsigned int & version){
		
		ar & node_id;
		//ar & leaf_id;
		ar & tree_id;
		
		ar & split_snp_id;
		ar & split_snp_left_val;
		
		ar & split_cont_id;
		ar & split_cont_thr;
		
		ar & split_cate_id;
		ar & split_cate_left_code;
		ar & split_cate_right_code;
		
		ar & split_by;
		ar & split_stat;
		ar & terminal;
		//ar & node_class;
		ar & node_risk_case;
		ar & node_risk_ctrl;
		//ar & node_prob[0];
		//ar & node_prob[1];
		
		ar & iparent;
		ar & ichild1;
		ar & ichild2;
		
	}
		
};

struct CATE_CELL{
	
	int code;
	int ncase;
	int nctrl;
	int N;
	double r;
	
	bool operator < (const CATE_CELL& cc) const{
		return (r > cc.r);
	}
	
	CATE_CELL(const int &c, const int &n1, const int &n0, const double &r0){
		code = c;
		ncase = n1;
		nctrl = n0;
		N = ncase + nctrl;
		r = r0;
	}
	
};

struct IMPORTANCE{

	float imp_gini;
	float imp_raw;
	float imp_breiman_cutler;
	float imp_liaw_wiener;
	
	float imp_sd_breiman_cutler;
	float imp_sd_liaw_wiener;
	float sort_by;
	int var_id;

	bool operator < (const IMPORTANCE& im) const{
		return (sort_by > im.sort_by);
	}

	IMPORTANCE(const float &ig, const float &ir, const float &ibc, const float &ilw, const int &sb){
		imp_gini = ig;
		imp_raw = ir;
		imp_breiman_cutler = ibc;
		imp_liaw_wiener = ilw;
		
		imp_sd_breiman_cutler = .0f;
		imp_sd_liaw_wiener = .0f;
		sort_by = .0f;
		var_id = sb;
	}

};

struct FITTED_PROB{

	int sample_id;
	int true_class;
	double post_prob1;
	double post_prob2;
	double sort_by;
	
	bool operator < (const FITTED_PROB& auc) const{
		return (sort_by > auc.sort_by);
	}
	
	FITTED_PROB(const int &sid, const int &tc, const double &pp1, const double &pp2, const double &sb){
		sample_id = sid;
		true_class = tc;
		post_prob1 = pp1;
		post_prob2 = pp2;
		sort_by = sb;
	}
	
};


inline int BitCount(uint64 x){

	x -= (x >> 1) & 0x5555555555555555ULL;
	x = (x & 0x3333333333333333ULL) + ((x >> 2) & 0x3333333333333333ULL);
	x = (x + (x >> 4)) & 0x0f0f0f0f0f0f0f0fULL;
	x += x >> 8;
	x += x >> 16;
	x += x >> 32;
	return x & 0x7f;

}

inline bitvec operator&(const bitvec& bv1, const bitvec& bv2){

	bitvec bv(bv1.size(), 0);

	for(int i = 0; i < bv1.size(); ++i){
		bv[i] = bv1[i] & bv2[i];
	}

	return bv;

}



class BAMBOO{
	
	private:
		
		char *path_bed;
		char *path_bim;
		char *path_fam;
		char *path_cont;
		char *path_cate;
		char *path_trainid;
		char *path_testid;
		char *path_snpid;
		char *path_samplewt;
		
		char *path_prox;
		char *path_imp;
		char *path_bam;
		char *path_err;
		char *path_cof;
		char *path_pred;
		char *path_auc;
		
		char *path_fam_test;
		char *path_bed_test;
		char *path_bim_test;
		char *path_cont_test;
		char *path_cate_test;
		
		string dir;
		
		bool snp_major;
		
		int mtry;
		int nthread;
		int ntree;
		int seed;
		int max_nleaf;
		int min_leaf_size;
		int imp_measure;
		double class_weight;//weight assigned to case group
		double cutoff;
		
		bool flip;
		bool output_prox;
		bool output_imp;
		bool output_bamboo;
		bool keepcovar;
		bool tidy;
		bool has_cont;
		bool has_cate;
		bool has_test;
		bool has_samplewt;
		bool balance;
		bool trace;
		
		bool pred_from_trained_model;
		bool pred_from_specified_model;
		
		unsigned char wordbits[65536];
		
		string version;
		vector<vector<MINI_NODE> > bamboo; //store the trained model
		vector<string> snp_used_in_forest; //store the names of SNPs involved in the trained model
		vector<string> cont_var_used_in_forest; //store the names of continuous variables involved in the trained model
		vector<string> cate_var_used_in_forest; //store the names of categorical variables involved in the trained model
		vector<int> snp_id_used_in_forest; //store SNPs' ids
		vector<int> cont_var_id_used_in_forest; //store continuous variables' ids
		vector<int> cate_var_id_used_in_forest; //store categorical variables' ids
		
		int nsub; //sample size
		int nsnp; //number of SNPs
		int nsnp_specified; //number of SNPs to be used in training
		int ncase; //number of cases
		int nctrl; //number of controls
		int ncont; //number of continuous covariates
		int ncate; //number of categorical covariates
		int nvar; //nsnp + ncont + ncate
		int nflip; //number of SNPs flipped in training data
		
		int nsub_test;//sample size in test dataset
		int nflip_test;//number of SNPs flipped in testing data
		
		vector<string> specified_trainid;
		vector<string> specified_testid;
		vector<string> specified_snpid;
		
		vector<int> inc_snp_id;//include snp in training
		
		vector<string> individual_id;//the individual IDs of all samples involved in the final analysis (without missing phenotype and covariates)
		vector<string> individual_id_test;//the individual IDs of all samples involved in prediction (without missing covariates)
		
		vector<const char*> bam_list;//when predicting testing data from several specified models in local files, bam_list store the names of BAM files in the directory specified by path_bam
		
		vector<string> snp_name;
		vector<string> snp_name_test;
		vector<int> y; //original outcome
		vector<vector<double> > xcont; //original continue covariates
		vector<vector<string> > xcate; //original categorical covariates
		vector<vector<int> > xcate_int; //categorical covariates recoded according to cate_code and cate_unique below
		
		vector<string> cont_var_name;
		vector<vector<double> > cont_cutpnt;//cut points of continuous covariates
		vector<vector<double> > cont_unique;//unique values in continuous covariates
		vector<vector<int> > cont_loc;//for each subject in each continuous covariate, cont_loc indicates its group defined by the cut points
		
		vector<string> cate_var_name;
		vector<int> cate_start;//for each level of categorical covariates, index of the start point of column in xcate64
		vector<int> cate_end;//for each level of categorical covariates, index of the end point of column in xcate64
		vector<vector<string> > cate_unique;//unique values in categorical covariate
		vector<vector<int> > cate_code;//recode the levels of categorical covariate by their uniuqe values
		
		bitvec y64; //original outcome in 64-bit format
		bitvec not_y64;
		bitmat geno64;
		bitmat not_geno64;
		bitmat xcate64;
		
		bitmat geno64_test;
		vector<vector<double> > xcont_test;
		vector<vector<string> > xcate_test;
		vector<vector<int> > xcate_int_test;
		
		vector<double> samplewt;
		
		bitmat oob_id64;//out-of-bag id, ntree x nblock
		bitmat ib_id64;//in-bag id, ntree x nblock
		
		vector<int> num_oob;//number of oob in each of trees. len = ntree
		
		int LEN_bed;
		uint8 PROBE1[4];
		uint8 PROBE2[4];
		
		vector<int> bitloc;
		bitvec MASK_offset;
		
		vector<int> bitloc_bed;
		bitvec MASK_offset_bed;
		
		vector<int> bitloc_test;
		bitvec MASK_offset_test;
		
		vector<int> bitloc_var;
		bitvec MASK_offset_var;
		
		vector<int> bitloc_tree;
		bitvec MASK_offset_tree;
		
		uint64 MASK;
		int LEN;
		int nblock;
		int nblock_test;
		int nblock_var;
		int nblock_tree;
		
		vector<bool> check_out;//indicate whether a tree has been created. used when estimate the progress
		vector<bool> check_in;
		//vector<vector<vector<double> > > model;
		bitmat var_id_used_in_tree_long;//indicate whether a var is used in a bamboo of forest, nvar x nblock_tree, each row is a snp/cont/cate
		bitmat var_id_used_in_tree_wide;//as above, but with size ntree x nblock_var, each row is a bamboo
		vector<int> num_var_in_tree;//number of used variables in each of trees. len = ntree
		
		
		vector<vector<uint16> > pred_node_id; //the node id of the leaf that a sample will fall in, predicted by a tree, ntree X nsub
		vector<double> oob_error;
		vector<vector<uint16> > prox;//nsub x nsub matrix for pair proximity
		
		
		vector<IMPORTANCE> importance;
		vector<int> num_vote_correct_class;//number of correct votes of each of trees. computed from oob only. len = ntree
		vector<vector<uint16> > num_vote_loss_permuted_correct_class;//each row corresponds to a tree. lengths of rows vary
		vector<int> oob_prediction;
		vector<vector<int> > confusion_matrix;
		
		vector<double> post_prob_case1_test;//posterior probability that a sample testing data is a case. Computed by definition 1
		vector<double> post_prob_case2_test;//posterior probability that a sample testing data is a case. Computed by definition 2
		
		vector<double> fitted_risk_ctrl;//risk as a control (oob) in training data, used in computing AUC
		vector<double> fitted_risk_case;//risk as a case (oob) in training data, used in computed AUC
		vector<FITTED_PROB> fitted_prob;//as above, but for training data (oob only)
		double auc;
		
		inline int PopCount(uint64 x){ return wordbits[x & 0xFFFF] + wordbits[(x >> 16) & 0xFFFF] + wordbits[(x >> 32) & 0xFFFF] + wordbits[x >> 48]; }
		
		void ComputeWordBits();
		void ParseModelPath(const char *const);
		void EncodeCont();
		void EncodeCate();
		void PrintPara();
		void LoadTrainingData();
		void CreateFile(char*, const char *const);
		void Walker_ProbSampleReplace(drand48_data&, const vector<double>&, const int, vector<int>&);
		void Bootstrap(const int, drand48_data&, bitmat&, bitmat&, bitmat&, vector<vector<int> >&, vector<int>&, vector<int>&, vector<int>&, int&);
		void ShuffleSNP(drand48_data&, vector<int>&);
		void ShuffleAllFeature(drand48_data&, vector<int>&, vector<int>&, vector<int>&);
		void ShuffleSNPKeepAllCovar(drand48_data&, vector<int>&, vector<int>&, vector<int>&);
		void Shuffle(drand48_data&, vector<int>&, vector<int>&, vector<int>&);
		bool SplitNode(NODE&, const bitmat&, const bitmat&, const vector<int>&, const vector<vector<int> >&, drand48_data&);
		void DestroyTree(vector<NODE*>&);
		void GrowTree(drand48_data&, const int);
		void PrintProgress(const time_t);
		void GrowForestMultiProc(const time_t);
		void SaveBamboo();
		void PutDownTrainingSampleToTree(const int, const int, int&);
		void PredictTrainingSample();
		void ComputeProximity();
		void WriteOOBError();
		void ComputeOOBError();
		void WriteOOBAUC();
		void ComputeOOBAUC();
		void WriteConfusionMatrix();
		void ComputeConfusionMatrix();
		void ComputeGiniImportance();
		void FisherYatesShuffle(drand48_data&, vector<int>&);
		void PutDownOOBSampleToTree(const int, const int, const int, const int, int&);
		void ComputePermutationImportance();
		void WriteImportance();
		void ComputeImportance();
		void LoadBamboo();
		void LoadTestingData();
		void PutDownTestingSampleToTree(const int, const int, int&);
		void PutDownTestingSampleToTree(const int, const vector<MINI_NODE>&, int&);
		void PredictTestingSampleFromSingleForest();
		void ScanMultipleForest();
		void PredictTestingSampleFromMultipleForest();
		void SavePrediction();
		
	public:
		
		BAMBOO(){};
		BAMBOO(const char *const, const int);
		BAMBOO(const char *const, const char *const, const char * const, const char * const, 
		const char * const, const char *const, const char *const, const char *const, const char *const, 
		const int, const int, const int, const int, const int, const int, const int, const double, 
		const double, const bool, const bool, const bool, const bool, const bool, const bool, 
		const bool, const bool);
		BAMBOO(const char *const, const char *const, const char *const, const char *const, const int, 
		const int);
		~BAMBOO();
		
		void PrintModelInfo();
		void GrowForest();
		void PredictTestingSample();
		
		friend bitvec operator&(const bitvec&, const bitvec&);
		friend ostream& operator<<(ostream&, const vector<double>&);
	
};


#endif

