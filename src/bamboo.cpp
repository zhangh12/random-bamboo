
#include "constant.h"
#include "bamboo.h"


CATE_CELL cate_cell0(-1, 0, 0, .0);

IMPORTANCE imp0(.0f, .0f, .0f, .0f, -1);

FITTED_PROB fitpr0(-1, -1, -1.0, -1.0, -1.0);

ostream& operator<<(ostream &os, const vector<double> &v){

	for(int i = 0; i < v.size(); ++i){
		os << v[i] << endl;
	}

	return os;

}


bool MINI_NODE::operator==(const MINI_NODE &mini_node) const {
	
	if(node_id != mini_node.node_id){
		return false;
	}
	
//	if(leaf_id != mini_node.leaf_id){
//		return false;
//	}
	
	if(tree_id != mini_node.tree_id){
		return false;
	}
	
	if(split_snp_id != mini_node.split_snp_id){
		return false;
	}
	
	if(split_snp_left_val != mini_node.split_snp_left_val){
		return false;
	}
	
	if(split_cont_id != mini_node.split_cont_id){
		return false;
	}
	
	if((split_cont_thr - mini_node.split_cont_thr) > 1e-12){
		return false;
	}
	
	if(split_cate_id != mini_node.split_cate_id){
		return false;
	}
	
	if(split_cate_left_code.size() != mini_node.split_cate_left_code.size()){
		return false;
	}
	if(split_cate_right_code.size() != mini_node.split_cate_right_code.size()){
		return false;
	}
	
	for(int i = 0; i < split_cate_left_code.size(); ++i){
		if(split_cate_left_code[i] != mini_node.split_cate_left_code[i]){
			return false;
		}
	}
	
	for(int i = 0; i < split_cate_right_code.size(); ++i){
		if(split_cate_right_code[i] != mini_node.split_cate_right_code[i]){
			return false;
		}
	}
	
	if(split_by != mini_node.split_by){
		return false;
	}
	
	if(split_stat != mini_node.split_stat){
		return false;
	}
	
	if(terminal != mini_node.terminal){
		return false;
	}
	
//	if(node_class != mini_node.node_class){
//		return false;
//	}
	
	if(abs(node_risk_case - mini_node.node_risk_case) > 1e-12){
		return false;
	}
	
	if(abs(node_risk_ctrl - mini_node.node_risk_ctrl) > 1e-12){
		return false;
	}
	
//	if(abs(node_prob[0] - mini_node.node_prob[0]) > 1e-12){
//		return false;
//	}
//	
//	if(abs(node_prob[1] - mini_node.node_prob[1]) > 1e-12){
//		return false;
//	}
	
	if(N != mini_node.N){
		return false;
	}
	
	if(iparent != mini_node.iparent){
		return false;
	}
	
	if(ichild1 != mini_node.ichild1){
		return false;
	}
	
	if(ichild2 != mini_node.ichild2){
		return false;
	}
	
	return true;
	
}


BAMBOO::~BAMBOO(){
	
	delete[] path_bed;
	delete[] path_bim;
	delete[] path_fam;
	
	if(path_prox){
		delete[] path_prox;
	}
	if(path_imp){
		delete[] path_imp;
	}
	if(path_bam){
		delete[] path_bam;
	}
	
	if(path_cont){
		delete[] path_cont;
	}
	if(path_cate){
		delete[] path_cate;
	}
	
	if(path_trainid){
		delete[] path_trainid;
	}
	
	if(path_testid){
		delete[] path_testid;
	}
	
	if(path_snpid){
		delete[] path_snpid;
	}
	
	if(path_samplewt){
		delete[] path_samplewt;
	}
	
	delete[] path_err;
	delete[] path_cof;
	delete[] path_auc;
	
	if(path_bed_test){
		delete[] path_bed_test;
	}
	if(path_bim_test){
		delete[] path_bim_test;
	}
	if(path_fam_test){
		delete[] path_fam_test;
	}
	if(path_cont_test){
		delete[] path_cont_test;
	}
	if(path_cate_test){
		delete[] path_cate_test;
	}
	
	if(path_pred){
		delete[] path_pred;
	}
	
}

void BAMBOO::ComputeWordBits(){
	
	for(int i = 0; i < 65536; ++i){
		wordbits[i] = BitCount(i);
	}
	
}

//find cut points of each continuous covariate
//find unique values of each continuous covariate
void BAMBOO::EncodeCont(){
	
	cont_cutpnt = vector<vector<double> > (ncont);
	cont_unique = vector<vector<double> > (ncont);
	cont_loc = vector<vector<int> > (ncont, vector<int> (nsub, 0));
	for(int j = 0; j < ncont; ++j){
		vector<double> v = xcont[j];
		sort(v.begin(), v.end());
		vector<double> u;
		u.push_back(v[0]);
		for(int i = 1; i < nsub; ++i){
			if(abs(v[i] - u[u.size()-1]) > 1e-6){
				u.push_back(v[i]);
			}
		}
		
		vector<double> c;
		for(int i = 0; i < u.size()-1; ++i){
			c.push_back((u[i] + u[i+1]) / 2.0);
		}
		cont_cutpnt[j] = c;
		cont_unique[j] = u;
		for(int k = 0; k < nsub; ++k){
			for(int i = 0; i < cont_cutpnt[j].size(); ++i){
				if(xcont[j][k] > cont_cutpnt[j][i]){
					++cont_loc[j][k];
				}
			}
		}
	}
	
	if(false){
		vector<double> sum_cont(ncont, .0);
		for(int j = 0; j < ncont; ++j){
			for(int i = 0; i < nsub; ++i){
				sum_cont[j] += xcont[j][i];
			}
			cout << "cont " << j+1 << " (" << cont_var_name[j] << "): " << sum_cont[j] << endl;
		}
	}
	
}

//find unique values of each categorical covariate
//recode categorical covariate into integers
//find locations of each categorical covariate in 64-bit 2D vectors
void BAMBOO::EncodeCate(){
	
	cate_start = vector<int> (ncate, -1);
	cate_end = vector<int> (ncate, -1);
	cate_unique = vector<vector<string> > (ncate);
	cate_code = vector<vector<int> > (ncate);
	
	for(int j = 0; j < ncate; ++j){
		vector<string> v = xcate[j];
		sort(v.begin(), v.end());
		vector<string> u;
		u.push_back(v[0]);
		for(int i = 1; i < nsub; ++i){
			if(!(v[i] == u[u.size()-1])){
				u.push_back(v[i]);
			}
		}
		
		vector<int> c;
		for(int i = 0; i < u.size(); ++i){
			c.push_back(i);
		}
		cate_code[j] = c;
		cate_unique[j] = u;
	}
	
	cate_start[0] = 0;
	cate_end[0] = cate_code[0].size() - 1;
	for(int j = 1; j < ncate; ++j){
		cate_start[j] = cate_end[j-1]+1;
		cate_end[j] = cate_start[j] + cate_code[j].size()-1;
	}
	
	xcate64 = bitmat (cate_end[ncate-1]+1, bitvec(nblock, (uint64) 0));
	xcate_int = vector<vector<int> > (ncate, vector<int> (nsub, -1));
	for(int j = 0; j < ncate; ++j){
		for(int k = 0; k < nsub; ++k){
			string s = xcate[j][k];
			int co = 0;
			for(int l = 0; l < cate_unique[j].size(); ++l){
				if(!(s == cate_unique[j][l])){
					++co;
				}else{
					break;
				}
			}
			if(co >= cate_unique[j].size()){
				cout << "Error: Failed in transforming categorical covariate " << j+1 << endl;
				exit(1);
			}
			xcate64[cate_start[j]+co][bitloc[k]] |= MASK_offset[k];
			xcate_int[j][k] = co;
		}
	}
	
	if(false){
		cout << "encoding rule" << endl;
		for(int j = 0; j < ncate; ++j){
			cout << "cate " << j+1 << " (" << cate_var_name[j] << "): " << endl;
			for(int k = 0; k < cate_unique[j].size(); ++k){
				int c = 0;
				for(int i = 0; i < nblock; ++i){
					c += PopCount(xcate64[cate_start[j] + cate_code[j][k]][i]);
				}
				cout << cate_unique[j][k] << " (" << cate_code[j][k] << "): " << c << "\t";
			}
			cout << endl;
		}
		cout << "-------" << endl;
	}
	
}


void BAMBOO::PrintPara(){
	
	cout << nsnp << " markers are loaded in SNP-major mode from [ " << path_bed << " ]" << endl;
	
	if(has_cont){
		cout << ncont << " continuous and/or ordinary covariates are loaded from [ " << path_cont << " ]" << endl;
	}
	
	if(has_cate){
		cout << ncate << " categorical covariates are loaded from [ " << path_cate << " ]" << endl;
	}
	
	cout << nsub << " individuals without missing phenotype and covariates are included" << endl;
	cout << ncase << " cases, " << nctrl << " controls and 0 missing" << endl;
	
	if(balance){
		cout << "The data is automatically balanced with weights 1 : " << class_weight << " (contol : case)" << endl;
	}else{
		if(class_weight > .0 && abs(class_weight - 1.0) > 1e-12){
			cout << "The data is adjusted with specified weights 1 : " << class_weight << " (control : case)" << endl;
		}else{
			cout << "No reweighting is applied even if the data is imbalanced" << endl;
		}
	}
	
	switch(imp_measure){
		case IMP_GINI:
			cout << "Gini importance will be calculated" << endl;
			break;
		case IMP_BREIMAN_CUTLER:
			cout << "Gini and permutation importances will be calculated" << endl;
			break;
		default:
			cout << "Error: Unknown value of option --imp" << endl;
			exit(1);
	}
	
	cout << ntree << " bamboo(s) are being grown with " << nthread << " thread(s)" << endl;
	cout << "Each bamboo is grown with mtry " << mtry << endl;
	
}


void BAMBOO::LoadTrainingData(){
	
	//loading individual ID from .fam
	//if .con and/or .cat files are available, loading their individual ID as well
	//only the subject without missing covariates/phenotype are included in the analysis
	//the following code determines the intersection set of individual ID
	
	cout << "Loading training data ..." << endl;
	
	if(!specified_trainid.empty()){
		specified_trainid.clear();
	}
	
	if(path_trainid){
		ifstream file_trainid(path_trainid);
		if(!file_trainid){
			cout << "Error: Cannot open IID file " << path_trainid << endl;
			exit(1);
		}
		
		for(string s; getline(file_trainid, s); ){
			istringstream sin(s);
			string iid;
			sin >> iid;
			specified_trainid.push_back(iid);
		}
		file_trainid.close();
		
		if(specified_trainid.size() == 0){
			cout << "Error: 0 individuals are specified in " << path_trainid << endl;
			exit(1);
		}else{
			cout << specified_trainid.size() << " individuals specified in [ " << path_trainid << " ] are used as training data" << endl;
		}
	}
	
	//loading individual ID with weight. Only the individual with a weight > 0 are used in training
	vector<string> trainid_with_weight;
	if(path_samplewt){
		ifstream file_samplewt1(path_samplewt);
		if(!file_samplewt1){
			cout << "Error: Cannot open SWT file " << path_samplewt << endl;
			exit(1);
		}
		
		for(string s; getline(file_samplewt1, s); ){
			istringstream sin(s);
			string ind_id;
			double swt;
			if(sin >> ind_id >> swt){
				if(swt > .0){
					trainid_with_weight.push_back(ind_id);
				}
			}else{
				cout << "Error: Invalid format in " << path_samplewt << endl;
				exit(1);
			}
		}
		file_samplewt1.close();
		if(trainid_with_weight.size() > 0){
			cout << "Unequal weights are loaded from [ " << path_samplewt << " ] for bootstrapping" << endl;
		}
	}
	
	//loading FAM file
	ifstream file_fam(path_fam);
	if(!file_fam){
		cout << "Error: Cannot open FAM file " << path_fam << endl;
		exit(1);
	}
	
	vector<string> individual_id_fam;//individual id without missing phenotype
	vector<int> index_individual_id_fam;
	int nrow = -1;
	int nsub_fam = 0;
	for(string s; getline(file_fam, s); ){//read by lines
		
		istringstream sin(s);
		string fam_id, ind_id, pat_id, mat_id, sex;
		int phen;
		++nrow;
		if(sin >> fam_id >> ind_id >> pat_id >> mat_id >> sex >> phen){//only six columns in .fam
			++nsub_fam;
			if(phen == 1 || phen == 2){//non-missing phenotype
				
				bool in_spec = true;
				if(path_trainid){
					in_spec = false;
					for(int i = 0; i < specified_trainid.size(); ++i){
						if(ind_id == specified_trainid[i]){
							in_spec = true;
							break;
						}
					}
				}
				
				if(in_spec && path_samplewt){
					in_spec = false;
					for(int i = 0; i < trainid_with_weight.size(); ++i){
						if(ind_id == trainid_with_weight[i]){
							in_spec = true;
							break;
						}
					}
				}
				
				if(in_spec){
					individual_id_fam.push_back(ind_id);
					index_individual_id_fam.push_back(nrow);
				}
				
			}
		}else{
			cout << "Error: Invalid format in " << path_fam << endl;
			exit(1);
		}
	}
	file_fam.close();
	
	//loading continuous covariates, if available
	vector<string> individual_id_cont;//individual id without missing continuous covariates
	vector<int> index_individual_id_cont;
	int nsub_cont = 0;
	
	if(path_cont){
		ifstream file_cont(path_cont);
		if(!file_cont){
			cout << "Error: Cannot open CON file " << path_cont << endl;
			exit(1);
		}
		
		int nrow = -2;//make sure that index_individual_id_cont starts from 0
		ncont = -1;//the first column should be individual id
		string missing_flag("NA");
		for(string s; getline(file_cont, s); ){
			istringstream sin(s);
			string header;
			++nrow;
			if(nrow == -1){//the first row should be header, we use it to count number of continuous covariates
				while(sin >> header){
					++ncont;
				}
				if(ncont <= 0){//empty file or only one column of individual id
					has_cont = false;
					ncont = 0;
					break;
				}
			}else{
				++nsub_cont;
				bool no_missing_this_row = true;
				string iid, str;
				sin >> iid;//store the individual id
				while(sin >> str){//if this row has missing entries
					if(str == missing_flag){
						no_missing_this_row = false;
						break;
					}
				}
				if(no_missing_this_row){
					individual_id_cont.push_back(iid);
					index_individual_id_cont.push_back(nrow);
				}
			}
		}
		file_cont.close();
	}else{
		has_cont = false;
		ncont = 0;
	}
	
	//load categorical covariates, if available
	vector<string> individual_id_cate;//individual id without missing categorical covariates
	vector<int> index_individual_id_cate;
	int nsub_cate = 0;
	
	if(path_cate){
		ifstream file_cate(path_cate);
		if(!file_cate){
			cout << "Error: Cannot open CAT file " << path_cate << endl;
			exit(1);
		}
		
		int nrow = -2;//make sure that index_individual_id_cate starts from 0
		ncate = -1;//the first column should be individual id
		string missing_flag("NA");
		for(string s; getline(file_cate, s); ){
			istringstream sin(s);
			string header;
			++nrow;
			if(nrow == -1){//the first row should be header, we use it to count number of categorical covariates
				while(sin >> header){
					++ncate;
				}
				if(ncate <= 0){//empty file or only one column of individual id
					has_cate = false;
					ncate = 0;
					break;
				}
			}else{
				++nsub_cate;
				bool no_missing_this_row = true;
				string iid, str;
				sin >> iid;//store the individual id
				while(sin >> str){
					if(str == missing_flag){
						no_missing_this_row = false;
						break;
					}
				}
				if(no_missing_this_row){
					individual_id_cate.push_back(iid);
					index_individual_id_cate.push_back(nrow);
				}
			}
		}
		file_cate.close();
	}else{
		has_cate = false;
		ncate = 0;
	}
	
	//determine the intersection of data files without missing entries
	vector<int> index_individual_id_fam_inc;
	vector<int> index_individual_id_cont_inc;
	vector<int> index_individual_id_cate_inc;
	
	if(ncont == 0 && ncate == 0){//only the genotypes are available
		individual_id = individual_id_fam;
		index_individual_id_fam_inc = index_individual_id_fam;
	}else{
		if(ncate == 0){//genotype and continuous covariates are available
			string str;
			for(int i = 0; i < individual_id_fam.size(); ++i){
				str = individual_id_fam[i];
				for(int j = 0; j < individual_id_cont.size(); ++j){
					if(str == individual_id_cont[j]){
						individual_id.push_back(str);
						index_individual_id_fam_inc.push_back(index_individual_id_fam[i]);
						index_individual_id_cont_inc.push_back(index_individual_id_cont[j]);
						break;
					}
				}
			}
		}else if(ncont == 0){//genotype and categorical covariates are available
			string str;
			for(int i = 0; i < individual_id_fam.size(); ++i){
				str = individual_id_fam[i];
				for(int j = 0; j < individual_id_cate.size(); ++j){
					if(str == individual_id_cate[j]){
						individual_id.push_back(str);
						index_individual_id_fam_inc.push_back(index_individual_id_fam[i]);
						index_individual_id_cate_inc.push_back(index_individual_id_cate[j]);
						break;
					}
				}
			}
		}else{//genotype, continuous and categorical covarates are available
			string str;
			for(int i = 0; i < individual_id_fam.size(); ++i){
				str = individual_id_fam[i];
				bool in_cont = false;
				int j1 = -1;
				for(int j = 0; j < individual_id_cont.size(); ++j){
					if(str == individual_id_cont[j]){
						in_cont = true;
						j1 = j;
						break;
					}
				}
				
				if(!in_cont){
					continue;
				}
				
				bool in_cate = false;
				int j2 = -1;
				for(int j = 0; j < individual_id_cate.size(); ++j){
					if(str == individual_id_cate[j]){
						in_cate = true;
						j2 = j;
						break;
					}
				}
				
				if(in_cont && in_cate){
					individual_id.push_back(str);
					index_individual_id_fam_inc.push_back(index_individual_id_fam[i]);
					index_individual_id_cont_inc.push_back(index_individual_id_cont[j1]);
					index_individual_id_cate_inc.push_back(index_individual_id_cate[j2]);
				}
			}
		}
		
	}
	
	//sample size and number of SNPs
	
	nsub = individual_id.size();
	if(nsub == 0){
		cout << "Sample size is 0. The program terminated" << endl;
		exit(1);
	}
	
	ifstream file_bim(path_bim);
	if(!file_bim){
		cout << "Error: Cannot open BIM file " << path_bim << endl;
		exit(1);
	}
	
//	char c;
//	nsnp = 0;
//	while(file_bim.get(c)){
//		if(c == '\n'){
//			++nsnp;
//		}
//	}
//	file_bim.close();
	
	nsnp = 0;
	for(string s; getline(file_bim, s); ){//read by lines
		istringstream sin(s);
		string chr, rs, allele_name1, allele_name2;
		int genetic_dist, base_pair_pos;
		if(sin >> chr >> rs >> genetic_dist >> base_pair_pos >> allele_name1 >> allele_name2){//only six columns in .bim
			snp_name.push_back(rs);
			++nsnp;
//			if(base_pair_pos >= 0){//keep this SNP
//				snp_name.push_back(rs);
//				++nsnp;
//			}else{//discard this SNP
//				//do something later
//			}
		}else{
			cout << "Error: Invalid format in " << path_bim << endl;
			exit(1);
		}
	}
	file_bim.close();
	
	nvar = nsnp + ncont + ncate;
	
	
	if(path_snpid){
		ifstream file_snpid(path_snpid);
		if(!file_snpid){
			cout << "Error: Cannot open SID file " << path_snpid << " to load SNP list" << endl;
			exit(1);
		}
		
		for(string s; getline(file_snpid, s); ){
			istringstream sin(s);
			string sid;
			sin >> sid;
			specified_snpid.push_back(sid);
		}
		file_snpid.close();
		
		if(specified_snpid.size() == 0){
			cout << "Warning: No SNP is specified in file " << path_snpid << ". All SNPs will be used in training" << endl;
			specified_snpid = snp_name;
			nsnp_specified = nsnp;
			inc_snp_id = vector<int> (nsnp, -1);
			for(int i = 0; i < nsnp; ++i){
				inc_snp_id[i] = i;
			}
		}else{
			nsnp_specified = specified_snpid.size();
			cout << nsnp_specified << " markers in [ " << path_snpid << " ] are used in training" << endl;
			
			vector<bool> tmp (nsnp, false);
			
			#pragma omp parallel num_threads(nthread)
			{
				#pragma omp for
				for(int i = 0; i < nsnp; ++i){
					string str = snp_name[i];
					for(int j = 0; j < nsnp_specified; ++j){
						if(str == specified_snpid[j]){
							tmp[i] = true;
							break;
						}
					}
				}
			}
			
			inc_snp_id.clear();
			for(int i = 0; i < nsnp; ++i){
				if(tmp[i]){
					inc_snp_id.push_back(i);
				}
			}
			
			if(inc_snp_id.size() == 0){
				cout << "Error: All markers in [ " << path_snpid << " ] are not found in the BED file" << endl;
				exit(1);
			}
			
			specified_snpid.clear();
			for(int i = 0; i < inc_snp_id.size(); ++i){
				specified_snpid.push_back(snp_name[inc_snp_id[i]]);
			}
			nsnp_specified = specified_snpid.size();
			
		}
		
	}else{
		specified_snpid = snp_name;
		nsnp_specified = nsnp;
		inc_snp_id = vector<int> (nsnp, -1);
		for(int i = 0; i < nsnp; ++i){
			inc_snp_id[i] = i;
		}
	}
	
	////
	
	//load genotypes by memory mapping
	stringstream file_bed;
	file_bed << path_bed;
	int fb = open(file_bed.str().c_str(), O_RDONLY);
	if(fb < 0){
		cout << "Error: Cannot open BED file " << file_bed.str().c_str() << endl;
		close(fb);
		exit(1);
	}
	
	LEN_bed = 8;
	int nblock_bed = (int) ceil((2 * (double) nsub_fam) / LEN_bed);
	
	uint8 *map;
	int file_size = (3 + nblock_bed * nsnp) * sizeof(uint8);
	map = (uint8*) mmap(0, file_size, PROT_READ, MAP_PRIVATE, fb, 0);
	if(map == MAP_FAILED){
		cout << "Error: Failed in memory mapping on the BED file " << file_bed.str().c_str() << endl;
		cout << strerror(errno) << " (" << errno << ")" << endl;
		close(fb);
		exit(1);
	}
	
	if(map[0] != 0x6c || map[1] != 0x1b){//test the magic number in BED file
		cout << "Error: " << file_bed.str().c_str() << " is not a plink BED file" << endl;
		munmap(map, file_size);
		close(fb);
		exit(1);
	}
	
	if(map[2] == 0x01){
		snp_major = true;
	}else if(map[2] == 0x00){
		snp_major = false;
		cout << "Error: bamboo only accepts BED file with SNP-major mode" << endl;
		munmap(map, file_size);
		close(fb);
		exit(1);
	}else{
		cout << "Error: Illegal BED flag" << endl;
		munmap(map, file_size);
		close(fb);
		exit(1);
	}
	
	PROBE1[0] = 0x01;
	PROBE2[0] = 0x02;
	PROBE1[1] = 0x04;
	PROBE2[1] = 0x08;
	PROBE1[2] = 0x10;
	PROBE2[2] = 0x20;
	PROBE1[3] = 0x40;
	PROBE2[3] = 0x80;
	
	
	MASK = 0x8000000000000000;
	LEN = 64;
	
	/////////
	
	nblock_var = (int) ceil (((double) nvar) / LEN);
	bitloc_var.reserve(nvar);
	MASK_offset_var.reserve(nvar);
	for(int k = 0; k < nvar; ++k){
		bitloc_var.push_back(k / LEN);
		MASK_offset_var.push_back(MASK >> (k % LEN));
	}
	
	nblock_tree = (int) ceil (((double) ntree) / LEN);
	bitloc_tree.reserve(ntree);
	MASK_offset_tree.reserve(ntree);
	for(int k = 0; k < ntree; ++k){
		bitloc_tree.push_back(k / LEN);
		MASK_offset_tree.push_back(MASK >> (k % LEN));
	}
	
	var_id_used_in_tree_long = bitmat(nvar, bitvec(nblock_tree, (uint64) 0));
	var_id_used_in_tree_wide = bitmat(ntree, bitvec(nblock_var, (uint64) 0));
	num_var_in_tree = vector<int> (ntree, 0);
	
	//////
	
	nblock = (int) ceil (((double) nsub) / LEN);
	bitloc.reserve(nsub);
	MASK_offset.reserve(nsub);
	
	bitloc_bed.reserve(LEN * nblock);
	MASK_offset_bed.reserve(LEN * nblock);
	for(int k = 0; k < nsub; ++k){
		bitloc.push_back(k / LEN); //indicate the block No. where the k-th sample locates in
		MASK_offset.push_back(MASK >> (k % LEN)); //make a mask for indicating the offset in a block with length LEN
		
		bitloc_bed.push_back(k / LEN);
		MASK_offset_bed.push_back(MASK >> (k % LEN));
	}
	
	for(int k = nsub; k < LEN * nblock; ++k){
		bitloc_bed.push_back(k / LEN);
		MASK_offset_bed.push_back(MASK >> (k % LEN));
	}
	
	geno64 = bitmat (nsnp * 3, bitvec(nblock, (uint64) 0));
	not_geno64 = bitmat (nsnp * 3, bitvec(nblock, (uint64) 0));
	
	#pragma omp parallel num_threads(nthread)
	{
		#pragma omp for
		for(int i = 0; i < nsnp; ++i){
			int k = -1;
			int t = 0;
			for(int j = 0; j < nblock_bed; ++j){
				uint8 b = map[3 + i * nblock_bed + j];
				for(int l = 0; l < 4; ++l){
					++k;
					if(k < nsub_fam){
						if(k == index_individual_id_fam_inc[t]){//this row doesn't have missing entries
							
							uint8 h1 = b & PROBE1[l];
							uint8 h2 = b & PROBE2[l];
							if(h1 && !h2){
								cout << "Error: bamboo does not allow missing genotypes" << endl;
								munmap(map, file_size);
								close(fb);
								exit(1);
							}
							
							if(!h1 && !h2){// 0/0
								geno64[i * 3][bitloc[t]] |= MASK_offset[t];
							}else{
								not_geno64[i * 3][bitloc[t]] |= MASK_offset[t];
							}
							
							if(!h1 && h2){// 0/1
								geno64[i * 3 + 1][bitloc[t]] |= MASK_offset[t];
							}else{
								not_geno64[i * 3 + 1][bitloc[t]] |= MASK_offset[t];
							}
							
							if(h1 && h2){// 1/1
								geno64[i * 3 + 2][bitloc[t]] |= MASK_offset[t];
							}else{
								not_geno64[i * 3 + 2][bitloc[t]] |= MASK_offset[t];
							}
							
							++t;
							if(t >= nsub){
								break;
							}
						}else{//row with missing, ignore
							continue;
						}
					}else{//all samples have been scanned
						break;
					}
				}
			}
			
		}
	}
	
	munmap(map, file_size);
	close(fb);
	
	nflip = 0;
	if(flip){//test if the genotypes have been correctly loaded
		
		vector<bool> flipped(nsnp, false);
		#pragma omp parallel num_threads(nthread)
		{
			#pragma omp for
			for(int i = 0; i < nsnp; ++i){
				int n0 = 0, n1 = 0, n2 = 0;
				for(int k = 0; k < nblock; ++k){
					n0 += PopCount(geno64[i*3][k]);
					n2 += PopCount(geno64[i*3+2][k]);
				}
				if(n0 < n2){// flip the genotype if necessary
					flipped[i] = true;
					for(int k = 0; k < nblock; ++k){
						geno64[i*3][k] = (~(geno64[i*3][k])) & (~(geno64[i*3+1][k]));
						geno64[i*3+2][k] = (~(geno64[i*3+2][k])) & (~(geno64[i*3+1][k]));
					}
					for(int k = nsub; k < LEN * nblock; ++k){
						geno64[i*3][bitloc_bed[k]] &= ~(MASK_offset_bed[k]);
						geno64[i*3+2][bitloc_bed[k]] &= ~(MASK_offset_bed[k]);
					}
				}
				
				if(false){
					n0 = 0, n1 = 0, n2 = 0;
					for(int k = 0; k < nblock; ++k){
						n0 += PopCount(geno64[i*3][k]);
						n1 += PopCount(geno64[i*3+1][k]);
						n2 += PopCount(geno64[i*3+2][k]);
					}
					cout << "SNP " << i+1 << " " << n0 << "\t" << n1 << "\t" << n2 << endl;
				}
			}
		}
		
		for(int i = 0; i < nsnp; ++i){
			if(flipped[i]){
				nflip++;
			}
		}
	}
	
	if(nflip > 0){
		cout << nflip << " markers with MAF > 0.5 are flipped" << endl;
	}
	
	//load y
	y64 = bitvec(nblock, (uint64) 0);
	not_y64 = bitvec(nblock, (uint64) 0);
	y = vector<int> (nsub, -1);
	ifstream file_phen(path_fam);
	if(!file_phen){
		cout << "Error: Cannot open FAM file " << path_fam << endl;
		exit(1);
	}
	
	int k = -1;
	int t = 0;
	ncase = 0;//phen == 2
	nctrl = 0;//phen == 1
	for(string s; getline(file_phen, s); ){
		istringstream sin(s);
		string fam_id, ind_id, pat_id, mat_id, sex;
		int phen;
		++k;
		if(k == index_individual_id_fam_inc[t]){
			if(sin >> fam_id >> ind_id >> pat_id >> mat_id >> sex >> phen){
				if(phen == 2){//case
					y64[bitloc[t]] |= MASK_offset[t];
					y[t] = 1;
					++ncase;
				}else if(phen == 1){//control
					not_y64[bitloc[t]] |= MASK_offset[t];
					y[t] = 0;
					++nctrl;
				}else{
					cout << "Error: Invalid coding in " << path_fam << endl;
					file_phen.close();
					exit(1);
				}
			}else{
				cout << "Error: Invalid coding in " << path_fam << endl;
				file_phen.close();
				exit(1);
			}
			
			++t;
		}else{
			continue;
		}
	}
	file_phen.close();
	
	//load continuous covariates
	if(has_cont){
		ifstream f_cont(path_cont);
		if(!f_cont){
			cout << "Error: Cannot open CON file " << path_cont << endl;
			exit(1);
		}
		
		vector<bool> index_inc (nsub_cont, false);
		for(int i = 0; i < index_individual_id_cont_inc.size(); ++i){
			index_inc[index_individual_id_cont_inc[i]] = true;
		}
		
		vector<vector<double> > xcont_all = vector<vector<double> > (ncont, vector<double> (nsub_cont, -999999.0));
		k = -2;
		t = 0;
		for(string s; getline(f_cont, s); ){
			istringstream sin(s);
			++k;
			string iid;
			if(k >= 0){
				if(index_inc[k]){//no missing in this row
					sin >> iid;
					for(int j = 0; j < ncont; ++j){
						sin >> xcont_all[j][k];
					}
				}
			}else if(k == -1){//skip the header
				sin >> iid;
				string header;
				while(sin >> header){
					cont_var_name.push_back(header);
				}
				continue;
			}else{//impossible
				cout << "Error: in loading continuous covarates" << endl;
				exit(1);
			}
		}
		f_cont.close();
		
		xcont = vector<vector<double> > (ncont, vector<double> (nsub, .0));
		for(int j = 0; j < ncont; ++j){
			for(int i = 0; i < nsub; ++i){
				xcont[j][i] = xcont_all[j][index_individual_id_cont_inc[i]];
			}
		}
		
		EncodeCont();
	}
	
	if(has_cate){
		ifstream f_cate(path_cate);
		if(!f_cate){
			cout << "Error: Cannot open CAT file " << path_cate << endl;
			exit(1);
		}
		
		vector<bool> index_inc (nsub_cate, false);
		for(int i = 0; i < index_individual_id_cate_inc.size(); ++i){
			index_inc[index_individual_id_cate_inc[i]] = true;
		}
		
		vector<vector<string> > xcate_all = vector<vector<string> > (ncate, vector<string> (nsub_cate, string("err")));
		k = -2;
		t = 0;
		for(string s; getline(f_cate, s); ){
			istringstream sin(s);
			++k;
			string iid;
			if(k >= 0){
				if(index_inc[k]){//no missing in this row
					sin >> iid;
					for(int j = 0; j < ncate; ++j){
						sin >> xcate_all[j][k];
					}
				}
			}else if(k == -1){//skip the header
				sin >> iid;
				string header;
				while(sin >> header){
					cate_var_name.push_back(header);
				}
				continue;
			}else{//impossible
				cout << "Error: in loading categorical covariates" << endl;
				exit(1);
			}
		}
		f_cate.close();
		
		xcate = vector<vector<string> > (ncate, vector<string> (nsub));
		for(int j = 0; j < ncate; ++j){
			for(int i = 0; i < nsub; ++i){
				xcate[j][i] = xcate_all[j][index_individual_id_cate_inc[i]];
			}
		}
		
		EncodeCate();
	}
	
	oob_id64 = bitmat (ntree, bitvec(nblock, (uint64) 0));
	ib_id64 = bitmat (ntree, bitvec(nblock, (uint64) 0));
	
	num_oob = vector<int> (ntree, 0);
	
	////////
	
	//loading sample weights, if available
	samplewt = vector<double>(nsub, -1.0);
	
	if(has_samplewt){
		ifstream file_samplewt(path_samplewt);
		if(!file_samplewt){
			cout << "Error: Cannot open SWT file " << path_samplewt << endl;
			exit(1);
		}
		
		for(string s; getline(file_samplewt, s); ){
			istringstream sin(s);
			string ind_id;
			double swt;
			if(sin >> ind_id >> swt){
				for(int i = 0; i < nsub; ++i){
					if(ind_id == individual_id[i]){
						//when I load fam file before, I have forced individual_id to include samples from trainid.iid (if any),
						// and samples from samplewt.swt (if any, with weights > 0)
						//thus, the following condition statement should never be met
						if(swt < .0){
							cout << "Error: Sample weights are supposed to be zero or positive. The program is needed to be debugged" << endl;
							exit(1);
						}
						samplewt[i] = swt;
					}
				}
			}else{
				cout << "Error: Invalid format in " << path_samplewt << endl;
				exit(1);
			}
		}
		file_samplewt.close();
		
		double sumwt = .0;
		for(int i = 0; i < nsub; ++i){
			if(samplewt[i] < .0){
				cout << "Error: All training data should have a sample weight" << endl;
				exit(1);
			}else{
				sumwt += samplewt[i];
			}
		}
		for(int i = 0; i < nsub; ++i){
			samplewt[i] /= sumwt;
		}
	}else{
		samplewt = vector<double>(nsub, 1.0/nsub);
	}
	
//	for(int i = 0; i < nsub; ++i){
//		cout << individual_id[i] << "\t" << samplewt[i] << endl;
//	}
//	exit(1);
	
}

void BAMBOO::PrintModelInfo(){
	
	if(bam_list.size() == 1){
		ifstream ifs(path_bam);
		if(!ifs){
			cout << "Error: Cannot open BAM file " << path_bam << endl;
			exit(1);
		}
		
		boost::archive::binary_iarchive ar(ifs);
			
		ar & version;
		ar & nsnp & ncont & ncate;
		ar & flip;
		ar & snp_used_in_forest & snp_id_used_in_forest;
		ar & cont_var_used_in_forest & cont_var_id_used_in_forest;
		ar & cate_var_used_in_forest & cate_var_id_used_in_forest;
		ar & cate_unique & cate_code;
		ar & bamboo;
		
		ntree = bamboo.size();
		
		cout << "The model in [ " << path_bam << " ] is built with Random Bamboo v " << RANDOM_BAMBOO_VERSION << endl;
		cout << "The model includes: " << endl;
		cout << ntree << " bamboos" << endl;
		cout << snp_used_in_forest.size() << " markers" << endl;
		if(cont_var_id_used_in_forest.size() > 0){
			cout << cont_var_used_in_forest.size() << " continuous covariate(s)" << endl;
		}
		if(cate_var_used_in_forest.size() > 0){
			cout << cate_var_used_in_forest.size() << " categorical covariate(s)" << endl;
		}
		if(flip){
			cout << "Markers with MAF > 0.5 are flipped" << endl;
		}
	}else{
		
		int nforest = bam_list.size();
		vector<string> ver(nforest);
		vector<int> ns(nforest);
		vector<int> nco(nforest);
		vector<int> nca(nforest);
		vector<bool> fl(nforest);
		
		for(int i = 0; i < nforest; ++i){
			ifstream ifs(bam_list[i]);
			if(!ifs){
				cout << "Error: Cannot open BAM file " << bam_list[i] << endl;
				exit(1);
			}
			
			boost::archive::binary_iarchive ar(ifs);
			
			string ver0;
			int ns0, nco0, nca0;
			bool fl0;
			ar & ver0;
			ar & ns0 & nco0 & nca0;
			ar & fl0;
			
			ver[i] = ver0;
			ns[i] = ns0;
			nco[i] = nco0;
			nca[i] = nca0;
			fl[i] = fl0;
			
		}
		
		for(int i = 1; i < nforest; ++i){
			if(ver[i] != ver[0]){
				cout << "Error: The specified models are built with Random Bamboo with different versions" << endl;
				exit(1);
			}
			
			if(ns[i] != ns[0] || nco[i] != nco[0] || nca[i] != nca[i]){
				cout << "Error: The specified models are built with different training data" << endl;
				exit(1);
			}
			
			if(fl[i] != fl[0]){
				cout << "Error: Some of the specified models flip the genotypes but others do not" << endl;
				exit(1);
			}
		}
		
		nsnp = ns[0];
		ncont = nco[0];
		ncate = nca[0];
		flip = fl[0];
		
		ComputeWordBits();//don't forget this otherwise the PopCount function fails
		
		LEN = 64;
		int nblock_snp = (int) ceil (((double) nsnp) / LEN);
		int nblock_cont = (int) ceil (((double) ncont) / LEN);
		int nblock_cate = (int) ceil (((double) ncate) / LEN);
		
		vector<int> bitloc_snp;
		vector<int> bitloc_cont;
		vector<int> bitloc_cate;
		
		bitvec MASK_offset_snp;
		bitvec MASK_offset_cont;
		bitvec MASK_offset_cate;
		
		bitloc_snp.reserve(nsnp);
		bitloc_cont.reserve(ncont);
		bitloc_cate.reserve(ncate);
		
		MASK_offset_snp.reserve(nsnp);
		MASK_offset_cont.reserve(ncont);
		MASK_offset_cate.reserve(ncate);
		
		MASK = 0x8000000000000000;
		
		for(int k = 0; k < nsnp; ++k){
			bitloc_snp.push_back(k / LEN);
			MASK_offset_snp.push_back(MASK >> (k % LEN));
		}
		
		for(int k = 0; k < ncont; ++k){
			bitloc_cont.push_back(k / LEN);
			MASK_offset_cont.push_back(MASK >> (k % LEN));
		}
		
		for(int k = 0; k < ncate; ++k){
			bitloc_cate.push_back(k / LEN);
			MASK_offset_cate.push_back(MASK >> (k % LEN));
		}
		
		bitmat snp_id_used_in_multiple_models64(nthread, bitvec(nblock_snp, (uint64) 0));
		bitmat cont_id_used_in_multiple_models64(nthread, bitvec(nblock_cont, (uint64) 0));
		bitmat cate_id_used_in_multiple_models64(nthread, bitvec(nblock_cate, (uint64) 0));
		
		vector<int> nt (nthread, 0);
		
		#pragma omp parallel num_threads(nthread)
		{
			#pragma omp for
			for(int i = 0; i < nforest; ++i){
				
				int tid = omp_get_thread_num();
				ifstream ifs(bam_list[i]);
				if(!ifs){
					cout << "Error: Cannot open BAM file " << bam_list[i] << endl;
					exit(1);
				}
				
				boost::archive::binary_iarchive ar(ifs);
				
				string version0;
				int nsnp0, ncont0, ncate0;
				bool flip0;
				vector<string> snp_used_in_forest0, cont_var_used_in_forest0, cate_var_used_in_forest0;
				vector<int> snp_id_used_in_forest0, cont_var_id_used_in_forest0, cate_var_id_used_in_forest0;
				vector<vector<string> > cate_unique0;
				vector<vector<int> > cate_code0;
				vector<vector<MINI_NODE> > bamboo0;
					
				ar & version0;
				ar & nsnp0 & ncont0 & ncate0;
				ar & flip0;
				ar & snp_used_in_forest0 & snp_id_used_in_forest0;
				ar & cont_var_used_in_forest0 & cont_var_id_used_in_forest0;
				ar & cate_var_used_in_forest0 & cate_var_id_used_in_forest0;
				ar & cate_unique0 & cate_code0;
				ar & bamboo0;
				
				nt[tid] += bamboo0.size();
				//swap is very slow but can save memory (to be confirmed)
//				vector<vector<MINI_NODE> >().swap(bamboo0);
//				vector<string>().swap(snp_used_in_forest0);
//				vector<string>().swap(cont_var_used_in_forest0);
//				vector<string>().swap(cate_var_used_in_forest0);
				
				for(int k = 0; k < snp_id_used_in_forest0.size(); ++k){
					int sid = snp_id_used_in_forest0[k];
					snp_id_used_in_multiple_models64[tid][bitloc_snp[sid]] |= MASK_offset_snp[sid];
				}
				
				for(int k = 0; k < cont_var_id_used_in_forest0.size(); ++k){
					int cid = cont_var_id_used_in_forest0[k];
					cont_id_used_in_multiple_models64[tid][bitloc_cont[cid]] |= MASK_offset_cont[cid];
				}
				
				for(int k = 0; k < cate_var_id_used_in_forest0.size(); ++k){
					int cid = cate_var_id_used_in_forest0[k];
					cate_id_used_in_multiple_models64[tid][bitloc_cate[cid]] |= MASK_offset_cate[cid];
				}
				
			}
			
		}
		
		int nsnp_used = 0;
		for(int j = 0; j < nblock_snp; ++j){
			uint64 u = (uint64) 0;
			for(int l = 0; l < nthread; ++l){
				u |= snp_id_used_in_multiple_models64[l][j];
			}
			nsnp_used += PopCount(u);
		}
		
		int ncont_used = 0;
		for(int j = 0; j < nblock_cont; ++j){
			uint64 u = (uint64) 0;
			for(int l = 0; l < nthread; ++l){
				u |= cont_id_used_in_multiple_models64[l][j];
			}
			ncont_used += PopCount(u);
		}
		
		int ncate_used = 0;
		for(int j = 0; j < nblock_cate; ++j){
			uint64 u = (uint64) 0;
			for(int l = 0; l < nthread; ++l){
				u |= cate_id_used_in_multiple_models64[l][j];
			}
			ncate_used += PopCount(u);
		}
		
		ntree = 0;
		for(int j = 0; j < nthread; ++j){
			ntree += nt[j];
		}
		
		cout << "The " << bam_list.size() << " models in [ " << dir << " ] include: " << endl;
		cout << ntree << " bamboos" << endl;
		cout << nsnp_used << " markers" << endl;
		if(ncont_used > 0){
			cout << ncont_used << " continuous covariate(s)" << endl;
		}
		if(ncate_used > 0){
			cout << ncate_used << " categorical covariate(s)" << endl;
		}
		if(flip){
			cout << "Markers with MAF > 0.5 are flipped" << endl;
		}
	}
	
}

void BAMBOO::ParseModelPath(const char *const input_path_bam){
	
	path p (input_path_bam);
	if(exists(p)){
		path_bam = NULL;
		if(is_directory(p)){
			
			dir = p.string();
			
			typedef vector<path> vec;
			vec v;
			
			vec bl;
			
			copy(directory_iterator(p), directory_iterator(), back_inserter(v));
			sort(v.begin(), v.end());
      for (vec::const_iterator it(v.begin()), it_end(v.end()); it != it_end; ++it){
        if((*it).extension() == string(".bam")){
        	bl.push_back((*it).native());
        }
      }
      
      if(bl.size() == 0){
      	cout << "Error: No BAM file is found in directory [ " << p.string() << " ]" << endl;
      	exit(1);
      }else if(bl.size() == 1){
      	path_bam = strdup(bl[0].string().c_str());
      	bam_list.push_back(path_bam);
      }else{
      	path_bam = NULL;//set it as NULL if more than one bam are used in prediction
      	for(int i = 0; i < bl.size(); ++i){
      		bam_list.push_back(strdup(bl[i].string().c_str()));
      	}
      	cout << "Loading " << bam_list.size() << " BAM files from [ " << p.string() << " ]" << endl;
      }
      
		}else{
			cout << "Error: Invalid input of --bam. It should be a BAM file's name (no extension) or a directory" << endl;
			exit(1);
		}
	}else{
		path_bam = new char[strlen(input_path_bam)+5];
		path_bam[0] = '\0';
		strcat(path_bam, input_path_bam);
		strcat(path_bam, ".bam");
		path file(path_bam);
		if(is_regular_file(file)){
			bam_list.push_back(path_bam);
			//cout << file << " is a file with size " << file_size(file) << endl;
		}else{
			cout << "Error: Invalid input of --bam. It should be a BAM file's name (no extension) or a directory" << endl;
			exit(1);
		}
	}
	
	
}

BAMBOO::BAMBOO(const char *const input_path_bam, const int input_nthread){
	
	nthread = input_nthread;
	ParseModelPath(input_path_bam);
	
	pred_from_trained_model = false;
	pred_from_specified_model = false;
		
}

BAMBOO::BAMBOO(const char *const input_path_out, const char *const input_path_test, 
	const char *const input_path_bam, const char *const input_path_testid, 
	const int input_nthread, const int input_ntree){
	
	time_t start_time;
	time(&start_time);
	cout << "Program started: " << ctime(&start_time);
	cout << "Random Bamboo is launched for predicting data ..." << endl;
	
	nthread = input_nthread;
	ntree = input_ntree;
	
	has_test = true;
	
	path_fam_test = new char[strlen(input_path_test)+5];
	path_bim_test = new char[strlen(input_path_test)+5];
	path_bed_test = new char[strlen(input_path_test)+5];
	path_pred = new char[strlen(input_path_out)+5];
	path_fam_test[0] = '\0';
	path_bim_test[0] = '\0';
	path_bed_test[0] = '\0';
	path_pred[0] = '\0';
	strcat(path_fam_test, input_path_test);
	strcat(path_bim_test, input_path_test);
	strcat(path_bed_test, input_path_test);
	strcat(path_pred, input_path_out);
	strcat(path_fam_test, ".fam");
	strcat(path_bim_test, ".bim");
	strcat(path_bed_test, ".bed");
	strcat(path_pred, ".prd");
	
	ParseModelPath(input_path_bam);
	
	path_cont_test = new char[strlen(input_path_test)+5];
	path_cont_test[0] = '\0';
	strcat(path_cont_test, input_path_test);
	strcat(path_cont_test, ".con");
	
	path_cate_test = new char[strlen(input_path_test)+5];
	path_cate_test[0] = '\0';
	strcat(path_cate_test, input_path_test);
	strcat(path_cate_test, ".cat");
	
	if(input_path_testid){
		path_testid = new char[strlen(input_path_testid)+5];
		path_testid[0] = '\0';
		strcat(path_testid, input_path_testid);
		strcat(path_testid, ".iid");
	}else{
		path_testid = NULL;
	}
	
	pred_from_trained_model = false;
	pred_from_specified_model = true;
}

//this constructor is used when --file is specified
//if --pred is also set, then the test dataset is predicted using the model trained from --file
BAMBOO::BAMBOO(const char *const input_path_plink, const char *const input_path_out, 
	const char * const input_path_cont, const char * const input_path_cate, 
	const char *const input_path_test, const char *const input_path_trainid, 
	const char *const input_path_testid, const char *const input_path_snpid, 
	const char *const input_path_samplewt, 
	const int input_ntree, const int input_mtry, const int input_max_nleaf, 
	const int input_min_leaf_size, const int input_imp_measure, const int input_seed, 
	const int input_nthread, const double input_class_weight, const double input_cutoff, 
	const bool input_flip, const bool input_output_prox, const bool input_output_imp, 
	const bool input_output_bamboo, const bool input_keepcovar, 
	const bool input_balance, const bool input_trace){
	
	time_t start_time;
	time(&start_time);
	cout << "Program started: " << ctime(&start_time);
	cout << "Random Bamboo is launched for training predictive model ..." << endl;
	
	//loading parameters
	ntree = input_ntree;
	mtry = input_mtry;
	seed = input_seed;
	imp_measure = input_imp_measure;
	nthread = input_nthread;
	cutoff = input_cutoff;
	
	flip = input_flip;
	output_prox = input_output_prox;
	output_imp = input_output_imp;
	output_bamboo = input_output_bamboo;
	has_cont = input_path_cont ? true : false;
	has_cate = input_path_cate ? true : false;
	has_test = input_path_test ? true : false;
	has_samplewt = input_path_samplewt ? true : false;
		
	keepcovar = (input_keepcovar && (has_cont || has_cate)) ? true : false;
	
	balance = input_balance;
	trace = input_trace;
	
	//////
	
	bamboo = vector<vector<MINI_NODE> > (ntree);
	
	check_in = vector<bool> (ntree, false);
	check_out = vector<bool> (ntree, false);
	
	//model = vector<vector<vector<double> > > (ntree);
	
	version = string(RANDOM_BAMBOO_VERSION);
	
	path_fam = new char[strlen(input_path_plink)+5];
	path_bim = new char[strlen(input_path_plink)+5];
	path_bed = new char[strlen(input_path_plink)+5];
	path_fam[0] = '\0';
	path_bim[0] = '\0';
	path_bed[0] = '\0';
	strcat(path_fam, input_path_plink);
	strcat(path_bim, input_path_plink);
	strcat(path_bed, input_path_plink);
	strcat(path_fam, ".fam");
	strcat(path_bim, ".bim");
	strcat(path_bed, ".bed");
	
	if(has_cont){
		path_cont = new char[strlen(input_path_cont)+5];
		path_cont[0] = '\0';
		strcat(path_cont, input_path_cont);
		strcat(path_cont, ".con");
	}else{
		path_cont = NULL;
		ncont = 0;
	}
	
	if(has_cate){
		path_cate = new char[strlen(input_path_cate)+5];
		path_cate[0] = '\0';
		strcat(path_cate, input_path_cate);
		strcat(path_cate, ".cat");
	}else{
		path_cate = NULL;
		ncate = 0;
	}
	
	if(input_path_trainid){
		path_trainid = new char[strlen(input_path_trainid)+5];
		path_trainid[0] = '\0';
		strcat(path_trainid, input_path_trainid);
		strcat(path_trainid, ".iid");
	}else{
		path_trainid = NULL;
	}
	
	if(input_path_snpid){
		path_snpid = new char[strlen(input_path_snpid)+5];
		path_snpid[0] = '\0';
		strcat(path_snpid, input_path_snpid);
		strcat(path_snpid, ".sid");
	}else{
		path_snpid = NULL;
	}
	
	if(input_path_samplewt){
		path_samplewt = new char[strlen(input_path_samplewt)+5];
		path_samplewt[0] = '\0';
		strcat(path_samplewt, input_path_samplewt);
		strcat(path_samplewt, ".swt");
	}else{
		path_samplewt = NULL;
	}
	
	ComputeWordBits();
	
	LoadTrainingData();
	
	
	if(output_imp){
		path_imp = new char[strlen(input_path_out)+5];
		path_imp[0] = '\0';
		strcat(path_imp, input_path_out);
		strcat(path_imp, ".imp");
	}else{
		path_imp = NULL;
	}
	
	if(output_prox){
		path_prox = new char[strlen(input_path_out)+5];
		path_prox[0] = '\0';
		strcat(path_prox, input_path_out);
		strcat(path_prox, ".prx");
	}else{
		path_prox = NULL;
	}
	
	if(output_bamboo){
		path_bam = new char[strlen(input_path_out)+5];
		path_bam[0] = '\0';
		strcat(path_bam, input_path_out);
		strcat(path_bam, ".bam");
	}else{
		path_bam = NULL;
	}
	
	path_err = new char[strlen(input_path_out)+5];
	path_err[0] = '\0';
	strcat(path_err, input_path_out);
	strcat(path_err, ".err");
	
	path_cof = new char[strlen(input_path_out)+5];
	path_cof[0] = '\0';
	strcat(path_cof, input_path_out);
	strcat(path_cof, ".cof");
	
	path_auc = new char[strlen(input_path_out)+5];
	path_auc[0] = '\0';
	strcat(path_auc, input_path_out);
	strcat(path_auc, ".auc");
	
	
	/////
	if(input_max_nleaf <= 0 || input_max_nleaf > nsub){
		max_nleaf = nsub;
	}else{
		max_nleaf = input_max_nleaf;
	}
	
	if(input_min_leaf_size <= 0){
		min_leaf_size = 1;
	}else if(input_min_leaf_size >= nsub){
		cout << "Error: The option --minleafsize is too large (" << input_min_leaf_size << "). Program terminates" << endl;
		exit(1);
	}else{
		min_leaf_size = input_min_leaf_size;
	}
	
	
	//other parameters
	
	
	if(balance){
		class_weight = nctrl * 1.0 / ncase;
	}else{
		class_weight = input_class_weight;
	}
	
	if(!mtry){// 0 is the default value of input_mtry, in that case, the "sqrt" rule is applied
		mtry = (int) sqrt(1.0 * (nsnp_specified + ncont + ncate));
	}else{
		if(mtry > nsnp_specified + ncont + ncate){
			cout << "Error: mtry is too large" << endl;
			exit(1);
		}
	}
	
	if(keepcovar && mtry < ncont + ncate){
		cout << "Error: mtry is too small. Try --selcovar" << endl;
		exit(1);
	}
	
	time_t data_loaded_time;
	time(&data_loaded_time);
	cout << "Data loaded: " << ctime(&data_loaded_time);
	
	PrintPara();
	
	/////prediction
	
	importance = vector<IMPORTANCE> (nvar, imp0);
	for(int i = 0; i < nvar; ++i){
		importance[i].var_id = i;
	}
	
	oob_error = vector<double> (ntree, .0);
	num_vote_correct_class = vector<int> (ntree, 0);
	if(imp_measure == IMP_BREIMAN_CUTLER || imp_measure == IMP_LIAW_WIENER){
		num_vote_loss_permuted_correct_class = vector<vector<uint16> > (ntree);
	}
	
	fitted_risk_ctrl = vector<double> (nsub, .0);
	fitted_risk_case = vector<double> (nsub, .0);
	auc = .0;
	
	/////
	
	if(has_test){
		path_fam_test = new char[strlen(input_path_test)+5];
		path_bim_test = new char[strlen(input_path_test)+5];
		path_bed_test = new char[strlen(input_path_test)+5];
		path_pred = new char[strlen(input_path_out)+5];
		path_fam_test[0] = '\0';
		path_bim_test[0] = '\0';
		path_bed_test[0] = '\0';
		path_pred[0] = '\0';
		strcat(path_fam_test, input_path_test);
		strcat(path_bim_test, input_path_test);
		strcat(path_bed_test, input_path_test);
		strcat(path_pred, input_path_out);
		strcat(path_fam_test, ".fam");
		strcat(path_bim_test, ".bim");
		strcat(path_bed_test, ".bed");
		strcat(path_pred, ".prd");
		
		if(has_cont){
			path_cont_test = new char[strlen(input_path_test)+5];
			path_cont_test[0] = '\0';
			strcat(path_cont_test, input_path_test);
			strcat(path_cont_test, ".con");
		}else{
			path_cont_test = NULL;
		}
		
		if(has_cate){
			path_cate_test = new char[strlen(input_path_test)+5];
			path_cate_test[0] = '\0';
			strcat(path_cate_test, input_path_test);
			strcat(path_cate_test, ".cat");
		}else{
			path_cate_test = NULL;
		}
		
		if(input_path_testid){
			path_testid = new char[strlen(input_path_testid)+5];
			path_testid[0] = '\0';
			strcat(path_testid, input_path_testid);
			strcat(path_testid, ".iid");
		}else{
			path_testid = NULL;
		}
		
		pred_from_trained_model = true;
		pred_from_specified_model = false;
	}else{
		pred_from_trained_model = false;
		pred_from_specified_model = false;
	}
	
}

void BAMBOO::CreateFile(char *file, const char *const bfile){
	
	if(!bfile){
		file = new char[strlen(bfile)+1];
		file[0] = '\0';
		strcat(file, bfile);
	}else{
		file = NULL;
	}
	
}

//Walker's alias method is modified from R's library by Yifan Yang
//this function is used for sampling from specified distribution with replacement
//p is the distribution to be sampled from
//nans is sample size
void BAMBOO::Walker_ProbSampleReplace(drand48_data &buf, const vector<double> &p, const int nans, vector<int> &ans){
	
	ans = vector<int>(nans, -1);
	
	int n = p.size();
	vector<int> a(n);
	
    double *q, rU;
    int i, j, k;
    int *HL, *H, *L;

    /* Create the alias tables.
       The idea is that for HL[0] ... L-1 label the entries with q < 1
       and L ... H[n-1] label those >= 1.
       By rounding error we could have q[i] < 1. or > 1. for all entries.
     */

	/* Slow enough anyway not to risk overflow */
	HL =(int *) calloc(n, sizeof(int));
	q = (double *)calloc(n, sizeof(double));

    H = HL - 1; L = HL + n;
    for (i = 0; i < n; i++) {
	q[i] = p[i] * n;
	if (q[i] < 1.) *++H = i; else *--L = i;
    }
    if (H >= HL && L < HL + n) { /* So some q[i] are >= 1 and some < 1 */
	for (k = 0; k < n - 1; k++) {
	    i = HL[k];
	    j = *L;
	    a[i] = j;
	    q[j] += q[i] - 1;
	    if (q[j] < 1.) L++;
	    if(L >= HL + n) break; /* now all are >= 1 */
	}
    }
    for (i = 0; i < n; i++) q[i] += i;

    /* generate sample */
    for (i = 0; i < nans; i++) {
    	double dr;
    	drand48_r(&buf, &dr);
	rU = dr * n;
	k = (int) rU;
	ans[i] = (rU < q[k]) ? k+1 : a[k]+1;
		ans[i] -= 1;//added by Han to make id starting from 0 rather than 1
    }

	free(HL);
	free(q);
}

//buf: used in random number generator
//y64_omp, not_y64_omp: storing y after bootstrapping, expanded as 64-bit matrix to represent replicated samples
//in_sample64_omp: indicate which element in y64_omp and not_y64_omp are included in the bootstrap samples, expanded and in bit-wise
//iter_omp: indicates which element in in_sample64_omp has at least one 1-bit
//oob_id_omp: the index of oob sample, subset of (0, ..., nsub-1)
//not_oob_id_omp: the index of bootstrap samples without replication, subset of (0, ..., nsub-1)
//id_cnt_omp: a vector of length nsub. Indicates every original sample is repeated by how many times in the bootstrap samples.
//ncase_omp: how many cases are included in the bootstrap samples
void BAMBOO::Bootstrap(const int tree_id, drand48_data &buf, bitmat &y64_omp, bitmat &not_y64_omp, 
	bitmat &in_sample64_omp, vector<vector<int> > &iter_omp, vector<int> &oob_id_omp, 
	vector<int> &not_oob_id_omp, vector<int> &id_cnt_omp, int &ncase_omp){
	
	vector<int> sel_id (nsub, -1);//storing ids of samples selected in bootstrap samples, with replication
	id_cnt_omp = vector<int> (nsub, 0);//count, initialized as 0
	
	if(!has_samplewt){
		for(int k = 0; k < nsub; ++k){
			double dr;
			drand48_r(&buf, &dr);
			sel_id[k] = (int) (dr * nsub);
			if(sel_id[k] < 0 || sel_id[k] >= nsub){
				cout << "Error: in function Bootstrap" << endl;
				exit(1);
			}
			
			id_cnt_omp[sel_id[k]]++;
		}
	}else{
		
		Walker_ProbSampleReplace(buf, samplewt, nsub, sel_id);
		for(int k = 0; k < nsub; ++k){
			if(sel_id[k] < 0 || sel_id[k] >= nsub){
				cout << sel_id[k] << "\t" << nsub << endl;
				cout << "Error: in function Bootstrap" << endl;
				exit(1);
			}
			id_cnt_omp[sel_id[k]]++;
		}
	}
	
	if(false){
		cout << "Bootstrap samples:\nsid\tcount" << endl;
		for(int k = 0; k < nsub; ++k){
			cout << k+1 << "\t" << id_cnt_omp[k] << endl;
		}
	}
	
	///////
	int max_cnt = 0;//find the max in count. max_cnt is used to determine the dim of y64_omp and not_y64_omp
	ncase_omp = 0;
	if(!oob_id_omp.empty()){
		oob_id_omp.clear();
	}
	if(!not_oob_id_omp.empty()){
		not_oob_id_omp.clear();
	}
	
	for(int k = 0; k < nsub; ++k){
		
		if(id_cnt_omp[k] > max_cnt){
			max_cnt = id_cnt_omp[k];
		}
		if(id_cnt_omp[k]){//the k-th sample is not an oob (included in the Bootstrap samples)
			not_oob_id_omp.push_back(k);
			if(y[k]){//y == 1
				ncase_omp += id_cnt_omp[k];
			}
		}else{//the k-th sample is an oob (not included in the Bootstrap samples)
			oob_id_omp.push_back(k);
			if(false) cout << k+1 << endl;//print index of oob
		}
	}
	
	///////
	y64_omp = bitmat (max_cnt, bitvec(nblock, (uint64) 0));
	not_y64_omp = bitmat (max_cnt, bitvec(nblock, (uint64) 0));
	in_sample64_omp = bitmat (max_cnt, bitvec(nblock, (uint64) 0));
	
	num_oob[tree_id] = oob_id_omp.size();
	for(int k = 0; k < nsub; ++k){
		if(id_cnt_omp[k] == 0){
			oob_id64[tree_id][bitloc[k]] |= MASK_offset[k];
		}else{
			ib_id64[tree_id][bitloc[k]] |= MASK_offset[k];
		}
		for(int j = 0; j < id_cnt_omp[k]; ++j){
			in_sample64_omp[j][bitloc[k]] |= MASK_offset[k];
			if(y[k]){//y = 1
				y64_omp[j][bitloc[k]] |= MASK_offset[k];
			}else{
				not_y64_omp[j][bitloc[k]] |= MASK_offset[k];
			}
		}
	}
	
	if(!iter_omp.empty()){
		iter_omp.clear();
	}
	iter_omp.reserve(max_cnt);
	for(int j = 0; j < max_cnt; ++j){
		vector<int> v;
		for(int l = 0; l < nblock; ++l){
			if(in_sample64_omp[j][l]){
				v.push_back(l);
			}
		}
		iter_omp.push_back(v);
	}
	
}

void BAMBOO::ShuffleSNP(drand48_data &buf, vector<int> &sel_snp_id_omp){
	
	vector<int> snp_id = inc_snp_id;//index of snp id, initialized to snp ids specified by the user or 1:nsnp
	
	for(int k = snp_id.size() - 1; k > 0; --k){
		long int li;
		lrand48_r(&buf, &li);
		int j = li % (k + 1);
		int s = snp_id[j];
		snp_id[j] = snp_id[k];
		snp_id[k] = s;
	}
	
	if(!sel_snp_id_omp.empty()){
		sel_snp_id_omp.clear();
	}
	
	if(mtry > inc_snp_id.size()){
		cout << "Error: debug ShuffleSNP" << endl;
	}
	
	for(int i = 0; i < mtry; ++i){//selecte the top mtry SNPs to split a node
		sel_snp_id_omp.push_back(snp_id[i]);
	}
	sort(sel_snp_id_omp.begin(), sel_snp_id_omp.end());
	
}

void BAMBOO::ShuffleAllFeature(drand48_data &buf, vector<int> &sel_snp_id_omp, 
	vector<int> &sel_cont_id_omp, vector<int> &sel_cate_id_omp){
	
	vector<int> fid = inc_snp_id;//feature id. Index of snp id and other covariates' id
	for(int i = nsnp; i < nvar; ++i){
		fid.push_back(i);
	}
	
	int nvar0 = fid.size();
	
	for(int k = nvar0 - 1; k > 0; --k){
		long int li;
		lrand48_r(&buf, &li);
		int j = li % (k + 1);
		int f = fid[j];
		fid[j] = fid[k];
		fid[k] = f;
	}
	
	if(!sel_snp_id_omp.empty()){
		sel_snp_id_omp.clear();
	}
	if(!sel_cont_id_omp.empty()){
		sel_cont_id_omp.clear();
	}
	if(!sel_cate_id_omp.empty()){
		sel_cate_id_omp.clear();
	}
	
	for(int i = 0; i < mtry; ++i){//select SNPs from the top mtry features to split a node
		if(fid[i] < nsnp){
			sel_snp_id_omp.push_back(fid[i]);
		}else if(fid[i] < nsnp + ncont){
			sel_cont_id_omp.push_back(fid[i]-nsnp);//select continuous features
		}else{
			sel_cate_id_omp.push_back(fid[i]-nsnp-ncont);//select catigorical features
		}
	}
	sort(sel_snp_id_omp.begin(), sel_snp_id_omp.end());
	sort(sel_cont_id_omp.begin(), sel_cont_id_omp.end());
	sort(sel_cate_id_omp.begin(), sel_cate_id_omp.end());
	
}

void BAMBOO::ShuffleSNPKeepAllCovar(drand48_data &buf, vector<int> &sel_snp_id_omp, 
	vector<int> &sel_cont_id_omp, vector<int> &sel_cate_id_omp){
	
	vector<int> snp_id = inc_snp_id;//index of snp id, initialized to snp ids specified by the user or 1:nsnp
	
	for(int k = snp_id.size() - 1; k > 0; --k){
		long int li;
		lrand48_r(&buf, &li);
		int j = li % (k + 1);
		int s = snp_id[j];
		snp_id[j] = snp_id[k];
		snp_id[k] = s;
	}
	
	if(!sel_snp_id_omp.empty()){
		sel_snp_id_omp.clear();
	}
	
	if(mtry > inc_snp_id.size() + ncate + ncont){
		cout << "Error: debug ShuffleSNP" << endl;
	}
	
	for(int i = 0; i < mtry - ncate - ncont; ++i){//selecte the top mtry-ncate-ncont SNPs to split a node
		sel_snp_id_omp.push_back(snp_id[i]);
	}
	for(int i = 0; i < ncate; ++i){//all categorical covariates included
		sel_cate_id_omp.push_back(i);
	}
	for(int i = 0; i < ncont; ++i){//all continuous covariates included
		sel_cont_id_omp.push_back(i);
	}
	
	sort(sel_snp_id_omp.begin(), sel_snp_id_omp.end());
	sort(sel_cont_id_omp.begin(), sel_cont_id_omp.end());
	sort(sel_cate_id_omp.begin(), sel_cate_id_omp.end());
	
}

void BAMBOO::Shuffle(drand48_data &buf, vector<int> &sel_snp_id_omp, 
	vector<int> &sel_cont_id_omp, vector<int> &sel_cate_id_omp){
	
	if((has_cont && ncont > 0) || (has_cate && ncate > 0)){//covar 
		if(keepcovar){
			ShuffleSNPKeepAllCovar(buf, sel_snp_id_omp, sel_cont_id_omp, sel_cate_id_omp);//all covariates are used as candidate split for all nodes
		}else{
			ShuffleAllFeature(buf, sel_snp_id_omp, sel_cont_id_omp, sel_cate_id_omp);//select snps and covariates as the candidate to be split
		}
	}else{
		ShuffleSNP(buf, sel_snp_id_omp);//select snps only
	}
	
}


bool BAMBOO::SplitNode(NODE &node, const bitmat &y64_omp, 
	const bitmat &not_y64_omp, const vector<int> &id_cnt_omp, 
	const vector<vector<int> > &iter_omp, drand48_data &buf){
		
	if(!node.ncase || !node.nctrl){//pure node, should be a leaf
		node.processed = true;
		node.terminal = true;
		return !node.terminal;
	}
	
	//if node.processed == true, then the node has been processed and useful information have been stored
	//if node.terminal == true, a terminal cannot be further split
	if(node.processed || node.terminal){//this node has been processed or is a terminal node
		return !node.terminal;
	}
	
	if(false) cout << "examining node " << node.node_id+1 << endl;
	
	vector<int> sel_snp_id_omp;
	vector<int> sel_cont_id_omp;
	vector<int> sel_cate_id_omp;
	//choose features as candidates for node splitting
	Shuffle(buf, sel_snp_id_omp, sel_cont_id_omp, sel_cate_id_omp);
	
	if(false){
		if(false){
			cout << "selected SNP id:" << endl;
			for(int i = 0; i < sel_snp_id_omp.size(); ++i){
				cout << sel_snp_id_omp[i]+1 << endl;
			}
		}
		if(false){
			cout << "selected continuous id:" << endl;
			for(int i = 0; i < sel_cont_id_omp.size(); ++i){
				cout << cont_var_name[sel_cont_id_omp[i]] << endl;
			}
		}
		if(false){
			cout << "selected categorical id:" << endl;
			for(int i = 0; i < sel_cate_id_omp.size(); ++i){
				cout << cate_var_name[sel_cate_id_omp[i]] << endl;
			}
		}
	}
	
	/////split the node
	
	//sample sizes in left and right children
	int ncase_left = 0, nctrl_left = 0, N_left = 0;
	int ncase_right = 0, nctrl_right = 0, N_right = 0;
	double max_split_stat = -1.0;
	bool split_by_cont = false;
	bool split_by_cate = false;
	bool split_by_snp = false;
	
	//find split by SNPs
	int split_snp_id = -1;
	double max_snp_stat = -1.0;
	int split_snp_left_val = -1;
	
	assert(sel_snp_id_omp.size() + sel_cont_id_omp.size() + sel_cate_id_omp.size() == mtry);
	for(int i = 0; i < sel_snp_id_omp.size(); ++i){//evaluate the SNP: sel_snp_id_omp[i]
		
		int snp_id = sel_snp_id_omp[i];//snp id
		bitvec &ref_geno64_g0 = geno64[snp_id * 3];
		bitvec &ref_geno64_g1 = geno64[snp_id * 3 + 1];
		
		//table is actually the contingency table of control/case x SNP genotypes
		int table[3][3] = {{0}};
		
		if(node.large_sample_size){//if sample size in the node is sufficiently large, using Boolean operations
			uint64 rsid;
			uint64 u0, u1;
			uint64 x0, x1;
			
			for(int j = 0; j < iter_omp.size(); ++j){
				bitvec &ref_sample_id64 = node.sample_id64[j];
				const bitvec &ref_not_y64 = not_y64_omp[j];
				const bitvec &ref_y64 = y64_omp[j];
				for(int l = 0; l < iter_omp[j].size(); ++l){
					int k = iter_omp[j][l];
					rsid = ref_sample_id64[k];
					if(rsid){
						u0 = rsid & ref_geno64_g0[k];//snp == 0
						u1 = rsid & ref_geno64_g1[k];//snp == 1
						
						x0 = ref_not_y64[k];//control
						x1 = ref_y64[k];//case
						
						table[0][0] += PopCount(u0 & x0);
						table[0][1] += PopCount(u1 & x0);
						table[1][0] += PopCount(u0 & x1);
						table[1][1] += PopCount(u1 & x1);
					}
				}
			}
		}else{//otherwise, using counting to create contingency table
			for(int k = 0; k < node.sample_id.size(); ++k){//using all samples not in oob to create contingency table
				int sample_id = node.sample_id[k];
				int l = bitloc[sample_id];
				uint64 u = MASK_offset[sample_id];
				
				if(ref_geno64_g0[l] & u){//SNP == 0
					int c = id_cnt_omp[sample_id];
					if(y[sample_id]){//y == 1
						table[1][0] += c;
					}else{//y == 0
						table[0][0] += c;
					}
				}else if(ref_geno64_g1[l] & u){//SNP == 1
					int c = id_cnt_omp[sample_id];
					if(y[sample_id]){//y == 1
						table[1][1] += c;
					}else{//y == 0
						table[0][1] += c;
					}
				}else{//SNP == 2
					;//do nothing. table[0][2] and table[1][2] can be obtained by simple algebra below
				}
			}
		}
		
		table[0][2] = node.nctrl - table[0][0] - table[0][1];
		table[1][2] = node.ncase - table[1][0] - table[1][1];
		table[2][0] = table[0][0] + table[1][0];
		table[2][1] = table[0][1] + table[1][1];
		table[2][2] = table[0][2] + table[1][2];
		
		//for array ct, rows represent control/case, columns represent left/right children (partitions by SNP: 0/12 or 01/2)
		int ct[3][3] = {{0}};
		ct[0][2] = node.nctrl;
		ct[1][2] = node.ncase;
		ct[2][2] = node.N;
		
		double snp_stat = -1.0;
		int cl = -1;
		int nl0 = 0, nl1 = 0, Nl = 0, nr0 = 0, nr1 = 0, Nr = 0;
		
		for(int t = 0; t < 2; ++t){
			ct[0][0] += table[0][t];
			ct[1][0] += table[1][t];
			ct[0][1] = ct[0][2] - ct[0][0];
			ct[1][1] = ct[1][2] - ct[1][0];
			ct[2][0] = ct[0][0] + ct[1][0];
			ct[2][1] = ct[0][1] + ct[1][1];
			
			double stat = -1.0;//can be other value < -1.0
			if(ct[2][0] == 0 || ct[2][1] == 0){//this SNP has been used before
				stat = -1.0;
				
				if(false){
					cout << "=======" << endl;
					cout << "node " << node.node_id+1 << endl;
					cout << "snp " << snp_id+1 << endl;
					for(int ii = 0; ii < 3; ++ii){
						for(int jj = 0; jj < 3; ++jj){
							cout << ct[ii][jj] << "\t";
						}
						cout << endl;
					}
					cout << "stat = " << stat << endl;
				}
				
				continue;
			}else{
				
				//different weights for cases and controls, then weighted gini impurity is computed
				double p_left, p_right, p_par, p_left0, p_left1, p_right0, p_right1, p_par0, p_par1;
				p_left = ct[0][0] * 1.0 / nctrl + ct[1][0] * class_weight / ncase;//Pr(A_L): probability of a sample falling into the left child
				p_right = ct[0][1] * 1.0 / nctrl + ct[1][1] * class_weight / ncase;//Pr(A_R): probability of a sample falling into the right child
				p_par = p_left + p_right;//Pr(A): probability of a sample falling into the parent node
				p_left0 = ct[0][0] * 1.0 / nctrl / p_left;//Pr(y=0 | A_L): probability of a sample in the left child is a control
				p_left1 = 1.0 - p_left0;//Pr(y=1 | A_L)
				p_right0 = ct[0][1] * 1.0 / nctrl / p_right;//Pr(y=0 | A_R)
				p_right1 = 1.0 - p_right0;//Pr(y=1 | A_R)
				p_par0 = ct[0][2] * 1.0 / nctrl / p_par;//Pr(y=0 | A):probability of a sample in the parent node is a control
				p_par1 = 1.0 - p_par0;//pr(y=1 | A)
				
				stat = p_par * (1.0-p_par0*p_par0-p_par1*p_par1)-p_left*(1.0-p_left0*p_left0-p_left1*p_left1)-p_right*(1.0-p_right0*p_right0-p_right1*p_right1);
				
				
				if(false){
					cout << "=======" << endl;
					cout << "node " << node.node_id+1 << endl;
					cout << "snp " << snp_id+1 << endl;
					for(int ii = 0; ii < 3; ++ii){
						for(int jj = 0; jj < 3; ++jj){
							cout << ct[ii][jj] << "\t";
						}
						cout << endl;
					}
					cout << "stat = " << stat << endl;
				}
				
				if(stat < .0){
					continue;
				}
			}
			
			if(stat > .0 && stat - snp_stat > MIN_INC){
				snp_stat = stat;
				cl = t;//cut point (included) for left child
				
				nl0 = ct[0][0];
				nl1 = ct[1][0];
				Nl = ct[2][0];
				nr0 = ct[0][1];
				nr1 = ct[1][1];
				Nr = ct[2][1];
			}

		}
		
		if(snp_stat > .0 && snp_stat - max_split_stat > MIN_INC){
			split_by_snp = true;//at least one SNP is found to split the node
			max_snp_stat = snp_stat;
			max_split_stat = snp_stat;
			split_snp_id = snp_id;
			split_snp_left_val = cl;
			
			nctrl_left = nl0;
			ncase_left = nl1;
			N_left = Nl;
			nctrl_right = nr0;
			ncase_right = nr1;
			N_right = Nr;
		}
		
	}
	
	//find split by continue covariates
	int split_cont_lower = -1;
	int split_cont_upper = -1;
	int split_cont_id = -1;
	double max_cont_stat = -1.0;
	
	for(int i = 0; i < sel_cont_id_omp.size(); ++i){//evaluate every selected continuous and/or ordinary covariates
		
		int cont_id = sel_cont_id_omp[i];
		int lo = node.cont_lower[cont_id];
		int up = node.cont_upper[cont_id];
		
		if(lo >= up){//the covariate cannot be used anymore
			continue;
		}
		
		vector<vector<int> > table = vector<vector<int> > (2, vector<int> (up-lo+1, 0));
		for(int k = 0; k < node.sample_id.size(); ++k){//every sample not in oob, create contingency table
			int sample_id = node.sample_id[k];
			if(y[sample_id]){
				table[1][cont_loc[cont_id][sample_id]-lo] += id_cnt_omp[sample_id];
			}else{
				table[0][cont_loc[cont_id][sample_id]-lo] += id_cnt_omp[sample_id];
			}
		}
		
		int ct[3][3] = {{0}};
		ct[0][2] = node.nctrl;
		ct[1][2] = node.ncase;
		ct[2][2] = node.N;
		
		double cont_stat = -1.0;
		int cl = -1;
		int nl0 = 0, nl1 = 0, Nl = 0, nr0 = 0, nr1 = 0, Nr = 0;
		
		for(int k = 0; k < up-lo+1 - (1); ++k){//examining all split on this covariate
			if(!table[0][k] && !table[1][k]){//no sample fall in this column
				continue;
			}
			ct[0][0] += table[0][k];
			ct[1][0] += table[1][k];
			ct[0][1] = ct[0][2] - ct[0][0];
			ct[1][1] = ct[1][2] - ct[1][0];
			ct[2][0] = ct[0][0] + ct[1][0];
			ct[2][1] = ct[0][1] + ct[1][1];
			
			double stat = -1.0;
			if(ct[2][0] == 0 || ct[2][1] == 0){
				stat = -1.0;
				continue;
			}else{
				//different weights for cases and controls, then weighted gini impurity is computed
				double p_left, p_right, p_par, p_left0, p_left1, p_right0, p_right1, p_par0, p_par1;
				p_left = ct[0][0] * 1.0 / nctrl + ct[1][0] * class_weight / ncase;//Pr(A_L): probability of a sample falling into the left child
				p_right = ct[0][1] * 1.0 / nctrl + ct[1][1] * class_weight / ncase;//Pr(A_R): probability of a sample falling into the right child
				p_par = p_left + p_right;//Pr(A): probability of a sample falling into the parent node
				p_left0 = ct[0][0] * 1.0 / nctrl / p_left;//Pr(y=0 | A_L): probability of a sample in the left child is a control
				p_left1 = 1.0 - p_left0;//Pr(y=1 | A_L)
				p_right0 = ct[0][1] * 1.0 / nctrl / p_right;//Pr(y=0 | A_R)
				p_right1 = 1.0 - p_right0;//Pr(y=1 | A_R)
				p_par0 = ct[0][2] * 1.0 / nctrl / p_par;//Pr(y=0 | A):probability of a sample in the parent node is a control
				p_par1 = 1.0 - p_par0;//pr(y=1 | A)
				
				stat = p_par * (1.0-p_par0*p_par0-p_par1*p_par1)-p_left*(1.0-p_left0*p_left0-p_left1*p_left1)-p_right*(1.0-p_right0*p_right0-p_right1*p_right1);
				
				if(stat < .0){
					continue;
				}
			}
			
			if(stat > .0 && stat - cont_stat > MIN_INC){
				cont_stat = stat;
				cl = k + lo;
				nl0 = ct[0][0];
				nl1 = ct[1][0];
				Nl = ct[2][0];
				nr0 = ct[0][1];
				nr1 = ct[1][1];
				Nr = ct[2][1];
			}
		}
		
		if(false && (node.node_id+1==11 || node.node_id+1==26)){
			cout << node.node_id+1 << "\t" << cont_stat << endl;
		}
		
		if(cont_stat > .0 && cont_stat - max_split_stat > MIN_INC){
			
			split_by_cont = true;
			split_by_snp = false;
			max_cont_stat = cont_stat;
			max_split_stat = cont_stat;
			split_cont_id = cont_id;
			split_cont_lower = cl;
			split_cont_upper = cl + 1;
			for(int k = split_cont_lower+1-lo; k < up-lo+1; ++k){
				if(table[0][k] || table[1][k]){
					split_cont_upper = k + lo;
					break;
				}
			}
			
			if(split_cont_upper == -1){
				cout << "Error: check split_cont_upper" << endl;
				exit(1);
			}
			
			nctrl_left = nl0;
			ncase_left = nl1;
			N_left = Nl;
			nctrl_right = nr0;
			ncase_right = nr1;
			N_right = Nr;
		}
	}
	
	
	//find split by categorical covariates
	int split_cate_id = -1;
	double max_cate_stat = -1.0;
	vector<int> split_cate_left_code, split_cate_right_code;
	
	for(int i = 0; i < sel_cate_id_omp.size(); ++i){
		int cate_id = sel_cate_id_omp[i];
		
		vector<CATE_CELL> table(cate_code[cate_id].size(), cate_cell0);
		
		for(int t = 0; t < table.size(); ++t){
			CATE_CELL &tab = table[t];
			tab.code = t;
			tab.ncase = 0;
			tab.nctrl = 0;
			tab.r = .0;
			
			if(node.large_sample_size){
				bitvec &ref_xcate64 = xcate64[cate_start[cate_id] + t];
				for(int j = 0; j < iter_omp.size(); ++j){
					const bitvec &ref_not_y64 = not_y64_omp[j];
					const bitvec &ref_y64 = y64_omp[j];
					const bitvec &ref_sample_id64 = node.sample_id64[j];
					for(int l = 0; l < iter_omp[j].size(); ++l){
						int k = iter_omp[j][l];
						uint64 u = ref_sample_id64[k] & ref_xcate64[k];
						tab.nctrl += PopCount(u & ref_not_y64[k]);
						tab.ncase += PopCount(u & ref_y64[k]);
					}
				}
			}else{
				for(int k = 0; k < node.sample_id.size(); ++k){
					int sample_id = node.sample_id[k];
					int c = xcate_int[cate_id][sample_id];
					if(c == tab.code){//fix this bug in v0.2.x when member large_sample_size is included
						if(y[sample_id]){//case
							tab.ncase += id_cnt_omp[sample_id];
						}else{
							tab.nctrl += id_cnt_omp[sample_id];
						}
					}
				}
			}
			
			tab.N = tab.ncase + tab.nctrl;
			tab.r = tab.N ? (tab.ncase * 1.0 / tab.N) : (-1.0);
		}
		
		sort(table.begin(), table.end());
		
		int ct[3][3] = {{0}};
		ct[0][2] = node.nctrl;
		ct[1][2] = node.ncase;
		ct[2][2] = node.N;
		
		double cate_stat = -1.0;
		int cl = -1;
		int nl0 = 0, nl1 = 0, Nl = 0, nr0 = 0, nr1 = 0, Nr = 0;
		
		for(int t = 0; t < table.size() - 1; ++t){
			if(table[t].r < .0){
				continue;
			}
			ct[0][0] += table[t].nctrl;
			ct[1][0] += table[t].ncase;
			ct[0][1] = ct[0][2] - ct[0][0];
			ct[1][1] = ct[1][2] - ct[1][0];
			ct[2][0] = ct[0][0] + ct[1][0];
			ct[2][1] = ct[0][1] + ct[1][1];
			
			double stat = -1.0;
			if(ct[2][0] == 0 || ct[2][1] == 0){
				stat = -1.0;
				continue;
			}else{
				//different weights for cases and controls, then weighted gini impurity is computed
				double p_left, p_right, p_par, p_left0, p_left1, p_right0, p_right1, p_par0, p_par1;
				p_left = ct[0][0] * 1.0 / nctrl + ct[1][0] * class_weight / ncase;//Pr(A_L): probability of a sample falling into the left child
				p_right = ct[0][1] * 1.0 / nctrl + ct[1][1] * class_weight / ncase;//Pr(A_R): probability of a sample falling into the right child
				p_par = p_left + p_right;//Pr(A): probability of a sample falling into the parent node
				p_left0 = ct[0][0] * 1.0 / nctrl / p_left;//Pr(y=0 | A_L): probability of a sample in the left child is a control
				p_left1 = 1.0 - p_left0;//Pr(y=1 | A_L)
				p_right0 = ct[0][1] * 1.0 / nctrl / p_right;//Pr(y=0 | A_R)
				p_right1 = 1.0 - p_right0;//Pr(y=1 | A_R)
				p_par0 = ct[0][2] * 1.0 / nctrl / p_par;//Pr(y=0 | A):probability of a sample in the parent node is a control
				p_par1 = 1.0 - p_par0;//pr(y=1 | A)
				
				stat = p_par * (1.0-p_par0*p_par0-p_par1*p_par1)-p_left*(1.0-p_left0*p_left0-p_left1*p_left1)-p_right*(1.0-p_right0*p_right0-p_right1*p_right1);
				
				if(stat < .0){
					continue;
				}
			}
			
			if(stat > .0 && stat - cate_stat > MIN_INC){
				cate_stat = stat;
				cl = t;
				nl0 = ct[0][0];
				nl1 = ct[1][0];
				Nl = ct[2][0];
				nr0 = ct[0][1];
				nr1 = ct[1][1];
				Nr = ct[2][1];
			}
		}
		
		if(cate_stat > .0 && cate_stat - max_split_stat > MIN_INC){
			split_by_cate = true;
			split_by_snp = false;
			split_by_cont = false;
			max_cate_stat = cate_stat;
			max_split_stat = cate_stat;
			split_cate_id = cate_id;
			split_cate_left_code.clear();
			split_cate_right_code.clear();
			
			for(int j = 0; j < table.size(); ++j){
				if(j <= cl){
					split_cate_left_code.push_back(table[j].code);
				}else{
					split_cate_right_code.push_back(table[j].code);
				}
			}
			
			nctrl_left = nl0;
			ncase_left = nl1;
			N_left = Nl;
			nctrl_right = nr0;
			ncase_right = nr1;
			N_right = Nr;
		}
		
	}
	
	if(split_by_snp || split_by_cont || split_by_cate){
		
		node.split_stat = max_split_stat;
		if(split_by_snp){//split node by a SNP
			split_by_cont = false;
			split_by_cate = false;
			node.split_by = CODE_SNP;
			
			//node.split_snp = snp_name[split_snp_id];
			node.split_snp_id = split_snp_id;
			node.split_snp_left_val = split_snp_left_val;
			
			node.split_cont_id = -1;
			node.split_cont_lower = -1;
			node.split_cont_upper = -1;
			node.split_cont_thr = -1e20;
			
			node.split_cate_id = -1;
			node.split_cate_left_code.clear();
			node.split_cate_right_code.clear();
		}else if(split_by_cate){//split node by a categorical covariate
			split_by_cont = false;
			split_by_snp = false;
			node.split_by = CODE_CATE;
			
			node.split_snp_id = -1;
			node.split_snp_left_val = -1;
			node.split_cont_id = -1;
			
			node.split_cont_lower = -1;
			node.split_cont_upper = -1;
			node.split_cont_thr = -1e20;
			
			//node.split_cate = cate_var_name[split_cate_id];
			node.split_cate_id = split_cate_id;
			node.split_cate_left_code = split_cate_left_code;
			node.split_cate_right_code = split_cate_right_code;
		}else if(split_by_cont){//split node by a continuous covariate
			split_by_snp = false;
			split_by_cate = false;
			node.split_by = CODE_CONT;
			
			node.split_snp_id = -1;
			node.split_snp_left_val = -1;
			
			//node.split_cont = cont_var_name[split_cont_id];
			node.split_cont_id = split_cont_id;
			node.split_cont_lower = split_cont_lower;
			node.split_cont_upper = split_cont_upper;
			node.split_cont_thr = (cont_unique[split_cont_id][node.split_cont_lower] + cont_unique[split_cont_id][node.split_cont_upper])/2.0;
			
			node.split_cate_id = -1;
			node.split_cate_left_code.clear();
			node.split_cate_right_code.clear();
		}else{//impossible
			cout << "Error" << endl;
			exit(1);
		}
		
		node.ncase_left = ncase_left;
		node.nctrl_left = nctrl_left;
		node.N_left = N_left;
		node.ncase_right = ncase_right;
		node.nctrl_right = nctrl_right;
		node.N_right = N_right;
		
		node.processed = true;
		node.terminal = false;
		return !node.terminal;
		
	}else{//cannot find any variable to split this node
		node.processed = true;
		node.terminal = true;
		return !node.terminal;
	}
	
	
}


void BAMBOO::DestroyTree(vector<NODE*> &pointer){

	for(int i = 0; i < pointer.size(); ++i){
		if(pointer[i]){
			delete pointer[i];
			pointer[i] = NULL;
		}
	}

}


void BAMBOO::GrowTree(drand48_data &buf, const int tree_id){
	
//	if(!model[tree_id].empty()){
//		model[tree_id].clear();
//	}
	
	check_in[tree_id] = true;
	
	vector<bitvec > y64_omp;
	vector<bitvec > not_y64_omp;
	vector<bitvec > in_sample64_omp;
	vector<vector<int> > iter_omp;
	vector<int> oob_id_omp;
	vector<int> not_oob_id_omp;
	vector<int> id_cnt_omp;
	int ncase_omp;
	
	//generate bootstrap samples for a tree
	//update class members oob_id64 and ib_id64, which will be used in prediction
	Bootstrap(tree_id, buf, y64_omp, not_y64_omp, in_sample64_omp, iter_omp, 
	oob_id_omp, not_oob_id_omp, id_cnt_omp, ncase_omp);
	
	//note that for binary tree, the sample size in a node will decrease dramatically
	//in that case, the contingency table should be created by counting directly rather than the Boolean operations
	//sample_threshold is a rough estimation of the computing complexity of Boolean operations given the sample size
	int sample_threshold = 0;
	for(int j = 0; j < iter_omp.size(); ++j){
		sample_threshold += iter_omp[j].size();
	}
	
	vector<NODE*> all_allocate_node;//store the pointers of all allocated node
	all_allocate_node.reserve(min(nsub, 500));
	vector<NODE*> pointer_leaf;//store the pointers of all leaves
	pointer_leaf.reserve(max_nleaf);
	vector<int> cont_lower(ncont, 0);//searching lower bound (index) of continuous covariates
	vector<int> cont_upper(ncont, -1);//searching upper bound (index) of continuous covariates
	for(int j = 0; j < ncont; ++j){
		cont_upper[j] = cont_cutpnt[j].size();
	}
	
	int node_id = 0;//root
	NODE *root = new NODE(ncase, nctrl, in_sample64_omp, not_oob_id_omp, 
	cont_lower, cont_upper, tree_id, node_id, false, false, class_weight, 
	ncase_omp, nsub-ncase_omp, sample_threshold, NULL, NULL, NULL);
	pointer_leaf.push_back(root);
	all_allocate_node.push_back(root);
	
	while(pointer_leaf.size() < max_nleaf){
		
		if(false) cout << "number of leaf: " << pointer_leaf.size() << endl;
		
		double max_stat = -1.0;
		bool valid = false;
		int sel_leaf_id = -1;
		
		for(int i = 0; i < pointer_leaf.size(); ++i){
			
			NODE &node = *(pointer_leaf[i]);
			
			if(SplitNode(node, y64_omp, not_y64_omp, id_cnt_omp, iter_omp, buf)){
				valid = true;
				if(node.split_stat > .0 && node.split_stat - max_stat > MIN_INC){
					max_stat = node.split_stat;
					sel_leaf_id = i;
				}
			}
			
		}
		
		if(!valid){
			if(trace){
				cout << "All the nodes cannot be further split" << endl;
			}
			break;
		}
		
		
		//a leaf has been selected to be split
		NODE &node = (*(pointer_leaf[sel_leaf_id]));
		
		//cout << "split node " << node.node_id+1 << " with sample size " << node.N << endl;
		
		//print details
		if(node.split_by == CODE_SNP){//split by SNP
			
			if(trace){
				cout << "|| Split node " << node.node_id + (1) << " to nodes (" << node_id +1 + (1) << ", " 
				<< node_id+2 + (1) << ") by SNP " << node.split_snp_id + (1) << " (";
				for(int ii = 0; ii < 3; ++ii){
					cout << ii;
					if(ii == node.split_snp_left_val){
						cout << " |";
					}
					if(ii < 2){
						cout << " ";
					}
				}
				cout << ") Max stat = " << max_stat << endl;
			}
			
		}else if(node.split_by == CODE_CONT){
			
			if(trace){
				cout << "|| Split node " << node.node_id + (1) << " to nodes (" << node_id +1 + (1) << ", " 
				<< node_id+2 + (1) << ") by continuous covariate " << node.split_cont_id + (1) << " (" 
				<< node.split_cont_thr << ") Max stat = " << max_stat << endl;
			}
			
		}else if(node.split_by == CODE_CATE){
			
			if(trace){
				cout << "|| Split node " << node.node_id + (1) << " to nodes (" << node_id +1 + (1) << ", " 
				<< node_id+2 + (1) << ") by categorical covariate " << node.split_cate_id + (1) << " (";
				for(int ii = 0; ii < node.split_cate_left_code.size(); ++ii){
					cout << cate_unique[node.split_cate_id][node.split_cate_left_code[ii]] << " ";
				}
				cout << "| ";
				for(int ii = 0; ii < node.split_cate_right_code.size(); ++ii){
					cout << cate_unique[node.split_cate_id][node.split_cate_right_code[ii]] << " ";
				}
				cout << ") Max stat = " << max_stat << endl;
			}
		}else{//impossible
			cout << "Error: Invalid splitting code" << endl;
			exit(1);
		}
		
		//updated the list of var used in this tree
		int used_var;
		if(node.split_by == CODE_SNP){
			used_var = node.split_snp_id;
		}else if(node.split_by == CODE_CONT){
			used_var = node.split_cont_id + nsnp;
		}else if(node.split_by == CODE_CATE){
			used_var = node.split_cate_id + nsnp + ncont;
		}else{
			cout << "Error: Invalid used var_id" << endl;
			exit(1);
		}
		
		var_id_used_in_tree_long[used_var][bitloc_tree[tree_id]] |= MASK_offset_tree[tree_id];
		var_id_used_in_tree_wide[tree_id][bitloc_var[used_var]] |= MASK_offset_var[used_var];
		
		////////
		
//		vector<double> md (21, -1);
//		md[0] = pointer_leaf.size();
//		md[1] = node.node_id+1;
//		md[2] = node_id+2;
//		md[3] = node_id+3;
//		md[4] = 0;
//		md[5] = (!node.ncase_left || !node.nctrl_left) ? 1 : 0;
//		md[6] = (!node.ncase_right || !node.nctrl_right) ? 1 : 0;
//		md[7] = node.ncase;
//		md[8] = node.ncase_left;
//		md[9] = node.ncase_right;
//		md[10] = node.nctrl;
//		md[11] = node.nctrl_left;
//		md[12] = node.nctrl_right;
//		md[13] = node.N;
//		md[14] = node.N_left;
//		md[15] = node.N_right;
//		md[16] = ((node.split_by == CODE_SNP) ? (node.split_snp_id + 1) : -1);
//		md[17] = ((node.split_by == CODE_SNP) ? node.split_snp_left_val : -1);
//		md[18] = ((node.split_by == CODE_CONT) ? (node.split_cont_id + 1) : -1);
//		md[19] = ((node.split_by == CODE_CONT) ? node.split_cont_thr : -1);
//		md[20] = ((node.split_by == CODE_CATE) ? (node.split_cate_id + 1) : -1);
//		model[tree_id].push_back(md);
		
		//prepare for the child nodes
		int max_cnt = 0;
		vector<bitvec > sample_id64_child1;
		vector<bitvec > sample_id64_child2;
		vector<int> sample_id_child1, sample_id_child2;
		sample_id_child1.reserve(node.N_left);
		sample_id_child2.reserve(node.N_right);
		
		if(node.large_sample_size){//save memory
			max_cnt = node.sample_id64.size();
			sample_id64_child1 = node.sample_id64;
			sample_id64_child2 = node.sample_id64;
		}
		
		if(node.split_by == CODE_SNP){//split by a SNP
			
			if(node.large_sample_size){
				
				if(node.split_snp_left_val == 0){
					for(int k = 0; k < max_cnt; ++k){
						for(int j = 0; j < nblock; ++j){
							uint64 u = node.sample_id64[k][j];
							if(u){
								sample_id64_child1[k][j] = u & geno64[node.split_snp_id * 3][j];
								sample_id64_child2[k][j] = u & (geno64[node.split_snp_id * 3 + 1][j] | geno64[node.split_snp_id * 3 + 2][j]);
							}
						}
					}
				}else if(node.split_snp_left_val == 1){
					for(int k = 0; k < max_cnt; ++k){
						for(int j = 0; j < nblock; ++j){
							uint64 u = node.sample_id64[k][j];
							if(u){
								sample_id64_child1[k][j] = u & (geno64[node.split_snp_id * 3][j] | geno64[node.split_snp_id * 3 + 1][j]);
								sample_id64_child2[k][j] = u & geno64[node.split_snp_id * 3 + 2][j];
							}
						}
					}
				}else{//impossible
					cout << "Error: children" << endl;
					exit(1);
				}
			}
				
			for(int k = 0; k < node.sample_id.size(); ++k){
				int sample_id = node.sample_id[k];
				if(node.split_snp_left_val == 0){
					if(geno64[node.split_snp_id * 3][bitloc[sample_id]] & MASK_offset[sample_id]){//snp == 0
						sample_id_child1.push_back(sample_id);
					}else{//snp == 1, 2
						sample_id_child2.push_back(sample_id);
					}
				}else if(node.split_snp_left_val == 1){
					if(geno64[node.split_snp_id * 3 + 2][bitloc[sample_id]] & MASK_offset[sample_id]){//snp == 2
						sample_id_child2.push_back(sample_id);
					}else{//snp == 0, 1
						sample_id_child1.push_back(sample_id);
					}
				}else{//impossible
					cout << "Error: children" << endl;
					exit(1);
				}
			}
			
		}else if(node.split_by == CODE_CONT){//split by continuous variable
			
			if(false) cout << "split by cont" << endl;
			
			for(int k = 0; k < node.sample_id.size(); ++k){
				int sample_id = node.sample_id[k];
				if(xcont[node.split_cont_id][sample_id] <= node.split_cont_thr){//child1, left
					sample_id_child1.push_back(sample_id);
					if(node.large_sample_size){
						for(int j = 0; j < max_cnt; ++j){
							sample_id64_child2[j][bitloc[sample_id]] &= ~(MASK_offset[sample_id]);
						}
					}
				}else{//child2, right
					sample_id_child2.push_back(sample_id);
					if(node.large_sample_size){
						for(int j = 0; j < max_cnt; ++j){
							sample_id64_child1[j][bitloc[sample_id]] &= ~(MASK_offset[sample_id]);
						}
					}
				}
			}
		}else if(node.split_by == CODE_CATE){//split by categorical variable
			
			if(node.large_sample_size){
				bitvec v64 (nblock, (uint64) 0);
				int cate_id = node.split_cate_id;
				for(int j = 0; j < nblock; ++j){
					for(int t = 0; t < node.split_cate_left_code.size(); ++t){
						v64[j] |= xcate64[cate_start[cate_id]+node.split_cate_left_code[t]][j];
					}
				}
				
				for(int k = 0; k < max_cnt; ++k){
					for(int j = 0; j < nblock; ++j){
						sample_id64_child1[k][j] = node.sample_id64[k][j] & v64[j];
						sample_id64_child2[k][j] = node.sample_id64[k][j] & (~(v64[j]));//could not be a problem
					}
				}
				
				for(int k = 0; k < node.sample_id.size(); ++k){
					int sample_id = node.sample_id[k];
					if(sample_id64_child1[0][bitloc[sample_id]] & MASK_offset[sample_id]){
						sample_id_child1.push_back(sample_id);
					}else{
						sample_id_child2.push_back(sample_id);
					}
				}
			}else{
				
				for(int k = 0; k < node.sample_id.size(); ++k){
					int sample_id = node.sample_id[k];
					int c = xcate_int[node.split_cate_id][sample_id];
					bool to_side = false;
					if(node.split_cate_left_code.size() <= node.split_cate_right_code.size()){
						for(int ii = 0; ii < node.split_cate_left_code.size(); ++ii){
							if(c == node.split_cate_left_code[ii]){
								to_side = true;
								sample_id_child1.push_back(sample_id);
								break;
							}
						}
						
						if(!to_side){
							for(int ii = 0; ii < node.split_cate_right_code.size(); ++ii){
								if(c == node.split_cate_right_code[ii]){
									to_side = true;
									sample_id_child2.push_back(sample_id);
									break;
								}
							}
						}
					}else{
						for(int ii = 0; ii < node.split_cate_right_code.size(); ++ii){
							if(c == node.split_cate_right_code[ii]){
								to_side = true;
								sample_id_child2.push_back(sample_id);
								break;
							}
						}
						
						if(!to_side){
							for(int ii = 0; ii < node.split_cate_left_code.size(); ++ii){
								if(c == node.split_cate_left_code[ii]){
									to_side = true;
									sample_id_child1.push_back(sample_id);
									break;
								}
							}
						}
					}
					
					if(!to_side){
						cout << "Error: Unknown level (" << c << ") of categorical variable " << node.split_cate_id + 1 << " detected" << endl;
						exit(1);
					}
				}
				
			}
		}else{//impossible
			cout << "Error: Invalid splitting code" << endl;
			exit(1);
		}
		
		node_id++;
		
		vector<int> cl1 = node.cont_lower;
		vector<int> cl2 = node.cont_lower;
		vector<int> cu1 = node.cont_upper;
		vector<int> cu2 = node.cont_upper;
		
		if(node.split_by == CODE_CONT){
			cu1[node.split_cont_id] = node.split_cont_lower;
			cl2[node.split_cont_id] = node.split_cont_upper;
		}
		
		NODE *pchild1 = new NODE(ncase, nctrl, sample_id64_child1, sample_id_child1, 
		cl1, cu1, tree_id, node_id, false, false, class_weight, node.ncase_left, 
		node.nctrl_left, sample_threshold, &node, NULL, NULL);
		
		node_id++;
		NODE *pchild2 = new NODE(ncase, nctrl, sample_id64_child2, sample_id_child2, 
		cl2, cu2, tree_id, node_id, false, false, class_weight, node.ncase_right, 
		node.nctrl_right, sample_threshold, &node, NULL, NULL);
		node.child1 = pchild1;
		node.child2 = pchild2;
		
		all_allocate_node.push_back(pchild1);
		all_allocate_node.push_back(pchild2);
		
		pointer_leaf.erase(pointer_leaf.begin()+sel_leaf_id);
		
		pointer_leaf.push_back(pchild1);
		pointer_leaf.push_back(pchild2);
		
	}
	
	
	//refine the model structure
	//make sure that all leaves are correctly tagged
	for(int i = 0; i < pointer_leaf.size(); ++i){
		NODE &node = *(pointer_leaf[i]);
		node.processed = true;
		node.terminal = true;
		node.leaf_id = i;
//		for(int j = 0; j < model[tree_id].size(); ++j){
//			if(model[tree_id][j][2] == node.node_id+1){
//				if(!model[tree_id][j][5]){
//					model[tree_id][j][5] = 1;
//					if(false){
//						cout << "Warning: Potential issue in the program (code: 005)" << endl;
//					}
//				}
//			}
//			if(model[tree_id][j][3] == node.node_id+1){
//				if(!model[tree_id][j][6]){
//					model[tree_id][j][6] = 1;
//					if(false){
//						cout << "Warning: Potential issue in the program (code: 005)" << endl;
//					}
//				}
//			}
//		}
	}
	
	//prepare to save model to local file
	//prevent using pointers in saved object
	//I will create a vector of MINI_NODE to save each bamboo
	//using integer indexes to specify the relationship between mini nodes
	for(int i = 0; i < all_allocate_node.size(); ++i){
		NODE &node = *(all_allocate_node[i]);
		node.iparent = node.parent ? (*(node.parent)).node_id : -1;
		node.ichild1 = node.child1 ? (*(node.child1)).node_id : -1;
		node.ichild2 = node.child2 ? (*(node.child2)).node_id : -1;
	}
	
	if(!bamboo[tree_id].empty()){
		bamboo[tree_id].clear();
	}
	
	bamboo[tree_id] = vector<MINI_NODE> (all_allocate_node.size());
	for(int i = 0; i < all_allocate_node.size(); ++i){
		MINI_NODE mini_node(*(all_allocate_node[i]));
		bamboo[tree_id][i] = mini_node;
		if(!( bamboo[tree_id][i] == mini_node )){
			cout << i+1 << " node not equal" << endl;
		}
	}
	
	/////
	
	DestroyTree(all_allocate_node);
	if(trace){
		cout << "The tree is grown and destroyed" << endl;
	}
	
	//check_out == true if a thead has create bamboo for the specified tree_id
	//with this tag, we can update the progress bar
	check_out[tree_id] = true;
	
	//count how many variables are used in this tree
	num_var_in_tree[tree_id] = 0;
	for(int i = 0; i < var_id_used_in_tree_wide[tree_id].size(); ++i){
		num_var_in_tree[tree_id] += PopCount(var_id_used_in_tree_wide[tree_id][i]);
	}
	
}


void BAMBOO::PrintProgress(const time_t elapsed_time){
	
	int completed_jobs = 0;
	for(int i = 0; i < check_out.size(); ++i){
		completed_jobs += check_out[i];
	}
	int prg = (int) (completed_jobs * 1.0 / ntree * 100);
	
	cout << "\r\033[K\033[?25l\033[0m|";
	for(int i = 1; i <= 100; ++i){
		if(i <= prg && i % 2 == 0){
			cout << "\033[?25l\033[47m\033[1m ";
		}
		if(i == prg){
			cout << "\033[0m " << prg << "% ";
		}
		if(i > prg && i % 2 == 0){
			cout << " ";
		}
	}
	cout << "|  ETA ~";
	
	if(completed_jobs > 0){
		time_t rem = (time_t) ceil(elapsed_time * 1.0 / completed_jobs * (ntree - completed_jobs));
		streamsize ss = cout.precision();
		cout << setprecision(2) << rem/3600.0 << "h\r" << flush;
		cout.precision(ss);
	}else{
		cout << " NA\r" << flush;
	}
	
}

void BAMBOO::GrowForestMultiProc(const time_t start){
	
	drand48_data buf;
	
	#pragma omp parallel num_threads(nthread) private(buf)
	{
		if(seed < 0){
			srand48_r(time(NULL) + omp_get_thread_num(), &buf);
		}else{
			srand48_r(seed + omp_get_thread_num(), &buf);
		}
		
		#pragma omp for
		for(int i = 0; i < ntree; ++i){
			
			int tid = omp_get_thread_num();
			
			GrowTree(buf, i);
			
			if(tid == 0){
				PrintProgress(time(NULL) - start);
			}
		}
	}
	
	cout << "\r\033[K\r\033[0m\033[?25h\r";
	
	time_t bamboo_fitted_time;
	time(&bamboo_fitted_time);
	cout << "Bamboo grown: " << ctime(&bamboo_fitted_time);
	
}


void BAMBOO::SaveBamboo(){
	
	if(!output_bamboo){
		return;
	}
	
	//figure out which variables are involved in the final forest
	//this list will be useful to loading test dataset for prediction
	//they are saved to local file as well
	snp_used_in_forest.clear();
	snp_id_used_in_forest.clear();
	cont_var_used_in_forest.clear();
	cont_var_id_used_in_forest.clear();
	cate_var_used_in_forest.clear();
	cate_var_id_used_in_forest.clear();
	
	for(int i = 0; i < nvar; ++i){
		for(int j = 0; j < nblock_tree; ++j){
			if(var_id_used_in_tree_long[i][j]){//var i is used in at least one bamboo
				if(i < nsnp){//a snp is used
					snp_id_used_in_forest.push_back(i);
					snp_used_in_forest.push_back(snp_name[i]);
				}else if(i < nsnp + ncont){//a continuous covariate is used
					cont_var_id_used_in_forest.push_back(i - nsnp);
					cont_var_used_in_forest.push_back(cont_var_name[i - nsnp]);
				}else if(i < nvar){
					cate_var_id_used_in_forest.push_back(i - nsnp - ncont);
					cate_var_used_in_forest.push_back(cate_var_name[i - nsnp - ncont]);
				}else{//impossible
					cout << "Error: Invalid var_id in SaveBamboo" << endl;
					exit(1);
				}
				break;
			}
		}
	}
	
	{//save forest to local file
		ofstream ofs(path_bam);
		boost::archive::binary_oarchive ar(ofs);
			
		ar & version;
		ar & nsnp & ncont & ncate;
		ar & flip;
		ar & snp_used_in_forest & snp_id_used_in_forest;
		ar & cont_var_used_in_forest & cont_var_id_used_in_forest;
		ar & cate_var_used_in_forest & cate_var_id_used_in_forest;
		ar & cate_unique & cate_code;
		ar & bamboo;
	}
	cout << "Bamboo with " << ntree << " bamboos is saved in [ " << path_bam << " ]" << endl;
	cout << snp_used_in_forest.size() << " markers";
	if(cont_var_id_used_in_forest.size() > 0){
		cout << ", " << cont_var_used_in_forest.size() << " continuous covariate(s)";
	}
	if(cate_var_used_in_forest.size() > 0){
		cout << ", " << cate_var_used_in_forest.size() << " categorical covariate(s)";
	}
	cout << " are used in fitted bamboo" << endl;
	
}


//this function test if the a sample will fall down to the leaf
void BAMBOO::PutDownTrainingSampleToTree(const int sample_id, const int tree_id, int &fall_into_node_id){
	
	if(sample_id < 0 || sample_id >= nsub){
		cout << "Error: Invalid sample ID in PutDownTrainingSampleToTree" << endl;
		exit(1);
	}
	
	vector<bool> track_path_to_left;//true if to the left
	
	const MINI_NODE *mini_node = &(bamboo[tree_id][0]);
	bool direction = true;//true if to the left, false if to the right
	while(!(*mini_node).terminal){
		if((*mini_node).split_by == CODE_SNP){
			int snp_id = (*mini_node).split_snp_id;
			int snp_left_val = (*mini_node).split_snp_left_val;
			if(snp_left_val == 0){//snp == 0 are going to the left
				direction = 
				(geno64[snp_id * 3][bitloc[sample_id]] & MASK_offset[sample_id]) ? true : false;//true if snp == 0
			}else if(snp_left_val == 1){//snp == 0 or are going to the left
				direction = (geno64[snp_id * 3 + 2][bitloc[sample_id]] & MASK_offset[sample_id]) ? false : true;//false if snp == 2
			}else{//impossible
				cout << "Error: Invalide snp_left_val when predicting samples" << endl;
				exit(1);
			}
		}else if((*mini_node).split_by == CODE_CONT){
			int cont_id = (*mini_node).split_cont_id;
			direction = (xcont[cont_id][sample_id] <= (*mini_node).split_cont_thr) ? true : false;//true if <= threshold
		}else if((*mini_node).split_by == CODE_CATE){
			int cate_id = (*mini_node).split_cate_id;
			int size_left = (*mini_node).split_cate_left_code.size();
			int size_right = (*mini_node).split_cate_right_code.size();
			int c = xcate_int[cate_id][sample_id];
			if(size_left < size_right){
				direction = false;
				for(int ii = 0; ii < size_left; ++ii){
					if(c == (*mini_node).split_cate_left_code[ii]){
						direction = true;
						break;
					}
				}
			}else{
				direction = true;
				for(int ii = 0; ii < (*mini_node).split_cate_right_code.size(); ++ii){
					if(c == (*mini_node).split_cate_right_code[ii]){
						direction = false;
						break;
					}
				}
			}
		}else{//impossible
			cout << "Error: Invalid splitting code in mini_node.split_by" << endl;
			exit(1);
		}
		
		if(direction){
			mini_node = &(bamboo[tree_id][(*mini_node).ichild1]);
			track_path_to_left.push_back(true);
		}else{
			mini_node = &(bamboo[tree_id][(*mini_node).ichild2]);
			track_path_to_left.push_back(false);
		}
	}
	
	fall_into_node_id = (*mini_node).node_id;
	
}

//fulfill pred_node_id if necessary (output_prox == true)
void BAMBOO::PredictTrainingSample(){
	
	//to save memory, we predict ALL training data only if proximity matrix is required
	if(!output_prox){
		return;
	}
	
	cout << "Predicting training data from fitted bamboo" << endl;
	
	pred_node_id = vector<vector<uint16> > (ntree, vector<uint16> (nsub, -1));
	
	if(ntree > nsub){
		#pragma omp parallel num_threads(nthread)
		{
			#pragma omp for
			for(int tree_id = 0; tree_id < ntree; ++tree_id){
				for(int sample_id = 0; sample_id < nsub; ++sample_id){
					int fall_into_node_id = -1;
					PutDownTrainingSampleToTree(sample_id, tree_id, fall_into_node_id);
					if(fall_into_node_id < 0 || fall_into_node_id >= bamboo[tree_id].size()){
						cout << "Error: Invalid node id in ComputeProximity" << endl;
						exit(1);
					}
					
					pred_node_id[tree_id][sample_id] = (uint16) fall_into_node_id;
				}
			}
		}
	}else{
		#pragma omp parallel num_threads(nthread)
		{
			#pragma omp for
			for(int sample_id = 0; sample_id < nsub; ++sample_id){
				for(int tree_id = 0; tree_id < ntree; ++tree_id){
					int fall_into_node_id = -1;
					PutDownTrainingSampleToTree(sample_id, tree_id, fall_into_node_id);
					if(fall_into_node_id < 0 || fall_into_node_id >= bamboo[tree_id].size()){
						cout << "Error: Invalid node id in ComputeProximity" << endl;
						exit(1);
					}
					
					pred_node_id[tree_id][sample_id] = (uint16) fall_into_node_id;
				}
			}
		}
	}
	
}

void BAMBOO::ComputeProximity(){
	
	if(!output_prox){
		return;
	}
	
	//pair proximity
	//nsub x nsub
	prox = vector<vector<uint16> > (nsub, vector<uint16> (nsub, (uint16) 0));//use uint16 to save memory
	
	if(nthread == 1){
		for(int i = 0; i < nsub; ++i){
			prox[i][i] = (uint16) ntree;
			for(int j = i+1; j < nsub; ++j){
				prox[i][j] = (uint16) 0;//set as 0. Not necessary but just in case
				for(int k = 0; k < ntree; ++k){
					if(pred_node_id[k][i] == pred_node_id[k][j]){//two samples fall into the same leaf
						prox[i][j] += (uint16) 1;
					}
				}
				prox[j][i] = prox[i][j];
			}
		}
	}else{
		#pragma omp parallel num_threads(nthread)
		{
			#pragma omp for
			for(int i = 0; i < nsub; ++i){
				for(int j = 0; j < nsub; ++j){
					if(i == j){
						prox[i][i] = (uint16) ntree;
						continue;
					}
					prox[i][j] = (uint16) 0;//set as 0. Not necessary but just in case
					for(int k = 0; k < ntree; ++k){
						if(pred_node_id[k][i] == pred_node_id[k][j]){//two samples fall into the same leaf
							prox[i][j] += (uint16) 1;
						}
					}
				}
			}
		}
	}
	
	ofstream file_prox;
	file_prox.open(path_prox);
	file_prox << "IND_ID\t";
	for(int j = 0; j < nsub; ++j){
		file_prox << individual_id[j];
		
		if(j == nsub-1){
			file_prox << "\n";
		}else{
			file_prox << "\t";
		}
	}
	
	for(int i = 0; i < nsub; ++i){
		file_prox << individual_id[i] << "\t";
		
		for(int j = 0; j < nsub; ++j){
			file_prox << prox[i][j];
			if(j == nsub-1){
				file_prox << "\n";
			}else{
				file_prox << "\t";
			}
		}
	}
	file_prox.close();
	
	cout << "The proximites are saved in " << path_prox << endl;
	
}


void BAMBOO::WriteOOBError(){
	
	ofstream file_err;
	file_err.open(path_err);
	file_err << "TREE_ID\tOOB_ERROR" << endl;
	for(int i = 0; i < oob_error.size(); ++i){
		file_err << i+1 << "\t" << oob_error[i] << endl;
	}
	file_err.close();
	
	cout << "The OOB errors are saved in [ " << path_err << " ]" << endl;
	
}

void BAMBOO::ComputeOOBError(){
	
	oob_error = vector<double> (ntree, .0);//cumulative oob error with increasing number of trees
	oob_prediction = vector<int> (nsub, -1);//final prediction for oob samples
	fitted_risk_ctrl = vector<double> (nsub, .0);
	fitted_risk_case = vector<double> (nsub, .0);
	
	for(int tree_id = 0; tree_id < ntree; ++tree_id){
		bool at_least_one_oob_in_forest = false;
		int eff_nsub = 0;//effective sample size (of oob) for top j trees
		
		num_vote_correct_class[tree_id] = 0;
		for(int sample_id = 0; sample_id < nsub; ++sample_id){
			if(oob_id64[tree_id][bitloc[sample_id]] & MASK_offset[sample_id]){//sample i is an oob of tree j
				int fall_into_node_id = -1;
				if(output_prox){
					fall_into_node_id = (int) pred_node_id[tree_id][sample_id];
				}else{
					PutDownTrainingSampleToTree(sample_id, tree_id, fall_into_node_id);
				}
				
				MINI_NODE &leaf = bamboo[tree_id][fall_into_node_id];
				fitted_risk_ctrl[sample_id] += leaf.node_risk_ctrl;
				fitted_risk_case[sample_id] += leaf.node_risk_case;
				
				int nc = (leaf.node_risk_case > leaf.node_risk_ctrl) ? 0 : 1;//node class predicted by risk
				if(nc == y[sample_id]){//a vote cast for the correct class
					num_vote_correct_class[tree_id] += 1;
				}
			}
			
			if(fitted_risk_ctrl[sample_id] > .0 || fitted_risk_case[sample_id] > .0){
				at_least_one_oob_in_forest = true;
				++eff_nsub;
				
				//?? is the condition here reasonable? should we used cut-off here?
				int mvc = (fitted_risk_case[sample_id] > fitted_risk_ctrl[sample_id]) ? 0 : 1;//predicted class by major vote
				if(tree_id == ntree - 1){
					oob_prediction[sample_id] = mvc;//final prediction using all tree with sample i as oob sample
				}
				
				if(mvc != y[sample_id]){//incorrect prediction
					oob_error[tree_id] += 1.0;
				}
			}else{//sample i is not oob of any of the trees in the forest
				//do something later
			}
		}
		
		oob_error[tree_id] /= eff_nsub;
		if(!at_least_one_oob_in_forest){//no oob for this whole forest
			oob_error[tree_id] = -1.0;
		}
	}
	
	WriteOOBError();
	
}

void BAMBOO::WriteOOBAUC(){
	
	ofstream file_auc;
	file_auc.open(path_auc);
	file_auc << "IND_ID\tOOB_SCORE\tTRUE_CLASS" << endl;
	for(int i = 0; i < nsub; ++i){
		file_auc << individual_id[fitted_prob[i].sample_id] << "\t" << fitted_prob[i].post_prob2 << "\t" << fitted_prob[i].true_class << endl;
	}
	file_auc.close();
	
	cout << "The posterior probabilities of OOB samples are saved in [ " << path_auc << " ]" << endl;
	cout << "The OOB AUC: " << auc << endl;
	
}

//this function is called only if ComputeOOBError has been called
void BAMBOO::ComputeOOBAUC(){
	
	fitted_prob = vector<FITTED_PROB> (nsub, fitpr0);
	for(int sample_id = 0; sample_id < nsub; ++sample_id){
		fitted_prob[sample_id].sample_id = sample_id;
		fitted_prob[sample_id].true_class = y[sample_id];
		fitted_prob[sample_id].post_prob1 = -1.0;
		fitted_prob[sample_id].post_prob2 = fitted_risk_ctrl[sample_id] / (fitted_risk_ctrl[sample_id] + fitted_risk_case[sample_id]);
		fitted_prob[sample_id].sort_by = fitted_prob[sample_id].post_prob2;
	}
	
	sort(fitted_prob.begin(), fitted_prob.end());
	
	
	vector<double> pred_rank(nsub);
	vector<double> co(nsub+2);
	vector<double> tp(nsub+2);
	vector<double> fp(nsub+2);
	
	double inc = 1.0/nsub;
	double tpinc = 1.0/ncase;
	double fpinc = 1.0/nctrl;
	
	pred_rank[0] = 1.0;
	for(int i = 2; i <= nsub; ++i){
		pred_rank[i-1] = pred_rank[i-2] - inc;
	}
	
	tp[0] = .0;
	for(int i = 1; i <= nsub; ++i){
		if(fitted_prob[i-1].true_class == 1){
			tp[i] = tp[i-1] + tpinc;
		}else{
			tp[i] = tp[i-1];
		}
	}
	tp[nsub+1] = 1.0;
	
	fp[0] = .0;
	for(int i = 1; i <= nsub; ++i){
		if(fitted_prob[i-1].true_class == 0){
			fp[i] = fp[i-1] + fpinc;
		}else{
			fp[i] = fp[i-1];
		}
	}
	fp[nsub+1] = 1.0;
	
	co[0] = 1.0;
	for(int i = 0; i < nsub; ++i){
		co[i+1] = pred_rank[i];
	}
	co[nsub+1] = .0;
	
	auc = .0;
	for(int i = 1; i <= nsub; ++i){
		auc += .5 * abs(fp[i-1]-fp[i]) * (tp[i] + tp[i-1]);
	}
	
	WriteOOBAUC();
	
}


void BAMBOO::WriteConfusionMatrix(){
	
	cout << "Confusion matrix of oob samples (rows / cols: true / pred classes)" << endl;
	cout << "       \tControl\tCase   \tError  " << endl;
	cout << "Control\t" << confusion_matrix[0][0] << "\t" << confusion_matrix[0][1] << "\t" << confusion_matrix[0][1] * 1.0 / nctrl << "\t(1-Specificity)" << endl;
	cout << "Case   \t" << confusion_matrix[1][0] << "\t" << confusion_matrix[1][1] << "\t" << confusion_matrix[1][0] * 1.0 / ncase << "\t(1-Sensitivity)" << endl;
	cout << "       \t       \t       \t";
	cout << (confusion_matrix[0][1] + confusion_matrix[1][0]) * 1.0 / (confusion_matrix[0][0] + confusion_matrix[0][1] + confusion_matrix[1][0] + confusion_matrix[1][1]) << endl;
	
	ofstream file_cof;
	file_cof.open(path_cof);
	file_cof << "OOB AUC: " << auc << endl;
	file_cof << "Confusion matrix of oob samples (rows / cols: true / pred classes)" << endl;
	file_cof << "       \tControl\tCase   \tError  " << endl;
	file_cof << "Control\t" << confusion_matrix[0][0] << "\t" << confusion_matrix[0][1] << "\t" << confusion_matrix[0][1] * 1.0 / nctrl << "\t(1-Specificity)" << endl;
	file_cof << "Case   \t" << confusion_matrix[1][0] << "\t" << confusion_matrix[1][1] << "\t" << confusion_matrix[1][0] * 1.0 / ncase << "\t(1-Sensitivity)" << endl;
	file_cof << "       \t       \t       \t";
	file_cof << (confusion_matrix[0][1] + confusion_matrix[1][0]) * 1.0 / (confusion_matrix[0][0] + confusion_matrix[0][1] + confusion_matrix[1][0] + confusion_matrix[1][1]) << endl;
	file_cof.close();
	
	cout << "The confusion matrix is saved in [ " << path_cof << " ]" << endl;
	
}


void BAMBOO::ComputeConfusionMatrix(){
	
	confusion_matrix = vector<vector<int> > (2, vector<int> (2, 0));
	for(int i = 0; i < nsub; ++i){
		if(y[i] == 0 && oob_prediction[i] == 0){
			confusion_matrix[0][0] += 1;
		}else if(y[i] == 0 && oob_prediction[i] == 1){
			confusion_matrix[0][1] += 1;
		}else if(y[i] == 1 && oob_prediction[i] == 0){
			confusion_matrix[1][0] += 1;
		}else if(y[i] == 1 && oob_prediction[i] == 1){
			confusion_matrix[1][1] += 1;
		}else{
			//if ntree is too small, there might be some oob_prediction = -1
			//don't worry about that, do nothing here
			//cout << "Warning: The confusion matrix might be inaccurate due to small --ntree" << endl;
		}
	}
	
	WriteConfusionMatrix();
	
}


void BAMBOO::ComputeGiniImportance(){
	
	for(int t = 0; t < ntree; ++t){
		for(int k = 0; k < bamboo[t].size(); ++k){
			MINI_NODE &mini_node = bamboo[t][k];
			if(mini_node.terminal){//leaf isn't counted in computing gini importance
				continue;
			}
			
			int var_id = -1;
			int split_by = mini_node.split_by;
			if(split_by == CODE_SNP){
				var_id = mini_node.split_snp_id;
			}else if(split_by == CODE_CONT){
				var_id = mini_node.split_cont_id + nsnp;
			}else if(split_by == CODE_CATE){
				var_id = mini_node.split_cate_id + nsnp + ncont;
			}else{//impossible
				cout << "Error: Invalid var_id in ComputeGiniImportance" << endl;
				exit(1);
			}
			importance[var_id].imp_gini += (float) mini_node.split_stat * mini_node.N;
		}
	}
	
	for(int i = 0; i < nvar; ++i){
		int var_id = importance[i].var_id;
		importance[var_id].imp_gini /= ntree;
		
		if(imp_measure == IMP_NULL || imp_measure == IMP_GINI){
			importance[var_id].sort_by = importance[var_id].imp_gini;
		}
	}
	
}


void BAMBOO::FisherYatesShuffle(drand48_data &buf, vector<int> &v){
	
	for(int k = v.size()-1; k > 0; --k){
		long int li;
		lrand48_r(&buf, &li);
		int j = li % (k+1);
		int u = v[j];
		v[j] = v[k];
		v[k] = u;
	}
	
}


void BAMBOO::PutDownOOBSampleToTree(const int sample_id, const int var_id, 
	const int tree_id, const int permuted_sample_id, int &fall_into_node_id){
	
	if(sample_id < 0 || sample_id >= nsub || permuted_sample_id < 0 || permuted_sample_id >= nsub){
		cout << "Error: Invalid sample ID or permuted sample ID" << endl;
		exit(1);
	}
	
	const MINI_NODE *mini_node = &(bamboo[tree_id][0]);
	bool direction = true;//true if to the left, false if to the right
	while(!(*mini_node).terminal){
		
		if((*mini_node).split_by == CODE_SNP){
			int snp_id = (*mini_node).split_snp_id;
			int snp_left_val = (*mini_node).split_snp_left_val;
			int sid = (snp_id == var_id) ? permuted_sample_id : sample_id;
			
			if(snp_left_val == 0){//snp == 0 are going to the left
				direction = 
				(geno64[snp_id * 3][bitloc[sid]] & MASK_offset[sid]) ? true : false;//true if snp == 0
			}else if(snp_left_val == 1){//snp == 0 or are going to the left
				direction = (geno64[snp_id * 3 + 2][bitloc[sid]] & MASK_offset[sid]) ? false : true;//false if snp == 2
			}else{//impossible
				cout << "Error: Invalide snp_left_val when predicting samples" << endl;
				exit(1);
			}
		}else if((*mini_node).split_by == CODE_CONT){
			int cont_id = (*mini_node).split_cont_id;
			int sid = (cont_id + nsnp == var_id) ? permuted_sample_id : sample_id;
			direction = (xcont[cont_id][sid] <= (*mini_node).split_cont_thr) ? true : false;//true if <= threshold
		}else if((*mini_node).split_by == CODE_CATE){
			int cate_id = (*mini_node).split_cate_id;
			int sid = (cate_id + nsnp + ncont == var_id) ? permuted_sample_id : sample_id;
			int size_left = (*mini_node).split_cate_left_code.size();
			int size_right = (*mini_node).split_cate_right_code.size();
			
			int c = xcate_int[cate_id][sid];
			if(size_left < size_right){
				direction = false;
				for(int ii = 0; ii < size_left; ++ii){
					if(c == (*mini_node).split_cate_left_code[ii]){
						direction = true;
						break;
					}
				}
			}else{
				direction = true;
				for(int ii = 0; ii < size_right; ++ii){
					if(c == (*mini_node).split_cate_right_code[ii]){
						direction = false;
						break;
					}
				}
			}
		}else{//impossible
			cout << "Error: Invalid splitting code in mini_node.split_by" << endl;
			exit(1);
		}
		
		if(direction){
			mini_node = &(bamboo[tree_id][(*mini_node).ichild1]);
		}else{
			mini_node = &(bamboo[tree_id][(*mini_node).ichild2]);
		}
		
	}
	
	fall_into_node_id = (*mini_node).node_id;
	
}


void BAMBOO::ComputePermutationImportance(){
	
	if(imp_measure == IMP_NULL || imp_measure == IMP_GINI){
		return;
	}
	
	cout << "Calculating permutation importance" << endl;
	
	drand48_data buf;
	
	#pragma omp parallel num_threads(nthread) private(buf)
	{
		
		if(seed < 0){
			srand48_r(time(NULL) + omp_get_thread_num(), &buf);
		}else{
			srand48_r(seed + omp_get_thread_num(), &buf);
		}
		
		#pragma omp for
		for(int t = 0; t < ntree; ++t){
			int tree_id = t;
			vector<int> original_oob_id;
			original_oob_id.reserve(num_oob[tree_id]);
			for(int j = 0; j < nsub; ++j){
				if(original_oob_id.size() > num_oob[tree_id]){//save computing time
					break;
				}
				if(oob_id64[tree_id][bitloc[j]] & MASK_offset[j]){//sample j is an oob of tree t
					original_oob_id.push_back(j);
				}
			}
			vector<int> permuted_oob_id = original_oob_id;
			FisherYatesShuffle(buf, permuted_oob_id);
			
			if(!num_vote_loss_permuted_correct_class[tree_id].empty()){
				num_vote_loss_permuted_correct_class[tree_id].clear();
			}
			int checked_nvar = 0;//checked number of variables used in the tree
			for(int i = 0; i < nvar; ++i){
				int var_id = i;
				if(checked_nvar > num_var_in_tree[tree_id]){//no more used variables
					break;
				}
				if(var_id_used_in_tree_wide[tree_id][bitloc_var[var_id]] & MASK_offset_var[var_id]){//var i is used in tree t
					++checked_nvar;
					
					uint16 delta = num_vote_correct_class[tree_id];
					for(int k = 0; k < original_oob_id.size(); ++k){//predict every oob of the tree
						int sample_id = original_oob_id[k];
						int permuted_sample_id = permuted_oob_id[k];
						
						bool change = true;
						if(var_id < nsnp && (!(geno64[var_id*3][bitloc[sample_id]] & MASK_offset[sample_id]) == !(geno64[var_id*3][bitloc[permuted_sample_id]] & MASK_offset[permuted_sample_id])) 
							&& (!(geno64[var_id*3+1][bitloc[sample_id]] & MASK_offset[sample_id]) == !(geno64[var_id*3+1][bitloc[permuted_sample_id]] & MASK_offset[permuted_sample_id])) 
							&& (!(geno64[var_id*3+2][bitloc[sample_id]] & MASK_offset[sample_id]) == !(geno64[var_id*3+2][bitloc[permuted_sample_id]] & MASK_offset[permuted_sample_id]))){
							change = false;
						}else{
							change = true;
						}
						
						int fall_into_node_id = -1;
						if(change){
							PutDownOOBSampleToTree(sample_id, var_id, tree_id, permuted_sample_id, fall_into_node_id);
						}else{
							if(output_prox){
								fall_into_node_id = (int) pred_node_id[tree_id][sample_id];
							}else{
								PutDownTrainingSampleToTree(sample_id, tree_id, fall_into_node_id);
							}
						}
						
						if(fall_into_node_id >= 0 && fall_into_node_id < bamboo[tree_id].size()){
							MINI_NODE &leaf = bamboo[tree_id][fall_into_node_id];
							int nc = (leaf.node_risk_case > leaf.node_risk_ctrl) ? 0 : 1;//node class predicted by risk
							if(nc == y[sample_id]){//Hey! It is NOT a bug using "==" rather than "!=" here
								delta -= (uint16) 1;
							}
						}else{
							cout << "Error: Cannot assign sample " << sample_id + 1 << " to any leaf of the tree " << tree_id + 1 << endl;
							exit(1);
						}
						
					}
					
					num_vote_loss_permuted_correct_class[tree_id].push_back(delta);
					
				}
			}
		}
		
	}
	
	//initialization
	
	#pragma omp parallel num_threads(nthread)
	{
		#pragma omp for
		for(int i = 0; i < nvar; ++i){
			importance[i].imp_raw = .0f;
			importance[i].imp_sd_breiman_cutler = .0f;
			importance[i].imp_sd_liaw_wiener = .0f;
		}
	}
	
	
	//cannot be parallelized
	for(int t = 0; t < ntree; ++t){
		int tree_id = t;
		int checked_nvar = 0;
		for(int i = 0; i < nvar; ++i){
			int var_id = i;
			if(checked_nvar > num_var_in_tree[tree_id]){//no more used variables
				break;
			}
			
			if(var_id_used_in_tree_wide[tree_id][bitloc_var[var_id]] & MASK_offset_var[var_id]){//var i is used in tree t
				++checked_nvar;
				
				float delta = num_vote_loss_permuted_correct_class[tree_id][checked_nvar-1] * 1.0 / num_oob[tree_id];
				importance[var_id].imp_raw += delta;
				importance[var_id].imp_sd_breiman_cutler += delta * delta;
				importance[var_id].imp_sd_liaw_wiener += delta * delta * num_oob[tree_id];
			}
		}
	}
	
	#pragma omp parallel num_threads(nthread)
	{
		#pragma omp for
		for(int i = 0; i < nvar; ++i){
			int var_id = i;
			
			importance[var_id].imp_raw /= ntree;
			importance[var_id].imp_sd_breiman_cutler = sqrt((importance[var_id].imp_sd_breiman_cutler / ntree - importance[var_id].imp_raw * importance[var_id].imp_raw) / ntree);
			importance[var_id].imp_sd_liaw_wiener = sqrt((importance[var_id].imp_sd_liaw_wiener / ntree - importance[var_id].imp_raw * importance[var_id].imp_raw) / ntree);
			importance[var_id].imp_breiman_cutler = (importance[var_id].imp_sd_breiman_cutler > 1e-12) ? (importance[var_id].imp_raw / importance[var_id].imp_sd_breiman_cutler) : (-999999.0);
			importance[var_id].imp_liaw_wiener = (importance[var_id].imp_sd_liaw_wiener > 1e-12) ? (importance[var_id].imp_raw / importance[var_id].imp_sd_liaw_wiener) : (-999999.0);
			
			if(imp_measure == IMP_BREIMAN_CUTLER){
				importance[var_id].sort_by = importance[var_id].imp_breiman_cutler;
			}else if(imp_measure == IMP_LIAW_WIENER){
				importance[var_id].sort_by = importance[var_id].imp_liaw_wiener;
			}else{
				cout << "Error: Invalid importance code" << endl;
				exit(1);
			}
		}
	}
	
}


void BAMBOO::WriteImportance(){
	
	if(!output_imp){
		return;
	}
	
	sort(importance.begin(), importance.end());
	
	ofstream file_imp;
	file_imp.open(path_imp);
	file_imp << "INDEX\tVAR_ID\tVAR_NAME\tGINI_IMP";
	if(imp_measure == IMP_BREIMAN_CUTLER){
		file_imp << "\tRAW_IMP\tBREIMAN_CUTLER_IMP\tLIAW_WIENER_IMP\tSD_BREIMAN_CUTLER\tSD_LIAW_WIENER";
	}
	file_imp << endl;
	for(int i = 0; i < nvar; ++i){
		int var_id = importance[i].var_id;
		if(var_id < nsnp){
			file_imp << i+1 << "\tSNP_" << var_id+1 << "\t" << snp_name[var_id] << "\t";
		}else if(var_id < nsnp + ncont){
			file_imp << i+1 << "\tCON_" << var_id+1-nsnp << "\t" << cont_var_name[var_id-nsnp] << "\t";
		}else if(var_id < nvar){
			file_imp << i+1 << "\tCAT_" << var_id+1-nsnp-ncont << "\t" << cate_var_name[var_id-nsnp-ncont] << "\t";
		}else{
			cout << "Error: Invalid var id in writing [ " << path_imp << " ]" << endl;
			exit(1);
		}
		file_imp << "\t" << importance[i].imp_gini;
		if(imp_measure == IMP_BREIMAN_CUTLER){
			file_imp << "\t" << importance[i].imp_raw << "\t";
			
			if(importance[i].imp_sd_breiman_cutler < -99999.0){
				file_imp << "NA" << "\t";
			}else{
				file_imp << importance[i].imp_breiman_cutler << "\t";
			}
			if(importance[i].imp_sd_liaw_wiener < -99999.0){
				file_imp << "NA" << "\t";
			}else{
				file_imp << importance[i].imp_liaw_wiener << "\t";
			}
			
			if(importance[i].imp_sd_breiman_cutler < -99999.0){
				file_imp << "NA" << "\t";
			}else{
				file_imp << importance[i].imp_sd_breiman_cutler << "\t";
			}
			
			if(importance[i].imp_sd_liaw_wiener < -99999.0){
				file_imp << "NA";
			}else{
				file_imp << importance[i].imp_sd_liaw_wiener;
			}
		}
		file_imp << endl;
	}
	file_imp.close();
	
	cout << "The variable importances are saved in [ " << path_imp << " ]" << endl;
	
}


void BAMBOO::ComputeImportance(){
	
	if(!output_imp){
		return;
	}
	
	ComputeGiniImportance();
	
	ComputePermutationImportance();
	
	WriteImportance();
	
}


void BAMBOO::GrowForest(){
	
	if(pred_from_specified_model){
		return;
	}
	
	time_t start, end;
	start = time(NULL);
	
	if(nthread > sysconf( _SC_NPROCESSORS_ONLN ) || nthread < 1){
		cout << "Error: Invalid number of threads (" << nthread << ")" << endl;
		exit(1);
	}
	
	GrowForestMultiProc(start);
	
	end = time(NULL);
	cout << "Elapsed time: " << (end - start) / 3600 << "h " << ((end - start) % 3600) / 60 << "m " << (end - start) % 60 << "s" << endl;
	
	SaveBamboo();
	
	//must call this function first once bamboo is grown and saved
	PredictTrainingSample();
	
	//and then
	ComputeProximity();
	
	ComputeOOBError();
	
	ComputeOOBAUC();
	
	ComputeConfusionMatrix();
	
	ComputeImportance();
	
	time_t complete_time;
	time(&complete_time);
	cout << "Training accomplished: " << ctime(&complete_time);
	
	end = time(NULL);
	cout << "Elapsed time: " << (end - start) / 3600 << "h " << ((end - start) % 3600) / 60 << "m " << (end - start) % 60 << "s" << endl;
	
}


void BAMBOO::LoadBamboo(){
	
	if(pred_from_trained_model){//bamboo should be already in the memory
		;
	}else{
		
		if(pred_from_specified_model){//load forest from local file
			ifstream ifs(path_bam);
			if(!ifs){
				cout << "Error: Cannot open BAM file " << path_bam << endl;
				exit(1);
			}
			
			boost::archive::binary_iarchive ar(ifs);
				
			ar & version;
			ar & nsnp & ncont & ncate;
			ar & flip;
			ar & snp_used_in_forest & snp_id_used_in_forest;
			ar & cont_var_used_in_forest & cont_var_id_used_in_forest;
			ar & cate_var_used_in_forest & cate_var_id_used_in_forest;
			ar & cate_unique & cate_code;
			ar & bamboo;
			
			has_cont = cont_var_used_in_forest.size() ? true : false;
			has_cate = cate_var_used_in_forest.size() ? true : false;
			
		}else{
			cout << "Error: Cannot load bamboo for prediction" << endl;
			exit(1);
		}
	}
	
	cout << "The model with " << bamboo.size() << " bamboos is engaged from [ " << path_bam << " ] for prediction" << endl;
	cout << snp_used_in_forest.size() << " markers";
	if(cont_var_id_used_in_forest.size() > 0){
		cout << ", " << cont_var_used_in_forest.size() << " continuous covariate(s)";
	}
	if(cate_var_used_in_forest.size() > 0){
		cout << ", " << cate_var_used_in_forest.size() << " categorical covariate(s)";
	}
	cout << " are used in prediction" << endl;
	
}


void BAMBOO::LoadTestingData(){
	
	if(!pred_from_specified_model && !pred_from_trained_model){
		return;
	}
	
	time_t start = time(NULL);
	
	if(!specified_testid.empty()){
		specified_testid.clear();
	}
	
	if(path_testid){
		ifstream file_testid(path_testid);
		if(!file_testid){
			cout << "Error: Cannot open IID file " << path_testid << endl;
			exit(1);
		}
		
		for(string s; getline(file_testid, s); ){
			istringstream sin(s);
			string iid;
			sin >> iid;
			specified_testid.push_back(iid);
		}
		file_testid.close();
		
		if(specified_testid.size() == 0){
			cout << "Error: 0 individuals are specified in " << path_testid << endl;
			exit(1);
		}else{
			cout << specified_testid.size() << " individuals specified in [ " << path_testid << " ] are used as testing data" << endl;
		}
	}
	
	
	cout << "Aligning individuals ..." << endl;
	
	ifstream file_fam(path_fam_test);
	if(!file_fam){
		cout << "Error: Cannot open FAM file " << path_fam_test << " for prediction" << endl;
		exit(1);
	}
	
	vector<string> individual_id_fam;
	vector<int> index_individual_id_fam;
	int nsub_fam = 0;
	for(string s; getline(file_fam, s); ){
		istringstream sin(s);
		string fam_id, ind_id, pat_id, mat_id, sex;
		int phen;
		if(sin >> fam_id >> ind_id >> pat_id >> mat_id >> sex >> phen){
			++nsub_fam;
			
			bool in_spec = true;
			if(path_testid){
				in_spec = false;
				for(int i = 0; i < specified_testid.size(); ++i){
					if(ind_id == specified_testid[i]){
						in_spec = true;
						break;
					}
				}
			}else{
				in_spec = true;
			}
			
			if(in_spec){
				individual_id_fam.push_back(ind_id);
				index_individual_id_fam.push_back(nsub_fam-1);
			}
		}else{
			cout << "Error" << endl;
			exit(1);
		}
	}
	file_fam.close();
	
	//load continuous covariates
	vector<string> individual_id_cont;
	vector<int> index_individual_id_cont;
	vector<bool> cont_used;
	vector<int> cont_col;
	int nsub_cont = 0;
	
	if(has_cont){
		ifstream file_cont(path_cont_test);
		if(!file_cont){
			cout << "Error: Cannot open CON file " << path_cont_test << " for predicton" << endl;
			exit(1);
		}
		
		string missing_flag("NA");
		vector<string> header;
		int nrow = -2;
		for(string s; getline(file_cont, s); ){
			istringstream sin(s);
			++nrow;
			if(nrow == -1){
				string ind_id;
				sin >> ind_id;
				string h;
				while(sin >> h){
					header.push_back(h);
				}
				
				if(header.size() < cont_var_used_in_forest.size()){
					cout << "Error: Cannot find all necessary covariates from " << path_cont_test << endl;
					exit(1);
				}
				
				cont_used = vector<bool> (header.size(), false);
				cont_col = vector<int> (header.size(), -1);
				
				for(int i = 0; i < cont_var_used_in_forest.size(); ++i){
					string str = cont_var_used_in_forest[i];
					bool b = false;
					for(int j = 0; j < header.size(); ++j){
						if(str == header[j]){
							b = true;
							cont_used[j] = true;
							cont_col[j] = cont_var_id_used_in_forest[i];
							break;
						}
					}
					if(!b){
						cout << "Error: Cannot find covariate " << str << " in " << path_cont_test << endl;
						file_cont.close();
						exit(1);
					}
				}
			}else{
				++nsub_cont;
				bool no_missing_this_row = true;
				string iid, str;
				sin >> iid;
				int j = -1;
				while(sin >> str){
					++j;
					if(cont_used[j] && str == missing_flag){
						no_missing_this_row = false;
						break;
					}
				}
				if(no_missing_this_row){
					individual_id_cont.push_back(iid);
					index_individual_id_cont.push_back(nrow);
				}
			}
		}
		file_cont.close();
	}
	
	
	//load categorical covariates
	vector<string> individual_id_cate;
	vector<int> index_individual_id_cate;
	vector<bool> cate_used;
	vector<int> cate_col;
	int nsub_cate = 0;
	
	if(has_cate){
		ifstream file_cate(path_cate_test);
		if(!file_cate){
			cout << "Error: Cannot find " << path_cate_test << endl;
			exit(1);
		}
		
		string missing_flag("NA");
		vector<string> header;
		int nrow = -2;
		for(string s; getline(file_cate, s); ){
			istringstream sin(s);
			++nrow;
			if(nrow == -1){
				string ind_id;
				sin >> ind_id;
				string h;
				while(sin >> h){
					header.push_back(h);
				}
				
				if(header.size() < cate_var_used_in_forest.size()){
					cout << "Error: Cannot find all necessary covariates from " << path_cate_test << endl;
					exit(1);
				}
				
				cate_used = vector<bool> (header.size(), false);
				cate_col = vector<int> (header.size(), -1);
				
				for(int i = 0; i < cate_var_used_in_forest.size(); ++i){
					string str = cate_var_used_in_forest[i];
					bool b = false;
					for(int j = 0; j < header.size(); ++j){
						if(str == header[j]){
							b = true;
							cate_used[j] = true;
							cate_col[j] = cate_var_id_used_in_forest[i];
							break;
						}
					}
					if(!b){
						cout << "Error: Cannot find covariate " << str << " in " << path_cate_test << endl;
						file_cate.close();
						exit(1);
					}
				}
			}else{
				++nsub_cate;
				bool no_missing_this_row = true;
				bool no_invalid_str_this_row = true;
				string iid, str;
				sin >> iid;
				int j = -1;
				while(sin >> str){
					++j;
					if(cate_used[j] && str == missing_flag){
						no_missing_this_row = false;
						break;
					}
					
					if(cate_used[j]){//check if this variable has valid level as in the training data
						int ccj = cate_col[j];
						bool valid_str = false;
						for(int ii = 0; ii < cate_unique[ccj].size(); ++ii){
							if(str == cate_unique[ccj][ii]){
								valid_str = true;
								break;
							}
						}
						if(!valid_str){
							no_invalid_str_this_row = false;
						}
					}
					
					if(!no_invalid_str_this_row){
						break;
					}
				}
				if(no_missing_this_row && no_invalid_str_this_row){
					individual_id_cate.push_back(iid);
					index_individual_id_cate.push_back(nrow);
				}
			}
		}
		file_cate.close();
	}
	
	//determine intersection
	vector<int> index_individual_id_fam_inc;
	vector<int> index_individual_id_cont_inc;
	vector<int> index_individual_id_cate_inc;
	individual_id_test.clear();
	
	if(!has_cont && !has_cate){
		individual_id_test = individual_id_fam;
		index_individual_id_fam_inc = index_individual_id_fam;
	}else{
		if(!has_cate){
			string str;
			for(int i = 0; i < individual_id_fam.size(); ++i){
				str = individual_id_fam[i];
				for(int j = 0; j < individual_id_cont.size(); ++j){
					if(str == individual_id_cont[j]){
						individual_id_test.push_back(str);
						index_individual_id_fam_inc.push_back(index_individual_id_fam[i]);
						index_individual_id_cont_inc.push_back(index_individual_id_cont[j]);
						break;
					}
				}
			}
		}else if(!has_cont){
			string str;
			for(int i = 0; i < individual_id_fam.size(); ++i){
				str = individual_id_fam[i];
				for(int j = 0; j < individual_id_cate.size(); ++j){
					if(str == individual_id_cate[j]){
						individual_id_test.push_back(str);
						index_individual_id_fam_inc.push_back(index_individual_id_fam[i]);
						index_individual_id_cate_inc.push_back(index_individual_id_cate[j]);
						break;
					}
				}
			}
		}else{
			string str;
			for(int i = 0; i < individual_id_fam.size(); ++i){
				str = individual_id_fam[i];
				bool in_cont = false;
				int j1 = -1;
				for(int j = 0; j < individual_id_cont.size(); ++j){
					if(str == individual_id_cont[j]){
						in_cont = true;
						j1 = j;
						break;
					}
				}
				
				if(!in_cont){
					continue;
				}
				
				bool in_cate = false;
				int j2 = -1;
				for(int j = 0; j < individual_id_cate.size(); ++j){
					if(str == individual_id_cate[j]){
						in_cate = true;
						j2 = j;
						break;
					}
				}
				
				if(in_cont && in_cate){
					individual_id_test.push_back(str);
					index_individual_id_fam_inc.push_back(index_individual_id_fam[i]);
					index_individual_id_cont_inc.push_back(index_individual_id_cont[j1]);
					index_individual_id_cate_inc.push_back(index_individual_id_cate[j2]);
				}
			}
		}
	}
	
	//usable sample size in testing data
	nsub_test = individual_id_test.size();
	if(nsub_test == 0){
		cout << "Sample size of test dataset is 0. The program terminated" << endl;
		exit(1);
	}
	
	ifstream file_bim(path_bim_test);
	if(!file_bim){
		cout << "Error: Cannot open BIM file " << path_bim_test << endl;
		exit(1);
	}
	
	
	int nsnp_test = 0;
	vector<string> snp_name_test;
	for(string s; getline(file_bim, s); ){
		istringstream sin(s);
		string chr, rs, allele_name1, allele_name2;
		int genetic_dist, base_pair_pos;
		if(sin >> chr >> rs >> genetic_dist >> base_pair_pos >> allele_name1 >> allele_name2){
			++nsnp_test;
			snp_name_test.push_back(rs);
		}else{
			cout << "Error: Invalid format in " << path_bim_test << endl;
			file_bim.close();
			exit(1);
		}
	}
	file_bim.close();
	
	cout << "Aligning markers ..." << endl;
	vector<bool> snp_used(snp_name_test.size(), false);
	vector<int> snp_col(snp_name_test.size(), -1);
	
	#pragma omp parallel num_threads(nthread)
	{
		#pragma omp for
		for(int i = 0; i < snp_used_in_forest.size(); ++i){
			string str = snp_used_in_forest[i];
			bool b = false;
			for(int j = 0; j < snp_name_test.size(); ++j){
				if(str == snp_name_test[j]){
					b = true;
					snp_used[j] = true;
					snp_col[j] = snp_id_used_in_forest[i];
					break;
				}
			}
			if(!b){
				cout << "Error: Cannot find SNP " << str << " in " << path_bim_test << endl;
				exit(1);
			}
		}
	}
	
	//load genotypes by memory mapping
	
	cout << "Loading markers from [ " << path_bed_test << " ]" << endl;
	
	stringstream file_bed;
	file_bed << path_bed_test;
	int fb = open(file_bed.str().c_str(), O_RDONLY);
	if(fb < 0){
		cout << "Error: Cannot open BED file " << path_bed_test << endl;
		close(fb);
		exit(1);
	}
	
	LEN_bed = 8;
	int nblock_bed = (int) ceil((2 * (double) nsub_fam) / LEN_bed);
	
	uint8 *map;
	int file_size = (3 + nblock_bed * nsnp_test) * sizeof(uint8);
	map = (uint8*) mmap(0, file_size, PROT_READ, MAP_PRIVATE, fb, 0);
	if(map == MAP_FAILED){
		cout << "Error: Failed in memory mapping on the BED file " << path_bed_test << endl;
		close(fb);
		exit(1);
	}
	
	if(map[0] != 0x6c || map[1] != 0x1b){
		cout << "Error: " << path_bed_test << " is not a plink BED file" << endl;
		munmap(map, file_size);
		close(fb);
		exit(1);
	}
	
	if(map[2] == 0x01){
		;
	}else if(map[2] == 0x00){
		cout << "Error: bamboo only accepts BED file with SNP-major mode" << endl;
		munmap(map, file_size);
		close(fb);
		exit(1);
	}else{
		cout << "Error: Illegal BED flag" << endl;
		munmap(map, file_size);
		close(fb);
		exit(1);
	}
	
	PROBE1[0] = 0x01;
	PROBE2[0] = 0x02;
	PROBE1[1] = 0x04;
	PROBE2[1] = 0x08;
	PROBE1[2] = 0x10;
	PROBE2[2] = 0x20;
	PROBE1[3] = 0x40;
	PROBE2[3] = 0x80;
	
	MASK = 0x8000000000000000;
	LEN = 64;
	
	nblock_test = (int) ceil (((double) nsub_test) / LEN);
	
	bitloc_test.reserve(nsub_test);
	MASK_offset_test.reserve(nsub_test);
	for(int k = 0; k < nsub_test; ++k){
		bitloc_test.push_back(k / LEN);
		MASK_offset_test.push_back(MASK >> (k % LEN));
	}
	
	geno64_test = bitmat(nsnp * 3);
	
	#pragma omp parallel num_threads(nthread)
	{
		#pragma omp for
		for(int i = 0; i < nsnp_test; ++i){
			if(!snp_used[i]){
				continue;
			}else{ 
				geno64_test[snp_col[i] * 3] = bitvec(nblock_test, (uint64) 0);
				geno64_test[snp_col[i] * 3 + 1] = bitvec(nblock_test, (uint64) 0);
				geno64_test[snp_col[i] * 3 + 2] = bitvec(nblock_test, (uint64) 0);
			}
			
			int k = -1;
			int t = 0;
			for(int j = 0; j < nblock_bed; ++j){
				uint8 b = map[3 + i * nblock_bed + j];
				for(int l = 0; l < 4; ++l){
					++k;
					if(k < nsub_fam){
						if(k == index_individual_id_fam_inc[t]){
							uint8 h1 = b & PROBE1[l];
							uint8 h2 = b & PROBE2[l];
							
							if(h1 && !h2){
								cout << "Error: bamboo does not allow missing genotypes" << endl;
								munmap(map, file_size);
								close(fb);
								exit(1);
							}else if(!h1 && !h2){
								geno64_test[snp_col[i] * 3][bitloc_test[t]] |= MASK_offset_test[t];
							}else if(!h1 && h2){
								geno64_test[snp_col[i] * 3 + 1][bitloc_test[t]] |= MASK_offset_test[t];
							}else{
								geno64_test[snp_col[i] * 3 + 2][bitloc_test[t]] |= MASK_offset_test[t];
							}
							
							++t;
							if(t >= nsub_test){
								break;
							}
						}else{
							continue;
						}
					}else{
						break;
					}
				}
			}
		}
	}
	
	munmap(map, file_size);
	close(fb);
	
	cout << snp_used_in_forest.size() << " markers of " << nsub_test << " individuals are loaded from [ " << path_bed_test << " ] for prediction" << endl;
	
	//flip genotypes with MAF > .5
	nflip_test = 0;
	if(flip){//test if the genotypes have been correctly loaded
		vector<bool> flipped (nsnp_test, false);
		#pragma omp parallel num_threads(nthread)
		{
			#pragma omp for
			for(int i = 0; i < nsnp_test; ++i){
				if(!snp_used[i]){
					continue;
				}
				int n0 = 0, n1 = 0, n2 = 0;
				for(int k = 0; k < nblock_test; ++k){
					n0 += PopCount(geno64_test[snp_col[i]*3][k]);
					n2 += PopCount(geno64_test[snp_col[i]*3+2][k]);
				}
				if(n0 < n2){// flip the genotype if necessary
					flipped[i] = true;
					for(int k = 0; k < nblock_test; ++k){
						uint64 u = geno64_test[snp_col[i]*3+2][k];
						geno64_test[snp_col[i]*3+2][k] = geno64_test[snp_col[i]*3][k];
						geno64_test[snp_col[i]*3][k] = u;
					}
				}
			}
		}
		
		for(int i = 0; i < nsnp_test; ++i){
			if(flipped[i]){
				nflip_test++;
			}
		}
	}
	
	if(nflip_test > 0){
		cout << nflip_test << " markers with MAF > 0.5 are flipped" << endl;
	};
	
	//load continuous covariates
	if(cont_var_used_in_forest.size() > 0){
		ifstream f_cont(path_cont_test);
		if(!f_cont){
			cout << "Error: Cannot open CON file " << path_cont_test << endl;
			exit(1);
		}
		
		vector<bool> index_inc(nsub_cont, false);
		for(int i = 0; i < index_individual_id_cont_inc.size(); ++i){
			index_inc[index_individual_id_cont_inc[i]] = true;
		}
		
		vector<vector<double> > xcont_all_test (ncont, vector<double> (nsub_cont, -999999.0));
		
		int nrow = -2;
		for(string s; getline(f_cont, s); ){
			istringstream sin(s);
			string iid;
			++nrow;
			if(nrow == -1){//skip header
				continue;
			}else if(nrow >= 0){
				if(index_inc[nrow]){//use this row
					sin >> iid;
					double d;
					int j = -1;
					while(sin >> d){
						++j;
						if(cont_col[j] == -1){
							continue;
						}else if(cont_col[j] >= 0 && cont_col[j] < ncont){
							xcont_all_test[cont_col[j]][nrow] = d;
						}else{
							cout << "Error: in loading " << path_cont_test << endl;
							f_cont.close();
							exit(1);
						}
					}
				}else{
					continue;
				}
			}else{//impossible
				cout << "Error: what?" << endl;
			}
		}
		f_cont.close();
		
		xcont_test = vector<vector<double> > (ncont, vector<double> (nsub_test, .0));
		for(int j = 0; j < ncont; ++j){
			for(int i = 0; i < nsub_test; ++i){
				xcont_test[j][i] = xcont_all_test[j][index_individual_id_cont_inc[i]];
			}
		}
		
		cout << ncont << " continuous covariate(s) are loaded from [ " << path_cont_test << " ] for prediction" << endl;
	}
	
	//load categorical covariates
	if(cate_var_used_in_forest.size() > 0){
		ifstream f_cate(path_cate_test);
		if(!f_cate){
			cout << "Error: Cannot open CAT file " << path_cate_test << endl;
			exit(1);
		}
		
		vector<bool> index_inc(nsub_cate, false);
		for(int i = 0; i < index_individual_id_cate_inc.size(); ++i){
			index_inc[index_individual_id_cate_inc[i]] = true;
		}
		
		vector<vector<string> > xcate_all_test (ncate, vector<string> (nsub_cate));
		vector<vector<int> > xcate_int_all_test (ncate, vector<int> (nsub_cate));
		
		int nrow = -2;
		for(string s; getline(f_cate, s); ){
			istringstream sin(s);
			string iid;
			++nrow;
			if(nrow == -1){//skip header
				continue;
			}else if(nrow >= 0){
				if(index_inc[nrow]){//use this row
					sin >> iid;
					string str;
					int j = -1;
					while(sin >> str){
						++j;
						int ccj = cate_col[j];
						if(ccj == -1){////this variable is not used in trained bamboo or not in the training data
							continue;
						}else if(ccj >= 0 && ccj < ncate){//this variable is useful in prediction
							xcate_all_test[ccj][nrow] = str;
							
							for(int ii = 0; ii < cate_unique[ccj].size(); ++ii){
								if(str == cate_unique[ccj][ii]){
									xcate_int_all_test[ccj][nrow] = cate_code[ccj][ii];
									break;
								}
							}
						}else{
							cout << "Error: in loading " << path_cate_test << endl;
							f_cate.close();
							exit(1);
						}
					}
				}
			}
		}
		
		xcate_test = vector<vector<string> > (ncate, vector<string> (nsub_test));
		xcate_int_test = vector<vector<int> > (ncate, vector<int> (nsub_test));
		
		for(int j = 0; j < ncate; ++j){
			for(int i = 0; i < nsub_test; ++i){
				xcate_test[j][i] = xcate_all_test[j][index_individual_id_cate_inc[i]];
				xcate_int_test[j][i] = xcate_int_all_test[j][index_individual_id_cate_inc[i]];
			}
		}
		f_cate.close();
		
		cout << ncate << " categorical covariate(s) are loaded from [ " << path_cate_test << " ] for prediction" << endl;
		
	}
	
	time_t current_time;
	time(&current_time);
	cout << "Data loaded: " << ctime(&current_time);
	
	time_t end = time(NULL);
	cout << "Elapsed time: " << (end - start) / 3600 << "h " << ((end - start) % 3600) / 60 << "m " << (end - start) % 60 << "s" << endl;
	
}


void BAMBOO::PutDownTestingSampleToTree(const int sample_id, const int tree_id, 
	int &fall_into_node_id){
	
	if(sample_id < 0 || sample_id >= nsub_test){
		cout << "Error: Invalid sample ID in PutDownTestingSampleToTree" << endl;
		exit(1);
	}
	
	vector<bool> track_path_to_left;//true if to the left
	
	const MINI_NODE *mini_node = &(bamboo[tree_id][0]);
	bool direction = true;//true if to the left, false if to the right
	while(!(*mini_node).terminal){
		if((*mini_node).split_by == CODE_SNP){
			int snp_id = (*mini_node).split_snp_id;
			int snp_left_val = (*mini_node).split_snp_left_val;
			if(snp_left_val == 0){//snp == 0 are going to the left
				direction = 
				(geno64_test[snp_id * 3][bitloc_test[sample_id]] & MASK_offset_test[sample_id]) ? true : false;//true if snp == 0
			}else if(snp_left_val == 1){//snp == 0 or are going to the left
				direction = (geno64_test[snp_id * 3 + 2][bitloc_test[sample_id]] & MASK_offset_test[sample_id]) ? false : true;//false if snp == 2
			}else{//impossible
				cout << "Error: Invalide snp_left_val when predicting samples" << endl;
				exit(1);
			}
		}else if((*mini_node).split_by == CODE_CONT){
			int cont_id = (*mini_node).split_cont_id;
			direction = (xcont_test[cont_id][sample_id] <= (*mini_node).split_cont_thr) ? true : false;//true if <= threshold
		}else if((*mini_node).split_by == CODE_CATE){
			int cate_id = (*mini_node).split_cate_id;
			int size_left = (*mini_node).split_cate_left_code.size();
			int size_right = (*mini_node).split_cate_right_code.size();
			int c = xcate_int_test[cate_id][sample_id];
			if(size_left < size_right){
				direction = false;
				for(int ii = 0; ii < size_left; ++ii){
					if(c == (*mini_node).split_cate_left_code[ii]){
						direction = true;
						break;
					}
				}
			}else{
				direction = true;
				for(int ii = 0; ii < (*mini_node).split_cate_right_code.size(); ++ii){
					if(c == (*mini_node).split_cate_right_code[ii]){
						direction = false;
						break;
					}
				}
			}
		}else{//impossible
			cout << "Error: Invalid splitting code in mini_node.split_by" << endl;
			exit(1);
		}
		
		if(direction){
			mini_node = &(bamboo[tree_id][(*mini_node).ichild1]);
			track_path_to_left.push_back(true);
		}else{
			mini_node = &(bamboo[tree_id][(*mini_node).ichild2]);
			track_path_to_left.push_back(false);
		}
	}
	
	fall_into_node_id = (*mini_node).node_id;

}


void BAMBOO::PutDownTestingSampleToTree(const int sample_id, const vector<MINI_NODE>& rb, 
	int &fall_into_node_id){
	
	if(sample_id < 0 || sample_id >= nsub_test){
		cout << "Error: Invalid sample ID in PutDownTestingSampleToTree" << endl;
		exit(1);
	}
	
	vector<bool> track_path_to_left;//true if to the left
	
	const MINI_NODE *mini_node = &(rb[0]);
	bool direction = true;//true if to the left, false if to the right
	while(!(*mini_node).terminal){
		if((*mini_node).split_by == CODE_SNP){
			int snp_id = (*mini_node).split_snp_id;
			int snp_left_val = (*mini_node).split_snp_left_val;
			if(snp_left_val == 0){//snp == 0 are going to the left
				direction = 
				(geno64_test[snp_id * 3][bitloc_test[sample_id]] & MASK_offset_test[sample_id]) ? true : false;//true if snp == 0
			}else if(snp_left_val == 1){//snp == 0 or are going to the left
				direction = (geno64_test[snp_id * 3 + 2][bitloc_test[sample_id]] & MASK_offset_test[sample_id]) ? false : true;//false if snp == 2
			}else{//impossible
				cout << "Error: Invalide snp_left_val when predicting samples" << endl;
				exit(1);
			}
		}else if((*mini_node).split_by == CODE_CONT){
			int cont_id = (*mini_node).split_cont_id;
			direction = (xcont_test[cont_id][sample_id] <= (*mini_node).split_cont_thr) ? true : false;//true if <= threshold
		}else if((*mini_node).split_by == CODE_CATE){
			int cate_id = (*mini_node).split_cate_id;
			int size_left = (*mini_node).split_cate_left_code.size();
			int size_right = (*mini_node).split_cate_right_code.size();
			int c = xcate_int_test[cate_id][sample_id];
			if(size_left < size_right){
				direction = false;
				for(int ii = 0; ii < size_left; ++ii){
					if(c == (*mini_node).split_cate_left_code[ii]){
						direction = true;
						break;
					}
				}
			}else{
				direction = true;
				for(int ii = 0; ii < (*mini_node).split_cate_right_code.size(); ++ii){
					if(c == (*mini_node).split_cate_right_code[ii]){
						direction = false;
						break;
					}
				}
			}
		}else{//impossible
			cout << "Error: Invalid splitting code in mini_node.split_by" << endl;
			exit(1);
		}
		
		if(direction){
			mini_node = &(rb[(*mini_node).ichild1]);
			track_path_to_left.push_back(true);
		}else{
			mini_node = &(rb[(*mini_node).ichild2]);
			track_path_to_left.push_back(false);
		}
	}
	
	fall_into_node_id = (*mini_node).node_id;

}


void BAMBOO::SavePrediction(){
	
	ofstream file_pred;
	file_pred.open(path_pred);
	file_pred << "IND_ID\tPOST_PROB_CASE1\tPRED_CLASS1\tPOST_PROB_CASE2\tPRED_CLASS2" << endl;
	for(int i = 0; i < nsub_test; ++i){
		
		file_pred << individual_id_test[i] << "\t" << post_prob_case1_test[i] << "\t";
		if(post_prob_case1_test[i] > .5){//?? is 0.5 a good cutoff?
			file_pred << "CASE";
		}else{
			file_pred << "CONTROL";
		}
		
		file_pred << "\t" << post_prob_case2_test[i] << "\t";
		
		if(post_prob_case2_test[i] > .5){//?? is 0.5 a good cutoff?
			file_pred << "CASE";
		}else{
			file_pred << "CONTROL";
		}
		file_pred << endl;
	}
	file_pred.close();
	
	cout << "The predictions of testing data are saved in [ " << path_pred << " ]" << endl;
	
}

void BAMBOO::ScanMultipleForest(){
	
	int nforest = bam_list.size();
	
	if(nforest == 0){
		cout << "Error: Cannot find trained models for prediction" << endl;
		exit(1);
	}else if(nforest == 1){
		cout << "Error: debug in ScanMultipleForest" << endl;
		exit(1);
	}else{//good
		;
	}
	
	vector<string> ver(nforest);
	vector<int> ns(nforest);
	vector<int> nco(nforest);
	vector<int> nca(nforest);
	vector<bool> fl(nforest);
	
	//scan all models to see if we can combine them into a single model
	for(int i = 0; i < nforest; ++i){
		ifstream ifs(bam_list[i]);
		if(!ifs){
			cout << "Error: Cannot open BAM file " << bam_list[i] << endl;
			exit(1);
		}
		
		boost::archive::binary_iarchive ar(ifs);
		
		string ver0;
		int ns0, nco0, nca0;
		bool fl0;
		ar & ver0;
		ar & ns0 & nco0 & nca0;
		ar & fl0;
		
		ver[i] = ver0;
		ns[i] = ns0;
		nco[i] = nco0;
		nca[i] = nca0;
		fl[i] = fl0;
		
	}
	
	for(int i = 1; i < nforest; ++i){
		if(ver[i] != ver[0]){
			cout << "Error: The specified models are built with Random Bamboo with different versions" << endl;
			exit(1);
		}
		
		if(ns[i] != ns[0] || nco[i] != nco[0] || nca[i] != nca[i]){
			cout << "Error: The specified models are built with different training data" << endl;
			exit(1);
		}
		
		if(fl[i] != fl[0]){
			cout << "Error: Some of the specified models flip the genotypes but others do not" << endl;
			exit(1);
		}
	}
	
	nsnp = ns[0];
	ncont = nco[0];
	ncate = nca[0];
	flip = fl[0];
	
	ComputeWordBits();//don't forget this otherwise the PopCount function fails
	
	LEN = 64;
	int nblock_snp = (int) ceil (((double) nsnp) / LEN);
	int nblock_cont = (int) ceil (((double) ncont) / LEN);
	int nblock_cate = (int) ceil (((double) ncate) / LEN);
	
	vector<int> bitloc_snp;
	vector<int> bitloc_cont;
	vector<int> bitloc_cate;
	
	bitvec MASK_offset_snp;
	bitvec MASK_offset_cont;
	bitvec MASK_offset_cate;
	
	bitloc_snp.reserve(nsnp);
	bitloc_cont.reserve(ncont);
	bitloc_cate.reserve(ncate);
	
	MASK_offset_snp.reserve(nsnp);
	MASK_offset_cont.reserve(ncont);
	MASK_offset_cate.reserve(ncate);
	
	MASK = 0x8000000000000000;
	
	for(int k = 0; k < nsnp; ++k){
		bitloc_snp.push_back(k / LEN);
		MASK_offset_snp.push_back(MASK >> (k % LEN));
	}
	
	for(int k = 0; k < ncont; ++k){
		bitloc_cont.push_back(k / LEN);
		MASK_offset_cont.push_back(MASK >> (k % LEN));
	}
	
	for(int k = 0; k < ncate; ++k){
		bitloc_cate.push_back(k / LEN);
		MASK_offset_cate.push_back(MASK >> (k % LEN));
	}
	
	bitmat snp_id_used_in_multiple_models64(nthread, bitvec(nblock_snp, (uint64) 0));
	bitmat cont_var_id_used_in_multiple_models64(nthread, bitvec(nblock_cont, (uint64) 0));
	bitmat cate_var_id_used_in_multiple_models64(nthread, bitvec(nblock_cate, (uint64) 0));
	
	vector<vector<string> > snp_used_in_multiple_model(nthread, vector<string> (nsnp, string("")));
	vector<vector<string> > cont_var_used_in_multiple_model(nthread, vector<string> (ncont, string("")));
	vector<vector<string> > cate_var_used_in_multiple_model(nthread, vector<string> (ncate, string("")));
	
	vector<int> nt (nthread, 0);
	
	#pragma omp parallel num_threads(nthread)
	{
		#pragma omp for
		for(int i = 0; i < nforest; ++i){
			
			int tid = omp_get_thread_num();
			
			ifstream ifs(bam_list[i]);
			if(!ifs){
				cout << "Error: Cannot open BAM file " << bam_list[i] << endl;
				exit(1);
			}
			
			boost::archive::binary_iarchive ar(ifs);
			
			string version0;
			int nsnp0, ncont0, ncate0;
			bool flip0;
			vector<string> snp_used_in_forest0, cont_var_used_in_forest0, cate_var_used_in_forest0;
			vector<int> snp_id_used_in_forest0, cont_var_id_used_in_forest0, cate_var_id_used_in_forest0;
			vector<vector<string> > cate_unique0;
			vector<vector<int> > cate_code0;
			vector<vector<MINI_NODE> > bamboo0;
				
			ar & version0;
			ar & nsnp0 & ncont0 & ncate0;
			ar & flip0;
			ar & snp_used_in_forest0 & snp_id_used_in_forest0;
			ar & cont_var_used_in_forest0 & cont_var_id_used_in_forest0;
			ar & cate_var_used_in_forest0 & cate_var_id_used_in_forest0;
			ar & cate_unique0 & cate_code0;
			ar & bamboo0;
			
			if(i == 0){
				cate_unique = cate_unique0;
				cate_code = cate_code0;
			}
			
			nt[tid] += bamboo0.size();
			//swap is very slow but can save memory (to be confirmed)
//				vector<vector<MINI_NODE> >().swap(bamboo0);
//				vector<string>().swap(snp_used_in_forest0);
//				vector<string>().swap(cont_var_used_in_forest0);
//				vector<string>().swap(cate_var_used_in_forest0);
			
			for(int k = 0; k < snp_id_used_in_forest0.size(); ++k){
				int sid = snp_id_used_in_forest0[k];
				snp_id_used_in_multiple_models64[tid][bitloc_snp[sid]] |= MASK_offset_snp[sid];
				snp_used_in_multiple_model[tid][sid] = snp_used_in_forest0[k];
			}
			
			for(int k = 0; k < cont_var_id_used_in_forest0.size(); ++k){
				int cid = cont_var_id_used_in_forest0[k];
				cont_var_id_used_in_multiple_models64[tid][bitloc_cont[cid]] |= MASK_offset_cont[cid];
				cont_var_used_in_multiple_model[tid][cid] = cont_var_used_in_forest0[k];
			}
			
			for(int k = 0; k < cate_var_id_used_in_forest0.size(); ++k){
				int cid = cate_var_id_used_in_forest0[k];
				cate_var_id_used_in_multiple_models64[tid][bitloc_cate[cid]] |= MASK_offset_cate[cid];
				cate_var_used_in_multiple_model[tid][cid] = cate_var_used_in_forest0[k];
			}
			
		}
		
	}
	
	
	ntree = 0;
	for(int j = 0; j < nthread; ++j){
		ntree += nt[j];
	}
	
	bitvec vec;
	vec = bitvec (nblock_snp, (uint64) 0);
	int nsnp_used = 0;
	for(int i = 0; i < nblock_snp; ++i){
		uint64 u = (uint64) 0;
		for(int j = 0; j < nthread; ++j){
			u |= snp_id_used_in_multiple_models64[j][i];
		}
		vec[i] = u;
		nsnp_used += PopCount(u);
	}
	
	if(!snp_used_in_forest.empty()){
		snp_used_in_forest.clear();
	}
	
	snp_used_in_forest.reserve(nsnp_used);
	
	if(!snp_id_used_in_forest.empty()){
		snp_id_used_in_forest.clear();
	}
	
	snp_id_used_in_forest.reserve(nsnp_used);
	
	for(int i = 0; i < nsnp; ++i){
		if(vec[bitloc_snp[i]] & MASK_offset_snp[i]){
			snp_id_used_in_forest.push_back(i);
			
			for(int j = 0; j < nthread; ++j){
				if(snp_used_in_multiple_model[j][i] != string("")){
					snp_used_in_forest.push_back(snp_used_in_multiple_model[j][i]);
					break;
				}
			}
		}
	}
	
	vec = bitvec (nblock_cont, (uint64) 0);
	int ncont_used = 0;
	for(int i = 0; i < nblock_cont; ++i){
		uint64 u = (uint64) 0;
		for(int j = 0; j < nthread; ++j){
			u |= cont_var_id_used_in_multiple_models64[j][i];
		}
		vec[i] = u;
		ncont_used += PopCount(u);
	}
	
	
	if(!cont_var_used_in_forest.empty()){
		cont_var_used_in_forest.clear();
	}
	
	cont_var_used_in_forest.reserve(ncont_used);
	
	if(!cont_var_id_used_in_forest.empty()){
		cont_var_id_used_in_forest.clear();
	}
	
	cont_var_id_used_in_forest.reserve(ncont_used);
	
	for(int i = 0; i < ncont; ++i){
		if(vec[bitloc_cont[i]] & MASK_offset_cont[i]){
			cont_var_id_used_in_forest.push_back(i);
			
			for(int j = 0; j < nthread; ++j){
				if(cont_var_used_in_multiple_model[j][i] != string("")){
					cont_var_used_in_forest.push_back(cont_var_used_in_multiple_model[j][i]);
					break;
				}
			}
		}
	}
	
	vec = bitvec (nblock_cate, (uint64) 0);
	int ncate_used = 0;
	for(int i = 0; i < nblock_cate; ++i){
		uint64 u = (uint64) 0;
		for(int j = 0; j < nthread; ++j){
			u |= cate_var_id_used_in_multiple_models64[j][i];
		}
		vec[i] = u;
		ncate_used += PopCount(u);
	}
	
	
	if(!cate_var_used_in_forest.empty()){
		cate_var_used_in_forest.clear();
	}
	
	cate_var_used_in_forest.reserve(ncate_used);
	
	if(!cate_var_id_used_in_forest.empty()){
		cate_var_id_used_in_forest.clear();
	}
	
	cate_var_id_used_in_forest.reserve(ncate_used);
	
	for(int i = 0; i < ncate; ++i){
		if(vec[bitloc_cate[i]] & MASK_offset_cate[i]){
			cate_var_id_used_in_forest.push_back(i);
			
			for(int j = 0; j < nthread; ++j){
				if(cate_var_used_in_multiple_model[j][i] != string("")){
					cate_var_used_in_forest.push_back(cate_var_used_in_multiple_model[j][i]);
					break;
				}
			}
		}
	}
	
	cout << "The " << bam_list.size() << " models include: " << endl;
	cout << ntree << " bamboos" << endl;
	cout << nsnp_used << " markers" << endl;
	if(ncont_used > 0){
		cout << ncont_used << " continuous covariate(s)" << endl;
	}
	if(ncate_used > 0){
		cout << ncate_used << " categorical covariate(s)" << endl;
	}
	if(flip){
		cout << "Markers with MAF > 0.5 are flipped" << endl;
	}
	
	has_cont = cont_var_used_in_forest.size() ? true : false;
	has_cate = cate_var_used_in_forest.size() ? true : false;
	
}

//this function is called only if --pred and --bam are specified, --bam is set as a directory containing more than one .bam files, 
//and --file is not specified. In this case, path_bam is set as NULL in ParseModelPath() and the names of .bam files are stored in bam_list
void BAMBOO::PredictTestingSampleFromMultipleForest(){
	
	if(!pred_from_specified_model && !pred_from_trained_model){
		return;
	}
	
	time_t start = time(NULL);
	
	ScanMultipleForest();//find those markers used in the whole model
	LoadTestingData();
	
	vector<vector<double> > post_prob_case1_test_omp (nthread, vector<double> (nsub_test, .0));
	vector<vector<double> > tmp0_omp (nthread, vector<double> (nsub_test, .0));
	vector<vector<double> > tmp1_omp (nthread, vector<double> (nsub_test, .0));
	
	int nforest = bam_list.size();
	
	#pragma omp parallel num_threads(nthread)
	{
		
		#pragma omp for
		for(int i = 0; i < nforest; ++i){
			
			int tid = omp_get_thread_num();
			ifstream ifs(bam_list[i]);
			if(!ifs){
				cout << "Error: Cannot open BAM file " << bam_list[i] << endl;
				exit(1);
			}
			boost::archive::binary_iarchive ar(ifs);
			
			string version0;
			int nsnp0, ncont0, ncate0;
			bool flip0;
			vector<string> snp_used_in_forest0, cont_var_used_in_forest0, cate_var_used_in_forest0;
			vector<int> snp_id_used_in_forest0, cont_var_id_used_in_forest0, cate_var_id_used_in_forest0;
			vector<vector<string> > cate_unique0;
			vector<vector<int> > cate_code0;
			vector<vector<MINI_NODE> > bamboo0;
			
			ar & version0;
			ar & nsnp0 & ncont0 & ncate0;
			ar & flip0;
			ar & snp_used_in_forest0 & snp_id_used_in_forest0;
			ar & cont_var_used_in_forest0 & cont_var_id_used_in_forest0;
			ar & cate_var_used_in_forest0 & cate_var_id_used_in_forest0;
			ar & cate_unique0 & cate_code0;
			ar & bamboo0;
			
			int nt = bamboo0.size();
			
			for(int sample_id = 0; sample_id < nsub_test; ++sample_id){
				
				for(int tree_id = 0; tree_id < nt; ++tree_id){
					int fall_into_node_id = -1;
					PutDownTestingSampleToTree(sample_id, bamboo0[tree_id], fall_into_node_id);
					if(fall_into_node_id >= 0 && fall_into_node_id < bamboo0[tree_id].size()){
						double r0 = bamboo0[tree_id][fall_into_node_id].node_risk_ctrl;
						double r1 = bamboo0[tree_id][fall_into_node_id].node_risk_case;
						post_prob_case1_test_omp[tid][sample_id] += r0 / (r0 + r1);
						tmp0_omp[tid][sample_id] += r0;
						tmp1_omp[tid][sample_id] += r1;
					}else{
						cout << "Error: Cannot assign sample " << sample_id + 1 << " in testing data to any leaf of the tree " << tree_id + 1 << " in " << bam_list[i] << endl;
						exit(1);
					}
				}
			}
			
		}
		
	}
	
	post_prob_case1_test = vector<double> (nsub_test, .0);
	post_prob_case2_test = vector<double> (nsub_test, .0);
	vector<double> tmp0 (nsub_test, .0);
	vector<double> tmp1 (nsub_test, .0);
	
	#pragma omp parallel num_threads(nthread)
	{
		#pragma omp for
		for(int i = 0; i < nsub_test; ++i){
			for(int j = 0; j < nthread; ++j){
				post_prob_case1_test[i] += post_prob_case1_test_omp[j][i];
				tmp0[i] += tmp0_omp[j][i];
				tmp1[i] += tmp1_omp[j][i];
			}
			
			post_prob_case1_test[i] /= ntree;
			post_prob_case2_test[i] = tmp0[i] / (tmp0[i] + tmp1[i]);
		}
	}
	
	SavePrediction();
	
	time_t complete_time;
	time(&complete_time);
	cout << "Prediction accomplished: " << ctime(&complete_time);
	
	
	time_t end = time(NULL);
	cout << "Elapsed time: " << (end - start) / 3600 << "h " << ((end - start) % 3600) / 60 << "m " << (end - start) % 60 << "s" << endl;
	
}


void BAMBOO::PredictTestingSampleFromSingleForest(){
	
	if(!pred_from_specified_model && !pred_from_trained_model){
		return;
	}
	
	time_t start = time(NULL);
	
	LoadBamboo();
	LoadTestingData();
	cout << "model engaged" << endl;
	
	if(ntree == -1){
		ntree = bamboo.size();
	}else{
		if(ntree > bamboo.size()){
			ntree = bamboo.size();
			cout << "ntree is reset as " << ntree << endl;
		}else if(ntree < bamboo.size()){
			cout << "Only " << ntree << " bamboos are used in prediction" << endl;
		}else{
			;
		}
	}
	
	post_prob_case1_test = vector<double> (nsub_test, .0);
	post_prob_case2_test = vector<double> (nsub_test, .0);
	
	#pragma omp parallel num_threads(nthread)
	{
		#pragma omp for
		for(int sample_id = 0; sample_id < nsub_test; ++sample_id){
			
			double tmp0 = .0;
			double tmp1 = .0;
			for(int tree_id = 0; tree_id < ntree; ++tree_id){
				
				int fall_into_node_id = -1;
				PutDownTestingSampleToTree(sample_id, tree_id, fall_into_node_id);
				
				if(fall_into_node_id >= 0 && fall_into_node_id < bamboo[tree_id].size()){
					double r0 = bamboo[tree_id][fall_into_node_id].node_risk_ctrl;
					double r1 = bamboo[tree_id][fall_into_node_id].node_risk_case;
					post_prob_case1_test[sample_id] += r0 / (r0 + r1);
					tmp0 += r0;
					tmp1 += r1;
				}else{
					cout << "Error: Cannot assign sample " << sample_id + 1 << " in testing data to any leaf of the tree " << tree_id + 1 << endl;
					exit(1);
				}
				
			}
			
			post_prob_case1_test[sample_id] /= ntree;
			post_prob_case2_test[sample_id] = tmp0 / (tmp0 + tmp1);
			
		}
	}
	
	SavePrediction();
	
	time_t complete_time;
	time(&complete_time);
	cout << "Prediction accomplished: " << ctime(&complete_time);
	
	
	time_t end = time(NULL);
	cout << "Elapsed time: " << (end - start) / 3600 << "h " << ((end - start) % 3600) / 60 << "m " << (end - start) % 60 << "s" << endl;
	
}

void BAMBOO::PredictTestingSample(){
	
	if(bam_list.size() > 1){
		PredictTestingSampleFromMultipleForest();
	}else{
		PredictTestingSampleFromSingleForest();
	}
	
}



