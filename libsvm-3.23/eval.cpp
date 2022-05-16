#include <stdio.h>
#include <ctype.h>
#include <vector>
#include <algorithm>
#include <errno.h>
#include <string.h>
#include "svm.h"
#include "eval.h"

#define Malloc(type,n) (type *)malloc((n)*sizeof(type))

// Prototypes of evaluation functions
double precision(const size_t, const double *dec_values, const int *ty);
double recall(const size_t, const double *dec_values, const int *ty);
double fscore(const size_t, const double *dec_values, const int *ty);
double bac(const size_t, const double *dec_values, const int *ty);
double auc(const size_t, const double *dec_values, const int *ty);
double accuracy(const size_t, const double *dec_values, const int *ty);

// Evaluation function pointer
// You can assign this pointer to any above prototype
double (*evaluation_function)(const size_t, const double *, const int *) = accuracy;

bool check_binary_beforeCV(const struct svm_problem *prob, int *ty);

// Evaluation functions below are the same for both LIBSVM and LIBLINEAR

double precision(const size_t size, const double *dec_values, const int *ty)
{
	size_t i;
	int    tp, fp;
	double precision;

	tp = fp = 0;

	for(i = 0; i < size; ++i) if(dec_values[i] >= 0){
		if(ty[i] == 1) ++tp;
		else           ++fp;
	}

	if(tp + fp == 0){
		fprintf(stderr, "warning: No postive predict label.\n");
		precision = 0;
	}else
		precision = tp / (double) (tp + fp);
	printf("Precision = %g%% (%d/%d)\n", 100.0 * precision, tp, tp + fp);
	
	return precision;
}

double recall(const size_t size, const double *dec_values, const int *ty)
{
	size_t i;
	int    tp, fn; // true_positive and false_negative
	double recall;

	tp = fn = 0;

	for(i = 0; i < size; ++i) if(ty[i] == 1){ // true label is 1
		if(dec_values[i] >= 0) ++tp; // predict label is 1
		else                   ++fn; // predict label is -1
	}

	if(tp + fn == 0){
		fprintf(stderr, "warning: No postive true label.\n");
		recall = 0;
	}else
		recall = tp / (double) (tp + fn);
	// print result in case of invocation in prediction
	printf("Recall = %g%% (%d/%d)\n", 100.0 * recall, tp, tp + fn);
	
	return recall; // return the evaluation value
}

double fscore(const size_t size, const double *dec_values, const int *ty)
{
	size_t i;
	int    tp, fp, fn;
	double precision, recall;
	double fscore;

	tp = fp = fn = 0;

	for(i = 0; i < size; ++i) 
		if(dec_values[i] >= 0 && ty[i] == 1) ++tp;
		else if(dec_values[i] >= 0 && ty[i] == -1) ++fp;
		else if(dec_values[i] <  0 && ty[i] == 1) ++fn;

	if(tp + fp == 0){
		fprintf(stderr, "warning: No postive predict label.\n");
		precision = 0;
	}else
		precision = tp / (double) (tp + fp);
	if(tp + fn == 0){
		fprintf(stderr, "warning: No postive true label.\n");
		recall = 0;
	}else
		recall = tp / (double) (tp + fn);

	
	if(precision + recall == 0){
		fprintf(stderr, "warning: precision + recall = 0.\n");
		fscore = 0;
	}else
		fscore = 2 * precision * recall / (precision + recall);

	printf("F-score = %g\n", fscore);
	
	return fscore;
}

double bac(const size_t size, const double *dec_values, const int *ty)
{
	size_t i;
	int    tp, fp, fn, tn;
	double specificity, recall;
	double bac;

	tp = fp = fn = tn = 0;

	for(i = 0; i < size; ++i) 
		if(dec_values[i] >= 0 && ty[i] == 1) ++tp;
		else if(dec_values[i] >= 0 && ty[i] == -1) ++fp;
		else if(dec_values[i] <  0 && ty[i] == 1)  ++fn;
		else ++tn;

	if(tn + fp == 0){
		fprintf(stderr, "warning: No negative true label.\n");
		specificity = 0;
	}else
		specificity = tn / (double)(tn + fp);
	if(tp + fn == 0){
		fprintf(stderr, "warning: No positive true label.\n");
		recall = 0;
	}else
		recall = tp / (double)(tp + fn);

	bac = (specificity + recall) / 2;
	printf("BAC = %g\n", bac);
	
	return bac;
}

// only for auc
class Comp{
	const double *dec_val;
	public:
	Comp(const double *ptr): dec_val(ptr){}
	bool operator()(int i, int j) const{
		return dec_val[i] > dec_val[j];
	}
};

double auc(const size_t size, const double *dec_values, const int *ty)
{
	double roc  = 0;
	size_t i;
	std::vector<size_t> indices(size);

	for(i = 0; i < size; ++i) indices[i] = i;

	std::sort(indices.begin(), indices.end(), Comp(&dec_values[0]));

	int tp = 0,fp = 0;
	for(i = 0; i < size; i++) {
		if(ty[indices[i]] == 1) tp++;
		else if(ty[indices[i]] == -1) {
			roc += tp;
			fp++;
		}
	}

	if(tp == 0 || fp == 0)
	{
		fprintf(stderr, "warning: Too few postive true labels or negative true labels\n");
		roc = 0;
	}
	else
		roc = roc / tp / fp;

	printf("AUC = %g\n", roc);

	return roc;
}

double accuracy(const size_t total, const double *dec_values, const int *ty)
{
	int    correct = 0;
	size_t i;

	for(i = 0; i < total; ++i)
		if(ty[i] == (dec_values[i] >= 0? 1: -1)) ++correct;

	printf("Accuracy = %g%% (%d/%d)\n",
	       (double)correct/(int) total*100, correct, (int) total);

	return (double) correct / (int) total;
}

double ap(const size_t size, const double *dec_values, const int *ty)
{
	size_t i;
	std::vector<size_t> indices(size);

	for(i = 0; i < size; ++i) indices[i] = i;
	std::sort(indices.begin(), indices.end(), Comp(&dec_values[0]));

	int p = 0, tp = 0;
	double prev_recall = 0, area = 0;

	for(i = 0; i < size; ++i) p += (ty[i] == 1);

	if(p == 0) {
		fprintf(stderr, "warning: Too few postive labels\n");
		return 0;
	}

	for(i = 0; i < size; ++i) {
		tp += (ty[indices[i]] == 1);
		
		if(i+1 < size && dec_values[indices[i]] == dec_values[indices[i+1]]) 
			continue;

		double recall = (double)tp/p;
		double precision = (double)tp/(double)(i+1);

		area += precision*(recall-prev_recall);
		prev_recall = recall;
	}

	printf("AP = %g\n", area);
	return area;
}

bool check_binary_beforeCV(const svm_problem *prob, int *ty)
{
	int label[3];
	int nr_class = 0;
	bool is_binary = false;

	for (int i=0; i<prob->l; i++)
	{
		int this_label = (int)prob->y[i];
		int j;
		for(j=0;j<nr_class;j++)
			if(this_label == label[j])
				break;
		if(j == nr_class)
		{
			if (nr_class > 2)
				break;
			label[nr_class] = this_label;
			++nr_class;
		}
		if (this_label == 1)
			ty[i] = 1;
		else if (this_label == -1)
			ty[i] = -1;
	}
	if (nr_class == 2)
	{
		if ((label[0] == 1 && label[1] == -1) ||
		    (label[0] == -1 && label[1] == 1))
			is_binary = true;
		else
		{
			fprintf(stderr,"ERROR: to use other evaluation criteria, labels must be +1/-1\n");
			exit(1);
		}
	}

	return is_binary;
}

double binary_class_cross_validation(const svm_problem *prob, const svm_parameter *param, int nr_fold)
{
	int i;
	int *fold_start = Malloc(int,nr_fold+1);
	int l = prob->l;
	int *perm = Malloc(int,l);
	int *labels;
	double *dec_values = Malloc(double,l);
	int *ty = Malloc(int,l);
	check_binary_beforeCV(prob, ty);
	
	for(i=0;i<l;i++) perm[i]=i;
	for(i=0;i<l;i++)
	{
		int j = i+rand()%(l-i);
		std::swap(perm[i],perm[j]);
	}
	for(i=0;i<=nr_fold;i++)
		fold_start[i]=i*l/nr_fold;

	for(i=0;i<nr_fold;i++)
	{
		int begin = fold_start[i];
		int end = fold_start[i+1];
		int j,k;
		struct svm_problem subprob;

		subprob.l = l-(end-begin);
		subprob.x = Malloc(struct svm_node*,subprob.l);
		subprob.y = Malloc(double,subprob.l);

		k=0;
		for(j=0;j<begin;j++)
		{
			subprob.x[k] = prob->x[perm[j]];
			subprob.y[k] = prob->y[perm[j]];
			++k;
		}
		for(j=end;j<l;j++)
		{
			subprob.x[k] = prob->x[perm[j]];
			subprob.y[k] = prob->y[perm[j]];
			++k;
		}
		struct svm_model *submodel = svm_train(&subprob,param);
		int svm_type = svm_get_svm_type(submodel);
	
		if(svm_type == NU_SVR || svm_type == EPSILON_SVR){
			fprintf(stderr, "wrong svm type");
			exit(1);
		}

		labels = Malloc(int, svm_get_nr_class(submodel));
		svm_get_labels(submodel, labels);

		if(svm_get_nr_class(submodel) > 2) 
		{
			fprintf(stderr,"Error: the number of class is not equal to 2\n");
			exit(-1);
		}

		for(j=begin;j<end;j++)
			svm_predict_values(submodel,prob->x[perm[j]], &dec_values[perm[j]]);

		svm_free_and_destroy_model(&submodel);
		free(subprob.x);
		free(subprob.y);
		free(labels);
	}		
	precision(l, dec_values, ty);
	recall(l, dec_values, ty);
	fscore(l, dec_values, ty);
	bac(l, dec_values, ty);
	accuracy(l, dec_values, ty);
	ap(l, dec_values, ty);

	double current_rate = evaluation_function(l, dec_values, ty);
	
	free(perm);
	free(fold_start);
	free(dec_values);
	free(ty);
	
	return current_rate;
}

static void exit_input_error(int line_num)
{
	fprintf(stderr,"Wrong input format at line %d\n", line_num);
	exit(1);
}

static char *line = NULL;
static int max_line_len;

static char* readline(FILE *input)
{
	int len;

	if(fgets(line,max_line_len,input) == NULL)
		return NULL;

	while(strrchr(line,'\n') == NULL)
	{
		max_line_len *= 2;
		line = (char *) realloc(line,max_line_len);
		len = (int) strlen(line);
		if(fgets(line+len,max_line_len-len,input) == NULL)
			break;
	}
	return line;
}

void binary_class_predict(FILE *input, FILE *output, const svm_model *model)
{
	int total = 0;
	int *labels;
	int max_nr_attr = 64;
	struct svm_node *x = Malloc(struct svm_node, max_nr_attr);
	int max_nr_data = 64;
	double *dec_values = Malloc(double, max_nr_data);
	int *true_labels = Malloc(int, max_nr_data);

	int svm_type=svm_get_svm_type(model);
	
	if (svm_type==NU_SVR || svm_type==EPSILON_SVR){
		fprintf(stderr, "wrong svm type.");
		exit(1);
	}

	int nr_class = svm_get_nr_class(model);
	labels = Malloc(int, nr_class);
	svm_get_labels(model, labels);

	if (nr_class == 2)
		for (int j=0; j<2; j++)
			if (labels[j] != 1 && labels[j] != -1)
			{
				fprintf(stderr,"ERROR: to use other evaluation criteria, labels must be +1/-1\n");
				exit(1);
			}
	
	max_line_len = 1024;
	line = (char *)malloc(max_line_len*sizeof(char));
	while(readline(input) != NULL)
	{
		int i = 0;
		double target_label, predict_label;
		char *idx, *val, *label, *endptr;
		int inst_max_index = -1; // strtol gives 0 if wrong format, and precomputed kernel has <index> start from 0

		label = strtok(line," \t\n");
		if(label == NULL) // empty line
			exit_input_error(total+1);

		target_label = strtod(label,&endptr);
		if(endptr == label || *endptr != '\0')
			exit_input_error(total+1);

		while(1)
		{
			if(i>=max_nr_attr-1)	// need one more for index = -1
			{
				max_nr_attr *= 2;
				x = (struct svm_node *) realloc(x,max_nr_attr*sizeof(struct svm_node));
			}

			idx = strtok(NULL,":");
			val = strtok(NULL," \t");

			if(val == NULL)
				break;
			errno = 0;
			x[i].index = (int) strtol(idx,&endptr,10);
			if(endptr == idx || errno != 0 || *endptr != '\0' || x[i].index <= inst_max_index)
				exit_input_error(total+1);
			else
				inst_max_index = x[i].index;

			errno = 0;
			x[i].value = strtod(val,&endptr);
			if(endptr == val || errno != 0 || (*endptr != '\0' && !isspace(*endptr)))
				exit_input_error(total+1);

			++i;
		}
		x[i].index = -1;

		predict_label = svm_predict(model,x);
		fprintf(output,"%g\n",predict_label);

		double dec_value;
		svm_predict_values(model, x, &dec_value);

		true_labels[total] = (target_label > 0)? 1: -1;
		dec_values[total] = dec_value;

		total++;
		if(total >= max_nr_data)
		{
			max_nr_data *= 2;
			dec_values = (double *) realloc(dec_values, max_nr_data*sizeof(double));
			true_labels = (int *) realloc(true_labels, max_nr_data*sizeof(int));
		}
	}	
	precision(total, dec_values, true_labels);
	recall(total, dec_values, true_labels);
	fscore(total, dec_values, true_labels);
	bac(total, dec_values, true_labels);
	accuracy(total, dec_values, true_labels);
	ap(total, dec_values, true_labels);

	evaluation_function(total, dec_values, true_labels);

	free(labels);
	free(x);
	free(dec_values);
	free(true_labels);
}
