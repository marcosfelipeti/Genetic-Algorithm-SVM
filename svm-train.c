#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>
#include "svm.h"
#include "time.h"
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>




#define Malloc(type,n) (type *)malloc((n)*sizeof(type))

//unsigned int myseed = (unsigned int)time(NULL);
unsigned int myseed = (unsigned int)10;

Individual best_one;
double mean;
double stdev;
int population_size = 10;
int generations_number = 10;


//void print_null(const char *s) {}

void exit_with_help()
{
	printf(
	"Usage: svm-train [options] training_set_file [model_file]\n"
	"options:\n"
	"-s svm_type : set type of SVM (default 0)\n"
	"	0 -- C-SVC		(multi-class classification)\n"
	"	1 -- nu-SVC		(multi-class classification)\n"
	"	2 -- one-class SVM\n"
	"	3 -- epsilon-SVR	(regression)\n"
	"	4 -- nu-SVR		(regression)\n"
	"-t kernel_type : set type of kernel function (default 2)\n"
	"	0 -- linear: u'*v\n"
	"	1 -- polynomial: (gamma*u'*v + coef0)^degree\n"
	"	2 -- radial basis function: exp(-gamma*|u-v|^2)\n"
	"	3 -- sigmoid: tanh(gamma*u'*v + coef0)\n"
	"	4 -- precomputed kernel (kernel values in training_set_file)\n"
	"-d degree : set degree in kernel function (default 3)\n"
	"-g gamma : set gamma in kernel function (default 1/num_features)\n"
	"-r coef0 : set coef0 in kernel function (default 0)\n"
	"-c cost : set the parameter C of C-SVC, epsilon-SVR, and nu-SVR (default 1)\n"
	"-n nu : set the parameter nu of nu-SVC, one-class SVM, and nu-SVR (default 0.5)\n"
	"-p epsilon : set the epsilon in loss function of epsilon-SVR (default 0.1)\n"
	"-m cachesize : set cache memory size in MB (default 100)\n"
	"-e epsilon : set tolerance of termination criterion (default 0.001)\n"
	"-h shrinking : whether to use the shrinking heuristics, 0 or 1 (default 1)\n"
	"-b probability_estimates : whether to train a SVC or SVR model for probability estimates, 0 or 1 (default 0)\n"
	"-wi weight : set the parameter C of class i to weight*C, for C-SVC (default 1)\n"
	"-v n: n-fold cross validation mode\n"
	"-q : quiet mode (no outputs)\n"
	);
	exit(1);
}

void exit_input_error(int line_num)
{
	fprintf(stderr,"Wrong input format at line %d\n", line_num);
	exit(1);
}

void parse_command_line(double v, double c);
void read_problem(const char *filename);
double do_cross_validation();

struct svm_parameter param;		// set by parse_command_line
struct svm_problem prob;		// set by read_problem
struct svm_model *model;
struct svm_node *x_space;
int cross_validation;
int nr_fold;

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

//check parameters error
void check_p_error()
{
	const char *error_msg;
	error_msg = svm_check_parameter(&prob,&param);
	if(error_msg)
	{
		fprintf(stderr,"ERROR: %s\n",error_msg);
		exit(1);
	}
}

//***********************************************************************************************************************************************************************************************


/*
* Genetic Algorithm starts here, honey
*/

//print one individual
void print_individual(Individual ind) {

	int j;

	printf("\nGenome: ");
	for (j=0; j<genome_size; j++) {
		printf("%d", ind.genome[j]);
	}

	printf("\nValue g: %f\n", ind.v_value); 
	printf("Value c: %f\n", ind.c_value); 
	printf("Fitness: %f\n", ind.fitness);

	printf("\n\n\n");
}



//print population
void print_pop(Individual *individuals) {

	int i, j;

	for (i=0; i<population_size; i++) {
		printf("\nIndex: %d\n", i);
		printf("Genome: ");
		for (j=0; j<genome_size; j++) {
			printf("%d", individuals[i].genome[j]);
		}
		printf("\nValue g: %f\n", individuals[i].v_value); 
		printf("VAlue c: %f\n", individuals[i].c_value); 

		printf("Fitness: %f\n", individuals[i].fitness);

		printf("\n\n\n");

	}
}


//creating random population
int random_pop (Individual* individuals)
{
	int i, j;

	for (i = 0; i < population_size; i++){
		for (j = 0; j < genome_size; j++){
			individuals[i].genome[j] = (rand_r(&myseed)% 2);
		}	
	}

	return 0;
}

//binary to decimal
int bin_to_dec(Individual* individuals, int index, int begining, int end) {

	int j, sum=0;
	int exp = 0;

	for(j=end; j >= begining; j--) {
		sum += (int)((pow(2, exp)) * individuals[index].genome[j]);

		exp++;
	}

  return sum;
}


// interpolation
void interpolation(Individual* individuals) {

	int i, b10;

	for(i=0; i < population_size; i++) {		

		b10 = bin_to_dec(individuals, i, 0, TAM_var-1);
		individuals[i].v_value = (double)MINg + ((b10/(pow(2, TAM_var)-1)) * (MAXg - MINg));


		b10 = bin_to_dec(individuals, i, TAM_C_parameter, genome_size-1);
		individuals[i].c_value = (double)MINc + ((b10/(pow(2, TAM_C_parameter)-1)) * (MAXc-MINc));	
		
		//printf("original= %d\n", b10);
		//printf("c value= %f\n",individuals[i].c_value);
		}
}

//How good are you?
void evaluate(Individual *individuals)
{
	int i;

	for(i=1; i < population_size; i++) {
		

		parse_command_line(individuals[i].v_value, individuals[i].c_value);
		check_p_error();
		
		individuals[i].fitness = do_cross_validation();
		
	}
}


//------------------------------------------------------------------quick sort--------------------------------------------------------
void swap(Individual* a, Individual* b) {

	Individual tmp;

	tmp = *a;
	*a = *b;
	*b = tmp;
}
 
int partition(Individual* vec, int left, int right) {

	int i, j;
 
	i = left;
	for (j = left + 1; j <= right; ++j) {
		if (vec[j].fitness < vec[left].fitness) {
      			++i;
			swap(&vec[i], &vec[j]);
    		}
  	}
  	swap(&vec[left], &vec[i]);
 
	return i;
}
 
void quick_sort(Individual* vec, int left, int right) {

	int r;
 
	if (right > left) {
		r = partition(vec, left, right);
		quick_sort(vec, left, r - 1);
		quick_sort(vec, r + 1, right);
  	}
}

// -------------------------------------------------------end of quick sort----------------------------------------------------------



/******************************************************************Roleta****************************************************************/
	double summation(Individual *individuals)
	{
		double sum = 0.0;
		int i;

		for(i = 0; i < population_size; i++)
		{
			sum += individuals[i].fitness;
		}

		return sum;
	}

	int jequiti(Individual *individuals)
	{
		double sum = summation(individuals);
		double r = (rand_r((&myseed))*1000)%((unsigned int)(sum*1000));
		double partial_sum = 0.0;
		int index = -1;

		do
		{
			index += 1;
			partial_sum = partial_sum + individuals[index].fitness*1000;
		} while(partial_sum <= r);

		return index;
	}
	
	
	//calculate the mean of the population
void calculate_mean(Individual *individuals) {
	double sum = summation(individuals);
	mean = sum/population_size;
	printf("\nMean = %f\n", mean);
}


//calculate stdev of the population
void calculate_stdv(Individual *individuals) {
	calculate_mean(individuals);
	double variance = 0;
	int i;
	
	for (i= 0; i < population_size; i++){
			double difference = pow ((individuals[i].fitness - mean), 2);
			variance += difference;
	}
	
	variance = variance/population_size;
	stdev = sqrt(variance);
	
		printf("\nStdev = %f\n\n", stdev);
}
	

//let's have sex, honey!
void do_crossover(Individual *individual1, Individual *individual2)
{
	int i, j;	


	for(j=0; j<population_size; j+=2)
	{
		int dad = jequiti(individual1);
		int mom = jequiti(individual1);	
		
		int cutting = (rand_r(&myseed)% genome_size);
		double check = (rand_r(&myseed)%10);

		if(cutting != 0 && check > crossover_percent)
		{
			for(i=0; i<genome_size; i++) {
				individual2[j].genome[i] = individual1[dad].genome[i];
				individual2[j+1].genome[i] = individual1[mom].genome[i];
			}
			for (i = cutting; i < genome_size; i++) {
				individual2[j].genome[i] = individual1[mom].genome[i];
				individual2[j+1].genome[i] = individual1[dad].genome[i];
			}
		}
	}
}


//let's have sex in a different way, honey!
void do_crossover1(Individual *individual1, Individual *individual2)
{
	int i, j;
	int dad = jequiti(individual1);
	int mom = jequiti(individual1);			


	for(j=0; j<population_size; j+=2)
	{
		double check = (rand_r(&myseed)%10);

		if(check > crossover_percent)
		{

			for(i=0; i<genome_size; i++) {
				if(individual1[dad].genome[i] == individual1[mom].genome[i]){
				
				individual2[j].genome[i] = 1;
				
				}
				else{
					individual2[j].genome[i] = 0;
				}
			}

		}
	}
}

//I love having sex with you, baby!
void do_crossover2(Individual *individual1, Individual *individual2)
{
	int i, j;


	for(j=0; j<population_size; j+=2)
	{
		int dad = jequiti(individual1);
		int mom = jequiti(individual1);	
		
		int locus1 = (rand_r(&myseed)% genome_size);
		int locus2 = (rand_r(&myseed)% genome_size);
		
		double check = (rand_r(&myseed)%1000);

		if((locus1 != locus2) && (locus1 != 0) && (locus2 != 0) && check > crossover_percent)
		{
			if ( locus1 > locus2 ){
				
				individual2[j] = individual1[dad];
				individual2[j+1] = individual1[mom];
				
				for(i = locus1; i <= locus2; i++) {
					individual2[j].genome[i] = individual1[mom].genome[i];
					
				}
				for(i = locus1; i <= locus2; i++) {
					individual2[j+1].genome[i] = individual1[dad].genome[i];
				}
			}
			
			else{
				
				
				individual2[j] = individual1[dad];
				individual2[j+1] = individual1[mom];
				
				for(i = locus2; i <= locus1; i++) {
					individual2[j].genome[i] = individual1[mom].genome[i];
					
				}
				
				for(i = locus2; i <= locus1; i++) {
					individual2[j+1].genome[i] = individual1[dad].genome[i];
				}
			}
		}
			
			
		}
		
}

//Ok, Honey. I don't think we are getting along so well... Let's try something else.
void do_crossover3(Individual *individual1, Individual *individual2)
{
	int i, j;	


	for(j=0; j<population_size; j+=2)
	{
		int dad = jequiti(individual1);
		int mom = jequiti(individual1);	
		
		double check = (rand_r(&myseed)%1000);

		if(check > crossover_percent)
		{
			for(i=0; i<genome_size; i++) {
				
				int cutting = (rand_r(&myseed)% 2);
				
				if(cutting == 0) {
					individual2[j].genome[i] = individual1[dad].genome[i];
				}
				else {
					individual2[j].genome[i] = individual1[mom].genome[i];
				}
				
				
			for (i=0; i<genome_size; i++) {
				
				int cutting = (rand_r(&myseed)% 2);
				
				if(cutting == 0) {
					individual2[j+1].genome[i] = individual1[dad].genome[i];
				}
				else {
					individual2[j+1].genome[i] = individual1[mom].genome[i];
					}
				}
			}
		}
	}
}

//mutation
void do_mutation(Individual* individuals){

	int child, bit, value, quant, i;
	
	child = (rand_r(&myseed)%population_size);
	quant = rand_r(&myseed)%4;
	
	for(i=0; i<quant; i++) {
	
		bit = (rand_r(&myseed)%genome_size);
		value = individuals[child].genome[bit];
		if(child != 0){
			if (value==0){
			individuals[child].genome[bit] = 1;
			}
			else {
				individuals[child].genome[bit] = 0;
			}
		}
	}

}


//elitism

void elitism(Individual* individual1, Individual* individual2) {
int i;

//pick the greatest fitness and select the greatest individual

	quick_sort(individual1, 0, population_size-1);
	individual1[0]=individual1[population_size-1];

	for (i=1; i < population_size-1; i++) {

		individual1[i]=individual2[i-1];
	}
	
	best_one = individual1[0];
}


//Write in file mean
void write_file_mean( ){

	FILE *m;
   m = fopen("results/mean.txt", "a");
   fprintf(m, "%f\n", mean);
   fclose(m);
}

//Write in file stdev
void write_file_stdev( ){

	FILE *s;
   s = fopen("results/stdev.txt", "a");
   fprintf(s, "%f\n", stdev);
   fclose(s);
}

//Write in file best
void write_file_best( ){

	FILE *b;
   b = fopen("results/best.txt", "a");
   fprintf(b, "%f\n", best_one.fitness);
   fclose(b);
}


//Write in file best individual complete
void write_file_best_complete( ){

	FILE *b;
   b = fopen("results/best_complete.txt", "a");
   fprintf(b, "%f\n", best_one.v_value);
	fprintf(b, "%f\n", best_one.c_value);
   fprintf(b, "%f\n", best_one.fitness);
   fprintf(b, "\n\n\n\n");
   fclose(b);
}

//--------------------------------------------------------------------------------------------------------------------------------
int main(int argc, char **argv)
{
	//Create results directory
	struct stat st = {0};
	if (stat("results/", &st) == -1) {
    mkdir("results/", 0700);
	}


	population_size = atoi(argv[1]);
	generations_number = atoi(argv[2]);

	clock_t tic = clock();
	//srand( (unsigned)time(NULL) );
	srand( 100000 );
	int mut;
	int n=0;

	Individual* parents = (Individual*)malloc(population_size*sizeof(Individual));
	Individual* children = (Individual*)malloc(population_size*sizeof(Individual));


	/*
	* You're supposed to put the file to train SVM here
	*/
	read_problem("enzymes.txt");
		//read_problem("frequenciaAbsolutaRelativaRegiao+Toda.txt");

//---------------------------------------------------------------------------------------First Population-------------------------------------------------------------------
	printf("Creating population...\n\n");

	random_pop(parents);
	
	interpolation(parents);
	
	print_pop(children);
	
	printf("Verifying adaptation\n\n");

	evaluate(parents);

	quick_sort(parents, 0, population_size-1);
	
//print_pop(parents);

	best_one = parents[population_size-1];
	print_individual(best_one);
		
	calculate_mean(parents);
	write_file_mean();
	calculate_stdv(parents);
	write_file_stdev();
	write_file_best();
	write_file_best_complete( );
	
	
	

//-------------------------------------------------------------------------repeat for each generation----------------------------------------------------------------------------
	while (n < generations_number) {
	//while (best_one.fitness < 80) {

	printf("\nRunning Generation number= %d\n\n", n+1);

		printf("\nDoing crossover\n\n");
		
		do_crossover2(parents, children);
		
		elitism(parents, children);
		
		mut = (rand_r(&myseed)% 10); 
			if (mut > mutation_percent) {
				
				printf("\nMutating...\n\n");
				do_mutation(children);
			}
	
				
		interpolation(parents);

		printf("Evaluating...----------------------------------------------------------------------------------------------------");
		evaluate(parents);
		
		print_individual(parents[population_size-1]);	
		
		calculate_mean(parents);
		write_file_mean();
		calculate_stdv(parents);
		write_file_stdev();
		write_file_best();
		write_file_best_complete( );
	

		printf("--------------------------------------------------------------------------------------------------------------------------------------------------\n\n\n");

		n++;
		svm_destroy_param(&param);
		quick_sort(parents, 0, population_size-1);

	}

	
	print_individual(best_one);
	
	 clock_t toc = clock();

    printf("\n\nElapsed: %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);
    
    // save time in file
	FILE *b;
	b = fopen("best.txt", "a");
	fprintf(b, "%f\n\n", (double)(toc - tic) / CLOCKS_PER_SEC);
	fclose(b);

	//svm_destroy_param(&param);
	free(prob.y);
	free(prob.x);
	free(x_space);
	free(line);

	return 0;
}

double do_cross_validation()
{
	int i;
	int total_correct = 0;
	double total_error = 0;
	double sumv = 0, sumy = 0, sumvv = 0, sumyy = 0, sumvy = 0;
	double *target = Malloc(double,prob.l);
	double acc;

	svm_cross_validation(&prob,&param,nr_fold,target);
	
	if(param.svm_type == EPSILON_SVR ||
	   param.svm_type == NU_SVR)
	{
		for(i=0;i<prob.l;i++)
		{
			double y = prob.y[i];
			double v = target[i];
			total_error += (v-y)*(v-y);
			sumv += v;
			sumy += y;
			sumvv += v*v;
			sumyy += y*y;
			sumvy += v*y;
		}
		printf("Cross Validation Mean squared error = %g\n",total_error/prob.l);
		printf("Cross Validation Squared correlation coefficient = %g\n",
			((prob.l*sumvy-sumv*sumy)*(prob.l*sumvy-sumv*sumy))/
			((prob.l*sumvv-sumv*sumv)*(prob.l*sumyy-sumy*sumy))
			);
			free(target);
			return(-1);
	}
	else
	{
		for(i=0;i<prob.l;i++)
			if(target[i] == prob.y[i])
				++total_correct;
		//printf("Cross Validation Accuracy = %g%%\n\n\n\n",100.0*total_correct/prob.l);
		acc = 100.0*total_correct/prob.l;
		free(target);
		return(acc);
		
	}

}

void parse_command_line(double v, double c)
{
	//void (*print_func)(const char*) = NULL;	// default printing to stdout

	// default values
	param.svm_type = C_SVC;
	param.kernel_type = RBF;
	param.degree = 3;
    param.gamma = 0.0;     
	param.gamma = v;
	param.coef0 = 0;
	param.nu = 0.5;
	param.cache_size = 100;
	param.C = c;
	param.eps = 1e-3;
	param.p = 0.1;
	param.shrinking = 1;
	param.probability = 0;
	param.nr_weight = 0;
	param.weight_label = NULL;
	param.weight = NULL;
	cross_validation = 1;
        nr_fold = 10;
	
	//svm_set_print_string_function(print_func);

	// determine filenames


}

// read in a problem (in svmlight format)

void read_problem(const char *filename)
{
	int max_index, inst_max_index, i;
	size_t elements, j;
	FILE *fp = fopen(filename,"r");
	char *endptr;
	char *idx, *val, *label;

	if(fp == NULL)
	{
		fprintf(stderr,"can't open input file %s\n",filename);
		exit(1);
	}

	prob.l = 0;
	elements = 0;

	max_line_len = 1024;
	line = Malloc(char,max_line_len);
	while(readline(fp)!=NULL)
	{
		char *p = strtok(line," \t"); // label

		// features
		while(1)
		{
			p = strtok(NULL," \t");
			if(p == NULL || *p == '\n') // check '\n' as ' ' may be after the last feature
				break;
			++elements;
		}
		++elements;
		++prob.l;
	}
	rewind(fp);

	prob.y = Malloc(double,prob.l);
	prob.x = Malloc(struct svm_node *,prob.l);
	x_space = Malloc(struct svm_node,elements);

	max_index = 0;
	j=0;
	for(i=0;i<prob.l;i++)
	{
		inst_max_index = -1; // strtol gives 0 if wrong format, and precomputed kernel has <index> start from 0
		readline(fp);
		prob.x[i] = &x_space[j];
		label = strtok(line," \t\n");
		if(label == NULL) // empty line
			exit_input_error(i+1);

		prob.y[i] = strtod(label,&endptr);
		if(endptr == label || *endptr != '\0')
			exit_input_error(i+1);

		while(1)
		{
			idx = strtok(NULL,":");
			val = strtok(NULL," \t");

			if(val == NULL)
				break;

			errno = 0;
			x_space[j].index = (int) strtol(idx,&endptr,10);
			if(endptr == idx || errno != 0 || *endptr != '\0' || x_space[j].index <= inst_max_index)
				exit_input_error(i+1);
			else
				inst_max_index = x_space[j].index;

			errno = 0;
			x_space[j].value = strtod(val,&endptr);
			if(endptr == val || errno != 0 || (*endptr != '\0' && !isspace(*endptr)))
				exit_input_error(i+1);

			++j;
		}

		if(inst_max_index > max_index)
			max_index = inst_max_index;
		x_space[j++].index = -1;
	}

	if(param.gamma == 0 && max_index > 0)
		param.gamma = 1.0/max_index;

	if(param.kernel_type == PRECOMPUTED)
		for(i=0;i<prob.l;i++)
		{
			if (prob.x[i][0].index != 0)
			{
				fprintf(stderr,"Wrong input format: first column must be 0:sample_serial_number\n");
				exit(1);
			}
			if ((int)prob.x[i][0].value <= 0 || (int)prob.x[i][0].value > max_index)
			{
				fprintf(stderr,"Wrong input format: sample_serial_number out of range\n");
				exit(1);
			}
		}

	fclose(fp);
}
