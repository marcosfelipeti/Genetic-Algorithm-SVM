#include "ag.h"
#include "svm_struct_ag.h"
#include "time.h"
#include <math.h>




void imprimirIndividuo(Individuo individuo) {

	int j;

	printf("Genoma: ");
	for (j=0; j<tamanhoGenoma; j++) {
		printf("%d", individuo.genoma[j]);
	}

	printf("\nValor t: %d\n", individuo.valor_t); 
	printf("Valor d: %.1f\n", individuo.valor_d); 
	printf("Valor c: %f\n", individuo.valor_c); 
	printf("Fitness: %f\n", individuo.fitness);
	printf("Posição 1: %d\n", individuo.pos1);
	printf("Posição 2: %d\n", individuo.pos2);
	printf("\n\n\n");
}

void imprimirIndividuos(Individuo* individuos) {

	int i, j;

	for (i=0; i<tamanhoPopulacao; i++) {
		printf("Índice: %d\n", i);
		printf("Genoma: ");
		for (j=0; j<tamanhoGenoma; j++) {
			printf("%d", individuos[i].genoma[j]);
		}
		printf("\nValor t: %d\n", individuos[i].valor_t); 
		printf("Valor d: %.1f\n", individuos[i].valor_d); 
		printf("Valor c: %f\n", individuos[i].valor_c); 

		printf("Fitness: %f\n", individuos[i].fitness);

		printf("Posição 1: %d\n", individuos[i].pos1);
		printf("Posição 2: %d\n", individuos[i].pos2);
		printf("\n\n\n");

	}
}

int geraAleatorio (Individuo* individuos){

	int i, j;

	for (i = 0; i < tamanhoPopulacao; i++){
		for (j = 0; j < tamanhoGenoma; j++){
			individuos[i].genoma[j] = (rand()% 2); //MAX ter. \a que ser definido como 2. Ai fica gerando 									numeros 0 e 1, ai nem precisa somar esse +1
		}	
	}
}

int calcularB10(Individuo* individuos, int indice, int inicio, int fim) {

	int j, soma=0;
	int exp = 0;

	for(j=inicio; j < fim; j++) {
		soma = soma + ((pow(2, exp)) * individuos[indice].genoma[j]);

		exp++;
	}
  return soma;
}

void tranformarBinDec(Individuo* individuos) {

	int i, b10;

	for(i=0; i < tamanhoPopulacao; i++) {
		
		b10 = calcularB10(individuos, i, 0, TAM_T);
		individuos[i].valor_t = round(MINt + (b10/pow(2, TAM_T)) * (MAXt - MINt));


		if (individuos[i].valor_t == 1) {
			b10 = calcularB10(individuos, i, TAM_T, TAM_C);
			individuos[i].valor_d = round(MINd + (b10/pow(2, TAM_D)) * (MAXd - MINd));
		}

		else if (individuos[i].valor_t == 2) {
			b10 = calcularB10(individuos, i, TAM_T, TAM_C);
			individuos[i].valor_d = MINg + (b10/pow(2, TAM_D)) * (MAXg - MINg);
		}
		else individuos[i].valor_d = 0;

		b10 = calcularB10(individuos, i, TAM_C, tamanhoGenoma);
		individuos[i].valor_c = MINc + (b10/pow(2, TAM_C)) * (MAXc-MINc);
	}
}

void calcularValory(char *filetrain, char *filetest, Individuo* individuos) { //calcula fitness


	int i;
	for (i=0; i < tamanhoPopulacao; i++)
	{
		//individuos[i].fitness = 1 - run_svm(individuos[i].valor_c, individuos[i].valor_t, individuos[i].valor_d);
		individuos[i].fitness = run_svm(individuos[i].valor_c, individuos[i].valor_t, individuos[i].valor_d);
		printf("\n %f %d %f", individuos[i].valor_c, individuos[i].valor_t, individuos[i].valor_d);
		printf(" => Fitness %f\n\n", individuos[i].fitness);
	}
}

void calcularFitness(Individuo* individuos) { //calcula percentual_roleta

	int i;
	int soma =0;

	for (i=0; i < tamanhoPopulacao; i++)
		soma = soma + individuos[i].fitness;
	for (i=0; i < tamanhoPopulacao; i++)
		individuos[i].percentual_roleta = soma/individuos[i].fitness;
}

void swap(Individuo* a, Individuo* b) {

	Individuo tmp;

	tmp = *a;
	*a = *b;
	*b = tmp;
}
 
int partition(Individuo* vec, int left, int right) {

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
 
void quickSort(Individuo* vec, int left, int right) {

	int r;
 
	if (right > left) {
		r = partition(vec, left, right);
		quickSort(vec, left, r - 1);
		quickSort(vec, r + 1, right);
  	}
}

int determinaQuantPos () {

	int i; int tam=0;

	//determinar quantas posicoes serão
	for (i = 0; i < tamanhoPopulacao+1; i++)
		tam = tam + i;
	printf("Número de posições = %d\n\n", tam);

	return tam;
}

void determinarPosicoes(Individuo* individuos, int tam) { //mudei porque é maximizacao

	int i;

	for (i = 0; i < tamanhoPopulacao; i++) {

		if (i == 0) {
			individuos[i].pos1 = 0;
			individuos[i].pos2 = tamanhoPopulacao - 1;
		}
		else if (i==tamanhoPopulacao-1) {
			individuos[i].pos1 = tam-1;
			individuos[i].pos2 = tam-1;
		}
		else {
			individuos[i].pos1 = individuos[i-1].pos2 + 1;
			individuos[i].pos2 = individuos[i].pos1 + (tamanhoPopulacao - (i+1));
		}
	}
}

void criarRank(Individuo* individuos, int tam) {

	quickSort(individuos, 0, tamanhoPopulacao-1);
	determinarPosicoes(individuos, tam);
}

int encontrarSorteado(Individuo* individuos, int n_sorteado) {

	int i;

	for(i=0; i < tamanhoPopulacao; i++) {
		if( (individuos[i].pos1 <= n_sorteado) && (individuos[i].pos2 >= n_sorteado) ) {
			return i;
		}
	}
}

void crossover(Individuo* individuos, Individuo* individuos2, int cont, int pai1, int pai2) {

	int pontodecorte;
	int i;

	pontodecorte = (rand()% tamanhoGenoma);

	for(i=0; i<tamanhoGenoma; i++) {
		individuos2[cont].genoma[i] = individuos[pai1].genoma[i];
		individuos2[cont+1].genoma[i] = individuos[pai2].genoma[i];
	}
	for (i = pontodecorte; i < tamanhoGenoma; i++) {
		individuos2[cont].genoma[i] = individuos[pai2].genoma[i];
		individuos2[cont+1].genoma[i] = individuos[pai1].genoma[i];
	}
}

void selecionarPais(Individuo* individuos, Individuo* individuos2, int tam) {

	int sorteio1, sorteio2;
	int cont=0;
	int pai1, pai2;
	int i, muta;

	for (i=0; i<(tamanhoPopulacao/2); i++) {
		sorteio1 = (rand()% tam);
		sorteio2 = (rand()% tam);
		pai1 = encontrarSorteado(individuos, sorteio1);
		pai2 = encontrarSorteado(individuos, sorteio2);
		crossover(individuos, individuos2, cont, pai1, pai2);
		cont=cont+2;	// contador para ir salvando os filhos no vetor de filhos[cont].		
	}
}

void mutacao(Individuo* individuos){

	int filho, bit, valor;
	
	filho = (rand()%tamanhoPopulacao);
	bit = (rand()%tamanhoGenoma);
	valor = individuos[filho].genoma[bit];
		if (valor==0){
			individuos[filho].genoma[bit] = 1;
		}
		else {
			individuos[filho].genoma[bit] = 0;
		}
}


//elitismo

void trocarFilhosPaisElitismo(Individuo* individuos, Individuo* individuos2) {
int i;
float maiorFitness1;
//pegar o maior fitness e selecionar o maior individuo

	maiorFitness1= individuos[0].fitness;

	for (i=0; i < tamanhoPopulacao; i++) {
		if (individuos[i+1].fitness > maiorFitness1){
			maiorFitness1=individuos[i+1].fitness;
			individuos[0]=individuos[i+1];
		}	
	//printf ("\n=====> FITNESS%f\n", maiorFitness1);

		individuos[i+1]=individuos2[i];
	}
}


void trocarFilhosPais(Individuo* individuos, Individuo* individuos2) {

	int i;

	for (i=0; i < tamanhoPopulacao; i++) {
		individuos[i] = individuos2[i];	//vetor de pais recebe o vetor de filhos - morrem todos os pais
	}
}

int verificarAccuracy(Individuo* individuos) {

	int i;

	for (i = 0; i < tamanhoPopulacao; i++) {
		if (individuos[i].fitness > 0.99) {
			return 1;
		}
	}
return 0;
}

int main(int argc, char *argv[]){

	srand(time(NULL));

	int i, ac;
	int tam, muta;
	int n=0;

	srand(time(NULL)); // inicializa o rand
	Individuo* pais = (Individuo*)malloc(tamanhoPopulacao*sizeof(Individuo));
	Individuo* filhos = (Individuo*)malloc(tamanhoPopulacao*sizeof(Individuo));

	printf("Gerando população\n\n");
	geraAleatorio(pais);
	printf("Gerando valor de x - Transformando Binário Decimal. \n\n");
	tranformarBinDec(pais);
	printf("Calculando valor de y\n\n");
	calcularValory(argv[1], argv[2], pais);



	//while ((n < nGeracoes) || (ac=verificarAccuracy(pais) == 0.0)) {
	while ((n < nGeracoes) || (ac=verificarAccuracy(pais) < 0.0)) {

		printf("Determinando Quantidade de posições\n\n");
		tam = determinaQuantPos ();

		printf("Criando rank\n\n");
		criarRank(pais, tam);

		imprimirIndividuo(pais[0]);

		imprimirIndividuo(pais[1]);

		printf("Selecionando Pais e Crossover+Mutação\n\n");
		selecionarPais(pais, filhos, tam);
                	
		muta = (rand()% porcentagemMutacao); //taxa de mutação = 10. Rand ira gerar numeros de 0 a 9. 
			if (muta == 1) {
				mutacao(filhos);
			}
		
		printf("vetor de filhos -> vetor de pai\n\n");
		//trocarFilhosPais(pais, filhos);
		 trocarFilhosPaisElitismo(pais, filhos);
					//printf ("=====> fitness%f\n", maiorFitness1);
		printf("Gerando valor de x - Transformando Binário Decimal. \n\n");
		tranformarBinDec(pais);

		printf("Calculando valor de y\n\n");
		calcularValory(argv[1], argv[2], pais);



		n++;
	}

	imprimirIndividuos(pais);
}
