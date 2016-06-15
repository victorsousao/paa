#include<stdio.h>
#include<stdlib.h>
#include <iostream>
#include <list>
#include <queue>
#include <ctime>
#define TAM 13
#define INFINITY 999
#define INF INT_MAX

using namespace std;

//Funções que criam grafos de teste
void construirGrafo1();
void construirGrafo2();
void construirGrafo3();
void construirGrafo4();
void construirGrafo5();

//Funções diversas
void conecta(int a, int b, int v, bool orientado);
void imprimeMatrizAdjacencia(int m[TAM][TAM]);
bool grafoEhCompleto();
void dfs(int m[TAM][TAM], int u, bool imprimir=false);
bool conexo(int m[TAM][TAM]);
void completarGrafo();
void complementoGrafo();
void verificarQuantidadeComponentes(int m[TAM][TAM]);
void buscaProfundidade(int m[TAM][TAM], int u);
void buscaLargura(int m[TAM][TAM], int u);
void descobrirGrauVertices(int m[TAM][TAM]);

//Funções para Verificar o ciclo de Hamilton
int cicloHamilton(int m[TAM][TAM]);
void imprimirSolucao(int path[]);
int cicloHamiltonAux(int m[TAM][TAM], int path[], int pos);
int verificaGrafo(int v, int m[TAM][TAM], int path[], int pos);

//Funções para TSP Força Bruta
int tsp(int c [TAM][TAM]);
int tsp_r (int c[][TAM], int tour[], int start, int n);


/*********************************************************************/
// Nós da arvore TSP BB
struct Node
{
    //Armazena as arestas da arvore
	//Ajuda no rastreamento do caminho quando encontrado
    vector<pair<int, int> > path;

    // armazena a matriz reduzida
    int matrizReduzida[TAM][TAM];

    // armazena o menor custo
    int cost;

	//armazena o vertice atual
    int vertex;

	// Armazena o nivel mais fundo
    int level;
};
/*********************************************************************/

// Funções para TSP Branch and Bound
Node* newNode(int parentMatrix[TAM][TAM], 
			vector<pair<int, int> > path, 
			int level, int i, int j);
int reduzirLinha(int matrizReduzida[TAM][TAM], int row[TAM]);
int reduzirColuna(int matrizReduzida[TAM][TAM], int col[TAM]);
int calcularCusto(int matrizReduzida[TAM][TAM]);
void printPath(vector<pair<int, int> > list);
int tsp_bb(int costMatrix[TAM][TAM]);


//Variaveis globais necessarias no programa
int visitado[TAM]; // vetor de visitação
int grau[TAM]; //Grau de cada vertice do grafo
int adj[TAM][TAM]; // matriz de adjacencia
int adj_c[TAM][TAM]; // cópia da matriz de adjacencia
int adj_complemento[TAM][TAM]; // armazenar o complemento do grafo
int adj_completo[TAM][TAM]; // armazenar o grafo completo
int n_componentes=0; //Numero de componentes do grafo

//Funções para a AGM
void ArvoreGeradoraMinima(int m[TAM][TAM]);
int find(int posicao[], int vertice);
void unionFind(int posicao[],int componente_01,int componente_02);
void ordena();
void imprimirArestaEPesoArvore();

// Struct de arestas para a AGM
typedef struct aresta
{
    int u,v,w;

} aresta;

typedef struct listaresta
{
    aresta data[TAM];
    int n;

} listaresta;

//Variaveis que serão utilizadas na AGM
listaresta e_lista;
listaresta dist_lista;



/**
* Função principal que chama os demais metodos, faz a inicialização básica de algumas váriaveis
*/
int main() 
{
	//construirGrafo1();
	//construirGrafo2();
	//construirGrafo3();
	//construirGrafo4();
	construirGrafo5();
	
	imprimeMatrizAdjacencia(adj);
	
	conexo(adj);
	
	grafoEhCompleto();
	
	verificarQuantidadeComponentes(adj);
	
	descobrirGrauVertices(adj);
	
	buscaLargura (adj, 0);
	
	buscaProfundidade(adj, 0);
	
	cicloHamilton(adj);
	
    ArvoreGeradoraMinima(adj);
    
   
	for (int i=0; i<TAM; i++){
		for (int j=0; j<TAM; j++){
			if(adj[i][j] == 0)
				adj_c[i][j]	 = 999;
			else
				adj_c[i][j] = adj[i][j];
		}
	}
	
   	time_t now = time(0);
	char* dt = ctime(&now);
	cout << "  Força Bruta Inicio: " << dt;
   	//tsp(adj_c);
	now = time(0);
	dt = ctime(&now);
	cout << "  Força Bruta Fim: " << dt << endl;
	
	
	for (int i=0; i<TAM; i++){
		for (int j=0; j<TAM; j++){
			if(adj[i][j] == 0)
				adj_c[i][j]	 = INF;
			else
				adj_c[i][j] = adj[i][j];
		}
	}
	
	
	now = time(0);
	dt = ctime(&now);
	cout << "  BB Inicio: " << dt << endl;
	tsp_bb(adj_c);
	now = time(0);
	dt = ctime(&now);
	cout << "  BB Fim: " << dt << endl;
	
	getchar();
	
	//complementoGrafo();
	//completarGrafo();
}



/*********************************************************************/
/**
 * Faz a busca em largura a partir de um determinado vertice
 */
void buscaLargura(int m[TAM][TAM], int u)
{
	// limpa todas as marcas de visitação
	for (int i=0; i < TAM; i++) visitado[i] = 0;
	
	printf("\n  Largura: \n");
	printf("  >> %d", u+1); //Imprime o vertie inicial
	
	queue<int> fila; //Cria a fila para percurso
	fila.push(u); //Enfileira o primeiro vertice
	visitado[u] = 1; //Marca o primeiro vertice como visitado
	
	while (!fila.empty()){
		u = fila.front();
		fila.pop();		
		for(int i=0; i< TAM; i++){
			//Verifica se tem peso na aresta e se o vertice
			//ainda não foi visitado
			if(m[u][i] >= 1 && !visitado[i]){
				printf("  >> %d ", i+1);
				fila.push(i);
				visitado[i] = 1;
			}
		}
	}
	
	printf("\n");
}
/*********************************************************************/



/*********************************************************************/
/**
 * Verifica o grau de cada vertice e armazena em um vetor
 */
void descobrirGrauVertices(int m[TAM][TAM])
{
	int d_grau;
	printf("\n  Grau de cada vertice: \n");
	for(int i=0; i<TAM; i++) {
		d_grau=0;
		for(int j=0; j<TAM; j++){
			if(m[i][j] >= 1){
				d_grau++;
			}
		}
		grau[i] = d_grau;
		printf("  Grau do vertice %d: %d \n", i+1, d_grau);
	}
}
/*********************************************************************/



/*********************************************************************/
/**
 * Faz a busca em profundidade a partir de um determinado vertice
 */
void buscaProfundidade(int m[TAM][TAM], int u)
{
	// limpa todas as marcas de visitação
    for (int i=0; i < TAM; i++) visitado[i] = 0;
    printf("\n  Profundidade: \n");
	    
    //Faz a busca em profundidade a partir do primeiro vertice
    for (int i=0; i < TAM; i++){
        if (!visitado[i]){
        	dfs(adj, i, true);
        	printf("\n");
		}
	}
}
/*********************************************************************/



/*********************************************************************/
/**
 * Percorre o grafo e se tiver alguma posição que o peso da aresta é 0
 * não sendo i = j então o grafo não é completo
 */
bool grafoEhCompleto()
{
	int i,j;
	for(i=0;i<TAM;i++){
		for(j=0;j<TAM;j++)
		{
			//Só completa se i != j e o peso da aresta for 0 para não ter loop
			if(i!=j && adj[i][j] == 0){
				printf("\n  Grafo nao e completo. \n");
				return false;
			}
		}
	}
	printf("\n  Grafo e completo. \n");
	return true;
}
/*********************************************************************/



/*********************************************************************/
/**
 * Percorre o grafo e onde estiver ligação 0 passa para um não percorre
 * quando i = j para não ter auto loop no vertice
 */
void completarGrafo()
{
	int i,j;
	for(i=0;i<TAM;i++){
		for(j=0;j<TAM;j++)
		{
			//Só completa se i != j e o peso da aresta for 0 para não ter loop
			if(i!=j && adj[i][j] == 0){
				adj[i][j] = 1;
			}
			
		}
	}
}
/*********************************************************************/



/*********************************************************************/
/**
 * Gera um grafo complementar e armazena em uma nova matriz adj_complemento
 */
void complementoGrafo()
{
	int i,j;
	for(i=0;i<TAM;i++){
		for(j=0;j<TAM;j++)
		{
			//Só completa se i != j para não ter loop
			if(i!=j){
				if(adj[i][j] == 0){
					adj_complemento[i][j] = 1;
				}
			}
		}
	}
}
/*********************************************************************/



/*********************************************************************/
/**
 * Conta quantos componentes tem um determinado grafo informado
 */
void verificarQuantidadeComponentes(int m[TAM][TAM]){
	
	//Inicia o numero de componentes com 0
	n_componentes=0;
	// limpa todas as marcas de visitação
    for (int i=0; i < TAM; i++) visitado[i] = 0;
 
    // inicia busca no nó de índice 0
    dfs(m, 0);
    
    // Já coloca o número de componentes sendo 1
	// porque já leu o primeiro vertice
    n_componentes++; 
 
    // testa se todos os nós foram visitados com uma única busca
    for (int i=0; i < TAM; i++){
        if (!visitado[i]){
        	dfs(m, i, false); //Faz a busca em profundidade do próximo vertice não visitado
        	n_componentes++; //Incrementa o numero de componentes
		}
	}
        	
	//Imprime o numero de componentes do grafo       
	printf("  Numero de componentes do grafo: %d \n", n_componentes);
}
/*********************************************************************/



/*********************************************************************/
/**
 * Faz a conexão de um vertice a outro colocando o valor da aresta
 */
void conecta(int a, int b, int v, bool orientado)
{
	adj[a-1][b-1]=v;
	if(!orientado)
		adj[b-1][a-1]=v;	
}
/*********************************************************************/



/*********************************************************************/
/**
 * Faz a busca por profundidade em todos os elementos da matriz
 * e armazena no vetor de visitado
 **/
void dfs(int m[TAM][TAM], int u, bool imprimir) {
    if (visitado[u])
        return;
    visitado[u] = 1;
    
    //Imprime o caminho que está sendo feito
    //Apenas para facilitar a validação
    //Não faz aprte da função realmente
    if(imprimir) printf("  >> %d", u+1); 
    
    for (int v=0; v < TAM; v++) if (m[u][v]) {
        dfs(m, v, imprimir);
    }
}
/*********************************************************************/



/*********************************************************************/
/**
 * Retorna verdadeiro ou falso para grafo conexo
 **/
bool conexo(int m[TAM][TAM]) {
    // limpa todas as marcas de visitação
    for (int i=0; i < TAM; i++) visitado[i] = 0;
 
    // inicia busca no nó de índice 0
    dfs(m, 0, false);
 
    // testa se todos os nós foram visitados com uma única busca
    for (int i=0; i < TAM; i++)
        if (!visitado[i]){
        	printf("GRAFO NAO CONEXO \n");
        	return false;
		}
            
	printf("GRAFO CONEXO \n");       
    return true;
}
/*********************************************************************/



/*********************************************************************/
/**
 * Imprime a matriz de adjacencia 
 */
void imprimeMatrizAdjacencia(int m[TAM][TAM])
{
	int i,j,k;
	//system("cls");
	i=j=0;
	k=1;
	printf("       ");
	while(k!=(TAM+1))
	{
	    printf("  %d  ",k); //
	    k++;
	}
	printf("  \n\n\n  ");
	k=1;
	while(i!=TAM)
	{
		printf("  %d  ",k);//
		while(j!=TAM)
		{
			printf("  %d  ",m[i][j]); //
			j++;
		}
		i++;
		j=0;
		k++;
		printf("  \n\n  ");
	}
}
/*********************************************************************/



/*********************************************************************/
/**
 * Função que testa se o vertice já não foi visitado, se tiver sido,
 * retorna 0 para ir a pilha de cima e começar novamente
 */
int verificaGrafo(int v, int m[TAM][TAM], int path[], int pos)
{
	int i;

	if(m[path[pos-1]][v] == 0)
	return 0;

	for(i = 0; i < pos; i++)
	   if(path[i] == v)
	   return 0;

	return 1;
}
/*********************************************************************/


/*********************************************************************/
/**
 * Função recursiva que passa por todos os vertices para encontrar o ciclo
 */
int cicloHamiltonAux(int m[TAM][TAM], int path[], int pos)
{
	int v;
	if (pos == TAM)
	   {
		   if (m[path[pos - 1]][path[0]] >= 1)
	       return 1;
	       else
	       return 0;
		}
	
	for(v = 0; v < TAM; v++)
	{
		if(verificaGrafo(v,m,path,pos))
		{
			path[pos] = v;
			if(cicloHamiltonAux(m,path,pos + 1) == 1)
			return 1;
	
			path[pos] = -1;
	
		}
	}
	return 0;
}
/*********************************************************************/


void imprimirSolucao(int path[])
{
	printf("  Ciclo: ");
	int i;
	for(i = 0; i < TAM; i++)
	   printf("  >> %d",path[i]+1);

	printf(" >> %d",path[0]+1);
	printf("\n");
}

/*********************************************************************/
/**
 * Função que verifica se o grafo é Hamiltoniano
 */
int cicloHamilton(int m[TAM][TAM])
{
	int path[TAM], i;

	printf("\n  Hamiltoniano: \n");

	for(i = 0; i < TAM; i++)
		path[i] = -1;

	path[0] = 0;

	if(cicloHamiltonAux(m,path,1) == 0)
	{
		printf("  Nao existe um ciclo Hamiltoniano no grafo informado!\n");
		return 0;
	}

	imprimirSolucao(path);
	return 1;
}
/*********************************************************************/



/*********************************************************************/
/**
 * Abaixo seguem todas as funções necessárias para AGM - Algoritmo Kruskal 
 */
void ArvoreGeradoraMinima(int m[TAM][TAM])
{
	
	printf("\n  Arvore Geradora Minima (AGM) - Kruskal\n");
    int posicao[TAM], i, j, cnp_01, cnp_02;
    e_lista.n = 0;

    for (i = 1; i < TAM; i++)
        for (j = 0; j < i; j++)
        {
            if (m[i][j] != 0)
            {
                e_lista.data[e_lista.n].u = i;
                e_lista.data[e_lista.n].v = j;
                e_lista.data[e_lista.n].w = m[i][j];
                e_lista.n++;
            }
        }

    ordena();

    for(i = 0; i < TAM; i++)
        posicao[i] = i;

    dist_lista.n = 0;

    for(i = 0; i < e_lista.n; i++)
    {
        cnp_01 = find(posicao,e_lista.data[i].u);
        cnp_02 = find(posicao,e_lista.data[i].v);

        if (cnp_01 != cnp_02)
        {
            dist_lista.data[dist_lista.n] = e_lista.data[i];
            dist_lista.n = dist_lista.n + 1;
            unionFind(posicao, cnp_01, cnp_02);
        }
    }
    
    imprimirArestaEPesoArvore();
}
/*********************************************************************/

int find(int posicao[], int vertice)
{
    return(posicao[vertice]);
}

void unionFind(int posicao[],int componente_01,int componente_02)
{
    int i;

    for(i = 0; i < TAM ; i++)
        if(posicao[i] == componente_02)
            posicao[i] = componente_01;
}

void ordena()
{
   int i,j;
   aresta temp;

   for(i = 1; i < e_lista.n; i++)
      for(j = 0; j < e_lista.n-1; j++)
         if(e_lista.data[j].w > e_lista.data[j + 1].w)
            {
               temp = e_lista.data[j];
               e_lista.data[j] = e_lista.data[j + 1];
               e_lista.data[j + 1] = temp;
            }
}


void imprimirArestaEPesoArvore()
{
   printf("  Aresta \tPeso\n");
   int i,cost = 0;

   for(i = 0; i < dist_lista.n; i++)
      {
        printf("\n  %d\t\t%d", dist_lista.data[i].u, dist_lista.data[i].w);
        cost = cost + dist_lista.data[i].w;
      }

      printf("\n\n  Custo da Arvore Geradora -> %d\n\n", cost);
}
/*********************************************************************/



/*********************************************************************/
/** 
 * Função de TSP força bruta que inicia a busca do melhor caminho
 * a partir do primeiro vertice O(n!)
 */
int tsp(int m[TAM][TAM])
{
	int i, j; 
	int caminho[TAM]; /* Matriz de caminho. */
	int custo; /* Custo Minimo. */
	 
	printf ("  Problema do caixeiro viajante - Forca Bruta: \n");
		 
	if(!grafoEhCompleto()){
		printf("  O grafo nao e completo, complete antes de prosseguir com o PCV. \n");
		return 0;
	}
	
		 
	for (i=0; i<TAM; i++)
		caminho[i] = i;
	 
	custo = tsp_r (m, caminho, 0, TAM);
	 
	printf ("  Custo minimo: %d.\n  Caminho: ", custo);
	for (i=0; i<TAM; i++)
		printf ("%d ", caminho[i]+1);
	printf ("1\n");
	return 1;
}
/*********************************************************************/


/*********************************************************************/
/**
 * Função que faz a permutação das posições para encontrar
 * o melhor caminho no caixeiro viajante
 */
int tsp_r (int c[][TAM], int tour[], int start, int n)
{
	int i, j, k; 
	int temp[TAM]; /* Vetor temporario para calculos. */
	int caminho_min[TAM]; /* Vetor de caminho minimo. */
	int custo_min; /* Custo minimo. */
	int custo_atual; /* Custo atual. */
	 
	/* Condição de parada da recursão. */
	if (start == n - 2)
		return c[tour[n-2]][tour[n-1]] + c[tour[n-1]][0];
	 
	/* Calcule o caminho a partir da cidade atual. */
	custo_min = INFINITY;
	for (i = start+1; i<n; i++)
	{
		for (j=0; j<n; j++)
			temp[j] = tour[j];
	 
		/* Ajusta as posições. */
		temp[start+1] = tour[i];
		temp[i] = tour[start+1];
		 
		/* Encontrado um ciclo melhor? (Recorrência derivável.) */
		if (c[tour[start]][tour[i]] + 
		(custo_atual = tsp_r (c, temp, start+1, n)) < custo_min) {
			custo_min = c[tour[start]][tour[i]] + custo_atual;
			for (k=0; k<n; k++)
				caminho_min[k] = temp[k];
		}
	}
	 
	/* Definir o vetor de caminho-min. */
	for (i=0; i<n; i++)
		tour[i] = caminho_min[i];
	 
	return custo_min;
}
/*********************************************************************/


//Funções para o TSP Branch and Bound

/*********************************************************************/
// Função para alocar um novo nó
// (i, j) corresponde a cidade visitada j apartir de i 
Node* newNode(int parentMatrix[TAM][TAM], 
			vector<pair<int, int> > path, 
			int level, int i, int j)
{
    Node* node = new Node;

	//armazena a aresta anterior
    node->path = path;
	// pula o nó principal
    if(level != 0 )
		//adiciona a aresta ao caminho
        node->path.push_back(make_pair(i, j));

    //copia dados do nó pai para o filho
    memcpy(node->matrizReduzida, parentMatrix, 
    		sizeof node->matrizReduzida);

	//Altera todas as entradas para INF 
	//Pula para o nó RAIZ
    for(int k = 0; level != 0 && k < TAM; k++)
		//definir as arestas de saida das cidade i para INF
        node->matrizReduzida[i][k] = INF,
		//definir as arestas de entrada da cidade j para INF
        node->matrizReduzida[k][j] = INF;

	// Definir (j, 0) para infinito
	// começar do nó 0
    node->matrizReduzida[j][0] = INF;

    //definir o numero de cidades visitadas
    node->level = level;
	
	// atribuir o número atual cidade
    node->vertex = j;

	// return node
    return node;
}
/*********************************************************************/


/*********************************************************************/
// Funções de reduzir cada linha, de tal maneira que
// deve haver pelo menos um zero em cada fileira
int reduzirLinha(int matrizReduzida[TAM][TAM], int row[TAM])
{
	// inicializar matriz de linha para INF
    fill_n(row, TAM, INF);
	
	// row[i] contem o minimo de i
    for(int i = 0; i < TAM; i++)
        for(int j = 0; j < TAM; j++)
            if( matrizReduzida[i][j] < row[i])
                row[i] = matrizReduzida[i][j];

	// reduzir o valor mínimo de cada elemento em cada linha.
    for(int i = 0; i < TAM; i++)
        for(int j = 0; j < TAM; j++)
            if(matrizReduzida[i][j] != INF && row[i] != INF)
                matrizReduzida[i][j] -= row[i];
}
/*********************************************************************/


/*********************************************************************/
// Funções de reduzir cada coluna de tal forma que
// deve haver pelo menos um zero em cada coluna
int reduzirColuna(int matrizReduzida[TAM][TAM], int col[TAM])
{
	// inicializar matriz de coluna para INF
    fill_n(col, TAM, INF);

	// col[j] contem o minimo de j
    for(int i = 0; i < TAM; i++)
        for(int j = 0; j < TAM; j++)
            if(matrizReduzida[i][j] < col[j])
                col[j] = matrizReduzida[i][j];

	// rreduzir o valor mínimo de cada elemento em cada coluna
    for(int i = 0; i < TAM; i++)
        for(int j = 0; j < TAM; j++)
            if(matrizReduzida[i][j] != INF && col[j] != INF)
                matrizReduzida[i][j] -= col[j];	
}
/*********************************************************************/

/*********************************************************************/
// Função para obter o limite inferior
// Sobre o caminho que começa no nó min atual
int calcularCusto(int matrizReduzida[TAM][TAM])
{
	// inicializar custo para 0
    int cost = 0;
	
	// Redução de linha
    int row[TAM];
	reduzirLinha(matrizReduzida, row);
	
	// Reduçao de coluna
    int col[TAM];
    reduzirColuna(matrizReduzida, col);

	// O custo total esperado
	// É a soma de todas as reduções
    for(int i = 0; i < TAM; i++)
        cost += (row[i] != INT_MAX) ? row[i] : 0,
        cost += (col[i] != INT_MAX) ? col[i] : 0;

    return cost;
}
/*********************************************************************/

/*********************************************************************/
// Lista de impressão de cidades visitadas segundo o menor custo
void printPath(vector<pair<int, int> > list)
{
    for(int i = 0; i < list.size(); i++)
        cout << list[i].first + 1<< " -> " <<
                list[i].second + 1<< endl;
}
/*********************************************************************/


/*********************************************************************/
// Objeto de comparação a ser usado para ordenar a pilha
struct comp
{
    bool operator()(const Node* lhs, const Node* rhs) const
    {
        return lhs->cost > rhs->cost;
    }
};
/*********************************************************************/


/*********************************************************************/
// Função para resolver TSP usando Branch and Bound.
int tsp_bb(int costMatrix[TAM][TAM])
{
	//Cria uma fila de prioridade para armazenar os nós
    priority_queue<Node*, std::vector<Node*>, comp> pq;
	
    vector<pair<int, int> > v;
    // Cria um nó raiz e calcula o seu custo
	// O TSP começa a partir de primeira cidade isto é, o nó 0
    Node* root = newNode(costMatrix, v, 0, -1, 0);
	// pega o limite inferior do caminho que começa no nó 0
    root->cost = calcularCusto(root->matrizReduzida);

    //Adiciona a raiz para lista de nós
    pq.push(root);

    
	// Encontra um nó com menor custo,
    // Adiciona seus filhos a lista de nós vivos e
    // Finalmente retira da lista.
    while(!pq.empty())
    {
        // Encontrar um nó com menor custo estimado
        Node* min = pq.top();

        // O nó encontrado é excluído da lista
        pq.pop();

		// i armazena o numero da cidade atual
        int i = min->vertex;

        // se todas as cidades foram visitadas
        if(min->level == TAM - 1)
        {
            // retorna para a cidade inicial
            min->path.push_back(make_pair(i, 0));
			// imprime o caminho realizado
            printPath(min->path);
			
			// retorna caminho otimizado
            return min->cost;
        }

		// For para cada filho min
		// (i, j) cria uma aresta
        for(int j = 0; j < TAM; j++)
        {
            if(min->matrizReduzida[i][j] != INF)
            {
                // Cria um nó filho e calcula o seu custo
                Node* child = newNode(min->matrizReduzida, min->path, 
										min->level+1, i, j);

				/* Custo do filho  =
					custo do pai +
					custo da matriz reduzida [i][j] +
					custo da aresta (i,j)
					limite inferior do caminho a partir do nó j
				*/
                child->cost = min->cost + min->matrizReduzida[i][j] + 
                			  calcularCusto(child->matrizReduzida);

                // Adiciona o filho a lista de nós
                pq.push(child);
            }
        }
		// Libera o nó do caminho
        free(min);
    }
}
/*********************************************************************/




/*********************************************************************/
/**
 * Cria um grafo para usar como teste Grafo 1
 */
void construirGrafo1() 
{
	conecta(1,2,1,false);
	conecta(1,4,1,false);
	conecta(2,3,1,false);
	conecta(3,4,1,false);
}
/*********************************************************************/


/*********************************************************************/
/**
 * Cria um grafo para usar como teste Grafo 2
 */
void construirGrafo2() 
{
	conecta(1,2,1,false);
	conecta(1,5,1,false);
	conecta(2,3,1,false);	
	conecta(2,6,1,false);
}
/*********************************************************************/


/*********************************************************************/
/**
 * Cria um grafo para usar como teste Grafo 3
 */
void construirGrafo3()
{
	conecta(1,2,1,false);
	conecta(1,3,1,false);
	conecta(1,4,1,false);	
	conecta(2,4,1,false);
	conecta(2,3,1,false);
	conecta(3,4,1,false);
}
/*********************************************************************/


/*********************************************************************/
/**
 * Cria um grafo para usar como teste Grafo 4
 */
void construirGrafo4()
{
	conecta(1,2,3,false);
	conecta(1,4,2,false);
	conecta(1,5,5,false);	
	conecta(2,3,8,false);
	conecta(2,5,5,false);
	conecta(3,4,5,false);
	conecta(3,5,4,false);
}
/*********************************************************************/


void construirGrafo5()
{
	conecta(1,2,1,false);
	conecta(1,3,6,false);
	conecta(1,4,4,false);
	conecta(1,5,10,false);
	conecta(1,6,1,false);
	conecta(1,7,5,false);
	conecta(1,8,9,false);
	conecta(1,9,6,false);
	conecta(1,10,2,false);
	conecta(1,11,10,false);
	conecta(1,12,7,false);
	conecta(1,13,1,false);
	conecta(2,3,6,false);
	conecta(2,4,1,false);
	conecta(2,5,9,false);
	conecta(2,6,1,false);
	conecta(2,7,5,false);
	conecta(2,8,9,false);
	conecta(2,9,9,false);
	conecta(2,10,7,false);
	conecta(2,11,7,false);
	conecta(2,12,1,false);
	conecta(2,13,1,false);
	conecta(3,4,2,false);
	conecta(3,5,6,false);
	conecta(3,6,10,false);
	conecta(3,7,10,false);
	conecta(3,8,3,false);
	conecta(3,9,8,false);
	conecta(3,10,2,false);
	conecta(3,11,10,false);
	conecta(3,12,7,false);
	conecta(3,13,5,false);
	conecta(4,5,4,false);
	conecta(4,6,2,false);
	conecta(4,7,3,false);
	conecta(4,8,9,false);
	conecta(4,9,10,false);
	conecta(4,10,10,false);
	conecta(4,11,8,false);
	conecta(4,12,10,false);
	conecta(4,13,7,false);
	conecta(5,6,4,false);
	conecta(5,7,2,false);
	conecta(5,8,2,false);
	conecta(5,9,4,false);
	conecta(5,10,10,false);
	conecta(5,11,3,false);
	conecta(5,12,2,false);
	conecta(5,13,7,false);
	conecta(6,7,7,false);
	conecta(6,8,7,false);
	conecta(6,9,6,false);
	conecta(6,10,10,false);
	conecta(6,11,3,false);
	conecta(6,12,8,false);
	conecta(6,13,3,false);
	conecta(7,8,6,false);
	conecta(7,9,7,false);
	conecta(7,10,4,false);
	conecta(7,11,9,false);
	conecta(7,12,6,false);
	conecta(7,13,3,false);
	conecta(8,9,2,false);
	conecta(8,10,5,false);
	conecta(8,11,9,false);
	conecta(8,12,7,false);
	conecta(8,13,10,false);
	conecta(9,10,10,false);
	conecta(9,11,2,false);
	conecta(9,12,9,false);
	conecta(9,13,4,false);
	conecta(10,11,8,false);
	conecta(10,12,3,false);
	conecta(10,13,7,false);
	conecta(11,12,7,false);
	conecta(11,13,9,false);
	conecta(12,13,5,false);
}
