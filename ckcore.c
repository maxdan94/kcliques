/*
Maximilien Danisch, Oana D. Balalau and Qinna Wang
February 2015
http://bit.ly/maxdan94
maximilien.danisch@telecom-paristech.fr

Info:
Feel free to use these lines as you wish. This program enumerates all k-cliques, computes the k-clique core decomposition and computes a k approximation of a k-cliques densest subgraph.

To compile:
"gcc ckcore.c -O9 -o ckcore".

To execute:
"./ckcore k edgelist.txt ckdeg.txt ckcore.txt ckdens.txt".
k of k-clique to enumerate.
"edgelist.txt" should contain the graph: one edge on each line separated by a space.
Will print the total number of k-cliques.
Will print the k-clique core number of the graph
Will print the k-clique density, the edge density and the size of the densest subgraph
Will write in ckdeg.txt the ID of each node followed by its k-clique degree
Will write in ckcore.txt the ID of each node followed by its k-clique core number in a k-clique core ordering
Will write in ckdens.txt the the k-clique density, the edge density and the size of the densest subgraph followed by the ID of each node in it

Note:
iterating over edges (iterating over nodes for the ckcore computation) and increasing core ordering
*/

#include <stdlib.h>
#include <stdio.h>

#define NLINKS 100000000 //maximum number of edges for memory allocation, will increase if needed

// heap data structure :

typedef struct {
	unsigned key;
	unsigned value;
} keyvalue;

typedef struct {
	unsigned n_max;	// max number of nodes.
	unsigned n;	// number of nodes.
	unsigned *pt;	// pointers to nodes.
	keyvalue *kv; // nodes.
} bheap;

bheap *construct(unsigned n_max){
	unsigned i;
	bheap *heap=malloc(sizeof(bheap));

	heap->n_max=n_max;
	heap->n=0;
	heap->pt=malloc(n_max*sizeof(unsigned));
	for (i=0;i<n_max;i++) heap->pt[i]=-1;
	heap->kv=malloc(n_max*sizeof(keyvalue));
	return heap;
}

inline void swap(bheap *heap,unsigned i, unsigned j) {
	keyvalue kv_tmp=heap->kv[i];
	unsigned pt_tmp=heap->pt[kv_tmp.key];
	heap->pt[heap->kv[i].key]=heap->pt[heap->kv[j].key];
	heap->kv[i]=heap->kv[j];
	heap->pt[heap->kv[j].key]=pt_tmp;
	heap->kv[j]=kv_tmp;
}

inline void bubble_up(bheap *heap,unsigned i) {
	unsigned j=(i-1)/2;
	while (i>0) {
		if (heap->kv[j].value>heap->kv[i].value) {
			swap(heap,i,j);
			i=j;
			j=(i-1)/2;
		}
		else break;
	}
}

inline void bubble_down(bheap *heap) {
	unsigned i=0,j1=1,j2=2,j;
	while (j1<heap->n) {
		j=( (j2<heap->n) && (heap->kv[j2].value<heap->kv[j1].value) ) ? j2 : j1 ;
		if (heap->kv[j].value < heap->kv[i].value) {
			swap(heap,i,j);
			i=j;
			j1=2*i+1;
			j2=j1+1;
			continue;
		}
		break;
	}
}

inline void insert(bheap *heap,keyvalue kv){
	heap->pt[kv.key]=(heap->n)++;
	heap->kv[heap->n-1]=kv;
	bubble_up(heap,heap->n-1);
}

inline void update(bheap *heap,unsigned key){
	unsigned i=heap->pt[key];
	if (i!=-1){
		((heap->kv[i]).value)--;
		bubble_up(heap,i);
	}
}

inline keyvalue popmin(bheap *heap){
	keyvalue min=heap->kv[0];
	heap->pt[min.key]=-1;
	heap->kv[0]=heap->kv[--(heap->n)];
	heap->pt[heap->kv[0].key]=0;
	bubble_down(heap);
	return min;
}

// graph datastructure:

typedef struct {
	unsigned s;
	unsigned t;
} edge;

typedef struct {
	unsigned n; //number of nodes
	unsigned e; //number of edges
	edge *edges;//list of edges

	unsigned *d0; //degrees
	unsigned *cd0; //cumulative degree: (start with 0) length=dim+1
	unsigned *adj0; //list of neighbors

	unsigned *rank; //degeneracy rankings of nodes

	unsigned *d; //truncated degrees
	unsigned *cd; //cumulative degree: (start with 0) length=dim+1
	unsigned *adj; //list of neighbors with higher rank

	unsigned core; //core number of the graph
} sparse;


//compute the maximum of three unsigned
inline unsigned max3(unsigned a,unsigned b,unsigned c){
	a=(a>b) ? a : b;
	return (a>c) ? a : c;
}

//reading the edgelist from file
sparse* readedgelist(char* edgelist){
	unsigned e1=NLINKS;
	sparse *g=malloc(sizeof(sparse));
	FILE *file;

	g->n=0;
	g->e=0;
	file=fopen(edgelist,"r");
	g->edges=malloc(e1*sizeof(edge));
	while (fscanf(file,"%u %u", &(g->edges[g->e].s), &(g->edges[g->e].t))==2) {
		g->n=max3(g->n,g->edges[g->e].s,g->edges[g->e].t);
		if (g->e++==e1) {
			e1+=NLINKS;
			g->edges=realloc(g->edges,e1*sizeof(edge));
		}
	}
	fclose(file);
	g->n++;

	g->edges=realloc(g->edges,g->e*sizeof(edge));

	return g;
}

//used in qsort
int cmpfunc (const void * a, const void * b){
	if (*(unsigned*)a>*(unsigned*)b){
		return 1;
	}
	return -1;
}

//Building the graph structure
void mkgraph(sparse *g){
	unsigned i;
	g->d0=calloc(g->n,sizeof(unsigned));

	for (i=0;i<g->e;i++) {
		g->d0[g->edges[i].s]++;
		g->d0[g->edges[i].t]++;
	}
	g->cd0=malloc((g->n+1)*sizeof(unsigned));
	g->cd0[0]=0;
	for (i=1;i<g->n+1;i++) {
		g->cd0[i]=g->cd0[i-1]+g->d0[i-1];
		g->d0[i-1]=0;
	}

	g->adj0=malloc(2*g->e*sizeof(unsigned));

	for (i=0;i<g->e;i++) {
		g->adj0[ g->cd0[g->edges[i].s] + g->d0[ g->edges[i].s ]++ ]=g->edges[i].t;
		g->adj0[ g->cd0[g->edges[i].t] + g->d0[ g->edges[i].t ]++ ]=g->edges[i].s;
	}

	for (i=0;i<g->n;i++) {
		qsort(&g->adj0[g->cd0[i]],g->d0[i],sizeof(unsigned),cmpfunc);
	}

}

//Building the heap structure with (key,value)=(node,degree) for each node
bheap* mkheap(sparse *g){
	unsigned i;
	keyvalue kv;
	bheap* heap=construct(g->n);
	for (i=0;i<g->n;i++){
		kv.key=i;
		kv.value=g->d0[i];
		insert(heap,kv);
	}
	return heap;
}

void freeheap(bheap *heap){
	free(heap->pt);
	free(heap->kv);
	free(heap);
}

//computing degeneracy ordering and core value
void kcore(sparse* g){
	unsigned i,j;
	keyvalue kv;
	unsigned c=0;//the core number
	bheap *heap=mkheap(g);
	g->rank=malloc(g->n*sizeof(unsigned));

	for (i=0;i<g->n;i++){
		kv=popmin(heap);
		g->rank[kv.key]=i;
		if (kv.value>c){
			c=kv.value;
		}
		for (j=g->cd0[kv.key];j<g->cd0[kv.key+1];j++){
			update(heap,g->adj0[j]);
		}
	}
	freeheap(heap);
	g->core=c;
}

//Add the special feature to the graph structure: a truncated neighborhood contains only nodes with higher rank
void mkspecial(sparse *g){
	unsigned i,j,k;
	g->d=calloc(g->n,sizeof(unsigned));
	g->cd=malloc((g->n+1)*sizeof(unsigned));
	g->adj=malloc((g->e)*sizeof(unsigned));
	g->cd[0]=0;
	for (i=0;i<g->n;i++) {
		g->cd[i+1]=g->cd[i];
		for (j=g->cd0[i];j<g->cd0[i+1];j++){
			k=g->adj0[j];
			if(g->rank[k]>g->rank[i]){
				g->d[i]++;
				g->adj[g->cd[i+1]++]=k;
			}
		}
	}
}

void freesparse(sparse *g){
	free(g->edges);
	free(g->d0);
	free(g->cd0);
	free(g->adj0);
	free(g->rank);
	free(g->d);
	free(g->cd);
	free(g->adj);
	free(g);
}

//store the intersection of list1 and list2 in list3 and return the size of list3 (the 3 lists are sorted)
inline unsigned merging(unsigned *list1, unsigned s1, unsigned *list2, unsigned s2,unsigned *list3){
	unsigned i=0,j=0,s3=0;
	unsigned x=list1[0],y=list2[0];
	while (i<s1 && j<s2){
		if(x<y){
			x=list1[++i];
			continue;
		}
		if(y<x){
			y=list2[++j];
			continue;
		}
		list3[s3++]=x;
		x=list1[++i];
		y=list2[++j];
	}
	return s3;
}

//the recursion to compute all possible intersections
void recursion(unsigned kmax, unsigned k, unsigned* merge, unsigned* size, sparse* g, unsigned* ck, unsigned long long* nck){
	unsigned t=(k-3)*g->core,t2=t+g->core;
	unsigned i, u;

	if (size[k-3]<kmax-k){//stop if we already know k-cliques cannot be formed
		return;
	}

	if (k==kmax){//increasing the k-clique degrees
		for (i=0;i<kmax-1;i++){
			nck[ck[i]]+=size[k-3];
		}
		for (i=0;i<size[k-3];i++){
			nck[merge[t+i]]++;
		}
		return;
	}

	for(i=0; i<size[k-3]; i++){
		ck[k-1]=merge[t+i];
		size[k-2]=merging(&g->adj[g->cd[ck[k-1]]],g->d[ck[k-1]],&merge[t],size[k-3],&merge[t2]);
		recursion(kmax, k+1, merge, size, g, ck, nck);
	}
}

//one pass over all k-cliques
unsigned long long *onepass(sparse *g,unsigned kmax){
	unsigned e,i;
	unsigned *merge,*size,*ck;
	unsigned long long *nck;

	merge=malloc((kmax-2)*g->core*sizeof(unsigned));
	size=malloc((kmax-2)*sizeof(unsigned));
	ck=malloc(kmax*sizeof(unsigned));
	nck=calloc(g->n,sizeof(unsigned long long));

	for(e=0; e<g->e; e++){
		ck[0]=g->edges[e].s;
		ck[1]=g->edges[e].t;
		size[0]=merging(&(g->adj[g->cd[ck[0]]]),g->d[ck[0]],&(g->adj[g->cd[ck[1]]]),g->d[ck[1]],merge);
		recursion(kmax,3,merge,size,g,ck,nck);
	}

	free(merge);
	free(size);
	free(ck);

	return nck;
}

//printing the k-clique degree of each node
unsigned long long printckdeg(unsigned long long *nck, unsigned n, char* ckdeg){
	unsigned i;
	unsigned long long tot=0;
	FILE* file=fopen(ckdeg,"w");
	for (i=0;i<n;i++){
		fprintf(file,"%u %llu\n",i,nck[i]);
		tot+=nck[i];
	}
	fclose(file);
	return tot;
}

// heap data structure :

typedef struct {
	unsigned key;
	unsigned long long value;
} keyvalueLLU;

typedef struct {
	unsigned n_max;// max number of nodes.
	unsigned n;// number of nodes.
	unsigned *pt;// pointers to nodes.
	keyvalueLLU *kv;// (node,nck)
} bheapLLU;

bheapLLU *constructLLU(unsigned n_max){
	unsigned i;
	bheapLLU *heap=malloc(sizeof(bheapLLU));

	heap->n_max=n_max;
	heap->n=0;
	heap->pt=malloc(n_max*sizeof(unsigned));
	for (i=0;i<n_max;i++) heap->pt[i]=-1;
	heap->kv=malloc(n_max*sizeof(keyvalueLLU));
	return heap;
}

inline void swapLLU(bheapLLU *heap,unsigned i, unsigned j) {
	keyvalueLLU kv_tmp=heap->kv[i];
	unsigned pt_tmp=heap->pt[kv_tmp.key];
	heap->pt[heap->kv[i].key]=heap->pt[heap->kv[j].key];
	heap->kv[i]=heap->kv[j];
	heap->pt[heap->kv[j].key]=pt_tmp;
	heap->kv[j]=kv_tmp;
}

inline void bubble_upLLU(bheapLLU *heap,unsigned i) {
	unsigned j=(i-1)/2;
	while (i>0) {
		if (heap->kv[j].value>heap->kv[i].value) {
			swapLLU(heap,i,j);
			i=j;
			j=(i-1)/2;
		}
		else break;
	}
}

inline void bubble_downLLU(bheapLLU *heap) {
	unsigned i=0,j1=1,j2=2,j;
	while (j1<heap->n) {
		j=( (j2<heap->n) && (heap->kv[j2].value<heap->kv[j1].value) ) ? j2 : j1 ;
		if (heap->kv[j].value < heap->kv[i].value) {
			swapLLU(heap,i,j);
			i=j;
			j1=2*i+1;
			j2=j1+1;
			continue;
		}
		break;
	}
}

inline void insertLLU(bheapLLU *heap,keyvalueLLU kv){
	heap->pt[kv.key]=(heap->n)++;
	heap->kv[heap->n-1]=kv;
	bubble_upLLU(heap,heap->n-1);
}

inline void updateLLU(bheapLLU *heap,unsigned key,unsigned long long delta){
	unsigned i=heap->pt[key];
	if (i!=-1){
		((heap->kv[i]).value)-=delta;
		bubble_upLLU(heap,i);
	}
}

inline keyvalueLLU popminLLU(bheapLLU *heap){
	keyvalueLLU min=heap->kv[0];
	heap->pt[min.key]=-1;
	heap->kv[0]=heap->kv[--(heap->n)];
	heap->pt[heap->kv[0].key]=0;
	bubble_downLLU(heap);
	return min;
}

//Building the heap structure with (key,value)=(node,k-clique degree) for each node
bheapLLU* mkheapLLU(unsigned long long* nck,unsigned n){
	unsigned i;
	keyvalueLLU kv;
	bheapLLU* heap=constructLLU(n);
	for (i=0;i<n;i++){
		kv.key=i;
		kv.value=nck[i];
		insertLLU(heap,kv);
	}
	return heap;
}

void freeheapLLU(bheapLLU *heap){
	free(heap->pt);
	free(heap->kv);
	free(heap);
}

typedef struct {
	unsigned n;//number of nodes;
	keyvalueLLU *ic;// (node,k-clique core) in k-clique core ordering
	unsigned long long ckcore;// k-clique core number of the graph
	unsigned size;// size of the k-clique dense subgraph
	double rho;// edge density of the denses subgraph found
	double ckrho;// k-clique density
} densest;

//one pass over all k-cliques
void oneshortpass(sparse *g, bheapLLU* heap, unsigned kmax, unsigned u, unsigned long long *nck){
	unsigned i,v,s=0;
	static unsigned *merge=NULL,*size,*ck,*adj;
	if (merge==NULL){
		merge=malloc((kmax-2)*g->core*sizeof(unsigned));
		size=malloc((kmax-2)*sizeof(unsigned));
		ck=malloc(kmax*sizeof(unsigned));
		adj=malloc(g->n*sizeof(unsigned));
	}

	for (i=g->cd0[u];i<g->cd0[u+1];i++){
		v=g->adj0[i];
		if (heap->pt[v]!=-1){
			adj[s++]=v;
		}
	}

	if (s>kmax-2){
		ck[0]=u;
		for (i=0;i<s;i++){
			ck[1]=adj[i];
			size[0]=merging(adj,s,&(g->adj[g->cd[ck[1]]]),g->d[ck[1]],merge);
			recursion(kmax,3,merge,size,g,ck,nck);
		}
	}

}

densest* kcliquecore(unsigned kmax,unsigned long long *nck, unsigned long long ncktot, sparse* g){
	unsigned i,j,u,v;
	keyvalueLLU kv;
	unsigned long long c=0;//the k-clique core number
	bheapLLU* heap=mkheapLLU(nck,g->n);
	unsigned e=g->e; //number of edges left
	densest* ds=malloc(sizeof(densest));
	ds->ic=calloc(g->n,sizeof(keyvalueLLU));
	ds->n=g->n;
	double ckrho_tmp;
	unsigned size_tmp;
	unsigned long long *nck2=calloc(g->n,sizeof(unsigned long long));
	unsigned long long ncktot2;

	ds->ckrho=((double)ncktot)/((double)g->n);
	ds->size=g->n;
	ds->rho=((double)(2*e))/((double)(g->n*(g->n-1)));
	for (i=0;i<g->n-1;i++){
		kv=popminLLU(heap);
		u=kv.key;
		if (kv.value>c){
			c=kv.value;
		}
		ds->ic[i].key=u;
		ds->ic[i].value=c;
		oneshortpass(g,heap,kmax,u,nck2);
		ncktot2=0;
		for (j=g->cd0[u];j<g->cd0[u+1];j++){
			v=g->adj0[j];
			if (heap->pt[v]!=-1){
				e--;
				updateLLU(heap,v,nck2[v]);
				ncktot2+=nck2[v];
				nck2[v]=0;
			}
		}
		ncktot-=ncktot2/((unsigned long long)(kmax-1));//kcliques are counted kmax-1 times
		size_tmp=g->n-i-1;
		ckrho_tmp=((double)ncktot)/((double)size_tmp);
		if (ckrho_tmp>ds->ckrho){
			ds->ckrho=ckrho_tmp;
			ds->size=size_tmp;
			ds->rho=((double)(2*e))/((double)((size_tmp)*(size_tmp-1)));
		}
	}
	kv=popminLLU(heap);
	ds->ic[i].key=kv.key;
	ds->ic[i].value=c;
	ds->ckcore=c;
	freeheapLLU(heap);
	return ds;
}

void printdensest(densest *ds, char* ckcore, char* ckdens){
	FILE *file;
	unsigned i;

	file=fopen(ckcore,"w");
	for(i=0;i<ds->n;i++){
		fprintf(file,"%u %llu\n",ds->ic[i].key,ds->ic[i].value);
	}
	fclose(file);

	file=fopen(ckdens,"w");
	fprintf(file,"%lf %lf %u",ds->ckrho,ds->rho,ds->size);
	for(i=ds->n-1;i>ds->n-ds->size-1;i--){
		fprintf(file," %u",ds->ic[i].key);
	}
	fprintf(file,"\n");
	fclose(file);
}

int main(int argc,char** argv){
	sparse* g;
	unsigned kmax=atoi(argv[1]);
	unsigned long long *nck;
	unsigned long long tot;
	densest *ds;

	printf("Reading edgelist from file: %s\n",argv[2]);
	g=readedgelist(argv[2]);
	printf("Number of nodes = %u\n",g->n);
	printf("Number of edges = %u\n",g->e);
	printf("Building the graph structure\n");
	mkgraph(g);
	printf("Computing degeneracy ordering\n");
	kcore(g);
	printf("Core number = %u\n",g->core);
	mkspecial(g);
	printf("Computing the %u-clique degree of each node\n",kmax);
	nck=onepass(g,kmax);
	printf("Writing %u-cliques degree of each node in file: %s\n",kmax, argv[3]);
	tot=printckdeg(nck,g->n,argv[3]);
	tot/=(unsigned long long)kmax;
	printf("Number of %u-cliques: %llu\n",kmax,tot);
	ds=kcliquecore(kmax,nck,tot,g);
	printf("%u-clique core number of the graph = %llu\n",kmax,ds->ckcore);
	printf("%u-clique density of the dense subgraph = %lf\n",kmax,ds->ckrho);
	printf("Edge density of the found dense subgraph = %lf\n",ds->rho);
	printf("Size of the found dense subgraph = %u\n",ds->size);
	printf("Writing %u-clique core number of each node in file: %s\n",kmax, argv[4]);
	printf("Writing the found dense subgraph in file: %s\n", argv[5]);
	printdensest(ds,argv[4],argv[5]);
	freesparse(g);
	return 0;
}
