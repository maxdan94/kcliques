/*
Maximilien Danisch and Qinna Wang
January 2015
http://bit.ly/maxdan94
maximilien.danisch@telecom-paristech.fr

Info:
Feel free to use these lines as you wish. This program computes the k-clique degrees of each node in the graph (that is the number of k-cliques the node belongs to).

To compile:
"gcc ckdeg.c -O9 -o ckdeg -fopenmp".

To execute:
"./ckdeg p k edgelist.txt kdeg.txt".
p is the number of processors to use.
k of k-clique to enumerate.
"edgelist.txt" should contain the graph: one edge on each line separated by a space.
Will print the total number of k-cliques.
"kdeg.txt" will contain the k-clique degrees (node ID followed and its k-clique degree separated by a space on each line).

Note:
parallelisation over edges and increasing core ordering
*/

#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <time.h>

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
	//edge list structure:
	unsigned n; //number of nodes
	unsigned e; //number of edges
	unsigned n2; //number of nodes with core value larger than one
	unsigned e2; //number of edges between nodes with core value larger than one
	edge *edges;//list of edges

	//to compute a degeneracy ordering:
	unsigned *d0; //degrees
	unsigned *cd0; //cumulative degree: (start with 0) length=dim+1
	unsigned *adj0; //list of neighbors
	unsigned *rank; //degeneracy rankings of nodes
	unsigned *map;//map[newlabel]=oldlabel
	unsigned core; //core number of the graph

	//truncated neighborhoods:
	unsigned *d; //truncated degrees
	unsigned *cd; //cumulative degree: (start with 0) length=dim+1
	unsigned *adj; //list of neighbors with higher rank
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
void kcore(sparse* g,unsigned kmax){
	unsigned i,j,r=0,n=g->n,k=kmax-1;
	keyvalue kv;
	unsigned c=0;//the core number
	bheap *heap=mkheap(g);
	g->rank=malloc(g->n*sizeof(unsigned));
	g->map=malloc(g->n*sizeof(unsigned));
	for (i=0;i<g->n;i++){
		kv=popmin(heap);
		if (kv.value>c){
			c=kv.value;
		}
		if (c<k){//remove node with core value less than kmax-1
			g->rank[kv.key]=-1;
			n--;
		}
		else{
			g->map[n-(++r)]=kv.key;
			g->rank[kv.key]=n-r;
		}
		for (j=g->cd0[kv.key];j<g->cd0[kv.key+1];j++){
			update(heap,g->adj0[j]);
		}
	}
	freeheap(heap);
	free(g->d0);
	free(g->cd0);
	free(g->adj0);
	g->core=c;
	g->n2=n;
}

void relabelnodes(sparse *g) {
	unsigned i,j,source,target;
	j=0;
	for (i=0;i<g->e;i++) {
		source=g->rank[g->edges[i].s];
		target=g->rank[g->edges[i].t];
		if (source==-1 || target==-1){
			continue;
		}
		if (source<target) {
			g->edges[j].s=target;
			g->edges[j++].t=source;
		}
		else {
			g->edges[j].s=source;
			g->edges[j++].t=target;
		}
	}
	g->e2=j;
	g->edges=realloc(g->edges,g->e2*sizeof(edge));
}

//for future use in qsort
int cmpfunc (const void * a, const void * b){
	if (*(unsigned*)a>*(unsigned*)b){
		return 1;
	}
	return -1;
}

//Building the special graph structure
void mkspecial(sparse *g){
	unsigned i;
	g->d=calloc(g->n2,sizeof(unsigned));

	for (i=0;i<g->e2;i++) {
		g->d[g->edges[i].s]++;
	}
	g->cd=malloc((g->n2+1)*sizeof(unsigned));
	g->cd[0]=0;
	for (i=1;i<g->n2+1;i++) {
		g->cd[i]=g->cd[i-1]+g->d[i-1];
		g->d[i-1]=0;
	}

	g->adj=malloc((g->e2)*sizeof(unsigned));

	for (i=0;i<g->e2;i++) {
		g->adj[g->cd[g->edges[i].s] + g->d[ g->edges[i].s ]++ ]=g->edges[i].t;
	}

	#pragma omp parallel for private(i)
	for (i=0;i<g->n2;i++) {
		qsort(&g->adj[g->cd[i]],g->d[i],sizeof(unsigned),cmpfunc);
	}

	//free(g->edges); Can be freed if node parallelisation is used instead of edge
}

void freesparse(sparse *g){
	free(g->edges);
	free(g->rank);
	free(g->map);
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
void recursion(unsigned kmax, unsigned k, sparse* g, unsigned* ck, unsigned* merge, unsigned* size, unsigned long long* ckdeg){
	unsigned t=(k-3)*g->core,t2=t+g->core;
	unsigned i, j, u;

	if (size[k-3]<kmax-k){//stop if we already know k-cliques cannot be formed
		return;
	}

	if (k==kmax){//increasing the k-clique degrees
		for (i=0;i<kmax-1;i++){
			ckdeg[ck[i]]+=size[k-3];
		}
		for (i=0;i<size[k-3];i++){
			ckdeg[merge[t+i]]++;
		}
		return;
	}

	for(i=0; i<size[k-3]; i++){
		ck[k-1]=merge[t+i];
		size[k-2]=merging(&g->adj[g->cd[ck[k-1]]],g->d[ck[k-1]],&merge[t],size[k-3],&merge[t2]);
		recursion(kmax, k+1, g, ck, merge, size, ckdeg);
	}
}

//one pass over all k-cliques
unsigned long long *onepass(sparse *g,unsigned kmax){
	unsigned e,i,u,v,k;
	unsigned *merge,*size,*ck;
	unsigned long long *ckdeg_p,*ckdeg;
	if (kmax>2){

		ckdeg=calloc(g->n2,sizeof(unsigned long long));
		#pragma omp parallel private(merge,size,ck,ckdeg_p,e,u,v,k) shared(ckdeg)
		{
			merge=malloc((kmax-2)*g->core*sizeof(unsigned));
			size=malloc((kmax-2)*sizeof(unsigned));
			ck=malloc(kmax*sizeof(unsigned));
			ckdeg_p=calloc(g->n2,sizeof(unsigned long long));

			#pragma omp for schedule(dynamic, 1) nowait
			for(e=0; e<g->e2; e++){
				ck[0]=g->edges[e].s;
				ck[1]=g->edges[e].t;
				size[0]=merging(&(g->adj[g->cd[ck[0]]]),g->d[ck[0]],&(g->adj[g->cd[ck[1]]]),g->d[ck[1]],merge);
				recursion(kmax,3,g,ck,merge,size,ckdeg_p);
			}
			#pragma omp critical
			{
			for (i=0;i<g->n2;i++){
				ckdeg[i]+=ckdeg_p[i];
			}
			}
			free(ckdeg_p);
			free(merge);
			free(size);
		}
	}

	return ckdeg;
}

unsigned long long printckdeg(unsigned long long *ckdeg,sparse *g, char* output){
	unsigned i;
	unsigned long long sum=0;
	FILE* file=fopen(output,"w");
	for (i=0;i<g->n2;i++){
		fprintf(file,"%u %llu\n",g->map[i],ckdeg[i]);
		sum+=ckdeg[i];
	}
	fclose(file);
	return sum;
}

int main(int argc,char** argv){
	sparse* g;
	unsigned i,	kmax=atoi(argv[2]);
	unsigned long long *ckdeg,nck;
	omp_set_num_threads(atoi(argv[1]));

	time_t t0,t1,t2;
	t1=time(NULL);
	t0=t1;

	printf("Reading edgelist from file %s\n",argv[3]);
	g=readedgelist(argv[3]);

	t2=time(NULL);
	printf("- Time = %ldh%ldm%lds\n",(t2-t1)/3600,((t2-t1)%3600)/60,((t2-t1)%60));
	t1=t2;

	printf("Number of nodes: %u\n",g->n);
	printf("Number of edges: %u\n",g->e);
	printf("Building the graph structure\n");
	mkgraph(g);
	printf("Computing degeneracy ordering\n");
	kcore(g,kmax);
	relabelnodes(g);
	printf("Number of nodes (with core value > %u): %u\n",kmax-2,g->n2);
	printf("Number of edges (between nodes with core value > %u): %u\n",kmax-2,g->e2);
	printf("Core number = %u\n",g->core);
	mkspecial(g);

	t2=time(NULL);
	printf("- Time = %ldh%ldm%lds\n",(t2-t1)/3600,((t2-t1)%3600)/60,((t2-t1)%60));
	t1=t2;

	printf("computing %u-clique degrees\n",kmax);
	ckdeg=onepass(g,kmax);

	nck=printckdeg(ckdeg,g,argv[4]);
	nck/=kmax;
	printf("Number of %u-cliques: %llu\n",kmax,nck);
	free(ckdeg);
	freesparse(g);
	t2=time(NULL);
	printf("- Time = %ldh%ldm%lds\n",(t2-t1)/3600,((t2-t1)%3600)/60,((t2-t1)%60));
	t1=t2;

	printf("- Overall time = %ldh%ldm%lds\n",(t2-t0)/3600,((t2-t0)%3600)/60,((t2-t0)%60));

	return 0;
}
