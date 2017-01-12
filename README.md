# README #

##Info##

ck.c counts all k-cliques.  
cklist.c list all k-cliques and write them in a file (note that as the number of k-cliques can be extremly large and thus it can use a lot of disk storage).  
ckcore.c computes a k-clique core decomosition of the graph. ckcorepar.c is a parallel version of the same code.

Both programs scale to real-world networks containing several billions of edges.  
For instance:  
- Friendster (http://snap.stanford.edu/data/com-Friendster.html) has exactly 487,090,833,092,739 10-cliques, that is a bit less than 0.5 quadrillions 10-cliques. 
- Twitter (https://an.kaist.ac.kr/traces/WWW2010.html) has exactly 3,388,795,307,518,264 5-cliques, that is a bit more than 3 quadrillions 5-cliques.

##To compile##

gcc ck.c -O3 -o ck -fopenmp  
gcc cklist.c -O3 -o cklist  
gcc ckcore.c -O3 -o ckcore  
gcc ckcorepar.c -O3 -o ckcorepar -fopenmp

##To execute##

./ck p k edgelist.txt

- "p" is the number of threads to use.
- "k" of k-clique to enumerate.
- "edgelist.txt" should contain the graph: one edge on each line separated by a space.
- Will print the number of l-cliques for l in [1,k].

./cklist k edgelist.txt kcliques.txt

- "k" of k-clique to enumerate.
- "edgelist.txt" should contain the graph: one edge on each line separated by a space.
- "kcliques.txt" will contain the kcliques (one kclique on each line that is k unsigned ints on each line).

./ckcore k edgelist.txt ckdeg.txt ckcore.txt ckdens.txt  
./ckcorepar p k edgelist.txt ckdeg.txt ckcore.txt ckdens.txt

- p is the number of threads to use.
- k of k-clique to enumerate.
- "edgelist.txt" should contain the graph: one edge on each line separated by a space.
- Will print the total number of k-cliques.
- Will print the k-clique core number of the graph
- Will print the k-clique density, the edge density and the size of the densest subgraph
- Will write in ckdeg.txt the ID of each node followed by its k-clique degree
- Will write in ckcore.txt the ID of each node followed by its k-clique core number in a k-clique core ordering
- Will write in ckdens.txt the the k-clique density, the edge density and the size of the densest subgraph followed by the ID of each node in it

##Initial contributors##

Maximilien Danisch, Qinna Wang, Oana Balalau and Mauro Sozio  
January 2016  
http://bit.ly/maxdan94  
maximilien.danisch@telecom-paristech.fr
