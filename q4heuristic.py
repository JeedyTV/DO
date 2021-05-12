import numpy as np
from math import factorial,exp,floor
from random import randint,uniform,shuffle
import time
from collections import Counter

#define the set of genes
graph = "Hsa/BGSE1456.txt" #G = (V,E)


#define the adjacency matrix as a vector of all the edges 
adj = set()

#define V as the set of all vertices involved in a edge.
V = set()

#define the number of try to create the neighbors.
global tent
tent = 8

#define the cooling of the temperature
global cooling
cooling = 0.99

#define the init temperature
global init_temp
init_temp = 50

#define the final temperature
global end_temp
end_temp = 1

#define the consecutive number of failure accepted before the cut off.
global over
over = 3
   
#read the file and instantiate G = (V,E)

with open(graph,'r') as file:
    
    for i in file:
        st = i.split('\t')
        if st[0] not in V:
            V.add(st[0])
        if st[1] not in V:
            V.add(st[1])
        adj.add((st[0],st[1]))

#convert V into a python set -> acces is O(1)
V = list(V)

def new_permutation(V,m,adj):
    """
    Defined a new permutation of the set V.
    Args:
        V (set): set of the vertices
        m (int): number of the wanted module
        adj (set): set of the edges

    Returns:
        set: new permutation of the set V
    """

    global tent

    perm = V.copy()
    
    """ try to select two vertices to swipe wisely. """
    
    #we select 1 vertex among the m first vertices
    p1 = randint(0,m-1)
    
    #we select 1 vertex among the vertices left
    p2 = randint(m,len(V)-1)

    def comp(p1,p2,adj,perm):
        """
        retrieve the degree of the 2 vertices and
        compare the degree
        Args:
            p1 (int): index of the vertex
            p2 (int): index of the vertex
            adj (set): set of the edges
            perm (int): current permutation of the vertices

        Returns:
            bool: true if degree of p2 is higher than the one of p1 ,false otherwise
        """
        #degree of p1
        f1 = 0
        #degree of p2
        f2 = 0
        
        #compute the degrees
        for i in range(m):
            if (V[p1],V[i]) in adj or (V[i],V[p1]) in adj:
                f1 += 1

        for i in range(m):
            if (V[p2],V[i]) in adj or (V[i],V[p2]) in adj:
                f2 += 1
       
        if f2 > f1:
            return True
        else:
            return False

    def check_prior(p1,p2,adj,perm,tent):
        """
        recursive function which try to swipe the 2 vertices 
        by comparing the degre of the vertexe.

        Args:
            p1 (int): index of the vertex
            p2 (int): index of the vertex
            d (set): set of the edges
            perm (set): new permutation
            tent (int): we fix the swipe process to tent try.

        Returns:
            int: the new neighbor aka permutation
        """
        
        #if the degree of the node p2 is higher or if the try is over we swipe 
        if comp(p1,p2,adj,perm) or tent == 0:
            temp = perm[p1]
            perm[p1] = perm[p2]
            perm[p2] = temp
            return perm
       
        tent -= 1
        
        #select a new vertex to swipe
        p2 = randint(m,len(V)-1)

        return check_prior(p1,p2,adj,perm,tent)
    
    return check_prior(p1,p2,adj,perm,tent)
    
def objective(V,m,adj):
    """
    Compute the nombre of missing edjes to obtain a module.

    Args:
        V (set): set of the vertices 
        m (int): number of the wanted module 
        adj (set): set of the edges 

    Returns:
        int: numbre of missing edges
    """
    #number of edges in G = (V,E)
    link = 0
    for i in range(m):
        for j in range(i,m):
            if i != j and (V[i],V[j]) in  adj or (V[j],V[i]) in  adj:
                    link += 1
    #number of expected edges to have a module
    edges = factorial(m)/(factorial(m-2)*2)
    obj =  edges - link

    return obj

def find_module(V,m,adj,temperature):
    """
    research a module of length m in the graph
    Args:
        V (set): set of the vertices 
        m (int): number of the wanted module 
        adj (set): set of the edges
        temperature (int): temperature as defined in the annealing algorithm

    Returns:
        set: the permutation which fit the condition of being a module
    """

    #compute the goal to approch a module as defined in the assigment
    f1 = objective(V, m, adj)
   
    #find a neighbor
    perm = new_permutation(V, m,adj)
    
    #compute the goal to approch a module as defined in the assigment
    f2 = objective(perm, m, adj)
    
    #compute the delta as describe in the annealing algorithm
    delta = f2-f1

    if delta < 0: #accept the new neihbor if we approch the goal
        
        return perm
     
    #to avoid to be block in a local minimum we accept wrong neighbours
    #with a certain probability
    
    else: 
        #compute the probability of accepting a wrong neighbours
        p = exp(-delta/temperature)
        if uniform(0,100) < p*100:
            return perm
    
    return V

def max_cardinality():
    """
    compute the maximum length of a module in a given graph

    Returns:
        int: the maximum module cardinality
    """
    #create a list containing the number of each vertex involvement.
    array = []
    for i in adj:
        array += [i[0],i[1]]

    #compute the degree by counting the involment
    degree = Counter(array).most_common()

    #retrieve the degree only
    degree_ = [ i[1] for i in degree]

    degree_ = np.array(degree_)
    
    max_m = None
    
    #check if m is valid
    for i in range(degree[0][1]+2)[2:]:
        
        #valid if there are at least m vertex with degree equals to at least m-1 
        if i < len(np.where(degree_>=i-1)[0]):
            max_m = i
        else:
            break
    
    
    return max_m+1
    
def main(V):
    global cooling
    global init_temp
    global end_temp
    global over
    
    #contain final solution
    solution = None

    #retrieve the max possible length
    max_m = max_cardinality()

    #nbre of failed
    nb = 0
    
    #research one module for each m
    while True:   
        nb = 0
        solution = None
        for m in range(max_m+1)[2:]:
            
            if m == len(V)-1:
                solution = V,m
                break
            
            if nb == over:
                break
            
            current_temp = init_temp
            succes = False
            
            vtemp = V.copy()
            shuffle(V) # initial solution

            #stop when it is cold
            while current_temp > end_temp:
                
                if objective(V, m, adj) <= 2 : 
                    succes = True
                    break

                V = find_module(V,m,adj,current_temp)
                
                current_temp *= cooling

            if succes:
                solution = V.copy(),m
                nb = 0
            else:
                nb += 1

            V = vtemp

        if solution != None:
            print(f'for m = {solution[1]} the module is: \n {solution[0][:solution[1]]}')
        
            for i in solution[0][:solution[1]]:
                
                V.remove(i)
            
        if len(V) <= 2:
            print(f'for m = {len(V)} the module is: \n {V}')
            break

if __name__ == '__main__':
    
    main(V)
    