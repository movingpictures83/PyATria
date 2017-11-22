import numpy
import math
import random
import sys
import networkx
############################################################
# Dijkstra's Algorithm
# Uses Modified Classes from the Python Data Structures Module
# One modification: heap is now max-heap
# Since this modifies distances IN PLACE, have to be careful
# Must prioritize the queues based on absolute value
# But must use the actual edge weights when computing distance
# Also: "start" is now an integer i
# I have to build a different graph ON THE FLY
from PageRank.pythonds.graphs import PriorityQueue, Graph, Vertex
numdiff = 0
import itertools

# Optimization: This function should be called *one* time.
# The graph will be built with all nodes duplicated.
# All subsequent graphs will use modifyGraph, which will
# change the starting vertex to have just one node, and
# reduplicate the old starting vertex
def buildATriaGraph(myfile):
    # Build graph with duplicate vertices, initially empty
    aGraph = Graph()
    
    # CSV   
    if (myfile[len(myfile)-3:] == "csv"):
     ###########################################################
     # Read the file
     # Put results in filestuff
     filestuff = open(myfile, 'r')
     firstline = filestuff.readline()
     bacteria = firstline.split(',')
     bacteria.remove('\"\"')
     n = len(bacteria)
     inf = float("infinity")
     ###########################################################

     # Add edges from everywhere to everywhere
     candidates = []
     added = numpy.zeros([n])
     i = 0
     for line in filestuff:
       contents = line.split(',')
       for j in range(n):
          value = float(contents[j+1])
          print value, " "
          if (i != j and value > 0):
              if (not added[i]):
                 aGraph.addVertex(bacteria[i].strip()+"+", i)
                 aGraph.addVertex(bacteria[i].strip()+"-", i)
                 added[i] = True
              if (not added[j]):
                 aGraph.addVertex(bacteria[j].strip()+"+", j)
                 aGraph.addVertex(bacteria[j].strip()+"-", j)
                 added[j] = True
              aGraph.addEdge(bacteria[i].strip()+"+", bacteria[j].strip()+"+", value)
              aGraph.addEdge(bacteria[i].strip()+"-", bacteria[j].strip()+"-", value)
          elif (i != j and value < 0):
              if (not added[i]):
                 aGraph.addVertex(bacteria[i].strip()+"+", i)
                 aGraph.addVertex(bacteria[i].strip()+"-", i)
                 added[i] = True
              if (not added[j]):
                 aGraph.addVertex(bacteria[j].strip()+"+", j)
                 aGraph.addVertex(bacteria[j].strip()+"-", j)
                 added[j] = True
              aGraph.addEdge(bacteria[i].strip()+"+", bacteria[j].strip()+"-", value)
              aGraph.addEdge(bacteria[i].strip()+"-", bacteria[j].strip()+"+", value)
       if (added[i]):
          candidates.append(i)
       i = i + 1
       print ""
     # Return the graph
    else:
      print "Making graph: "
      G = networkx.read_gml(myfile)
      bacteria = G.nodes()
      n = len(bacteria)
      first = True
      weighted = False
      candidates = []
      eps = 1 / (len(G.edges()))
      added = numpy.zeros([n])
      for name in bacteria:
         i = bacteria.index(name)
         for key in G.adj[name]:
            j = bacteria.index(key)
            if (first):
               value = G.get_edge_data(name, key)['weight']
               if (value != 0):
                  weighted = True
               else:
                  value = 0.5# - eps
                  #value = 1. - eps
            elif (weighted):
               value = G.get_edge_data(name, key)['weight']
            # I hate code duplication, but leaving it for now...
            if (i != j and value > 0):
              if (not added[i]):
                 candidates.append(i)
                 aGraph.addVertex(bacteria[i].strip()+"+", i)
                 aGraph.addVertex(bacteria[i].strip()+"-", i)
                 added[i] = True
              if (not added[j]):
                 candidates.append(j)
                 aGraph.addVertex(bacteria[j].strip()+"+", j)
                 aGraph.addVertex(bacteria[j].strip()+"-", j)
                 added[j] = True
              aGraph.addEdge(bacteria[i].strip()+"+", bacteria[j].strip()+"+", value)
              aGraph.addEdge(bacteria[i].strip()+"-", bacteria[j].strip()+"-", value)
            elif (i != j and value < 0):
              if (not added[i]):
                 candidates.append(i)
                 aGraph.addVertex(bacteria[i].strip()+"+", i)
                 aGraph.addVertex(bacteria[i].strip()+"-", i)
                 added[i] = True
              if (not added[j]):
                 candidates.append(j)
                 aGraph.addVertex(bacteria[j].strip()+"+", j)
                 aGraph.addVertex(bacteria[j].strip()+"-", j)
                 added[j] = True
              aGraph.addEdge(bacteria[i].strip()+"+", bacteria[j].strip()+"-", value)
              aGraph.addEdge(bacteria[i].strip()+"-", bacteria[j].strip()+"+", value)
      candidates.sort()  
    print "Done."
    return bacteria, aGraph, candidates


G = None
L = None
u = None
losses = None
gains = None
nolosses = None
nogains = None
numgains = None
numlosses = None
edgedepends = dict()
nopluspath = None
nominuspath = None
def initGL(n):
 global G, L, u, depends, pathdepends, edgedepends, losses, gains, nopluspath, nominuspath, nolosses, nogains, numlosses, numgains
 G = numpy.zeros([n, n])
 L = numpy.zeros([n, n])
 u = numpy.zeros([n])
 losses = []
 gains = []
 numlosses = numpy.zeros([n])
 numgains = numpy.zeros([n])
 nolosses = []
 nogains = []
 #nopluspath = []
 #nominuspath = []
 for i in range(n):
    losses.append([])
    gains.append([])
    nolosses.append(range(i+1, n))
    nogains.append(range(i+1, n))

 #losses = numpy.zeros([n])
 #gains = numpy.zeros([n])
 nopluspath = numpy.zeros([n, n])
 nominuspath = numpy.zeros([n, n])
 print "Initializing edge depends..."
 #  edgedepends.append([])
 #edgedepends = numpy.ndarray([n])
 #  for j in range(n):
 ##     if (j <= i):
 #        edgedepends[i].append(None)
 #     else:
 #        edgedepends[i].append(set())
 print "Done"


# www.peterbe.com/plog/uniqifiers-benchmark
def uniqueify(seq, idfun=None):
   # order preserving
   if idfun is None:
       def idfun(x): return x
   seen = {}
   result = []
   for item in seq:
       marker = idfun(item)
       # in old Python versions:
       # if seen.has_key(marker)
       # but in new ones:
       if marker in seen: continue
       seen[marker] = 1
       result.append(item)
   return result

def addUnique(seq, value):
   if (seq.count(value) == 0):
      seq.append(value)
   return seq

def clearEdgeDependency(i, j):
   edgedepends[i][j] = []

def addEdgeDependency(i, j, value):
   if (not edgedepends.has_key(i)):
      edgedepends[i] = dict()
   if (not edgedepends[i].has_key(j)):
      edgedepends[i][j] = []
   edgedepends[i][j] = addUnique(edgedepends[i][j], value)

def mergeEdgeDependency(i, j, otherdep):
   if (not edgedepends.has_key(i)):
      edgedepends[i] = dict()
   if (not edgedepends[i].has_key(j)):
      edgedepends[i][j] = []
   edgedepends[i][j] = uniqueify(edgedepends[i][j] + otherdep)

def getEdgeDependency(i, j):
   if (not edgedepends.has_key(i)):
      return []
   if (not edgedepends[i].has_key(j)):
      return []
   else:
      return edgedepends[i][j]


def dijkstra(start, bacteria, aGraph, eliminated, first):
    # Number of nodes
    n = len(bacteria)
    # Change starting node
    #reactiv = set()
    nodes = []
    for v in aGraph:
       sign = v.id[len(v.id)-1]
       if (v.index == start):
          pass
       elif (v.index < start ):
          if (sign == '+'):
             value = G[v.index][start]
          else:
             value = L[v.index][start]
          if (value != 0):
             v.setDistance(value)
             #v.setDistance(L[v.index][start])
             v.setPred(None)
          #if (v.getDistance() != 0 and len(v.active) != 0):
             nodes.append((abs(value), v))
       elif (not first and eliminated[start] and eliminated[v.index]):
             pass
             #v.setDistance(0)
             #v.setPred(None)
             #aGraph.deactivate(v.id)
             #reactiv.add(v.id)
       elif (sign == '+' and G[start][v.index] != 0):
             v.setDistance(G[start][v.index])
             v.setPred(None)
             #aGraph.deactivate(v.id)
             #reactiv.add(v.id)
             #if (len(v.active) != 0):
             nodes.append((G[start][v.index], v))
       elif (sign == '-' and L[start][v.index] != 0):
             v.setDistance(L[start][v.index])
             v.setPred(None)
             #aGraph.deactivate(v.id)
             #reactiv.add(v.id)
             #if (len(v.active) != 0):
             nodes.append((-L[start][v.index], v))
       else:
          if (sign == '+' and nopluspath[start][v.index] == 0  or
              sign == '-' and nominuspath[start][v.index] == 0 ):
             nodes.append((0, v))
             v.setDistance(0)
             aGraph.reactivateVertex(v.id)
          #else:
          #   aGraph.deactivate(v.id)
          #   reactiv.add(v.id)
          #   print "Not appending ", v.id, " for: ", bacteria[start].strip(), " ", start
       #nodes.append((abs(v.getDistance()), v))
    startV = aGraph.getVertex(bacteria[start].strip()+"+")
    startV.setDistance(1)
    startVNeg = aGraph.getVertex(bacteria[start].strip()+"-")
    startVNeg.setDistance(2)  # Artificial, will be removed
    nodes.append((1, startV))
    nodes.append((2, startVNeg))
    # Create Priority Queue
    # Underlying max-heap should be populated with absolute values of distances
    # When the algorithm starts, startV's distance is 1 and the rest are 0
    # So, startV will be the first vertex out of the max-heap
    #print "Running Dijkstra"
    pq = PriorityQueue()
    pq.buildHeap(nodes)
    #pq.buildHeap([(abs(v.getDistance()),v) for v in aGraph])
    pq.delMax()  # Negative starting vertex, unused
    

    
    #print "Queue..."
    while not pq.isEmpty():
        currentVert = pq.delMax()  # Get the vertex with the maximum distance
        j = currentVert.index
        dist = currentVert.getDistance()
        #print "Running: ", start, " ", j, " ", dist, " ", G[start][j], " ", L[start][j]
        if (dist == 0):
           #for v in reactiv:
           #   aGraph.reactivateVertex(v)
           return
        sign = currentVert.sign()
        #print "Path..."
        if (start < j and ((G[start][j] == 0 and dist > 0) or (L[start][j] == 0 and dist < 0))):
        #if ((G[start][j] == 0 and dist > 0) or (L[start][j] == 0 and dist < 0)):
              u[start] += dist
              u[j] += dist
              if (sign == '+'):
                 G[start][j] = dist
                 G[j][start] = dist
                 if (first):
                   numgains[start] += 1
                   numgains[j] += 1
                   gains[start].append(j)
              else:
                 L[start][j] = dist
                 L[j][start] = dist
                 if (first):
                   numlosses[start] += 1
                   numlosses[j] += 1
                   losses[start].append(j)

              # Also pay nodes along the path
              curr = currentVert
              prev = currentVert.getPred()
              d = 1
              # prev could be None for a deactivated vertex
              #print "Loop..."
              #tmpedge = set()
              tmpedge = []
              previd = -1
              while (prev != None and prev != startV):
                 previd = prev.index
                 #print "PREVID: ", prev.index, " START: ", start
                 addEdgeDependency(start, j, tuple(sorted((previd, curr.index))))
                 #edgedepends[start][j].add(tuple(sorted((previd, curr.index))))
                 if (start < previd):
                    d *= prev.getWeight(curr)
                    flag = False
                    if (d < 0 and L[previd][j] == 0):
                       char = "-"
                       L[previd][j] = d
		       L[j][previd] = d
                       if (first):
                          numlosses[previd] += 1
                          numlosses[j] += 1
                          if (previd < j):
                             losses[previd].append(j)
                          else:
                             losses[j].append(previd)
                       flag = True
                    elif (d > 0 and G[previd][j] == 0):
                       char = "+"
                       G[previd][j] = d
		       G[j][previd] = d
                       if (first):
                         numgains[previd] += 1
                         numgains[j] += 1
                         if (previd < j):
                             gains[previd].append(j)
                         else:
                             gains[j].append(previd)
                       flag = True
                    
                    if (flag):
                       u[j] += d
                       u[previd] += d
                    tmpedge = addUnique(tmpedge, tuple(sorted((previd, curr.index))))
                    #tmpedge.add(tuple(sorted((previd, curr.index))))
                    if (previd < j):
                       mergeEdgeDependency(previd, j, tmpedge)
		       #edgedepends[previd][j] = edgedepends[previd][j].union(tmpedge)
                    else:
                       mergeEdgeDependency(j, previd, tmpedge)
		       #edgedepends[j][previd] = edgedepends[j][previd].union(tmpedge)
                    newprev = prev.getPred()
                    curr = prev
                    prev = newprev
                 # Break if we find a node less than start, we know everything about it already.
                 else:
                    mergeEdgeDependency(start, j, getEdgeDependency(previd, start))
                    #edgedepends[start][j] = edgedepends[start][j].union(edgedepends[previd][start])
                    break
              # This would be true if we hit a vertex for which we already know its value, and 
              if (previd > start and prev != startV):
                       mergeEdgeDependency(start, j, getEdgeDependency(start, previd))
                       #edgedepends[start][j] = edgedepends[start][j].union(edgedepends[start][previd])
              if (prev == startV):
                       addEdgeDependency(start, j, tuple(sorted((start, curr.index))))
                       #edgedepends[start][j].add(tuple(sorted((start, curr.index))))
             
        # Loop over all of the neighbors of this vertex
        #print "Active: ", currentVert.index
        for nextVert in currentVert.active:
        #for nextVert in currentVert.getConnections():
            # The newly computed distance is the current vertex's distance
            # from the starting vertex, times the edge weight of this neighbor
            newDist = currentVert.getDistance() * currentVert.getWeight(nextVert)
  
	    # We use this edge IF:
            # This new distance is positive, larger than the current distance
            # of this neighbor, and the neighbor is the "+" one, OR: 
            # This new distance is negative, smaller than the current distance
            # of this neighbor, and the neighbor is the "-" one.          
	    if ( (newDist > 0 and newDist > nextVert.getDistance() and nextVert.sign() == '+') or
                 (newDist < 0 and newDist < nextVert.getDistance() and nextVert.sign() == '-') ):
                     # Update the distance for the neighboring vertex
                     nextVert.setDistance ( newDist )
                     # Set its predecessor (Don't use this currently, but may eventually)
                     nextVert.setPred (currentVert)
                     # Set its new value in the table (max-heap)
                     pq.increaseKey(nextVert, abs(newDist))
     
     # End Dijkstra's Algorithm
     ###############################################################################################
    #for v in reactiv:
    #   aGraph.reactivateVertex(v)

####################################################################################################

####################################################################################################
# ATria Centrality Algorithm
def atria_centrality(bacteria, myGraph, candidates):
 global nopluspath, nominuspath, nogains, nolosses

 n = len(bacteria)
 initGL(n)
 # u will be a temporary vector of payoffs
 # U will be the final vector of payoffs
 U = numpy.zeros([n])
 eps = 1. / 10**6  
 # Table header
 print "Rank\tOTU\tPAY\n"

 # Build initial graph for Dijkstra
 
 print "Initializing Graph..."
 #myGraph = buildGraphAndCandidates(bacteria, ADJ, candidates, currentcandidates)
 currentcandidates = []
 for i in range(len(candidates)):
    currentcandidates.append(candidates[i])
 print "Done." 
 maxindex = -1
 # Will run n times, if there are n nodes
 eliminated = numpy.zeros([n]) 
 for i in range(n):
   # Here we determine the one to mark
   maxm = -1
   skip = []
   numcan = len(currentcandidates)
   for j in candidates:
      # Only rerun Dijkstra for unmarked nodes
        # Update the graph, have only one copy of the start node
        # Duplicate the former start node
        if (numcan == 0):
            print "Breaking early"
            break
        #  if (eliminated[j]):
        #     myGraph.deactivate(bacteria[k].strip())
        #if (maxm > max(len(losses[j]), len(gains[j])) and j not in skip):
        #   print "We would skip: ", bacteria[j].strip()
        if (j not in skip):
           #print "Running Dijkstra for node: ", j
           dijkstra(j, bacteria, myGraph, eliminated, (i==0))
           #print "Done.  Updating graph..."
           g = len(gains[j])
           l = len(losses[j])
           for x in range(min(g, l)):
              gain = gains[j][x]
              loss = losses[j][x]
              if (gain in nogains[j]):
                 nogains[j].remove(gain)
              if (loss in nolosses[j]):
                 nolosses[j].remove(loss)
           if (g < l):
              for x in range(g, l):
                 loss = losses[j][x]
                 if (loss in nolosses[j]):
                    nolosses[j].remove(loss)
           else:
              for x in range(l, g):
                 gain = gains[j][x]
                 if (gain in nogains[j]):
                    nogains[j].remove(gain)
  
           nG = len(nogains[j])
           nL = len(nolosses[j])           

           for x in range(min(g, l)):
              gain = gains[j][x]
              loss = losses[j][x]
              for y in range(min(nG, nL)):
                 nogain = nogains[j][y] 
                 noloss = nolosses[j][y]
                 nopluspath[gain][nogain] = True
                 if (nogain in nogains[gain]):
                    nogains[gain].remove(nogain)
                 nopluspath[loss][noloss] = True
                 if (noloss in nogains[loss]):
                    nogains[loss].remove(noloss)
                 nominuspath[gain][noloss] = True
                 if (noloss in nolosses[gain]):
                    nolosses[gain].remove(noloss)
                 nominuspath[loss][nogain] = True
                 if (nogain in nolosses[loss]):
                    nolosses[loss].remove(nogain)
              if (nG < nL):
               for y in range(nG, nL):
                 noloss = nolosses[j][y]
                 nopluspath[loss][noloss] = True
                 if (noloss in nogains[loss]):
                    nogains[loss].remove(noloss)
                 nominuspath[gain][noloss] = True
                 if (noloss in nolosses[gain]):
                    nolosses[gain].remove(noloss)
              else:
               for y in range(nL, nG):
                 nogain = nogains[j][y]
                 nopluspath[gain][nogain] = True
                 if (nogain in nogains[gain]):
                    nogains[gain].remove(nogain)
                 nominuspath[loss][nogain] = True
                 if (nogain in nolosses[loss]):
                    nolosses[loss].remove(nogain)

           if (g < l):
            for x in range(g+1, l):
              loss = losses[j][x]
              for y in range(min(nG, nL)):
                 nogain = nogains[j][y] 
                 noloss = nolosses[j][y]
                 nopluspath[loss][noloss] = True
                 if (noloss in nogains[loss]):
                   nogains[loss].remove(noloss)
                 nominuspath[gain][noloss] = True
                 nominuspath[loss][nogain] = True
                 if (nogain in nolosses[loss]):
                    nolosses[loss].remove(nogain)
              if (nG < nL):
               for y in range(nG+1, nL):
                 noloss = nolosses[j][y]
                 nopluspath[loss][noloss] = True
                 if (noloss in nogains[loss]):
                    nogains[loss].remove(noloss)
              else:
               for y in range(nL+1, nG):
                 nogain = nogains[j][y]
                 nominuspath[loss][nogain] = True
                 if (nogain in nolosses[loss]):
                    nolosses[loss].remove(nogain)
           
           else:
            for x in range(l+1, g):
              gain = gains[j][x]
              for y in range(min(nG, nL)):
                 nogain = nogains[j][y] 
                 noloss = nolosses[j][y]
                 nopluspath[gain][nogain] = True
                 if (nogain in nogains[gain]):
                    nogains[gain].remove(nogain)
                 nominuspath[gain][noloss] = True
                 if (noloss in nolosses[gain]):
                    nolosses[gain].remove(noloss)
              if (nG < nL):
               for y in range(nG, nL):
                 noloss = nolosses[j][y]
                 nominuspath[gain][noloss] = True
                 if (noloss in nolosses[gain]):
                    nolosses[gain].remove(noloss)
              else:
               for y in range(nL, nG):
                 nogain = nogains[j][y]
                 nopluspath[gain][nogain] = True
                 if (nogain in nogains[gain]):
                    nogains[gain].remove(nogain)


           #for x in range(len(gains[j])):
           #   index1 = gains[j][x]
           #   if (index1 in nogains[j]): 
           #         nogains[j].remove(index1)
           #for x in range(len(losses[j])):
           #   index1 = losses[j][x]
           #   if (index1 in nolosses[j]): 
           #         nolosses[j].remove(index1)
           #for index1 in gains[j]:
           #   for index2 in nogains[j]:
           #      nopluspath[index1][index2] = True
           #      if (index2 in nogains[index1]):
           #         nogains[index1].remove(index2)
           #   for index2 in nolosses[j]:
           #      nominuspath[index1][index2] = True
           #      if (index2 in nolosses[index1]):
           #         nolosses[index1].remove(index2)
           #for index1 in losses[j]:
           #   for index2 in nogains[j]:
           #      nominuspath[index1][index2] = True
           #      if (index2 in nolosses[index1]):
           #         nolosses[index1].remove(index2)
           #   for index2 in nolosses[j]:
           #      nopluspath[index1][index2] = True
           #      if (index2 in nogains[index1]):
           #         nogains[index1].remove(index2)
        #print "PAY FOR ", bacteria[j].strip(), ": ", u[j]
        myGraph.deactivateAll()
        #myGraph.deactivate(bacteria[j].strip())
        if (not eliminated[j]):
           currentcandidates.remove(j)
           numcan -= 1
        if (abs(u[j]) > maxm):
             index = j
             maxm = abs(u[j])
             #if (i != 0):
             result = []
             for k in currentcandidates:
                if ((i != 0 and maxm > max(numlosses[k], numgains[k])) or
                    (i == 0 and (n - 1) - j + max(numgains[k], numlosses[k]) < maxm)):
                   #print "Eliminating: ", bacteria[k].strip()
                   eliminated[k] = True
                   numcan -= 1
                else:
                   result.append(k)
             currentcandidates = result
             val = u[j]
   candidates.remove(index)
   # As soon as we find a gain or loss of zero, we can break.  The rest must also be zero.
   if (maxm < eps):
      break

   # The node that is the most central receives its final pay, and we print that node
   # and reset its pay to zero for the next round
   U[index] = val
   maxindex = index
   print "SELECTED: ", bacteria[index].strip(), " WITH PAY: ", U[index], ".  REMOVING TRIADS..."
   # Remove the triads from the one to mark
   # Trying to optimize triad removal TMC
   vertices = [myGraph.vertices[bacteria[index].strip()+"+"], myGraph.vertices[bacteria[index].strip()+"-"]]
   edgesremoved = set()
   for vertex in vertices:
    for vertex2 in vertex.connectedTo:
      toRemove = []
      for vertex3 in vertex2.connectedTo:
         if (vertex3 != vertex): # Dictionary, so will check redundantly unfortunately...
            if (vertex.connectedTo.has_key(vertex3)):
               #numneg = 0
               #e1 = vertex.connectedTo[vertex2]
               #e2 = vertex.connectedTo[vertex3]
               #e3 = vertex2.connectedTo[vertex3]
               #if (e1 < 0):
               #   numneg += 1
               #if (e2 < 0):
               #   numneg += 1
               #if (e3 < 0):
               #   numneg += 1
               #if (numneg % 2 != 1 and abs(e3) <= abs(e2) and abs(e3) <= abs(e1)):
                 toRemove.append(vertex3.id)
                 edgesremoved.add(tuple(sorted((vertex2.index, vertex3.index))))
      for v3id in toRemove:
        myGraph.removeEdge(vertex2.id, v3id)
   # Take out the edges.  Note we had to wait for this, because there could be multiple triads involving
   # the same two edges
   myGraph.removeVertex(bacteria[index].strip()+"+")
   myGraph.removeVertex(bacteria[index].strip()+"-")

   #myGraph.reactivate()
   currentcandidates = []
   for a in range(len(candidates)):
         candidate = candidates[a]
         currentcandidates.append(candidate)
         eliminated[candidate] = False
         u[candidate] -= (G[candidate][maxindex] + L[candidate][maxindex])
         skipflag = True  # If completely independent of maxindex and removed edge, totally skip it
         for b in range(a+1, len(candidates)):
            candidate2 = candidates[b]
            flag = False
            for edge in getEdgeDependency(candidate, candidate2):
               if (edge.count(maxindex) != 0 or edge in edgesremoved):
                  flag = True
                  break
            if (flag):
               skipflag = False
               u[candidate] -= (G[candidate][candidate2] + L[candidate][candidate2])
               u[candidate2] -= (G[candidate][candidate2] + L[candidate][candidate2])
               G[candidate][candidate2] = 0
               L[candidate][candidate2] = 0
               G[candidate2][candidate] = 0
               L[candidate2][candidate] = 0
               clearEdgeDependency(candidate, candidate2)
               #edgedepends[candidate][candidate2].clear()
               #edgedepends[candidate2][candidate].clear()
         if (skipflag):
            skip.append(candidate)

 # Return the pay vector
 return U
####################################################################################################

class PyATriaPlugin:
   def input(self, file):
      self.bacteria, self.graph, self.candidates = buildATriaGraph(file)
   def run(self):
      self.U = atria_centrality(self.bacteria, self.graph, self.candidates)
   def output(self, file):
     UG = []
     for i in range(len(self.U)):
      UG.append((abs(self.U[i]), self.bacteria[i].strip()))
     UG.sort()
     UG.reverse()
     # Formatted for Cytoscape
     outfile = open(file, 'w')
     outfile.write("Name\tCentrality\tRank\n")
     #data = [['Name', 'Centrality']]
     centvals = numpy.zeros([len(UG)])
     for i in range(len(UG)):
       print (UG[i][1], UG[i][0])
       bac = UG[i][1]
       if (bac[0] == '\"'):
          bac = bac[1:len(bac)-1]
       if (UG[i][0] != UG[len(UG)-1][0]):
         outfile.write(bac+"\t"+str(abs(UG[i][0]))+"\t"+str(len(UG)-i)+"\n")
       else:
         outfile.write(bac+"\t"+str(abs(UG[i][0]))+"\t"+"0\n")
       #if (i > 2):
       centvals[i] = abs(UG[i][0])

     print "Wrote file: ", file
     print "Min centrality: ", numpy.min(centvals)
     print "Max centrality: ", numpy.max(centvals)
     mymean = numpy.mean(centvals)
     stddev = numpy.std(centvals)
     print "Standard Deviation: ", stddev
     print "Two STDs back: ", mymean - 2*stddev
     print "One STD back: ", mymean - stddev
     print "One STD forward: ", mymean + stddev
     print "Two STDs forward: ", mymean + 2*stddev


