import numpy as np
from scipy import linalg
import matplotlib.pyplot as plt
import random

# FINDS COMMUNITIES OR MODULES IN DIRECTED NETWORKS AND EVALUATES THEIR MODULARITY INDEX. THE ALGORITHM WHICH OPTIMIZES A MODULARITY FUNCTION AND MAKES CONTINUOUS BISECTIONS UNTIL NO FURTHER IMPROVEMENT ON THE MODULARITY FUNCTION IS POSSIBLE, IS FROM THE FOLLOWING PAPER:
# Leicht, Elizabeth A., and Mark EJ Newman. "Community structure in directed networks." Physical review letters 100.11 (2008): 118703.

# I HAVE IMPLEMENTED THE UNDIRECTED VERSION (https://github.com/rentzi/netRewireAnalyze ). THE CORRESPONDING PAPER IS IN
# Newman, M. E. (2006). Modularity and community structure in networks. Proceedings of the national academy of sciences, 103(23), 8577-8582.

#Note that connections are represented in nXn adjacency matrices A, where Aij is the connection from j to i, i.e. j->i

# MAKEMODULARITYMATRIX takes the adjacency matrix and returns the modularity matrix (Bij = Aij - (din(i)*dout(j))/(m)) that will be used for finding communities within the graph
# INPUT:
# A: the adjacency matrix for the directed graph: Aij -> the connection from j to i, j -> i
# OUTPUT:
# B: the modularity matrix
# m: number of in-degrees (should be the same as the number of outdegrees)
def makeModularityMatrix(A):

    degOut = np.sum(A, axis=0)
    degIn = np.sum(A, axis=1)

    vertices = degIn.size
    m = np.sum(degIn)

    expConnMatrix = np.zeros((vertices, vertices))
    for i in range(vertices):
        for j in range(vertices):
            expConnMatrix[i, j] = (degIn[i] * degOut[j]) / m

    B = A - expConnMatrix

    return B, m


# MAKEMODULARITYMATRIXFROMPARTITION estimates the modularity matrix of a partition or subgraph. It is equation [6] from (Newman,PNAS, 2006). Used for the iterative partitioning
# INPUT:
# B: modularity matrix of the full graph
# partitionInd: indeces of the nodes of the subgraph
# OUTPUT:
# Bpart: the modularity matrix of the partition. Similarly to the modularity matrix of the full graph, the sum of its rows or columns sum to zero
def makeModularityMatrixFromPartition(B, partitionInd):

    Btemp = B.copy()

    # kind of like an outer product to get the right indices, it does not work like matlab (dahhh)
    Bpart = Btemp[np.ix_(partitionInd, partitionInd)]

    for i in np.arange(Bpart.shape[0]):
        Bpart[i, i] -= (0.5 * np.sum(Bpart[i, :]) + 0.5 * np.sum(Bpart[i, :]))

    return Bpart


# DIV2COM divides the graph or subgraph into two communities and gives the modularity index.
# will be used in PREORDERPARTIONING
# INPUT:
# B: the modularity matrix
# OUTPUT:
# s: vector which can take two values (-1,+1) depending on the community that a node is. For example s[j] = -1 indicates that node j is in community 1, s[k] = 1 that node k is in community 2.
# Q: modularity index; up to a multiplicative constant, the number of edges falling within groups minus the expected number in an equivalent network with edges placed at random

#B, totalConnections = makeModularityMatrix(A)
# or
#Bpart = makeModularityMatrixFromPartition(B, partitionInd)
def div2Com(B, totalConnections):

    BSym = B + np.transpose(B)

    # ascending order of eigenvalues
    lambdasAsc, vAsc = np.linalg.eigh(BSym)
    # reverse them: start from largest going to smallest
    lambdas = lambdasAsc[::-1]
    v = np.flip(vAsc, 1)

    # pick the eigenvector corresponding to the largest eigenvalue
    # v1 = v[:, 0:1]  # it has rank 2 by doing 0:1 instead of just 0
    v1 = v[:, 0]  # this has rank 1 (n,)

    # find the positive and negative elements of the eigenvector
    indPos = np.where(v1 > 0)[0]
    indNeg = np.where(v1 <= 0)[0]

    # makes a dictionary with the two partitions,each containing the indices of the nodes belonging to them
    partitionsInd = {}
    partitionsInd[-1] = indNeg
    partitionsInd[1] = indPos

    # when positive make element of s +1, when negative make it -1
    s = np.zeros((v1.size, 1))

    s[indPos] = 1
    s[indNeg] = -1

    # calculating the modularity index
    Qtemp = (s.T @ BSym @ s) / (4 * totalConnections)
    Q = float(np.squeeze(Qtemp))

    return partitionsInd, v1, Q


# The structure from which you make repeated bisection using Newman,PNAS, 2006 paper.
# community is a class representing a cluster of nodes. The variables left and right,
# if the cluster of nodes can be partitioned, point to the partitioned clusters, otherwise
# they do not point anywhere. communityInd is the indices of the nodes in the partition
# Q is the modularity index
class community:
    def __init__(self, communityInd=None):
        self.left = None  # left child
        self.right = None  # right child
        self.communityInd = communityInd
        self.Q = None  # this scalar is positive if the nodes in communityInd can be further split. Add them


# the partitionBinaryTree helps build the tree. It initiates its root (self.root), the indices of all the nodes.
# It also gets the modularity matrix B of the original adjacency matrix, and the total connections of the adjacency matrix. These variables will be used throughout.
class partitionBinaryTree:
    def __init__(self, B, totalConnections):
        communityInd = np.arange(B.shape[0])
        self.root = community(communityInd)  # the indices of the root are 0,1,2....numofVertices
        self.B = B
        self.totalConnections = totalConnections

    # PREORDERPARTIONING iteratively partitions a network of nodes using Newman's method
    # INPUT:
    # startNode: the starting point (partitionBinaryTree.root) of the network containing all the indices of the nodes
    # OUTPUT:
    # Qlist: The list of all the modularity indices from the partitions. Add them all together to get the total modularity index
    # communitiesDict: dictionary with the indices of each community. The keys are positive integers starting from 1 to the number of communities
    def preorderPartitioning(self, startNode, Qlist=[], communitiesDict={}):
        # Root ->Left->Right
        if startNode is not None:
            partB = makeModularityMatrixFromPartition(self.B, startNode.communityInd)
            communitiesInd, v1, startNode.Q = div2Com(partB, self.totalConnections)
            #print('Q is ' ,startNode.Q)
            #print('Community1 size is %d, and community 2 size is %d'%(communitiesInd[-1].size,communitiesInd[1].size))
            if startNode.Q > 0 and communitiesInd[-1].size > 0 and communitiesInd[1].size > 0:
                Qlist.append(startNode.Q)
                startNode.left = community(startNode.communityInd[communitiesInd[-1]])
                startNode.right = community(startNode.communityInd[communitiesInd[1]])
                self.preorderPartitioning(startNode.left, Qlist, communitiesDict)
                self.preorderPartitioning(startNode.right, Qlist, communitiesDict)
            else:
                if not communitiesDict:  # if the dictionary is empty, first timer
                    communitiesDict[1] = startNode.communityInd
                else:
                    maxKey = np.max(list(communitiesDict.keys()))
                    newKey = maxKey + 1
                    communitiesDict[newKey] = startNode.communityInd

        return Qlist, communitiesDict


# GETMODULARITYINDEX uses the functions above to get the modularity index value Q of the adjacency matrix A
# INPUT:
# A: the adjacency matrix
# OUTPUT:
# Q: the modularity index
# communitiesDict: dictionary with the indices of each community. The keys are positive integers starting from 1 to the number of communities
def getModularityIndex(A):

    B, totalConnections = makeModularityMatrix(A)
    graph = partitionBinaryTree(B, totalConnections)
    Qlist, communitiesDict = graph.preorderPartitioning(graph.root, Qlist=[], communitiesDict={})
    #print('for probability = %f and time = %f '%(p,t))
    # print(Qlist)
    Q = np.sum(Qlist)

    return Q, communitiesDict
