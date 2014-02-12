# treeoptimizer.py
# -------------------------
# Jul 2, 2013; Alex Safatli
# -------------------------
# Using a tree structure, perform
# a breadth-first tree-based optimization
# of the reference amongst a set of pairwise
# alignments in order to acquire the best
# possible multiple structure alignment.
# Works solely on the basis of indices.

# Constants

WARNING_LEVEL = 6

# Auxiliary Structures

class treenode:
    def __init__(self,data):
        self.data = data
        self.branches = []
    def add(self,data):
        node = treenode(data)
        self.branches.append(node)
        return node
    def __iter__(self):
        for b in self.branches: yield b

class tree:
    def __init__(self):
        self.root = treenode(None)
    def setRoot(self,data):
        self.root.data = data
    def __iter__(self):
        for b in self.root: yield b

class alignment_profile:
    ''' Indicates a cluster of structures that possess
    Hamming distance 0 between them in terms of landmark
    presence. '''
    def __init__(self,key,landlist):
        self.ambassador = key
        self.profile = [key]
        self.landlist = landlist
        self.states = [landlist[key]]
    def isInProfile(self,key):
        li = self.landlist[key]
        return (self.states[0] == li)
    def addToProfile(self,key):
        li = self.landlist[key]
        self.profile.append(key)
        self.states.append(li)

class alignment_state:
    def __init__(self,landmarks,gaplist,omit,omitted,n):
        self.score = 0
        self.level = n
        self.landmarks = landmarks
        self.gaplist = gaplist
        self.omit = omit
        self.gaps = 0
        if omit:
            noneInOmit = True
            for key in omit.profile:
                if key in omitted:
                    noneInOmit = False
                    break
            if noneInOmit:
                self.gaps = gaplist[omit.ambassador]
            self.omitted = omitted + omit.profile
        else: self.omitted = omitted
        self.benefit     = self.gaps*n
    def scoreState(self):
        landlist         = self.landmarks
        omit             = self.omitted
        keys             = landlist.keys()
        locations = range(len(landlist[keys[0]]))
        for key in keys:
            if not key in omit:
                vector = landlist[key]
                for l in locations:
                    if vector[l] == 0:
                        locations.remove(l)
        self.score = len(locations)*(len(landlist)-len(omit)+1)
    
class treeoptimizer:
    def __init__(self,landmarks,omitted,n,logf):
        gaplist = {}
        for key in landmarks: gaplist[key] = landmarks[key].count(0)
        root = alignment_state(landmarks,gaplist,None,omitted,n)
        root.scoreState()
        self.tree = tree()
        self.tree.setRoot(root)
        self.gaplist = gaplist
        self.landmarks = landmarks
        self.toplevel = n
        self.logf = logf
        self.warned = 0
        self.optimum = None
        self.omitted = omitted
        self.__defineprofiles__()
    def __defineprofiles__(self):
        self.sortedLandmarks = list(
            reversed(sorted(self.landmarks,key=lambda d: self.landmarks[d].count(0))))
        self.profiles = []
        delegated = []
        for delegate in self.sortedLandmarks:
            if delegate in delegated: continue
            profile = alignment_profile(delegate,self.landmarks)
            delegated.append(delegate)
            for it in self.sortedLandmarks:
                if it in delegated: continue
                elif profile.isInProfile(it):
                    delegated.append(it)
                    profile.addToProfile(it)
            self.profiles.append(profile)
    def __df__(self,vertex):
        '''
        Recursive depth-first addition and analysis of the optimizing
        tree.
        '''
        node = vertex.data
        n    = node.level
        # Check all profiles.
        if n == 2: return
        elif n <= (self.toplevel-WARNING_LEVEL) and not self.warned == n:
            self.warned = n
            self.logf.write('WARNING; Optimizing search hit a high number of nodes (%d).' \
                  % (WARNING_LEVEL))
        for profile in self.profiles:
            # Check profile.
            state = alignment_state(node.landmarks,self.gaplist,profile,node.omitted,
                                    n-len(profile.profile))
            if state.benefit > 0: state.scoreState()
            else: continue
            if state.score > node.score:
                branch = vertex.add(state)
                self.__df__(branch)
            # Cleanup (remove all but highest branch).
            benefits = [x.data.benefit for x in vertex]
            if len(benefits) == 0: continue
            maxBenefit = max(benefits)
            ind = benefits.index(maxBenefit)
            vertex.branches = [vertex.branches[ind]]
    def __getPath__(self,vertex):
        '''
        Recursive search for (only) path of optimizing tree.
        '''
        branches = [x for x in vertex]
        if len(branches) > 1: return None
        elif len(branches) == 0: return vertex
        return self.__getPath__(vertex.branches[-1])
    def optimize(self):
        '''
        Perform an optimization.
        '''
        # Perform recursive depth-first benefit calculation. 
        self.__df__(self.tree.root)
        # Determine the best path.
        opt = self.__getPath__(self.tree.root).data
        self.optimum = opt
    def getNumLandmarks(self):
        if self.optimum: return self.optimum.score
        else: return self.tree.root.data.score
    def getOmitList(self):
        if self.optimum: return self.optimum.omitted
        else: return self.omitted
            
            
        