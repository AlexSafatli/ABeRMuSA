''' Quick reference method index determination. '''

import random

class quickref:
    
    ''' Randomize what reference to choose from given possibilities. '''
    
    def __init__(self,filelist,numiters=2,random=True):
        self.filelist = filelist
        self.indices = range(len(filelist))
        self.numiters = numiters
        self.random = random
    
    def get(self):
        indices, numiters = self.indices, self.numiters
        if numiters > len(indices) or numiters < 1: numiters = len(indices)        
        if self.random:
            random.shuffle(indices)
            picked = []
            for i in xrange(numiters):
                picked_row = indices.pop(0)
                picked.append(self.filelist[picked_row])
            return picked
        else: return []
    