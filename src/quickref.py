''' This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. 

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

E-mail: asafatli@dal.ca ++

quickref.py
May 6th, 2013; Alex Safatli

Quick reference method index determination. '''

# Imports

import random

class quickref:
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
    