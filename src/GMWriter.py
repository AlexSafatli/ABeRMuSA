''' This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. 

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

E-mail: asafatli@dal.ca ++

GMWriter.py
May 7th, 2013; Alex Safatli; Based on original code by Kyle Nguyen

Take a whole set of pairwise alignments and creates a final GM and landmarks file. '''

# Imports

import glob, os, operator, cPickle, re, math
from scipy import mean
from utils import PDBnet, IO, FASTAnet
from utils.PDBnet import aa_names
from utils.logfile import logfile, XMLfile
from treeoptimizer import treeoptimizer
from scoring import scoreFile

# Internals

class gmdata:
    def __init__(self,fastas,bad,alnlen,ldata,best,carb,outpre,log):
        self.fastas = fastas
        self.bad    = bad
        self.alnlen = alnlen
        self.ldata  = ldata
        self.best   = best
        self.carbon = carb
        self.logf   = log
        self.outpre = outpre
    def write(self):
        ldata  = self.ldata
        fastas = self.fastas
        bad    = self.bad   
        outfiles = []
        # Get length of a chain label to write in GM and landmarks file
        lenlabel = 1
        while (26**lenlabel < (len(fastas) - len(bad) + 1)): lenlabel += 1

        # Get a list of valid alignment positions.
        locations = [x for x in range(self.alnlen)]
        for chain in ldata:
            if chain in bad: continue
            for l in range(len(ldata[chain])):
                if ldata[chain][l] == 0 and l in locations: locations.remove(l)  
        valid = locations
        if len(valid) == 0: return None

        # Create final gm and landmarks files
        fgm = open(self.outpre + '.gm', 'w')
        outfiles.append(self.outpre + '.gm')
        fland = open(self.outpre + '.landmarks', 'w')
        outfiles.append(self.outpre + '.landmarks')
        refwritten = False # Check if the reference structure is written into GM.        

        index = 0 # index to keep track of chain name

        # Loop over fastas.

        for name in fastas:

            # Get relevant data.
            f = IO.getFileName(name)
            seqs = FASTAnet.FASTAstructure(name,uniqueOnly=False)
            if f in bad or len(seqs.orderedSequences) != 2: continue
            m = re.search('(?<=-).*',f)
            f = m.group(0) # Remove aligner extension.
            f = f[len(self.best)+1:] # Remove best reference name.
            teststruc = seqs.orderedSequences[0].sequence
            refstruc  = seqs.orderedSequences[1].sequence            
            pdb = name[:name.rfind('.')].strip('.') + '.pdb'
            pdbdata = PDBnet.PDBstructure(pdb)

            # read from pdb and write to gm and landmarks here
            fgm.write('>' + f + ':' + f[:4] + ';')
            fland.write('>' + f + '\n')
            fastaindex = 0 # char index in test structure's fasta chain
            landindex = 0 # index in landmarks file
            refindex = -1 # char index in reference structure's fasta chain
            frstch = pdbdata.orderofchains[0]
            lstch = pdbdata.orderofchains[-1]
            if (lstch == frstch):
                self.logf.write('ERROR; First and last chain names in %s identical.'%(f))
                fgm.close()
                fland.close()
                return None            
            for i in range(len(refstruc)):
                if teststruc[i] != '-': fastaindex += 1
                if refstruc[i] != '-':
                    refindex += 1
                    if refindex in valid:
                        fland.write(str(landindex) + '\t')
                        landindex += 1
                        if len(pdbdata.chains[frstch]) <= fastaindex-1:
                            self.logf.write('ERROR; PDB data incomplete for %s. Stopping.' % (f))
                            fgm.close()
                            fland.close()
                            return None
                        resid = pdbdata.chains[frstch][pdbdata.chains[frstch].GetIndices()[fastaindex-1]]
                        # Check for sequence mismatch
                        if teststruc[i] == '-' or resid.name != aa_names[teststruc[i]]:
                            self.logf.write('ERROR; Sequence mismatch in %s at %d {%s,%s}. Stopping.' 
                                            % (f,i,resid.name,teststruc[i]))
                            fgm.close()
                            fland.close()                            
                            return None
                        fland.write('%s\t%s\n' % (resid.index,resid.name))
                        if self.carbon=='CA':
                            alpha = resid.GetCA()
                            x,y,z = (alpha.x,alpha.y,alpha.z)
                        else:
                            resid.Centroid()
                            centr = resid.centroid
                            x,y,z = centr.x,centr.y,centr.z
                        fgm.write('%f;%f;%f;' % (x,y,z))

            fgm.write('\n')
            index += 1
            if not refwritten:
                fastaindex = 0
                landindex = 0
                fgm.write('>' + self.best + ':' + self.best[:4] + ';')
                fland.write('>' + self.best + '\n')
                for i in range(len(refstruc)):
                    if refstruc[i] != '-':
                        fastaindex += 1
                        if (fastaindex-1) in valid:
                            if len(pdbdata.chains[lstch]) <= fastaindex-1:
                                self.logf.write('ERROR; PDB data incomplete for %s. Stopping.' % (f))
                                fgm.close()
                                fland.close()
                                return None                            
                            resid = pdbdata.chains[lstch][pdbdata.chains[lstch].GetIndices()[fastaindex-1]]
                            # Check for sequence mismatch
                            if resid.name != aa_names[refstruc[i]]:
                                self.logf.write('ERROR; Sequence mismatch in %s at %d {%s,%s}. Stopping.' 
                                                % (f,i,resid.name,refstruc[i]))
                                fgm.close()
                                fland.close()                                
                                return None	
                            fland.write('%s\t%s\t%s\n' % (landindex,resid.index,resid.name))
                            landindex += 1
                            if self.carbon=='CA':
                                alpha = resid.GetCA()
                                x,y,z = (alpha.x,alpha.y,alpha.z)
                            else:
                                resid.Centroid()
                                centr = resid.centroid
                                x,y,z = centr.x,centr.y,centr.z
                            fgm.write('%f;%f;%f;' % (x,y,z))
                refwritten = True
                index += 1
                fgm.write('\n')

        fgm.close()
        fland.close()      
        return outfiles

class gmwriter(object):

    def __init__(self,args,prefix,bestref,reffldr,logf,exe,alphaC=False,
                 optimize=False,touchup=False):
        self.args    = args
        self.prefix  = prefix
        self.outpre  = prefix + '_final'
        self.optpre  = prefix + '_optimized'
        self.bestref = bestref
        self.reffldr = reffldr
        self.logf    = logf
        self.exe     = exe
        self.alphaC  = alphaC
        self.optim   = optimize
        self.touchup = touchup

    ''' Read a GM and return a dictionary of landmarks for all structures/chains. '''

    def GetLandmarks(self):
        landmarks = {}
        f = open(self.outpre + '.gm')
        for line in f:
            name = line[line.find(':')+1:line.find(';')]
            if line.strip()[-1] == ';': temp = line.strip().split(';')[1:-1]
            else: temp = line.strip().split(';')[1:]
            landmarks[name] = [[], [], []]
            for p in range(len(temp)): landmarks[name][p%3].append(float(temp[p]))
        return landmarks

    ''' Write the GM file. '''

    def write(self):       

        # Setup alpha carbon detection.
        if self.alphaC: carbon = 'CA'
        else: carbon = 'C'

        # Alert about alpha carbons.
        if self.alphaC: self.logf.write('NOTE; GM will contain Alpha Carbon coordinates.')        

        # Prepare a list of all files created in this process.
        outfiles = []

        # Get the reference structure filename.
        best = self.reffldr
        self.logf.write('Reference folder <%s>' % (self.reffldr))
        if self.bestref: self.logf.write('Structure file <%s>' % (self.bestref[0]))        

        # Get all fastas found in the reference folder (legacy support; uses
        # arg list to populate this list).

        fastas = []
        if self.exe: f_e = self.exe.plugin.fasta_ext            
        else:        f_e = 'fasta'
        if self.args:      

            for f in self.args:
                key = IO.getFileName(f)
                if key == best: continue
                cmd = self.exe.plugin.default_exe
                fasta = '%s/%s-%s.%s.%s' % (best,cmd,best,key,f_e)
                if os.path.isfile(fasta): fastas.append(fasta)
                else: self.logf.write('WARNING; Missing file %s. ' % (fasta) +\
                                      'Will not be included in GM file.')

        else:

            # If no arg list defined, use the XML manifest found in ref folder.
            manifest = '%s/manifest.xml' % (best)
            if not os.path.isfile(manifest):
                self.logf.write('WARNING; No manifest file %s ' % (manifest) +\
                                'to populate file list for GM writing.')
                self.logf.write('Resorting to legacy support of reference pickle file.'+\
                                ' Failed files may not be catalogued.')
                o = open('%s/ref.pickl' % (best))
                scores,_,_,_,_ = cPickle.load(o)
                o.close()
                for item in scores:
                    fasta = '%s/%s.%s' % (best,item,f_e)
                    if os.path.isfile(fasta): fastas.append(fasta)
                    else: self.logf.write('WARNING; Missing file %s. ' % (fasta) +\
                                          'Will not be included in GM file.')
                
            else:
                xml = XMLfile(manifest)
                xml.read() # Read the XML file.
                if xml.root.tag != 'reference':
                    self.logf.write('ERROR; Invalid manifest file %s.' % (manifest))
                    return None
                for child in xml.root:
                    if child.tag == 'succeeded':
                        fasta = '%s/%s.%s' % (best,child.text,f_e)
                        if os.path.isfile(fasta): fastas.append(fasta)
                        else: self.logf.write('WARNING; Missing file %s. ' % (fasta) +\
                                              'Will not be included in GM file.')
                    else: self.logf.write('NOTE; Failed pairwise %s disregarded in GM file.' \
                                          % (child.text))


        # Setup landmark lists.
        ldata = {} # landmark lists
        ldata[best] = []
        alnlen = -1

        # Read fasta files.
        self.logf.write('Reading fasta files for final alignment...')
        for fi in fastas:
            seqs = FASTAnet.FASTAstructure(fi,uniqueOnly=False)
            if len(seqs.sequences) != 2:
                self.logf.write('WARNING; <%s> does not contain exactly 2 sequences.' % (fi))
            if len(seqs.sequences) < 2: continue
            teststruc = seqs.orderedSequences[0].sequence
            refstruc  = seqs.orderedSequences[1].sequence
            if (teststruc == refstruc):
                self.logf.write('NOTE; <%s> contains identical sequences.' % (fi))
            elif len(refstruc) != len(teststruc):
                self.logf.write('ERROR; <%s> has non-matching lengths.' % (fi))
                return None
            f = IO.getFileName(fi) # Get name of PDB.
            ldata[f] = []
            if not ldata[best]:
                for i in range(len(refstruc)):
                    if refstruc[i] != '-': ldata[best].append(1) # 1 for landmark, 0 for no landmark   
            for i in range(len(refstruc)):
                if refstruc[i] != '-':
                    if teststruc[i] == '-': ldata[f].append(0)
                    else: ldata[f].append(1)
            if alnlen == -1: alnlen = len(ldata[f])

        # List containing how many "None"s are in each structure's landmarks.
        # First item in each item (tuple) is chain number, second item is number of empty landmarks.
        landmarkList = []
        for key in ldata: landmarkList.append((key, ldata[key].count(0)))
        landmarkList = sorted(landmarkList, key=operator.itemgetter(1), reverse=True)

        # Optimize the number of landmarks.
        bad = []

        # Filter out non-significant P-value alignments.
        self.logf.write('Filtering non-significant P-value alignments (RRMSD)...')
        for key, count in landmarkList:
            if key == best: continue
            fname = os.path.join(self.reffldr,'%s.pickl' % (key))
            if not os.path.isfile(fname):
                bad.append(key)
                self.logf.write('Ignoring because not scored; <%s>.' % (key))
                continue
            scf = scoreFile(fname)
            scs = scf.getScores()
            for sn,sc,pv in scs:
                if sn == 'RRMSD' and pv < 0.05:
                    bad.append(key)
                    self.logf.write('Omitting for P-value %f; <%s>.' % (pv,key))

        # Optimize total landmarks.
        #if self.optim:
            #print '' # Newline
            #self.logf.write('Optimizing for landmarks...')
            #n = len(self.args)
            #totalNum = NumLandmarks(ldata, bad)
            #denom = n-len(bad)+1
            #if (denom == 0):
                #totalNum = 0
                #denom = 1
            #self.logf.write('Number of landmarks each before: %d.' % (totalNum/(denom)))
            #self.logf.write('Number of landmarks in total before: %d.' % (totalNum))
            #for key, count in landmarkList:
                #if key in bad: continue
                #if count > 0:
                    #if NumLandmarks(ldata, bad + [key]) > totalNum:
                        #self.logf.write('Omitting for empty landmarks; <%s>.' % (key))
                        #bad.append(key)
                        #totalNum = NumLandmarks(ldata, bad)
            #denom = n-len(bad)+1
            #if (denom == 0):
                #totalNum = 0
                #denom = 1          
            #self.logf.write('Number of landmarks after: %d.' % (totalNum/(denom)))
            #self.logf.write('Number of landmarks in total after: %d.' % (totalNum))
        #else: self.logf.write('No optimization was done.')
        if self.optim:
            self.logf.write('Optimizing for landmarks...')
            if self.args: n = len(self.args)
            else:         n = len(ldata.keys())
            if len(ldata.keys())-1 <= 0:
                self.logf.write('ERROR; Landmark list too small.')
                return None
            optli = {}
            for key in ldata.keys():
                if key == best: continue
                else: optli[key] = ldata[key]
            t = treeoptimizer(optli,bad,n-1,self.logf)
            totalNum = t.getNumLandmarks()
            denom = n-len(bad)+1
            if (denom == 0):
                totalNum = 0
                denom = 1
            self.logf.write('Number of landmarks each before: %d.' % (totalNum/(denom)))
            self.logf.write('Number of landmarks in total before: %d.' % (totalNum))
            t.optimize()
            totalNum = t.getNumLandmarks()
            o_bad = t.getOmitList()
            detected = list(set(o_bad).difference(set(bad)))
            for d in detected: self.logf.write('Omitting for bad landmarks; <%s>.' % (d))
            denom = n-len(o_bad)+1
            if (denom == 0):
                totalNum = 0
                denom = 1          
            self.logf.write('Number of landmarks after: %d.' % (totalNum/(denom)))
            self.logf.write('Number of landmarks in total after: %d.' % (totalNum))
            self.logf.write('Number of landmarks removed by optimization: %d.' 
                            % (len(detected)))
            # Write the optimized GM file.
            data = gmdata(fastas,o_bad,alnlen,ldata,best,carbon,self.optpre,self.logf)
            d = data.write()
            if d:
                outfiles.extend(d)
                self.logf.write('Optimized GM data successfully written to file.')
            else: self.logf.write('ERROR; No common columns for optimized data.')
        else: self.logf.write('No optimization was done.')        

        ########## Write the final GM and landmarks files ###########

        data = gmdata(fastas,bad,alnlen,ldata,best,carbon,self.outpre,self.logf)
        d = data.write()
        if not d:
            self.logf.write('ERROR; No common columns for non-optimized data.')
            return outfiles
        else:
            outfiles.extend(d)
            self.logf.write('Non-optimized GM data successfully written to file.')

        if self.touchup:
            # Get a dictionary of landmark coordinates from the GM file.
            landmarks = self.GetLandmarks()
            # Go to first landmarks item.
            first = landmarks.keys()[0]
            nLand = len(landmarks[first][0]) # the number of landmarks per structure
            idealstruc = [[], [], []]
            for i in range(nLand):
                idealstruc[0].append(median([landmarks[x][0][i] for x in landmarks.keys()]))
                idealstruc[1].append(median([landmarks[x][1][i] for x in landmarks.keys()]))
                idealstruc[2].append(median([landmarks[x][2][i] for x in landmarks.keys()]))
            idealcenter = [mean(idealstruc[0]), mean(idealstruc[1]), mean(idealstruc[2])]

            # New coordinates for landmark points.
            newlandmarks = {}
            for struc in landmarks:
                xV = idealcenter[0] - mean(landmarks[struc][0])
                yV = idealcenter[1] - mean(landmarks[struc][1])
                zV = idealcenter[2] - mean(landmarks[struc][2])
                translationV = [xV, yV, zV] ## translation vector
                newlandmarks[struc] = [[], [], []]
                for i in range(nLand):
                    newlandmarks[struc][0].append(landmarks[struc][0][i] + xV)
                    newlandmarks[struc][1].append(landmarks[struc][1][i] + yV)
                    newlandmarks[struc][2].append(landmarks[struc][2][i] + zV)

            # Write the new GM file.
            f = open(self.prefix + '_new.gm', 'w')
            outfiles.append(self.prefix + '_new.gm')
            for struc in newlandmarks:
                f.write('>%s_new:%s;'%(self.prefix, struc))
                for i in range(nLand):
                    f.write('%s;%s;%s;'%(str(round(newlandmarks[struc][0][i],3)), \
                                         str(round(newlandmarks[struc][1][i],3)), \
                                         str(round(newlandmarks[struc][2][i],3))))
                f.write('\n')
            f.close()

        return outfiles

# For purposes of generating a multiple str. alignment RMSD.

def GMtoMatrix(prefix):
    f = open(prefix + '.gm')
    s = {} ## s for structures
    for line in f:
        if line.strip().endswith(';'):
            line = line.strip()[:-1]
        temp = [[], [], []]
        bline = line[line.find(';')+1:].strip().split(';')
        for i in range(len(bline)):
            temp[i%3].append(float(bline[i]))
        s[line[:line.find(';')]] = temp
    return s

def RMSDfromMEAN(p,m):
    # Calculates the RMSD for the meanshape (m) and the problem shape (p)
    sum_dist=0.0
    count = 0
    for j in range(len(m[0])):
        d =(m[0][j]-p[0][j])**2 + (m[1][j]-p[1][j])**2 + (m[2][j]-p[2][j])**2
        sum_dist += d
        count += 1
    # calculate the sqrt of the mean deviation
    RMSD = math.sqrt(sum_dist/count)

    return RMSD

class pseudoalignment(gmwriter):

    def __init__(self,refldr):
        super(pseudoalignment,self).__init__(None,'temp',None,refldr,
                                             logfile('',False,True),None,True)

    def getRMSD(self):
        matrix   = []
        outfiles = super(pseudoalignment,self).write()
        if outfiles: structures = GMtoMatrix(self.outpre)
        else:        return -1.0
        names = sorted(structures.keys())
        for i in names:
            for j in names:
                if i != j:
                    matrix.append(RMSDfromMEAN(structures[j],structures[i]))
        for f in outfiles: os.remove(f) # Clean outfiles.
        return mean(matrix)
