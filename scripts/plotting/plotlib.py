# plotlib.py
# -------------------------
# Dec 16, 2013; Alex Safatli
# -------------------------
# Asset functions and classes
# to allow plotting functionality
# of ABeRMuSA output. Used by plot.py.

# Imports.

import os, pickle, utils, utils.io, numpy, cPickle, \
       pylab, scipy.stats, glob, operator
import ABeRMuSA.GMWriter as gmw
from utils.logfile import parseTime, calculateTime
from ABeRMuSA.scripts.output.summarizer import damastesRun as aRun
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import pyplot as plt

# Constants.

RMSD_THRES_VALUE = 15

# Visualization-related.

class plotter:
    
    ''' Handles the abstraction of matplotlib plotting
    functions and interfacing with a symbolManager to
    create uniformly formatted plots. '''
    
    def __init__(self,datafolder,plotfolder,symbMang):
        self.datafolder    = datafolder
        self.plotfolder    = plotfolder
        self.symbolManager = symbMang
        self.figures       = {}
    def isFigure(self,name):
        return (name in self.figures)
    def __makenewfig__(self,name):
        fol = os.path.join(self.plotfolder,name)
        if not os.path.isdir(fol): os.mkdir(fol)
        self.figures[name] = fol
    def new(self,name,suffix=None):
        plt.cla()
        if not self.isFigure(name): self.__makenewfig__(name)
        figpath = self.figures[name]
        if suffix: outfile = os.path.join(figpath,'%s.%s.pdf' % (name,suffix))
        else:      outfile = os.path.join(figpath,'%s.pdf' % (name))
        fi = plotfile(name,outfile)
        return fi         
    def writeData(self,name,suffix=None,li=None):
        figpath = self.figures[name]
        if suffix: outfile = os.path.join(figpath,'%s.%s.plotdata' % (name,suffix))
        else:      outfile = os.path.join(figpath,'%s.plotdata' % (name))
        fi = plotdata(name,outfile)
        if li: fi.write(li)
        return fi
    def plot(self,xdata,ydata,fi=None,**kwargs):
        colrs   = self.symbolManager.colors
        symbs   = self.symbolManager.symbols
        if fi:
            if not 'color' in kwargs and not 'marker' in kwargs:
                plt.plot(xdata,ydata,color=colrs[fi.name],marker=symbs[fi.name],**kwargs)
            elif 'marker' in kwargs:
                plt.plot(xdata,ydata,color=colrs[fi.name],**kwargs)
            elif 'color' in kwargs:
                plt.plot(xdata,ydata,marker=symbs[fi.name],**kwargs)
        else:  plt.plot(xdata,ydata,**kwargs)

class plotdata:
        
    ''' An abstraction of a tab-delimited file containing
    the contents of a list of tuples as single lines. Meant
    to be data coming along with a plot of that data. '''
    
    def __init__(self,figname,fpath):
        self.figure = figname
        self.fpath  = fpath
        if not os.path.isfile(fpath):
            fh = open(fpath,'w')
            fh.close()
    def write(self,li):
        o = open(self.fpath,'w')
        for l in li: o.write('\t'.join([str(x) for x in l])+'\n') 
        o.close()

class plotfile:
    
    ''' An abstraction of the actual plot file containing the
    matplotlib figure contents as defined by the current figure.
    When save function is called, saves to a pdf file. '''
    
    def __init__(self,figname,fpath):
        self.figure = figname
        self.fpath  = fpath
        if not os.path.isfile(fpath):
            fh = open(fpath,'w')
            fh.close()
    def save(self):
        pp = PdfPages(self.fpath)
        plt.savefig(pp,format='pdf')
        pp.close()

class symbolManager:
    
    ''' An object that handles the classification of given
    database/aligner pairs into symbols and colors for the use
    of visualization onto a series of plots with a uniform color
    and symbol scheme.'''
    
    def __init__(self):
        # Available symbols.
        self.savail = ['o','*','d','s','D','p','h']
        # Available colors.
        self.cavail = ['b','w','0.8','r','k','g']
        # Indices/dicts for personal use.
        self.cused = []
        self.sind = -1
        self.cind = -1
        self.__dbs = {}
        self.__als = {}
        # Color and symbol maps.
        self.colors = {}
        self.symbols = {}
    def __clr__(self,ali):
        if ali in self.__als:
            return self.__als[ali]
        self.cind += 1
        c = self.cavail[self.cind]
        self.__als[ali] = c
        return c
    def __symb__(self,db):
        if db in self.__dbs:
            return self.__dbs[db]
        self.sind += 1
        symb = self.savail[self.sind]
        self.__dbs[db] = symb
        return symb
    def writeSymbolTable(self,fipath):
        o = open(fipath,'w')
        for key in self.__dbs:
            o.write('%s\t%s\n' % (key,self.__dbs[key]))
        o.close()
    def writeColorTable(self,fipath):
        o = open(fipath,'w')
        for key in self.__als:
            o.write('%s\t%s\n' % (key,self.__als[key]))
        o.close()
    def add(self,key,db,align):
        colr = self.__clr__(align)
        symb = self.__symb__(db)
        self.colors[key]  = colr
        self.symbols[key] = symb

# Fundamental data structures abstracting ABeRMuSA raw output.

class dataset:
    
    ''' Abstracts an entire database/aligner pair. '''
    
    def __init__(self):
        self.name = ''
        self.path = ''
        self.name_header = ''
        self.database = ''
        self.names = []
        self.headers = []
        self.folders = []
        self.data = {}
    def getNumStrForMSTA(self,name):
        mstafi = os.path.join(self.path,'matt-msta*.sh')
        findsh = glob.glob(mstafi)
        if len(findsh) < 1: return (-1)
        o     = open(findsh[0])
        lines = o.readlines()
        k     = '%s/' % (name)
        s     = '%s.elapsed' % (name)
        i     = 0
        o.close()
        while (not s in lines[i] and not k in lines[i]): i += 1
        line  = lines[i]
        return line.count('.atm')+line.count('.pdb')+line.count('.ent')

class refset:
    
    ''' Abstracts a set of pairwise alignments to a single reference. '''
    
    def __init__(self):
        self.name    = ''
        self.refFile = ''
        self.scores = {}
        self.rmsds  = {}
        self.files  = {}
        self.mrmsd  = 0 # Multiple Str. Alignment RMSD
        self.avgscr = 0
        self.highpv = 0
        self.rank   = 0

class folder:
    
    ''' Abstracts a single group within a database/aligner pair. '''
    
    def __init__(self):
        self.name  = ''
        self.best  = None
        self.worst = None
        self.avg   = None
        self.refs  = []

# Reading and parsing of entire collections of ABeRMuSA data.

class databaseList:

    ''' Reads and processes a database list file exclusive to
    the plotting system used here. '''

    def __init__(self,fi,quiet=False):
        self.filename = fi
        self.args     = [] # Argument List
        self.dbs      = {} # Database List
        self.als      = {} # Aligner  List
        self.names    = []
        self.setnames = []
        self.filehand = open(fi)
        self.quiet    = quiet
        
    def close(self):
        self.filehand.close()
        
    def read(self):

        # Initialize structures.

        names = []
        setnames = []        

        args = [] # Argument List
        dbs = {}  # Database List
        als = {}  # Aligner List

        # Read.

        dbl = self.filehand.readlines()
        for line in dbl:

            # Parse line.

            if line.startswith('//'): continue
            crede = line.split('\t')
            if len(crede) != 4: continue

            # Get three fields.

            dbname, dbpath, aliname, setname = crede

            # Determine database and aligner type.

            dbs[dbpath] = dbname.lower()
            als[dbpath] = aliname.lower()
            args.append(dbpath)
            setnames.append(setname.strip())            

        self.args     = args
        self.dbs      = dbs
        self.als      = als
        self.names    = names
        self.setnames = setnames

class datafiles:
    
    ''' Contains different MStA and ABeRMuSA datasets and 
    appropriate symbol manager for plotting purposes. Wraps
    a databaseList instance. '''
    
    def __init__(self,dbList,symbManager,quiet=False):
        self.datalist = dbList
        self.symbmang = symbManager
        self.args     = dbList.args
        self.als      = dbList.als
        self.dbs      = dbList.dbs
        self.setnames = dbList.setnames
        self.names    = dbList.names
        self.sumf     = [] # Summary file list.
        self.elaf     = [] # Elapsed file list.
        self.ssets    = [] # Summary set list.
        self.esets    = [] # Elapsed set list.
        self.quiet    = quiet
        
    def flatten(self,datafolder):
        
        ''' Collects all data associated with the raw output
        of ABeRMuSA files and does computation to put it into 
        a readable format. '''
        
        # Flatten all summary files.
        
        for fi in self.sumf:
            
            # Establish flat file.
            
            name = fi.name
            n = os.path.join(datafolder,name+'.dataset')
            
            if os.path.isfile(n):
                
                # Data has already been collected.
                
                if not self.quiet:
                    print '\tAcquiring metadata for %s from file...' % (name)
                o = open(n)
                loaded = pickle.load(o)    
                self.ssets.append(loaded)
                continue
            
            if not self.quiet:
                print '\tCompiling and flattening %s to file. This will take time.' % (name)
            
            gcnt = 0
            for group in fi.names:
                
                # Set up a folder abstraction of this group.
                
                fold = folder()
                
                # Go to every group in this d/a pair; acquire XML/log file.
                
                x = os.path.join(fi.path,group+'.xml') # XML File.
                l = os.path.join(fi.path,group+'.log') # Log File.
                
                if not os.path.isfile(x):
                    print '\t\t(%d/%d) %s:%s File does not exist.' \
                          % (gcnt+1,len(fi.names),fi.name,x)
                    continue
                
                # Parse run data from log and XML file using an appropriate interface.
                
                run = aRun(x)
                 
                # Get reference raw data.
                
                refs = run.getReferences(l)
                for ref in refs:
                    
                    # Get OS path details for all files.
                    refName  = utils.io.getFileName(ref)
                    refPath  = os.path.join(fi.path,refName)
                    refPickl = os.path.join(refPath,'ref.pickl')
                    
                    # See if files actually exist.
                    if not os.path.isfile(refPickl):
                        if not self.quiet:
                            print '\t\t(%d/%d) %s --> %s Reference does not exist.' \
                                  % (gcnt+1,len(fi.names),fi.name,refName)
                        continue
                    
                    # Set up a refset object abstraction of this reference.
                
                    refObj         = refset()
                    refObj.name    = refName
                    refObj.refFile = refPickl
                    o              = open(refPickl)
                    load           = cPickle.load(o)
                    refObj.scores  = load[0]
                    refObj.avgscr  = load[1]
                    refObj.highpv  = load[2]
                    refObj.rank    = load[3]
                    o.close()
                    
                    # Individual RMSD scores.
                    for score in load[0]:
                        f = os.path.join(refPath,score+'.pickl')
                        o = open(f)
                        _,r,_ = cPickle.load(o)
                        o.close()
                        refObj.files[score] = f                        
                        refObj.rmsds[score] = r                        
                    
                    # Multiple RMSD score for entire reference.
                    ali          = gmw.pseudoalignment(refPath)
                    refObj.mrmsd = ali.getRMSD()
                    
                    # Append this object to the folder.
                    fold.refs.append(refObj)
                        
                # Determine best and worst references.
             
                fold.refs = sorted(fold.refs,key=lambda d: d.rank)
             
                if run.reference:
                    for ref in fold.refs:
                        if ref.name == run.reffldr:
                            best = ref
                            break
                else:
                    m, best = None, None
                    for ref in fold.refs:
                        if m == None and ref.rank >= 0:
                            m = ref.rank
                            best = ref
                        if ref.rank < m and ref.rank >= 0:
                            m = ref.rank
                            best = ref
                
                m, worst = None, None
                for ref in fold.refs:
                    if m == None and ref.rank >= 0:
                        m = ref.rank
                        worst = ref
                    if ref.rank > m and ref.rank >= 0:
                        m = ref.rank
                        worst = ref
                        
                filtered = [x for x in fold.refs if x.rank >= 0] 
                mid       = len(filtered)/2
                if mid < len(filtered): avg = filtered[mid]
                else:                   avg = None
                avgscc   = -1
                
                if best: bestna, bestra, avgscc = best.name, best.rank, best.avgscr
                else:    bestna, bestra = 'None', -1.0
                
                if not self.quiet: 
                    print '\t\t(%d/%d) Registered %s (best: %s, avgscr: %f, rank: %f).' % \
                          (gcnt+1,len(fi.names),group,bestna,avgscc,bestra)
                
                fold.best, fold.worst, fold.avg, fold.name = best, worst, avg, group
                fi.folders.append(fold)
                gcnt += 1
                
            fi.omitlist = filterSset(fi,datafolder)
            print '\tSaving %s to file...' % (name)
            f = open(n,'w')
            pickle.dump(fi,f)
            f.close()
            self.ssets.append(fi)
            
        # Flatten all elapsed files.
    
        for fi in self.elaf:
        
            # Reformat name in case of overwrite purposes.
            name = fi.name
            while (name in self.names): name + '_'
            self.names.append(name)
            fi.name = name
            
            # Establish flat file.
            n = os.path.join(datafolder,name+'.mstaset')
            if os.path.isfile(n):
                o = open(n)
                loaded = pickle.load(o)
                self.esets.append(loaded)
                o.close()
                if not self.quiet: print '\tAcquiring metadata for %s from file...' % (name)
                continue
            if not self.quiet: print '\tCompiling and flattening %s to file...' % (name)
            f = open(n,'w')
            pickle.dump(fi,f)
            f.close()
            self.esets.append(fi)        
            
    def read(self):

        # For every argument given.

        for arg in self.args:

            f = open(arg)
            lines = f.readlines()
            name_header = ''
            dnames = []
            headers = []
            data = {}

            # For every line in argument.

            for line in lines:

                if line == lines[0]:
                    # If first line, header line.
                    fields = line.strip().split('\t')
                    for field in fields[1:]: headers.append(field)
                    name_header = fields[0]
                else:
                    # Is a line of data.
                    fields = line.strip().split('\t')
                    name = fields[0]
                    if len(fields[1:]) != len(headers):
                        if not self.quiet:
                            print '\t%s; Line %d does not have a valid number of columns.' % \
                                  (arg,lines.index(line))
                        continue
                    if name in dnames:
                        if not self.quiet:
                            print '\t%s; Line %d had a duplicate name in its first column.' % \
                                  (arg,lines.index(line))
                        continue
                    dnames.append(name)
                    data[name] = {}
                    i = 0
                    for field in fields[1:]:
                        data[name][headers[i]] = field
                        i += 1
                        
            d = dataset()
            d.data        = data
            d.headers     = headers
            d.names       = dnames
            d.name_header = name_header
            d.aligner     = self.als[arg]
            d.database    = self.dbs[arg]
            d.path        = os.path.dirname(arg)
            d.name        = self.setnames[self.args.index(arg)]
            self.symbmang.add(d.name,d.database,d.aligner)
            if 'Num.Str.' in headers: self.elaf.append(d)
            else: self.sumf.append(d)
            f.close()
            
def filterSset(sset,datafolder):
    filterlist = []
    filterfile = os.path.join(datafolder,'%s.filterlist' % (sset.name))
    if os.path.isfile(filterfile):
        #print '\tOpening filter list for %s...' % (sset.name)
        o = open(filterfile)
        for line in o: filterlist.append(line.split('\t')[0])
        o.close()
    else:
        o = open(filterfile,'w')
        for fo in sset.folders:
            for ref in fo.refs:
                if fo.name in filterlist: continue
                for key in ref.rmsds:
                    rmsd = ref.rmsds[key]
                    if rmsd > RMSD_THRES_VALUE:
                        # Omit this and add to filterlist.
                        filterlist.append(fo.name)
                        o.write('%s\t%s\t%s\t%f\n' % (fo.name,ref.refFile,key,rmsd))
                        break
        o.close()
    return filterlist

def parseUnixTime(time):
    splitted = time.split(':')
    vector = [0,0,0,0,0,0]
    i = -1
    for split in reversed(splitted):
        vector[i] = int(float(split))
        i -= 1
        if (i == -7): break
    return vector