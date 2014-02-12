#!/usr/bin/env python

# plot.py
# -------------------------
# Dec 16, 2013; Alex Safatli
# -------------------------
# Flatten data needed for plots
# for the ABeRMuSA paper and save
# to PDF files dynamically. 

'''
Paper plan for plots, presentation of data since Nov. 2013.
--------------------------------------------------------------
1. There are limitations to current MStA implementations:
    - Demonstrate that current MStA fail (%)
    - What can we say about these that fail? Size? 
    - Our method will not fail, and if it does, it simply exclude the offending structure.

2. We can tell significant alignment from others using RRMSD and an empirical distribution of score.
   - Empirical distributions.
   - Normal distribution, test.

3. Not all structures, used as reference, give equally good results.
   - What is the gain from an average to the best reference?
   - What is the difference between picking the best and the worst reference?

4. Heuristic for identifying a reasonable reference
   - What is the RRMSD improvement with respect to exploring % of the dataset

// Test case, something big.

5. It works well with other pairwise alignment implementations.  
-------------------------------------------------------------
'''

# Imports :::::::::::::::::::::::::::::::::::::::::::::

import sys, os, pickle, cPickle, numpy, pylab, scipy.stats, glob
from random import randint
from utils.logfile import parseTime, calculateTime
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import pyplot as plt
from plotlib import *

# Constants :::::::::::::::::::::::::::::::::::::::::::

DATA_FOLDER    = 'plot_data'
PLOT_FOLDER    = 'plot_figs'
TEST_ALIGNER   = 'matt'
QUICK_SCRIPT   = 'quickModeTester.py'
NUM_PDBS_THRES = 12

# Initialization ::::::::::::::::::::::::::::::::::::::

# Make folders.
if not os.path.isdir(DATA_FOLDER): os.mkdir(DATA_FOLDER)
if not os.path.isdir(PLOT_FOLDER): os.mkdir(PLOT_FOLDER)

# Set up color/symbol manager.
colrs = symbolManager()

# Get database file.
if len(sys.argv) != 2: print 'usage: %s database_file' % (sys.argv[0])
databaseFile = sys.argv[1]

# Read and process database file.
print 'Reading database file (%s)...' % (databaseFile)
dbFile = databaseList(databaseFile)
dbFile.read()
dbFile.close()

# Process and flatten all output data from ABeRMuSA.
data = datafiles(dbFile,colrs)
data.read()
data.flatten(DATA_FOLDER)

# Set up plotting manager.
pltr = plotter(DATA_FOLDER,PLOT_FOLDER,colrs)

# Plotting ::::::::::::::::::::::::::::::::::::::::::::

''' [0] Write symbol table and color table. '''

symbTable = os.path.join(DATA_FOLDER,'symboltable')
colrTable = os.path.join(DATA_FOLDER,'colortable')
print '>>> Writing symbol table to %s.' % (DATA_FOLDER)
colrs.writeSymbolTable(symbTable)
print '>>> Writing color table to %s.' % (DATA_FOLDER)
colrs.writeColorTable(colrTable)

''' [1] Demonstrate that current MStA implementations fail. '''

figName = '1_msta_implem_fails'
plot    = pltr.new(figName,'scop')
datf    = pltr.writeData(figName,'scop.fails')
print '>>> Aligner to be utilized for testing is set as <<%s>>.' % (TEST_ALIGNER)
print '>>> Gathering data and plotting figure [%s].' % (figName)

f_x, f_y, fails = [], [], []
for fi in data.esets:

    # Set up graph data.
    x,y,fx,fy = [0],[0],[],[]
    
    # Go over each group.
    for name in fi.names:
        if not name in fi.data: continue
        val     = fi.data[name]['Elapsed']  # Get Elapsed Time data.
        elapsed = parseUnixTime(val)        # Parse it.
        time    = calculateTime(elapsed)    # Get in seconds.
        failed  = False 
        numstr  = fi.data[name]['Num.Str.'] # Number of structures.
        num     = int(numstr)
        if (num == 0):
            # A zero-value in the elapsed file implies failure of alignment.
            failed = True
            # Get correct number of structures by interpolation.
            num = fi.getNumStrForMSTA(name)
            # Record failure into a list for later writing.
            fails.append((name,num))
        if failed:
            fx.append(num)
            fy.append(time)
        else:
            x.append(num)
            y.append(time)
        f_x.append(num)
        f_y.append(time)
    
    # Put acquired data up on plot.
    pltr.plot(x,y,fi,color='k',ls='None',ms=10) # Not failed.
    pltr.plot(fx,fy,fi,color='w',ls='None',ms=10) # Failed.
    
plt.yscale('log')
plt.ylim(0,max(f_y)+10)
plt.xlim(0,max(f_x)+10)
plt.ylabel('Time (s)',fontsize=20)
plt.xlabel('Number of Structures',fontsize=20)
datf.write(fails) # Write records of failures.
plot.save() # Save figure.

plot = pltr.new(figName,'abermusa')

for fi in [x for x in data.ssets if x.aligner.lower() == TEST_ALIGNER]:
    
    # Set up graph data.
    grd,x,y,yerr = {},[0],[0],[0]
    
    # Go over each group.
    for fld in fi.folders:
        name    = fld.name
        if not name in fi.data or fi.data[name]['Output'] != 'Y':
            continue
        val     = fi.data[name]['Elapsed']
        elapsed = parseTime(val)
        time    = calculateTime(elapsed)
        if (time == 0): continue
        try: nums = int(fi.data[name]['Num. Structs'])
        except:         continue
        if nums in grd: grd[nums].append(time)
        else:           grd[nums] = [time]
    
    # Go over acquired data and calculate error.
    for nums in grd:
        x.append(nums)
        y.append(numpy.average(grd[nums]))
        stderr = numpy.std(grd[nums])/((len(grd[nums]))**(1/2.0))
        yerr.append(stderr)
        
    # Plot errorbars.
    plt.errorbar(x,y,yerr=yerr,ls='None',ecolor='0.4')
    
    # Plot actual points.
    pltr.plot(x,y,fi,ls='None',label=fi.name,ms=10)

plt.yscale('log')
plt.ylim(0,max(f_y)+10)
plt.xlim(0,max(f_x)+10)
plt.ylabel('Time (s)',fontsize=20)
plt.xlabel('Number of Structures',fontsize=20)
# Without legend.
plot.save()
# With legend.
plot.fpath = plot.fpath[:-4] + '.legend.pdf'
plt.legend()
plot.save()

''' [2] We can tell significant alignment from others using RRMSD 
        and an empirical distribution of score. '''

figName = '2_stat_distr'
plot    = pltr.new(figName)
datf    = pltr.writeData(figName)
print '>>> Gathering data and plotting figure [%s].' % (figName)

maxval, graphset, stats, full = 0, {}, [], []
for fi in [x for x in data.ssets if x.aligner.lower() == TEST_ALIGNER]:
    
    # Gather all RRMSD data.
    graphset[fi.name] = []
    for fld in fi.folders:
        for ref in fld.refs:
            r = [x for x in ref.scores.values() if x]
            graphset[fi.name].extend(r)
    if len(graphset[fi.name]) <= 1:
        print '\tWARNING (%s): Distribution too small!' % (fi.name)
        continue
    
    # Perform normality test.
    dat   = graphset[fi.name]
    me,st = scipy.stats.norm.fit(dat)
    _,pva = scipy.stats.normaltest(dat)
    print '\t(%s) had mean %.3f, stdv %.3f, pval %.3f' % (fi.name,me,st,pva)
    stats.append((fi.name,me,st,pva))
    full.extend(dat)
    
# Project all data onto histogram.
n,_,_ = plt.hist(full,bins=50,color='0.8')
maxval = max(maxval,max(n))

# Plot found normal distribution used in ABeRMuSA.
rang = numpy.arange(0,1,0.001)
print '\tABeRMuSA has fixed distribution; mean %.3f, stdv %.3f' % (0.177,0.083)
pdfr = scipy.stats.norm.pdf(rang,0.177,0.083)
pltr.plot(rang,pdfr*(maxval/max(pdfr)),None,color='red',
          ls='solid',label='Derived Distribution')
plt.xlim(0.0,1.0)
plt.ylabel('Frequency',fontsize=20)
plt.xlabel('RRMSD',fontsize=20)
plt.legend()
plot.save() # Save figure.
datf.write(stats) # Write record of found statistics.

''' [3] Reference matters. '''

figName = '3_ref_matters'
plotr   = pltr.new(figName,'rmsds')
rmsdf   = pltr.writeData(figName,'rmsds')
rrmsdf  = pltr.writeData(figName,'rrmsds')
print '>>> Gathering data and plotting figure [%s].' % (figName)

# Get data.
graphsetr,graphsetrr= [],[]
bestr, worstr, avgr = [],[],[]
bestrr,worstrr,avgrr= [],[],[]

for fi in [x for x in data.ssets if x.aligner.lower() == TEST_ALIGNER]:
    
    # Go over each group.
    for fld in fi.folders:
        name    = fld.name
        if not name in fi.data or fi.data[name]['Output'] != 'Y':
            continue
        best    = fld.best
        worst   = fld.worst
        avg     = fld.avg
        if best and worst and avg:
            # Ignore outlying or failing sets.
            if best.mrmsd >= 0 and worst.mrmsd >= 0 and avg.mrmsd >= 0:
                bestr.append(best.mrmsd)
                worstr.append(worst.mrmsd)
                avgr.append(avg.mrmsd)
                graphsetr.append((fi.name,fld.name,best.mrmsd,
                                  worst.mrmsd,avg.mrmsd))
            if best.avgscr >= 0 and worst.avgscr >= 0 and avg.avgscr >= 0:
                bestrr.append(best.avgscr)
                worstrr.append(worst.avgscr)
                avgrr.append(avg.avgscr)
                graphsetrr.append((fi.name,fld.name,best.avgscr,
                                   worst.avgscr,avg.avgscr))
        

# Set up histograms for RMSDs.
fig    = plt.figure()
ax     = fig.add_subplot(111)
ax.hist(bestrr,bins=50,color='0.8')
ax.hist(avgrr,bins=50,color='0.5')
ax.text(meanBest-0.007,121,'B',fontsize=12)
ax.text(meanAvg-0.007,121,'A',fontsize=12)
ax.text(meanWor-0.007,121,'W',fontsize=12)
ax.hist(worstrr,bins=50,color='0.2')
plt.ylabel('Frequency',fontsize=20)
plt.xlabel('RMSD',fontsize=20)
ax.spines["right"].set_visible(False)
ax.spines["top"  ].set_visible(False)
ax.tick_params(axis="both",direction="out")
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()
plotr.save()

# Set up histograms for RRMSDs.
plotrr = pltr.new(figName,'rrmsds')
fig    = plt.figure()
ax     = fig.add_subplot(111)
meanBest = numpy.mean(bestrr)
meanAvg  = numpy.mean(avgrr )
meanWor  = numpy.mean(worstrr)
ax.axvline(x=meanBest,color='0.8')
ax.axvline(x=meanAvg,color='0.5')
ax.axvline(x=meanWor,color='0.2')
ax.annotate('',(meanBest,115),(meanAvg,115),
             arrowprops={'arrowstyle':
                         '|-|,widthA=0.25,widthB=0.25'})
ax.annotate('',(meanBest,110),(meanWor,110),
             arrowprops={'arrowstyle':
                         '|-|,widthA=0.25,widthB=0.25'})
ax.text    (meanBest-0.0275,116,'%.3f' % (meanAvg-meanBest),
            fontsize=8)
ax.text    (meanWor+0.005,111,'%.3f' % (meanWor-meanBest),
            fontsize=8)
ax.text    (meanBest-0.007,121,'B',fontsize=12)
ax.text    (meanAvg-0.007,121,'A',fontsize=12)
ax.text    (meanWor-0.007,121,'W',fontsize=12)
ax.hist(bestrr,bins=50,color='0.8')
ax.hist(avgrr,bins=50,color='0.5')
ax.hist(worstrr,bins=50,color='0.2')
plt.ylabel('Frequency',fontsize=20)
plt.xlabel('RRMSD',fontsize=20)
ax.spines["right"].set_visible(False)
ax.spines["top"  ].set_visible(False)
ax.tick_params(axis="both",direction="out")
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()
plotrr.save()

# Save to datafiles.
rmsdf.write(graphsetr)
rrmsdf.write(graphsetrr)

''' [4] The heuristic works. '''

figName = '4_ref_heuristic'
plotr   = pltr.new(figName,'perc')
figf    = pltr.writeData(figName,'perc')
print '>>> Gathering data and plotting figure [%s].' % (figName)

# Get data.
## percdata: % of database which found best reference.
## diffdata: deltaRRMSD data.
## numsdata: Number of structures on average per.
## realdata: Dataset actual instances by key.
## refsdata: Actual reference data by key.
## partdata: Get the value for 3/numstructs on average.
percdata, diffdata, numsdata = {}, {}, {}
realdata, refsdata, partdata = {}, {}, {}

for fi in [x for x in data.ssets if x.aligner.lower() == TEST_ALIGNER]:
    ## Following code is extracted from abPlots; could be optimized further.
    
    # Extract found_in data created by quick search analyzer script.
    found_in = os.path.join(fi.path,'found_in.txt')
    
    # If found_in data does not exist, refer to the quick search script.
    if not os.path.isfile(found_in):
        print '\tINFO (%s): No quick search analysis found.' % (fi.name)
        print '\tINFO (%s): Please run %s on all XML files in path <%s>.' \
              % (fi.name,QUICK_SCRIPT,fi.path)
        continue

    # Set up data structures.
    percdata[fi.name] = {}
    diffdata[fi.name] = {}
    refsdata[fi.name] = {}
    numsdata[fi.name] = {}
    realdata[fi.name] = fi
    
    # Get data from found_in.
    foundinf = open(found_in)
    for li in foundinf:
        spl         = li.split('\t')
        if len(spl) != 3: continue
        n,perc,diff = spl
        percdata[fi.name][n] = int(perc)
        diffdata[fi.name][n] = float(diff.strip('\n'))
    foundinf.close()
    
    # Get partial portion / dataset; assessed on avg.
    parts = []
    for foldr in fi.folders:
        if not 'homstrad' in fi.name:
            numstr = len(foldr.refs)
            if numstr != 0: parts.append(min(1,float(3)/numstr))
            else:           parts.append(0)
        if len(foldr.refs) == 0:
            print '\tINFO (%s): %s contained 0 items.' % (fi.name,foldr.name)
            continue
        aref       = foldr.refs[0]
        numInGroup = len(aref.scores)
        numsdata[fi.name][foldr.name] = numInGroup
        if foldr.best: refsdata[fi.name][foldr.name] = foldr.best.name
        else:          refsdata[fi.name][foldr.name] = None
    if not 'homstrad' in fi.name and len(parts) > 0:
        partdata[fi.name] = sum(parts)/float(len(parts))
    
    # Plot this specific dataset.
    x = sorted([int(_) for _ in percdata[fi.name].keys()])
    y = [percdata[fi.name][str(_)] for _ in x]
    pltr.plot(x,y,fi,ls='solid',label=fi.name,ms=10)

# Finish off plot.
plt.ylabel('% Actual Best Reference',fontsize=20)
plt.xlabel('Number of Quick Search Iterations',fontsize=20)
plotr.save()
figf.write(percdata)

# Plot portion of dataset vs. average deltaRRMSD.
plotrr = pltr.new(figName,'delta')
figf   = pltr.writeData(figName,'delta')
plt.ylabel('$\Delta$ RRMSD',fontsize=20)
plt.xlabel('Average Portion of Dataset Assessed',fontsize=20)
for key in partdata:
    x,y  = [],[]
    rang = sorted([int(_) for _ in diffdata[key].keys()])
    for n in rang:
        _  = abs(diffdata[key][str(n)])
        xv = min(1,n*partdata[key])
        if xv in x:
            i  = x.index(xv)
            yv = y[i]
            if _ < yv:
                y[i] = _
                continue
        if xv == 1: _ = 0
        x.append(xv)
        y.append(_ )
    pltr.plot(x,y,realdata[key],ls='solid',label=key,ms=10)
plotrr.save()
figf.write(diffdata)

# Plot #/iterations vs. #/structures.
plotrr = pltr.new(figName,'nums')
figf   = pltr.writeData(figName,'nums')
plt.ylabel('Avg. Number of Structures in Dataset')
plt.xlabel('Number of Quick Search Iterations')
for key in numsdata:
    x,y = [],[]
    for n in numsdata[key]:
        rang = sorted([int(_) for _ in percdata[key]])
        for i in rang:
            _   = '%s_q%s.xml' % (n,str(i))
            _   = os.path.join(realdata[key].path,_)
            if not os.path.isfile(_):
                print '\tWARNING (%s): File %s was not found.' % (key,_)
                continue
            dat = aRun(_)
            if dat.reference and refsdata[key][n] in dat.reference:
                x.append(i)
                y.append(numsdata[key][n])
                break
    pltr.plot(x,y,realdata[key],ls='none',ms=10)
plotrr.save()
figf.write(numsdata)

''' [5] Other alignment implementations also work. '''

figName = '5_other_aligners'
plot    = pltr.new(figName)
datf    = pltr.writeData(figName)
print '>>> Gathering data and plotting figure [%s].' % (figName)

graphset = {}
for fi in data.ssets:
    
    # Gather all RRMSD data.
    if not fi.aligner in graphset:
        graphset[fi.aligner] = []
    g = graphset[fi.aligner]
    for fld in fi.folders:
        for ref in fld.refs:
            r = [x for x in ref.scores.values() if x]
            g.extend(r)
    
# Project data onto histogram (per aligner).
g = graphset
_ = reversed(sorted(g,key=lambda d: len(g[d])))
for ali in _:
    weights = numpy.ones_like(g[ali])/len(g[ali])
    plt.hist(g[ali],bins=numpy.arange(0,1.,0.025),
             color=colrs.__clr__(ali),
             weights=weights)

# Save file(s).
plt.ylabel('Density')
plt.xlabel('RRMSD')
plot.save()
datf.write(graphset)
