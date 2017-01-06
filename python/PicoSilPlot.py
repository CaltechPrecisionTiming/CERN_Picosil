from optparse import OptionParser
import ROOT as rt
from array import array
import numpy as np
from math import sqrt
import sys
import json

if __name__ == '__main__':
    parser = OptionParser()

    parser.add_option('-c','--cut',dest="cut",type="string",default='1',
                  help="cut to apply before averaging, e.g. to look only at event 101, do event==101; ring0int = 1 central hexagon, ring1int = 6 central hexagons, ring2int = 18 central hexagons")
    parser.add_option('-w','--wire-mapping',dest="wireMapping",type="string",default='wireMapping_06182016_1230.txt',
                  help="wire mapping text file")
    parser.add_option('-p','--plot',dest="plot",type="choice",action='store',default='amp',choices=['amp','int','intfull'],
                  help="plot varible on the z-axis: amp or int or intfull")
    parser.add_option('-d','--out-dir',dest="outDir",type="string",default='./',
                  help="output directory to store ")
    
    (options,args) = parser.parse_args()

    f = rt.TFile.Open(args[0])

    tree = f.Get("pulse")

    empty = rt.TH2D("empty","empty",100, -6, 6, 100, -6, 6)
    
    h2p = rt.TH2Poly()
    h2p.SetName("picosil")
    h2p.SetTitle("picosil")

    xVertex = array('d',[1, 0.5, -0.5, -1, -0.5, 0.5])
    yVertex = array('d',[0, -sqrt(3.)/2.,-sqrt(3.)/2., 0, sqrt(3.)/2., sqrt(3.)/2.])

    centerPositions = [[-3, -sqrt(3)], [-3, 0], [-3, sqrt(3)],
                       [-1.5, -5*sqrt(3)/2], [-1.5, -3*sqrt(3)/2], [-1.5, -sqrt(3)/2], [-1.5, sqrt(3)/2], [-1.5, 3*sqrt(3)/2], [-1.5, 5*sqrt(3)/2],
                       [0, -3*sqrt(3)], [0, -2*sqrt(3)], [0, -1*sqrt(3)],
                       [0, 0],
                       [0, sqrt(3)], [0, 2*sqrt(3)], [0, 3*sqrt(3)],
                       [1.5, -5*sqrt(3)/2], [1.5, -3*sqrt(3)/2], [1.5, -sqrt(3)/2], [1.5, sqrt(3)/2], [1.5, 3*sqrt(3)/2], [1.5, 5*sqrt(3)/2],
                       [3, -sqrt(3)], [3, 0], [3, sqrt(3)]
                       ]

    wireMappingFile = open(options.wireMapping,'r')
    indexValues = []
    attValues = {}
    for l in wireMappingFile.readlines():
        if l[0]!='#':
            l = l.replace('\n','')
            lArray = l.split(' ')
            indexValues.append(lArray[0])
            if len(lArray)>1:
                level = float(lArray[1].replace('db',''))
                att = pow(10.,-level/20.)
                attValues[lArray[0]] = 1.0/att
            else:
                attValues[lArray[0]] = 1.0

    mapArrayToCenterPos = {
        indexValues[0]: [-3, sqrt(3)],      
        indexValues[1]: [-3, 0],   
        indexValues[2]: [-3, -sqrt(3)], 
        indexValues[3]: [-1.5, 5*sqrt(3)/2],   
        indexValues[4]: [-1.5, 3*sqrt(3)/2],   
        indexValues[5]: [-1.5, sqrt(3)/2],
        indexValues[6]: [-1.5, -sqrt(3)/2],  
        indexValues[7]: [-1.5, -3*sqrt(3)/2],       
        indexValues[8]: [-1.5, -5*sqrt(3)/2],  
        indexValues[9]: [0, 3*sqrt(3)],        
        indexValues[10]: [0, 2*sqrt(3)],   
        indexValues[11]: [0, sqrt(3)],      
        indexValues[12]: [0, 0],        
        indexValues[13]: [0, -sqrt(3)],      
        indexValues[14]: [0, -2*sqrt(3)],   
        indexValues[15]: [0, -3*sqrt(3)],        
        indexValues[16]: [1.5, 5*sqrt(3)/2],   
        indexValues[17]: [1.5, 3*sqrt(3)/2], 
        indexValues[18]: [1.5, sqrt(3)/2],       
        indexValues[19]: [1.5, -sqrt(3)/2],     
        indexValues[20]: [1.5, -3*sqrt(3)/2], 
        indexValues[21]: [1.5, -5*sqrt(3)/2],       
        indexValues[22]: [3, sqrt(3)],   
        indexValues[23]: [3, 0],
        indexValues[24]: [3, -sqrt(3)]
        }

    ring0Index = []
    ring1Index = []
    ring2Index = []
    for key, val in mapArrayToCenterPos.iteritems():
        if key == '': continue
        if np.linalg.norm(np.array(val)) <= 0: ring0Index.append(key)
        if np.linalg.norm(np.array(val)) <= sqrt(3): ring1Index.append(key)
        if np.linalg.norm(np.array(val)) <= 2*sqrt(3): ring2Index.append(key)


    for centerPos in centerPositions:
        x2 = array('d',[xi+centerPos[0] for xi in xVertex])
        y2 = array('d',[yi+centerPos[1] for yi in yVertex])
        h2p.AddBin(6,x2,y2)

                        
    options.cut = options.cut.replace('ring0intfull','+'.join(['%f*intfull[%s]'%(attValues[index],index) for index in ring0Index]))
    options.cut = options.cut.replace('ring1intfull','+'.join(['%f*intfull[%s]'%(attValues[index],index) for index in ring1Index]))
    options.cut = options.cut.replace('ring2intfull','+'.join(['%f*intfull[%s]'%(attValues[index],index) for index in ring2Index]))
    
    options.cut = options.cut.replace('ring0int','+'.join(['%f*int[%s]'%(attValues[index],index) for index in ring0Index]))
    options.cut = options.cut.replace('ring1int','+'.join(['%f*int[%s]'%(attValues[index],index) for index in ring1Index]))
    options.cut = options.cut.replace('ring2int','+'.join(['%f*int[%s]'%(attValues[index],index) for index in ring2Index]))
                                        
    options.cut = options.cut.replace('ring0amp','+'.join(['%f*amp[%s]'%(attValues[index],index) for index in ring0Index]))
    options.cut = options.cut.replace('ring1amp','+'.join(['%f*amp[%s]'%(attValues[index],index) for index in ring1Index]))
    options.cut = options.cut.replace('ring2amp','+'.join(['%f*amp[%s]'%(attValues[index],index) for index in ring2Index]))
                                    
    nevents = tree.Draw('>>elist',options.cut,'entrylist')
        
    elist = rt.gDirectory.Get('elist')
    
    entry = -1

    ring0sum = 0
    ring1sum = 0
    ring2sum = 0
    while True:
        entry = elist.Next()
        if entry == -1: break
        tree.GetEntry(entry)
        for key, val in mapArrayToCenterPos.iteritems():
            if key == '': continue
            h2p.Fill(val[0],val[1], eval('%f*tree.%s[%s]'%(attValues[key],options.plot,key)))
            
        ring0sum += eval('+'.join(['max(%f*tree.%s[%s],0)'%(attValues[index],options.plot,index) for index in ring0Index]))
        ring1sum += eval('+'.join(['max(%f*tree.%s[%s],0)'%(attValues[index],options.plot,index) for index in ring1Index]))
        ring2sum += eval('+'.join(['max(%f*tree.%s[%s],0)'%(attValues[index],options.plot,index) for index in ring2Index]))
            
        #print tree.event
        
    h2p.Scale(1.0/nevents)
    ring0sum/=nevents
    ring1sum/=nevents
    ring2sum/=nevents
    
    rt.gStyle.SetOptStat(0)
    rt.gStyle.SetOptTitle(0)
    rt.gStyle.SetPaintTextFormat("+3.3f")

    c = rt.TCanvas('c','c',500,500)
    c.SetRightMargin(0.2)
        
    empty.Draw("")
    if options.plot=='amp':
        h2p.GetZaxis().SetTitle("Amplitude [V]")
    elif options.plot=='int' or options.plot=='intfull':
        h2p.GetZaxis().SetTitle("Charge [pC]")
    h2p.GetZaxis().SetTitleOffset(2)
    if options.plot=='amp':
        h2p.SetMaximum(0.1)
    else:
        h2p.SetMaximum(4)
    h2p.Draw("colztextsame")
    
    l = rt.TLatex()
    l.SetTextAlign(22)
    l.SetTextSize(0.03)
    l.SetTextFont(42)
    l.SetNDC()
    l.DrawLatex(0.5,0.94,"cut: %s #rightarrow %i events"%(options.cut,nevents))
    
    l.SetTextAlign(12)
    l.DrawLatex(0.6,0.85,"ring0 sum: %.3f"%(ring0sum))
    l.DrawLatex(0.6,0.80,"ring1 sum: %.3f"%(ring1sum))
    l.DrawLatex(0.6,0.75,"ring2 sum: %.3f"%(ring2sum))

    print "ring0 sum: %.3f"%(ring0sum)
    print "ring1 sum: %.3f"%(ring1sum)
    print "ring2 sum: %.3f"%(ring2sum)
    
    c.Print('%s/picosil_%s.pdf'%(options.outDir,args[0].split('/')[-1].replace('.root','')))
    c.Print('%s/picosil_%s.png'%(options.outDir,args[0].split('/')[-1].replace('.root','')))
    c.Print('%s/picosil_%s.C'%(options.outDir,args[0].split('/')[-1].replace('.root','')))
    
    f.Close()
