#!/bin/python
from array import array
from collections import defaultdict
from math import sqrt
import glob as gb
import ROOT 

##### FUNCTIONS #####

def makeHist(filename,binEdges):
  Hist = ROOT.TH1D(filename,filename,len(binEdges)-1,array('f',binEdges)) 
  return Hist

def normaliseHist(h):
	normalisationFactor = h.Integral()
	h.Scale(1.,"width")
	h.Scale(1./normalisationFactor)

def fillHist(h):
	for counter, theLine in enumerate(ifile.readlines()):
		strippedLine = theLine.strip()
		tokenizedLine = strippedLine.split(" ")
		binEdge = float(tokenizedLine[0])
		binContent = float(tokenizedLine[1])
		binNumber = h.FindFixBin(binEdge)
		extBinContent = h.GetBinContent(binNumber)
		h.SetBinContent(binNumber,extBinContent+binContent)
		h.SetBinError(binNumber,sqrt(extBinContent+binContent))

def fillWeightedHist(filename,h,chiLogBins,dKey1):
	chiVal,minvVal,weightVal = ([] for i in range(3))
	with open(filename,"r") as inputfile:	
		#lineNo = 0
		for line in inputfile:
			#lineNo += 1
			columns = [float(x) for x in line.split()]
			chi,minv,weight = columns[2],columns[3],columns[4]
			chiVal.append(chi)
			minvVal.append(minv)
			weightVal.append(weight)
			#print lineNo,chi,minv,weight
	#h.Sumw2()
	sumOfWeights = 0
	lowerBound = float(dKey1.split("_")[2])
	upperBound = float(dKey1.split("_")[3])
	print "lower bound: %d , upper bound: %d" %(lowerBound,upperBound)
	#print "BinNo, lowerChiBound, upperChiBound, lineNo, mInv, chi, weight, sum of weights:"
	for k in range(0,len(chiLogBins)-1):
		for j in range(0,len(weightVal)):
			if (minvVal[j] >= lowerBound and minvVal[j] < upperBound and chiVal[j] >= chiLogBins[k] and chiVal[j] < chiLogBins[k+1]):
				#h.Fill(chiVal[j],weightVal[j])
				sumOfWeights += weightVal[j]
				#print k,chiLogBins[k-1],chiLogBins[k],j,minvVal[j],chiVal[j],weightVal[j],sumOfWeights
		print k,sumOfWeights
		h.SetBinContent(k+1,sumOfWeights)
		h.SetBinError(k+1,sqrt(sumOfWeights))
		sumOfWeights = 0
	print "Integral before normalisation: %f" %h.Integral()
	print "Bin contents before normalisation:"
	for k in range(0,len(chiLogBins)):
		print k,h.GetBinContent(k)
	normaliseHist(h)
	print "Integral after normalisation: %f" %h.Integral()
	print "Bin contents after normalisation:"
	for k in range(0,len(chiLogBins)):
		print k,h.GetBinContent(k)

def fillDict(filetype,binEdges,wValuesFile):
	d,dw,dr = (defaultdict(dict) for i in range(3))
	for elmt,filename in enumerate(sorted(gb.iglob(filetype)),start=1):
		ifile = open(filename)
		splitFilename = filename.split("_")
		dKey1 = splitFilename[0]+"_"+splitFilename[1]+"_"+splitFilename[2]+"_"+splitFilename[3]
		dKey2 = splitFilename[4]
		print elmt,filename,dKey1,dKey2
		h = makeHist(filename,binEdges)
		print h
		hw = h.Clone()
		hr = h.Clone()
		for counter, theLine in enumerate(ifile.readlines()):
			strippedLine = theLine.strip()
			tokenizedLine = strippedLine.split(" ")
			binEdge = float(tokenizedLine[0])
			binContent = float(tokenizedLine[1])
			binNumber = h.FindFixBin(binEdge)
			extBinContent = h.GetBinContent(binNumber)
			h.SetBinContent(binNumber,extBinContent+binContent)
			h.SetBinError(binNumber,sqrt(extBinContent+binContent))
		print "Unweighted bin contents before normalisation:"
		for k in range(0,len(binEdges)):
			print k,h.GetBinContent(k)
		normaliseHist(h)
		print "Unweighted bin contents after normalisation:"
		for k in range(0,len(binEdges)):
			print k,h.GetBinContent(k)
		print "Printout for unweighted histogram %s:" %dKey1
		h.Print("all")
		d[dKey1][dKey2] = h
		fillWeightedHist(wValuesFile,hw,binEdges,dKey1)
		print "Printout for weighted histogram %s:" %dKey1
		hw.Print("all")
		dw[dKey1][dKey2] = hw
		for k in range(0,len(binEdges)):
			if (hw.GetBinContent(k) > 0):
				hr.SetBinContent(k,(hw.GetBinContent(k)/h.GetBinContent(k)))
			else:
				hr.SetBinContent(k,0)
		print "Printout for ratio histogram (weighted/unweighted):"
		hr.Print("all")
		dr[dKey1][dKey2] = hr
	return d,dw,dr

def drawHist(alphaOrder,c,d,colourN,styleN,L,label):
	for keyN, dKey in enumerate(sorted(d),start=1):
		print keyN,dKey
		h = d[dKey][alphaOrder[0]]
		print h
		h.SetMarkerColor(colourN)
		h.SetLineColor(colourN)
		h.SetMarkerStyle(styleN)
		h.SetMarkerSize(0.5)
		h.GetXaxis().SetLabelSize(0.06)
		h.GetYaxis().SetLabelSize(0.06)
		h.SetTitleSize(0.05)
		c.cd(keyN)
		if ("ratio" not in label):
			ROOT.gPad.cd(1)
			ROOT.gPad.SetPad(0.,0.4,1.,1.)
			ROOT.gPad.SetLogx()
			if (keyN == 1):
				h.SetTitle("#bf{2500 GeV < m_{jj} < 2800 GeV}")
				h.SetMaximum(0.045)
				h.SetMinimum(0.025)
			if (keyN == 2):
				h.SetTitle("#bf{2800 GeV < m_{jj} < 3100 GeV}")
				h.SetMaximum(0.045)
				h.SetMinimum(0.020)
			if (keyN == 3):
				h.SetTitle("#bf{3100 GeV < m_{jj} < 3400 GeV}")
				h.SetMaximum(0.045)
				h.SetMinimum(0.020)
			if (keyN == 4):
				h.SetTitle("#bf{3400 GeV < m_{jj} < 4000 GeV}")
				h.SetMaximum(0.045)
				h.SetMinimum(0.020)
			if (keyN == 5):
				h.SetTitle("#bf{4000 GeV < m_{jj} < 4600 GeV}")
				h.SetMaximum(0.045)
				h.SetMinimum(0.020)
			if (keyN == 6):
				h.SetTitle("#bf{4600 GeV < m_{jj} < 5400 GeV}")
				h.SetMaximum(0.050)
				h.SetMinimum(0.010)
			if (keyN == 7):
				h.SetTitle("#bf{5400 GeV < m_{jj}}")
				h.SetMaximum(0.055)
				h.SetMinimum(0.005)
			h.GetXaxis().SetTitle("#chi")
			h.GetYaxis().SetTitle("#frac{1}{N} #frac{dN}{d#chi}")
			h.Draw("same")
		else:
			ROOT.gPad.cd(2)
			ROOT.gPad.SetPad(0.,0.,1.,0.4)
			ROOT.gPad.SetLogx()
			h.SetTitle("Ratio of trijet alpha_s (order) = 1, PhaseSpace:RsepMin = 1.0 (weighted/unweighted)")
			h.GetXaxis().SetTitle("#chi")
			h.GetYaxis().SetTitle("R")
			h.Draw("same")
			neutralityLine.Draw("same")
	L.AddEntry(h,label)

##### MAIN CODE #####

neutralityLine = ROOT.TF1("neutral line","1",0.,30.)
neutralityLine.SetLineColor(16)
neutralityLine.SetLineStyle(2)
neutralityLine.SetLineWidth(1)

kFacFile = ROOT.TFile("kfacOriginal.root")
if kFacFile.IsOpen():
	print("File %s opened successfully.\n") %(kFacFile)

h_2500_2800 = ROOT.TH1D(kFacFile.Get("kfacSmooth_2500Mjj2800"))
h_2800_3100 = ROOT.TH1D(kFacFile.Get("kfacSmooth_2800Mjj3100"))
h_3100_3400 = ROOT.TH1D(kFacFile.Get("kfacSmooth_3100Mjj3400"))
h_3400_3700 = ROOT.TH1D(kFacFile.Get("kfacSmooth_3400Mjj3700"))
#h_3700_4000 = ROOT.TH1D(kFacFile.Get("kfacSmooth_3700Mjj4000"))
h_4000_4300 = ROOT.TH1D(kFacFile.Get("kfacSmooth_4000Mjj4300"))
#h_4300_4600 = ROOT.TH1D(kFacFile.Get("kfacSmooth_4300Mjj4600"))
h_4600_4900 = ROOT.TH1D(kFacFile.Get("kfacSmooth_4600Mjj4900"))
#h_4900_5400 = ROOT.TH1D(kFacFile.Get("kfacSmooth_4900Mjj5400"))
h_5400_6500 = ROOT.TH1D(kFacFile.Get("kfacSmooth_5400Mjj6500"))

khist = [h_2500_2800,h_2800_3100,h_3100_3400,h_3400_3700,h_4000_4300,h_4600_4900,h_5400_6500]
#khist = [h_2500_2800,h_2800_3100,h_3100_3400,h_3400_3700,h_4000_4300,h_4300_4600,h_4600_4900,h_4900_5400,h_5400_6500]

alphaOrder = ["alpha1.txt"]

canv = ROOT.TCanvas("","",3200,2400)
canv.Divide(2,4)
for i in range(1,8):
	canv.cd(i)
	ROOT.gPad.Divide(1,2)

Ddi,DdiNLO = (defaultdict(dict) for i in range(2))
DtriR01,DtriR01w,DtriR01r,DtriR02,DtriR02w,DtriR02r,DtriR03,DtriR03w,DtriR03r,DtriR04,DtriR04w,DtriR04r,DtriR05,DtriR05w,DtriR05r,DtriR06,DtriR06w,DtriR06r,DtriR07,DtriR07w,DtriR07r,DtriR08,DtriR08w,DtriR08r,DtriR09,DtriR09w,DtriR09r,DtriR10,DtriR10w,DtriR10r = (defaultdict(dict) for i in range(30))
chiLogBins = [1.,1.34986,1.82212,2.4596,3.32012,4.48169,6.04965,8.16617,11.0232,14.8797,20.0855,30.]

for elmt, filename in enumerate(sorted(gb.iglob("Dijet/dijet_chi*alpha1.txt")), start=1):
	ifile = open(filename)
	splitFilename = filename.split("_")
	chidictKey1 = splitFilename[0]+"_"+splitFilename[1]+"_"+splitFilename[2]+"_"+splitFilename[3]
	chidictKey2 = splitFilename[4]
	print(elmt, filename, chidictKey1, chidictKey2)
	h = makeHist(filename,chiLogBins)
	fillHist(h)
	hNLO = h.Clone()
	normaliseHist(h)
	Ddi[chidictKey1][chidictKey2] = h
	if (chidictKey1 == "Dijet/dijet_chi_2500_2800"):
		for i in range(0,len(chiLogBins)):
			hNLO.SetBinContent(i,(hNLO.GetBinContent(i))*(khist[0].GetBinContent(i)))
	if (chidictKey1 == "Dijet/dijet_chi_2800_3100"):
		for i in range(0,len(chiLogBins)):
			hNLO.SetBinContent(i,(hNLO.GetBinContent(i))*(khist[1].GetBinContent(i)))
	if (chidictKey1 == "Dijet/dijet_chi_3100_3400"):
		for i in range(0,len(chiLogBins)):
			hNLO.SetBinContent(i,(hNLO.GetBinContent(i))*(khist[2].GetBinContent(i)))
	if (chidictKey1 == "Dijet/dijet_chi_3400_4000"):
		for i in range(0,len(chiLogBins)):
			hNLO.SetBinContent(i,(hNLO.GetBinContent(i))*(khist[3].GetBinContent(i)))
	if (chidictKey1 == "Dijet/dijet_chi_4000_4600"):
		for i in range(0,len(chiLogBins)):
			hNLO.SetBinContent(i,(hNLO.GetBinContent(i))*(khist[4].GetBinContent(i)))
	if (chidictKey1 == "Dijet/dijet_chi_4600_5400"):
		for i in range(0,len(chiLogBins)):
			hNLO.SetBinContent(i,(hNLO.GetBinContent(i))*(khist[5].GetBinContent(i)))
	if (chidictKey1 == "Dijet/dijet_chi_5400_13000"):
		for i in range(0,len(chiLogBins)):
			hNLO.SetBinContent(i,(hNLO.GetBinContent(i))*(khist[6].GetBinContent(i)))
	normaliseHist(hNLO)
	DdiNLO[chidictKey1][chidictKey2] = hNLO

#DtriR01,DtriR01w = fillDict("RsepMin01m2/trijet_chi*alpha1.txt",chiLogBins,"RsepMin01m2/trijet01R01m2rV.txt")
#DtriR02,DtriR02w = fillDict("RsepMin02m2/trijet_chi*alpha1.txt",chiLogBins,"RsepMin02m2/trijet01R02m2rV.txt")
#DtriR03,DtriR03w = fillDict("RsepMin03m2/trijet_chi*alpha1.txt",chiLogBins,"RsepMin03m2/trijet01R03m2rV.txt")
#DtriR04,DtriR04w = fillDict("RsepMin04m2/trijet_chi*alpha1.txt",chiLogBins,"RsepMin04m2/trijet01R04m2rV.txt")
#DtriR05,DtriR05w = fillDict("RsepMin05m2/trijet_chi*alpha1.txt",chiLogBins,"RsepMin05m2/trijet01R05m2rV.txt")
#DtriR06,DtriR06w = fillDict("RsepMin06m2/trijet_chi*alpha1.txt",chiLogBins,"RsepMin06m2/trijet01R06m2rV.txt")
#DtriR07,DtriR07w = fillDict("RsepMin07m2/trijet_chi*alpha1.txt",chiLogBins,"RsepMin07m2/trijet01R07m2rV.txt")
#DtriR08,DtriR08w = fillDict("RsepMin08m2/trijet_chi*alpha1.txt",chiLogBins,"RsepMin08m2/trijet01R08m2rV.txt")
#DtriR09,DtriR09w = fillDict("RsepMin09m2/trijet_chi*alpha1.txt",chiLogBins,"RsepMin09m2/trijet01R09m2rV.txt")
DtriR10,DtriR10w,DtriR10r = fillDict("RsepMin10m2/trijet_chi*alpha1.txt",chiLogBins,"RsepMin10m2/trijet01R10m2rV.txt")

L = ROOT.TLegend(0.0,0.0,1.0,1.0)
drawHist(alphaOrder,canv,Ddi,4,33,L,"dijet alpha_s (order) = 1 (@LO)")
drawHist(alphaOrder,canv,DdiNLO,38,27,L,"dijet alpha_s (order) = 1 (@LO + NLO)")
#drawHist(alphaOrder,canv,DtriR01,2,20,L,"trijet alpha_s (order) = 1 (@LO), PhaseSpace:RsepMin = 0.1")
drawHist(alphaOrder,canv,DtriR01w,46,24,L,"trijet alpha_s (order) = 1 (@LO, weighted), PhaseSpace:RsepMin = 0.1")
#drawHist(alphaOrder,canv,DtriR02w,46,24,L,"trijet alpha_s (order) = 1 (@LO, weighted), PhaseSpace:RsepMin = 0.2")
#drawHist(alphaOrder,canv,DtriR04,3,21,L,"trijet alpha_s (order) = 1 (@LO), PhaseSpace:RsepMin = 0.4")
#drawHist(alphaOrder,canv,DtriR04w,30,25,L,"trijet alpha_s (order) = 1 (@LO, weighted), PhaseSpace:RsepMin = 0.4")
#drawHist(alphaOrder,canv,DtriR06w,13,27,L,"trijet alpha_s (order) = 1 (@LO, weighted), PhaseSpace:RsepMin = 0.6")
#drawHist(alphaOrder,canv,DtriR07,4,22,L,"trijet alpha_s (order) = 1 (@LO), PhaseSpace:RsepMin = 0.7")
#drawHist(alphaOrder,canv,DtriR07w,38,26,L,"trijet alpha_s (order) = 1 (@LO, weighted), PhaseSpace:RsepMin = 0.7")
#drawHist(alphaOrder,canv,DtriR08w,38,26,L,"trijet alpha_s (order) = 1 (@LO, weighted), PhaseSpace:RsepMin = 0.8")
#drawHist(alphaOrder,canv,DtriR10,28,23,L,"trijet alpha_s (order) = 1 (@LO), PhaseSpace:RsepMin = 1.0")
drawHist(alphaOrder,canv,DtriR10w,26,32,L,"trijet alpha_s (order) = 1 (@LO, weighted), PhaseSpace:RsepMin = 1.0")
#drawHist(alphaOrder,canv,DtriR10r,1,20,L,"ratio of trijet alpha_s (order) = 1, PhaseSpace:RsepMin = 1.0 (weighted/unweighted)")
canv.cd(8)
L.Draw("same")
canv.SaveAs("Results/ChiPlots/5.5.3.1chiN_diNLO_tri10RW_ratiodiNLOvstriW.pdf")

