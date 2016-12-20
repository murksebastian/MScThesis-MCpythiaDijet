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

def fillDict(filetype):
	DinvM,Dphi,DpT,Dy,DyB,DyS = (defaultdict(dict) for i in range(6))
	for elmt,filename in enumerate(sorted(gb.iglob(filetype)),start=1):
		ifile = open(filename)
		splitFilename = filename.split("_")
		if ("invMass" in filename) or ("yBoost" in filename) or ("yStar" in filename):
			dKey1 = splitFilename[0]+"_"+splitFilename[1]
			dKey2 = splitFilename[2]
		else:
			dKey1 = splitFilename[0]+"_"+splitFilename[1]+"_"+splitFilename[2]
			dKey2 = splitFilename[3]
		binEdges,binContents = ([] for i in range(2))
		for counter, theLine in enumerate(ifile.readlines()):
			strippedLine = theLine.strip()
			tokenizedLine = strippedLine.split(" ")
			binEdges.append(float(tokenizedLine[0]))
			binContents.append(float(tokenizedLine[1]))
		h = makeHist(filename,binEdges)
		print elmt,filename,dKey1,dKey2
		print h
		ifile.seek(0)
		for counter, theLine in enumerate(ifile.readlines()):
			strippedLine = theLine.strip()
			tokenizedLine = strippedLine.split(" ")
			binEdge = float(tokenizedLine[0])
			binContent = float(tokenizedLine[1])
			binNumber = h.FindFixBin(binEdge)
			#extBinContent = h.GetBinContent(binNumber)
			h.SetBinContent(binNumber,binContent)
			h.SetBinError(binNumber,sqrt(binContent))
		if (splitFilename[1] == "invMass"):
			DinvM[dKey1][dKey2] = h
		if (splitFilename[2] == "phi"):
			Dphi[dKey1][dKey2] = h
		if (splitFilename[2] == "pT"):
			DpT[dKey1][dKey2] = h
		if (splitFilename[2] == "y"):
			Dy[dKey1][dKey2] = h
		if (splitFilename[1] == "yBoost"):
			DyB[dKey1][dKey2] = h
		if (splitFilename[1] == "yStar"):
			DyS[dKey1][dKey2] = h
		h.Print("all")
	return DinvM,Dphi,DpT,Dy,DyB,DyS

def drawHist(alphaOrder,c,D,colourN,styleN,L,label,ofilename):
	for keyN,dKey in enumerate(sorted(D),start=1):
		print keyN,dKey
		c.cd(keyN)
		if ("pT" in dKey):
			ROOT.gPad.SetLogx()
		if ("= 0" in label):
			h = D[dKey][alphaOrder[0]]
			h.SetLineStyle(2)
		if ("= 1" in label):
			h = D[dKey][alphaOrder[1]]
			h.SetLineStyle(1)
		if ("= 2" in label):
			h = D[dKey][alphaOrder[2]]
			h.SetLineStyle(3)
		print h
		h.SetMarkerColor(colourN)
		h.SetLineColor(colourN)
		h.SetMarkerStyle(styleN)
		h.SetMarkerSize(0.5)
		h.SetLineWidth(1)
		#h.GetXaxis().SetLabelSize(0.02)
		#h.GetYaxis().SetLabelSize(0.02)
		h.Draw("same")
		if ("invM" in dKey):
			h.SetTitle("dijet vs. trijet invariant mass distribution")
			xAxisLabel = ROOT.TLatex()
			xAxisLabel.SetTextAlign(22)
			xAxisLabel.SetTextFont(43)
			xAxisLabel.SetTextSize(15)
			xAxisLabel.SetNDC()
			xAxisLabel.DrawLatex(0.92,0.06,"m_{jj} [GeV]")
			yAxisLabel = ROOT.TLatex()
			yAxisLabel.SetTextAlign(22)
			yAxisLabel.SetTextFont(43)
			yAxisLabel.SetTextSize(15)
			yAxisLabel.SetNDC()
			yAxisLabel.DrawLatex(0.028,0.90,"N")
		if ("phi" in dKey):
			h.GetXaxis().SetTitle("#phi")
			h.GetYaxis().SetTitle("N")
		if ("pT" in dKey):
			h.GetXaxis().SetTitle("p_{T}")
			h.GetYaxis().SetTitle("N")
		if ("ing_y" in dKey):
			h.GetXaxis().SetTitle("y")
			h.GetYaxis().SetTitle("N")
		if ("yBoost" in dKey):
			h.SetTitle("dijet vs. trijet y_{B}")
			xAxisLabel = ROOT.TLatex()
			xAxisLabel.SetTextAlign(22)
			xAxisLabel.SetTextFont(43)
			xAxisLabel.SetTextSize(15)
			xAxisLabel.SetNDC()
			xAxisLabel.DrawLatex(0.93,0.06,"y_{B}")
			yAxisLabel = ROOT.TLatex()
			yAxisLabel.SetTextAlign(22)
			yAxisLabel.SetTextFont(43)
			yAxisLabel.SetTextSize(15)
			yAxisLabel.SetNDC()
			yAxisLabel.DrawLatex(0.028,0.90,"N")
		if ("yStar" in dKey):
			h.SetTitle("dijet vs. trijet y^{*}")
			xAxisLabel = ROOT.TLatex()
			xAxisLabel.SetTextAlign(22)
			xAxisLabel.SetTextFont(43)
			xAxisLabel.SetTextSize(15)
			xAxisLabel.SetNDC()
			xAxisLabel.DrawLatex(0.93,0.06,"y^{*}")
			yAxisLabel = ROOT.TLatex()
			yAxisLabel.SetTextAlign(22)
			yAxisLabel.SetTextFont(43)
			yAxisLabel.SetTextSize(15)
			yAxisLabel.SetNDC()
			yAxisLabel.DrawLatex(0.028,0.90,"N")
	L.AddEntry(h,label)
	if ("invM" in dKey or "yBoost" in dKey or "yStar" in dKey):
		c.cd(keyN)
		L.Draw("same")
	else:
		c.cd(keyN+1)
		L.Draw("same")
	c.SaveAs(ofilename)

##### MAIN CODE #####

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

alphaOrder = ["alpha0.txt","alpha1.txt","alpha2.txt"]

canvinvM,canvphi,canvpT,canvy,canvyB,canvyS = (ROOT.TCanvas() for i in range(6))
canvinvM.Divide(2,1)
canvphi.Divide(2,2)
canvpT.Divide(2,2)
canvy.Divide(2,2)
canvyB.Divide(2,1)
canvyS.Divide(2,1)

#di0DinvM,di0Dphi,di0DpT,di0Dy,di0DyB,di0DyS = fillDict("Dijet/dijet_[ilsy]*alpha0.txt")
di1DinvM,di1Dphi,di1DpT,di1Dy,di1DyB,di1DyS = fillDict("Dijet/dijet_[ilsy]*alpha1.txt")
#di2DinvM,di2Dphi,di2DpT,di2Dy,di2DyB,di2DyS = fillDict("Dijet/dijet_[ilsy]*alpha2.txt")
tri1m2DinvM,tri1m2Dphi,tri1m2DpT,tri1m2Dy,tri1m2DyB,tri1m2DyS = fillDict("RsepMin10m2/trijet_[ilsy]*alpha1.txt")
tri1m3DinvM,tri1m3Dphi,tri1m3DpT,tri1m3Dy,tri1m3DyB,tri1m3DyS = fillDict("RsepMin10m3/trijet_[ilsy]*alpha1.txt")

LinvM,Lphi,LpT,Ly,LyB,LyS = (ROOT.TLegend(0.63,0.60,0.98,0.75) for i in range(6))

#drawHist(alphaOrder,canvinvM,di0DinvM,4,21,LinvM,"dijet #Omicron(#alpha_{s}) = 0 (@LO)","Results/KinVarPlots/invM_di0.pdf")
#drawHist(alphaOrder,canvphi,di0Dphi,4,21,Lphi,"dijet #Omicron(#alpha_{s}) = 0 (@LO)","Results/KinVarPlots/phi_di0.pdf")
#drawHist(alphaOrder,canvpT,di0DpT,4,21,LpT,"dijet #Omicron(#alpha_{s}) = 0 (@LO)","Results/KinVarPlots/pT_di0.pdf")
#drawHist(alphaOrder,canvy,di0Dy,4,21,Ly,"dijet #Omicron(#alpha_{s}) = 0 (@LO)","Results/KinVarPlots/y_di0.pdf")
#drawHist(alphaOrder,canvyB,di0DyB,4,21,LyB,"dijet #Omicron(#alpha_{s}) = 0 (@LO)","Results/KinVarPlots/yB_di0.pdf")
#drawHist(alphaOrder,canvyS,di0DyS,4,21,LyS,"dijet #Omicron(#alpha_{s}) = 0 (@LO)","Results/KinVarPlots/yS_di0.pdf")

drawHist(alphaOrder,canvinvM,di1DinvM,1,20,LinvM,"dijet #Omicron(#alpha_{s}) = 1 (@LO)","Results/KinVarPlots/invM_di1.pdf")
drawHist(alphaOrder,canvphi,di1Dphi,1,20,Lphi,"dijet #Omicron(#alpha_{s}) = 1 (@LO)","Results/KinVarPlots/phi_di1.pdf")
drawHist(alphaOrder,canvpT,di1DpT,1,20,LpT,"dijet #Omicron(#alpha_{s}) = 1 (@LO)","Results/KinVarPlots/pT_di1.pdf")
drawHist(alphaOrder,canvy,di1Dy,1,20,Ly,"dijet #Omicron(#alpha_{s}) = 1 (@LO)","Results/KinVarPlots/y_di1.pdf")
drawHist(alphaOrder,canvyB,di1DyB,1,20,LyB,"dijet #Omicron(#alpha_{s}) = 1 (@LO)","Results/KinVarPlots/yB_di1.pdf")
drawHist(alphaOrder,canvyS,di1DyS,1,20,LyS,"dijet #Omicron(#alpha_{s}) = 1 (@LO)","Results/KinVarPlots/yS_di1.pdf")

#drawHist(alphaOrder,canvinvM,di2DinvM,2,22,LinvM,"dijet #Omicron(#alpha_{s}) = 2 (@LO)","Results/KinVarPlots/invM_di0_di1_di2.pdf")
#drawHist(alphaOrder,canvphi,di2Dphi,2,22,Lphi,"dijet #Omicron(#alpha_{s}) = 2 (@LO)","Results/KinVarPlots/phi_di0_di1_di2.pdf")
#drawHist(alphaOrder,canvpT,di2DpT,2,22,LpT,"dijet #Omicron(#alpha_{s}) = 2 (@LO)","Results/KinVarPlots/pT_di0_di1_di2.pdf")
#drawHist(alphaOrder,canvy,di2Dy,2,22,Ly,"dijet #Omicron(#alpha_{s}) = 2 (@LO)","Results/KinVarPlots/y_di0_di1_di2.pdf")
#drawHist(alphaOrder,canvyB,di2DyB,2,22,LyB,"dijet #Omicron(#alpha_{s}) (order) = 2 (@LO)","Results/KinVarPlots/yB_di0_di1_di2.pdf")
#drawHist(alphaOrder,canvyS,di2DyS,2,22,LyS,"dijet #Omicron(#alpha_{s}) = 2 (@LO)","Results/KinVarPlots/yS_di0_di1_di2.pdf")

drawHist(alphaOrder,canvinvM,tri1m2DinvM,38,24,LinvM,"trijet_{m2} #Omicron(#alpha_{s}) = 1 (@LO)","Results/KinVarPlots/invM_di1_tri1m2.pdf")
drawHist(alphaOrder,canvphi,tri1m2Dphi,38,24,Lphi,"trijet_{m2} #Omicron(#alpha_{s}) = 1 (@LO)","Results/KinVarPlots/phi_di1_tri1m2.pdf")
drawHist(alphaOrder,canvpT,tri1m2DpT,38,24,LpT,"trijet_{m2} #Omicron(#alpha_{s}) = 1 (@LO)","Results/KinVarPlots/pT_di1_tri1m2.pdf")
drawHist(alphaOrder,canvy,tri1m2Dy,38,24,Ly,"trijet_{m2} #Omicron(#alpha_{s}) = 1 (@LO)","Results/KinVarPlots/y_di1_tri1m2.pdf")
drawHist(alphaOrder,canvyB,tri1m2DyB,38,24,LyB,"trijet_{m2} #Omicron(#alpha_{s}) = 1 (@LO)","Results/KinVarPlots/yB_di1_tri1m2.pdf")
drawHist(alphaOrder,canvyS,tri1m2DyS,38,24,LyS,"trijet_{m2} #Omicron(#alpha_{s}) = 1 (@LO)","Results/KinVarPlots/yS_di1_tri1m2.pdf")

drawHist(alphaOrder,canvinvM,tri1m3DinvM,27,21,LinvM,"trijet_{m3} #Omicron(#alpha_{s}) = 1 (@LO)","Results/KinVarPlots/invM_di1_trim2_tri1m3.pdf")
drawHist(alphaOrder,canvphi,tri1m3Dphi,27,21,Lphi,"trijet_{m3} #Omicron(#alpha_{s}) = 1 (@LO)","Results/KinVarPlots/phi_di1_trim2_tri1m3.pdf")
drawHist(alphaOrder,canvpT,tri1m3DpT,27,21,LpT,"trijet_{m3} #Omicron(#alpha_{s}) = 1 (@LO)","Results/KinVarPlots/pT_di1_trim2_tri1m3.pdf")
drawHist(alphaOrder,canvy,tri1m3Dy,27,21,Ly,"trijet_{m3} #Omicron(#alpha_{s}) = 1 (@LO)","Results/KinVarPlots/y_di1_trim2_tri1m3.pdf")
drawHist(alphaOrder,canvyB,tri1m3DyB,27,21,LyB,"trijet_{m3} #Omicron(#alpha_{s}) = 1 (@LO)","Results/KinVarPlots/yB_di1_trim2_tri1m3.pdf")
drawHist(alphaOrder,canvyS,tri1m3DyS,27,21,LyS,"trijet_{m3} #Omicron(#alpha_{s}) = 1 (@LO)","Results/KinVarPlots/yS_di1_trim2_tri1m3.pdf")

