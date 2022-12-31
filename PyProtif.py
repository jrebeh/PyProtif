# Authors: Gilbert El Khoury, Wael Azzam, Joseph Rebehmed
# Supervised by: Joseph Rebehmed
# Last modified: 31/12/2022

######################### functions to install non-installed packages #########################################

import os
import re
import sys
import urllib.request
import subprocess
from pymol import cmd
from pymol.cmd import _cmd, DEFAULT_ERROR, DEFAULT_SUCCESS, _raising, is_ok, is_error, is_list, space_sc, safe_list_eval, is_string, loadable
import time
import platform
import json
import os.path
from os import path
import webbrowser
from datetime import date
import base64
from collections import OrderedDict

def install_package(name):
	subprocess.check_call(["pip", "install", name])

if sys.version_info[0] < 3:
    import Tkinter as tk
    from Tkinter import *
    from tkinter import messagebox
    from tkinter import ttk
    #webdrivermanager
else:
    import tkinter as tk
    from tkinter import *
    from tkinter import messagebox
    from tkinter import ttk

try:
	from selenium import webdriver
	from selenium.webdriver.common.by import By
	import chromedriver_autoinstaller
	from selenium.webdriver.support.ui import WebDriverWait
	from selenium.webdriver.support import expected_conditions as EC
	from selenium.webdriver.common.action_chains import ActionChains
	from selenium.webdriver.common.keys import Keys
	from selenium.webdriver.chrome.service import Service
	from subprocess import CREATE_NO_WINDOW
	chromedriver_autoinstaller.install()
except:
	install_package("selenium")
	install_package("chromedriver-autoinstaller")
	from selenium import webdriver
	from selenium.webdriver.common.by import By
	from selenium.webdriver.support.ui import WebDriverWait
	from selenium.webdriver.support import expected_conditions as EC
	from selenium.webdriver.common.action_chains import ActionChains
	from selenium.webdriver.common.keys import Keys
	from selenium.webdriver.chrome.service import Service
	from subprocess import CREATE_NO_WINDOW
	import chromedriver_autoinstaller
	chromedriver_autoinstaller.install()

try:
	import plotly.express as px
except:
	install_package("pandas")
	install_package("plotly")
	import plotly.express as px
	#webdriver-manager
########################## Motifs definition based on text processing ############################
 
def sheets(line):
	sheets_Info={}
	stripped =  line.split("::")
	stripped = [elem.strip() for elem in stripped]
	sheets_Info["Chain"] = stripped[0]
	sheets_Info["Sheet ID"] = stripped[1]
	sheets_Info["Number of strands"] = stripped[2]
	sheets_Info["Strands"]=[]
	return sheets_Info
def representSheets(result,pdb_id):
	global motifDicti
	global visitedChains
	for sheets in result["sheets"]:
		chain =  sheets["Chain"]
		for SimilarChains in similarityChain[chain]:
			visitedChains.append(SimilarChains)
			sheetId= sheets["Sheet ID"]
			strands = sheets["Strands"]
			selector=pdb_id+"."+SimilarChains + ".sheet."+sheetId
			selectStrandsMotifs(selector,strands,pdb_id,SimilarChains)
			group(pdb_id+"."+SimilarChains+".sheets",selector)
			motifDicti["sheets"][SimilarChains].append(selector)
			group(pdb_id + "."+SimilarChains+".structural_motifs",pdb_id+"."+SimilarChains+".sheets")
def bab_motifs(line):
	bab_motifs_Info={}
	stripped =  line.split("::")
	stripped = [elem.strip() for elem in stripped]
	bab_motifs_Info["Chain"] = stripped[0]
	bab_motifs_Info["Length of the loop"] = stripped[11]
	bab_motifs_Info["Residues"] = {"0":[stripped[1],stripped[3]],"1":[stripped[5],stripped[7]]}
	bab_motifs_Info["Number of helices"] = stripped[12]
	bab_motifs_Info["Helix length"] = stripped[13]
	bab_motifs_Info["Helices"] = {"0":[stripped[14],stripped[15]],"1":[stripped[16],stripped[17]],"2":[stripped[18],stripped[19]],"3":[stripped[20],stripped[21]],"4":[stripped[22],stripped[23]]}
	return bab_motifs_Info
def representBab_Motifs(result,pdb_id):
	global visitedChains
	global motifDicti
	for bab_motif in result["bab_motifs"]:
		chain =  bab_motif["Chain"]
		for SimilarChains in similarityChain[chain]:
			visitedChains.append(SimilarChains)
			selector  = pdb_id+"."+SimilarChains + "." + "bab_motif." 
			selector+=bab_motif["Residues"]["0"][0] +"_" + bab_motif["Residues"]["0"][1]  + "."  +bab_motif["Residues"]["1"][0] +"_" + bab_motif["Residues"]["1"][1] 
			arrayToUse=[]
			for residue in bab_motif["Residues"]:
				bab_motif["Residues"][residue][0]=str(bab_motif["Residues"][residue][0])
				bab_motif["Residues"][residue][1]=str(bab_motif["Residues"][residue][1])
				arrayToUse.append(bab_motif["Residues"][residue])
			for helix in bab_motif["Helices"]:
				if bab_motif["Helices"][helix][0]!='00000':
					bab_motif["Helices"][helix][0]=str(bab_motif["Helices"][helix][0])
					bab_motif["Helices"][helix][1]=str(bab_motif["Helices"][helix][1])
					arrayToUse.append(bab_motif["Helices"][helix])
			selectMultiResiduesMotifs(selector,arrayToUse,SimilarChains,pdb_id)
			group(pdb_id+"."+SimilarChains+".bab_motifs",selector)
			motifDicti["bab_motifs"][SimilarChains].append(selector)
			group(pdb_id + "."+SimilarChains+".structural_motifs",pdb_id+"."+SimilarChains+".bab_motifs")
def beta_hairpins(line):
	beta_hairpins_Info={}
	stripped =  line.split("::")
	stripped = [elem.strip() for elem in stripped]
	beta_hairpins_Info["Chain"] = stripped[0]
	beta_hairpins_Info["Type"] = stripped[11]
	beta_hairpins_Info["Residues"] = {"0":[stripped[1],stripped[3]],"1":[stripped[5],stripped[7]]}
	beta_hairpins_Info["Number of residues"] = stripped[9]
	beta_hairpins_Info["Number of bridging residues"] = stripped[10]
	beta_hairpins_Info["Sequence"] = stripped[10]
	beta_hairpins_Info["Sequence"] = stripped[15]
	return beta_hairpins_Info
def representBeta_Hairpins(result,pdb_id):
	global visitedChains
	global motifDicti
	for hairpin in result["beta_hairpins"]:
		chain =  hairpin["Chain"]
		for SimilarChains in similarityChain[chain]:
			visitedChains.append(SimilarChains)
			selector  = pdb_id+"."+SimilarChains + "." + "beta_hairpin." 
			selector+=hairpin["Residues"]["0"][0] +"_" + hairpin["Residues"]["0"][1]  + "."  +hairpin["Residues"]["1"][0] +"_" + hairpin["Residues"]["1"][1] 
			selectTwoResiduesMotifs(selector,str(int(hairpin["Residues"]["0"][0])),str(int(hairpin["Residues"]["1"][1])),SimilarChains,pdb_id)
			group(pdb_id+"."+SimilarChains+".beta_hairpins",selector)
			motifDicti["beta_hairpins"][SimilarChains].append(selector)
			group(pdb_id + "."+SimilarChains+".structural_motifs",pdb_id+"."+SimilarChains+".beta_hairpins")
def psi_loops(line):
	psi_loops_Info={}
	stripped =  line.split("::")
	stripped = [elem.strip() for elem in stripped]
	psi_loops_Info["Chain"] = stripped[0]
	psi_loops_Info["Length of the loop"] = stripped[11]
	psi_loops_Info["Residues"] = {"0":[stripped[1],stripped[3]],"1":[stripped[5],stripped[7]]}
	return psi_loops_Info
def representPsi_Loops(result,pdb_id):
	global visitedChains
	global motifDicti
	for psi_loop in result["psi_loops"]:
		chain =  psi_loop["Chain"]
		for SimilarChain in similarityChain[chain]:
			grouppingA=pdb_id+"."+SimilarChain+".psi_loops"
			visitedChains.append(SimilarChain)
			selector = pdb_id + "." + SimilarChain+ ".psi_loop" 
			line = ""
			for ind in psi_loop["Residues"].keys():
				elem = psi_loop["Residues"][ind]
				selector +=   "." + str(elem[-1])
				if line=="":
					line +=  " ( resi " + str(elem[-1])+ " and chain " + str(SimilarChain) + " and "  + pdb_id + " ) "
				else:
					line +=  " or ( resi " + str(elem[-1])+ " and chain " + str(SimilarChain) + " and "  + pdb_id + " ) "
			line = "select "+selector+ " , " + line
			cmd.do(line)
			group(grouppingA,selector)
			motifDicti["psi_loops"][SimilarChain].append(selector)
			group(pdb_id + "."+SimilarChain+".structural_motifs",grouppingA)
def beta_bulges(line):
	beta_bulge_Info={}
	hasResidues=[False,False,False,False,False]
	residues=[[],[]]
	strands=[[],[],[]]
	stripped =  line.split("::")
	stripped = [elem.strip() for elem in stripped]
	for elem in stripped:
		if not elem.isalnum():
			stripped.remove(elem)
	indexOfFirstTrue = stripped.index("TRUE")
	indexOfFirstFalse = stripped.index("FALSE")
	startOfExistenceIndicator = min(indexOfFirstFalse,indexOfFirstTrue)
	psiPhiToSkip=0

	for ind in range(startOfExistenceIndicator,startOfExistenceIndicator+5):
		if "T" in stripped[ind]:
			hasResidues[ind-startOfExistenceIndicator]=True
		psiPhiToSkip+=2
	for ind in range(1,startOfExistenceIndicator-2):
		if stripped[ind].isalpha() and len(stripped[ind])==1:
			residues[0].append(stripped[ind])
		else:
			if stripped[ind].isnumeric():
				residues[1].append(stripped[ind])
	strands[0].append(stripped[startOfExistenceIndicator+psiPhiToSkip])
	strands[1].append(stripped[startOfExistenceIndicator+psiPhiToSkip+1])
	strands[2].append(stripped[startOfExistenceIndicator+psiPhiToSkip+2])
	strands[0].append(stripped[startOfExistenceIndicator+psiPhiToSkip+3])
	strands[1].append(stripped[startOfExistenceIndicator+psiPhiToSkip+4])
	strands[2].append(stripped[startOfExistenceIndicator+psiPhiToSkip+5])
	beta_bulge_Info["Chain"] = stripped[0]
	beta_bulge_Info["Type"] = stripped[startOfExistenceIndicator-2]
	beta_bulge_Info["Class"] = stripped[startOfExistenceIndicator-1]
	beta_bulge_Info["Residues"] = {}
	count=0
	for ind in range(0,len(hasResidues)):
		if hasResidues[ind]:
			beta_bulge_Info["Residues"][str(ind)]=[residues[0][count],residues[1][count]]
			count+=1
	beta_bulge_Info["Strands"] = {"0":[strands[0][0],strands[1][0],strands[2][0]],"1":[strands[0][1],strands[1][1],strands[2][1]]}
	return beta_bulge_Info
def representBeta_Bulges(result,pdb_id):
	global visitedChains
	global motifDicti
	for beta_bulge in result["beta_bulges"]:
		chain =  beta_bulge["Chain"]
		for SimilarChain in similarityChain[chain]:
			grouppingA=pdb_id+"."+SimilarChain+".beta_bulges"
			visitedChains.append(SimilarChain)
			selector = pdb_id + "." + SimilarChain+ ".beta_bulge" 
			line = ""
			for ind in beta_bulge["Residues"].keys():
				elem = beta_bulge["Residues"][ind]
				selector +=   "." + str(elem[-1])
				if line=="":
					line +=  " ( resi " + str(elem[-1])+ " and chain " + str(SimilarChain) + " and "  + pdb_id + " ) "
				else:
					line +=  " or ( resi " + str(elem[-1])+ " and chain " + str(SimilarChain) + " and "  + pdb_id + " ) "
			line = "select "+selector+ " , " + line
			cmd.do(line)
			group(grouppingA,selector)
			motifDicti["beta_bulges"][SimilarChain].append(selector)
			group(pdb_id + "."+SimilarChain+".structural_motifs",grouppingA)
def strands(line):
	strands_Info={}
	stripped =  line.split("::")
	stripped = [elem.strip() for elem in stripped]
	strands_Info["Chain"] = stripped[0]
	strands_Info["Strand ID"] = stripped[1]
	strands_Info["Residues"] = {"0":[stripped[2],stripped[4]]}
	strands_Info["Sheets ID"]=[stripped[6],stripped[7]]
	strands_Info["Sequences"]=stripped[9]
	return strands_Info                           
def representStrands(result,pdb_id):
	global visitedChains
	global motifDicti
	for strand in result["strands"]:
		chain =  strand["Chain"]
		for SimilarChains in similarityChain[chain]:
			visitedChains.append(SimilarChains)
			selector  = pdb_id+"."+SimilarChains + "." + "strand." 
			selector+=strand["Residues"]["0"][0] +"_" + strand["Residues"]["0"][1]  
			selectTwoResiduesMotifs(selector,str(int(strand["Residues"]["0"][0])),str(int(strand["Residues"]["0"][1])),SimilarChains,pdb_id)
			group(pdb_id+"."+SimilarChains+".strands",selector)
			motifDicti["strands"][SimilarChains].append(selector)
			group(pdb_id + "."+SimilarChains+".structural_motifs",pdb_id+"."+SimilarChains+".strands")
def helices(line):
	helices_Info={}
	stripped =  line.split("::")
	stripped = [elem.strip() for elem in stripped]
	helices_Info["Chain"] = stripped[0]
	helices_Info["Helix ID"] = stripped[1]
	helices_Info["Residues"] = {"0":[stripped[2],stripped[4]]}
	helices_Info["Type"]=stripped[6]
	return helices_Info                         
def representHelices(result,pdb_id):
	global visitedChains
	global motifDicti
	for helice in result["helices"]:
		chain =  helice["Chain"]
		for SimilarChains in similarityChain[chain]:
			visitedChains.append(SimilarChains)
			selector  = pdb_id+"."+SimilarChains + "." + "helice." 
			selector+=helice["Residues"]["0"][0] +"_" + helice["Residues"]["0"][1]  
			selectTwoResiduesMotifs(selector,str(int(helice["Residues"]["0"][0])),str(int(helice["Residues"]["0"][1])),SimilarChains,pdb_id)
			group(pdb_id+"."+SimilarChains+".helices",selector)
			motifDicti["helices"][SimilarChains].append(selector)
			group(pdb_id + "."+SimilarChains+".structural_motifs",pdb_id+"."+SimilarChains+".helices")
def helix_interactions(line):
	helix_interactions_Info={}
	stripped =  line.split("::")
	stripped = [elem.strip() for elem in stripped]
	helix_interactions_Info["Chain"] = stripped[0]
	helix_interactions_Info["Interacting Helices"] = [[stripped[1],stripped[2]],[stripped[4],stripped[5]]]
	helix_interactions_Info["Distance"] = stripped[7]
	helix_interactions_Info["Helices"] = []
	return helix_interactions_Info
def representHelix_Interactions(result,pdb_id):
	global motifDicti
	global visitedChains
	for helix_interaction in result["helix_interactions"]:
		chain =  helix_interaction["Chain"]
		helix_interactionID= helix_interaction["Interacting Helices"][0][0]+"."+helix_interaction["Interacting Helices"][0][1]+"."+helix_interaction["Interacting Helices"][1][0]+"."+helix_interaction["Interacting Helices"][1][1]
		helices = helix_interaction["Helices"]
		selector=pdb_id+"."+chain + ".helix_interaction."+helix_interactionID
		selectHelixInteractionMotifs(selector,helices,pdb_id,chain)
		group(pdb_id+"."+chain+".helix_interactions",selector)
		motifDicti["helix_interactions"][chain].append(selector)
		group(pdb_id + "."+chain+".structural_motifs",pdb_id+"."+chain+".helix_interactions")
def beta_turns(line):
	beta_turns_Info={}
	stripped =  line.split("::")
	stripped = [elem.strip() for elem in stripped]
	beta_turns_Info["Chain"] = stripped[0]
	beta_turns_Info["Residues"] = [stripped[1],stripped[3],stripped[5],stripped[7]]
	beta_turns_Info["Sequence"]=stripped[9]
	beta_turns_Info["Type"]=stripped[10]
	return beta_turns_Info
def representBeta_Turns(result,pdb_id):
	global motifDicti
	global visitedChains
	for beta_turn in result["beta_turns"]:
		chain =  beta_turn["Chain"]
		for SimilarChains in similarityChain[chain]:
			visitedChains.append(SimilarChains)
			beta_turnsId= str(beta_turn["Residues"][0])+"."+str(beta_turn["Residues"][1])+"."+str(beta_turn["Residues"][2])+"."+str(beta_turn["Residues"][3])
			residues = beta_turn["Residues"]
			selector=pdb_id+"."+SimilarChains + ".beta_turn."+beta_turnsId
			selectAcrossMultiResiduesMotifs(selector,residues,SimilarChains,pdb_id)
			group(pdb_id+"."+SimilarChains+".beta_turns",selector)
			motifDicti["beta_turns"][SimilarChains].append(selector)
			group(pdb_id + "."+SimilarChains+".structural_motifs",pdb_id+"."+SimilarChains+".beta_turns")
def gamma_turns(line):
	gamma_turns_Info={}
	stripped =  line.split("::")
	stripped = [elem.strip() for elem in stripped]
	gamma_turns_Info["Chain"] = stripped[0]
	gamma_turns_Info["Residues"] = [stripped[1],stripped[3],stripped[5]]
	gamma_turns_Info["Sequence"]=stripped[7]
	gamma_turns_Info["Type"]=stripped[8]
	return gamma_turns_Info
def representGamma_Turns(result,pdb_id):
	global motifDicti
	global visitedChains
	for gamma_turn in result["gamma_turns"]:
		chain =  gamma_turn["Chain"]
		for SimilarChains in similarityChain[chain]:
			visitedChains.append(SimilarChains)
			gamma_turnID= str(gamma_turn["Residues"][0])+"."+str(gamma_turn["Residues"][1])+"."+str(gamma_turn["Residues"][2])
			residues = gamma_turn["Residues"]
			selector=pdb_id+"."+SimilarChains + ".gamma_turn."+gamma_turnID
			selectAcrossMultiResiduesMotifs(selector,residues,SimilarChains,pdb_id)
			group(pdb_id+"."+SimilarChains+".gamma_turns",selector)
			motifDicti["gamma_turns"][SimilarChains].append(selector)
			group(pdb_id + "."+SimilarChains+".structural_motifs",pdb_id+"."+SimilarChains+".gamma_turns")
def disulphides(line):
	disulphides_Info={}
	stripped =  line.split("::")
	stripped = [elem.strip() for elem in stripped]
	disulphides_Info["Chain"] = stripped[0]
	disulphides_Info["Residues"] = [[stripped[1],stripped[2]],[stripped[3],stripped[4]]]
	return disulphides_Info
def representDisulphides(result,pdb_id):
	global motifDicti
	global visitedChains
	for disulphide in result["disulphides"]:
		chain =  disulphide["Chain"]
		disulphideID= disulphide["Residues"][0][0]+"."+disulphide["Residues"][0][1]+"."+disulphide["Residues"][1][0]+"."+disulphide["Residues"][1][1]
		residues = disulphide["Residues"]
		selector=pdb_id+"."+chain + ".disulphide."+disulphideID
		selectDisulphidesMotifs(selector,residues[0][1],residues[0][0],residues[1][1],residues[1][0],pdb_id)
		group(pdb_id+"."+chain+".disulphides",selector)
		motifDicti["disulphides"][chain].append(selector)
		group(pdb_id + "."+chain+".structural_motifs",pdb_id+"."+chain+".disulphides")

########################################################

#####################Analyzing the Promotif file #####################

def make_sheets(dictionary):
	for sheet in dictionary["sheets"]:
		sheetID=sheet["Sheet ID"]
		for strand in dictionary["strands"]:
			if str(sheetID) in strand["Sheets ID"]:
				sheet["Strands"].append(strand)
def make_helix_interactions(dictionary):
	for helices_interaction in dictionary["helix_interactions"]:
		helices_interactionID=helices_interaction["Interacting Helices"]
		for helix in dictionary["helices"]:
			if str(helix["Helix ID"]) in helices_interactionID[0] or str(helix["Helix ID"]) in helices_interactionID[1] :
				helices_interaction["Helices"].append(helix)
def make_dicto():
	dico = { 
			'sheets': [],
			'bab_motifs' : [],
			'beta_hairpins' : [],
			'psi_loops' : [],
			'beta_bulges' : [],
			'strands'  :  [],
			'helices' : [],
			'helix_interactions' : [],
			'beta_turns' : [],
			'gamma_turns' : [],
			'disulphides' : []
			}
	return dico
 
def make_dictio():
	dico = { 
			'sheets': {},
			'bab_motifs' : {},
			'beta_hairpins' : {},
			'psi_loops' : {},
			'beta_bulges' : {},
			'strands'  :  {},
			'helices' : {},
			'helix_interactions' : {},
			'beta_turns' : {},
			'gamma_turns' : {},
			'disulphides' : {}
			}
	return dico

def AnalyzePromotifData(code,pdbId=""):
	global similarityChain
	#function to initialize a dictionary that has the information related to the motifs data from the promotif files.
	dico = make_dicto()
	raw_data  = {'sheets': ""    ,
	'bab_motifs' : "",
	'beta_hairpins' : "",
	'psi_loops' : "",
	'beta_bulges' : "",
	'strands'  :  "",
	'helices' : "",
	'helix_interactions' : "",
	'beta_turns' : "",
	'gamma_turns' : "",
	'disulphides' : ""
	}
	array =  code.split("START:")
	array = array[1: (len(array)-1)]
	for  i in array:
		identifier  = i.split("\n")[0].lower()
		raw_data[identifier]  = i
	for i in raw_data:
		lines =  list(filter ( lambda x: x!="" , raw_data[i].split("\n") ) )
		lines = lines[2:]
		for j in lines:
			if i== "sheets": 
				dico[i].append(sheets(j))
			if i== "bab_motifs": 
				dico[i].append(bab_motifs(j))
			if i== "beta_hairpins":
				dico[i].append(beta_hairpins(j))
			if i== "psi_loops":
				dico[i].append(psi_loops(j))
			if i== "strands":
				dico[i].append(strands(j))
			if i== "helices":
				dico[i].append(helices(j))
			if i== "helix_interactions":
				dico[i].append(helix_interactions(j))
			if i== "disulphides":
				dico[i].append(disulphides(j))
			if i == "gamma_turns":
				dico[i].append(gamma_turns(j))
			if i == "beta_bulges":
				dico[i].append(beta_bulges(j))
			if i == "beta_turns":
				dico[i].append(beta_turns(j))
	make_sheets(dico)
	make_helix_interactions(dico)
	json_string = json.dumps(dico)
	with open("PyProtif/"+pdbId.lower()+"/Structural_motifs/"+pdbId.lower()+"_promotif_Analyzed.json", 'w') as outfile:
		outfile.write(json_string)
	return dico
def getFasta(pdb_id):
	Adjustment=adjustForResidueDifference(pdb_id)
	dictionary=findResiduesToRemoveFromFasta(pdb_id,Adjustment)
	chains= findChains(pdb_id)[1]
	fastaDict=OrderedDict()
	for chain in chains:
		fastaDict[str(chain)]=""
		if chain in dictionary.keys():
			fastaDict[str(chain)]=dictionary[str(chain)][-1]["AdjustedChain"]
		else:
			fastaDict[str(chain)]=getFastaChain(pdb_id,str(chain))
	fastaFileAdjusted=""
	f = open("PyProtif/"+pdb_id.lower()+"/"+pdb_id.lower()+"_adjusted.fasta",'w')
	for key in fastaDict.keys():
		if fastaDict[key] and fastaDict[key]!= "not found" and len(fastaDict[key])>0:
			f.write(">"+str(pdb_id)+"_"+str(key)+"\n"+fastaDict[key].strip()+"\n")
			fastaFileAdjusted+=">"+str(pdb_id)+"_"+str(key)+"\n"+fastaDict[key].strip()+"\n"
	f.close()
	return fastaFileAdjusted
def getFastaTemp(pdb_id):
	url = "https://www.rcsb.org/fasta/entry/{:s}/display".format(pdb_id)
	response = urllib.request.urlopen(url)
	fastaFile = response.read().decode('utf-8')
	with open("PyProtif/"+pdb_id.lower()+"/"+pdb_id.lower()+".fasta",'w') as f:
		f.write(fastaFile)
		f.close()
	fastaFileAdjusted=""
	with open("PyProtif/"+pdb_id.lower()+"/"+pdb_id.lower()+".fasta") as f:
		f2 = open("PyProtif/"+pdb_id.lower()+"/"+pdb_id.lower()+"_adjusted_temp.fasta",'w')
		for line in f:
			if ">" not in line:
				f2.write(line)
				fastaFileAdjusted+=line
			else:
				fastaHead = line.split("|")[0].split("_")[0]+"_"
				chains=line.split("|")[1].split(" ")[1:]
				for elem in range(0,len(chains)):
					if "auth" in chains[elem]:
						chains[elem]=chains[elem].replace("[auth","")
						chains[elem+1]=chains[elem+1].replace("]","")
						chains[elem]=chains[elem].replace(chains[elem],chains[elem+1])
				chains = list(set(chains))
				for elem in chains:
					if elem==chains[-1]:
						fastaHead+=elem.replace(",","")
					else:
						fastaHead+=elem.replace(",","")+"_"
				fastaHead+="\n"
				f2.write(fastaHead)
				fastaFileAdjusted+=fastaHead
		f.close()
		f2.close()
	return fastaFileAdjusted

def findChains(pdb_id):
	#function to analyze fasta file and get the chain numbers as well as chain labels
    dicti=adjustForResidueDifference(pdb_id)
    with open("PyProtif/"+pdb_id.lower()+"/"+pdb_id.lower()+"_chains.txt",'w') as f:
    	strin=str(len(dicti.keys()))+" "
    	for chain in dicti.keys():
    		if chain == list(dicti.keys())[-1]:
    			strin+=chain
    		else:
    			strin+=chain+" "
    	f.write(strin)
    	f.close()
    return [len(dicti.keys()), list(dicti.keys())]

def readChains(pdb_id):
	with open("PyProtif/"+pdb_id.lower()+"/"+pdb_id.lower()+"_chains.txt") as f:
		for line in f:
			return [line.split(" ")[0],line.split(" ")[1:]]
def getPromotifData(pdb_id):
	try:
		pdb =  pdb_id.lower()
		url  = "http://www.ebi.ac.uk/thornton-srv/databases/PDBsum/{:s}{:s}/{:s}/promotif.dat".format(pdb[1],pdb[2],pdb)
		result_page =  urllib.request.urlopen(url)
		text  = result_page.read().decode('utf-8')
		with open("PyProtif/"+pdb_id.lower()+"/Structural_motifs/"+pdb_id.lower()+"_promotif_data.txt",'w') as f:
			f.write(text)
			f.close()
		return text
	except:
		messagebox.showwarning("Warning","Couldn't get structural motif data for the pdb: "+pdb)
		return ""

def motifCounter( chains , promotif):
    #function that counts the existence of each motif given a motif dictionary for each chain
    motifCounter  = {}
    for i  in range(len(chains)):
        #initialize the motif count to 0 for each chain
        motifCounter[chains[i]] = { 'sheets': 0,
                     'bab_motifs' : 0,
                     'beta_hairpins' : 0,
                     'psi_loops' : 0,
                     'beta_bulges' : 0,
                     'strands'  :  0,
                     'helices' : 0,
                     'helix_interactions' : 0,
                    'beta_turns' : 0,
                    'gamma_turns' : 0, 
                    'disulphides' : 0
                    }
    for identifier in promotif:
        listr = promotif[identifier]
        for j in listr:
            chain =  j[0]
            motifCounter[chain][identifier] +=1    
    return motifCounter

######################################################################################

#Structure CLuster finds other proteins with similar 3D-structure via PDB ( this has been automatically adjusted by PDB on december 9 2020 check this link 
#https://www.rcsb.org/news/5fc9176809ae2a096d081e28
    
#doubles adjustment if needed

##############################Initialization methods #####################################

def fetchPDB(pdb_id):
	global loadedPDBs
	global validPDB
	validPDB=False
	filename=pdb_id+".pdb"
	filename2=pdb_id+".cif"
	loadedPDBs=cmd.get_object_list('solvent')
	if pdb_id.lower() in loadedPDBs:
		if pdb_id.casefold() in os.listdir("PyProtif"):
			if filename.casefold() in os.listdir("PyProtif"+fileSeparator+pdb_id) or filename2.casefold() in os.listdir("PyProtif"+fileSeparator+pdb_id):
				validPDB=True
				return 1
		
	#check if the pdb exist
	returnCode = testExistence(pdb_id)
	pdb_id = pdb_id.lower()
	#check the existence of the local files (.pdb) first:
	if pdb_id.casefold() in os.listdir("PyProtif"):
		if filename.casefold() in os.listdir("PyProtif"+fileSeparator+pdb_id):
			loadedPDBs.append(pdb_id)
			return fetchPDBOffline(pdb_id)
	#check the existence of the local files (.cif) first:
	if pdb_id.casefold() in os.listdir("PyProtif"):
		if filename2.casefold() in os.listdir("PyProtif"+fileSeparator+pdb_id):
			pdbformat  = "load PyProtif"+fileSeparator+pdb_id+fileSeparator+pdb_id+".cif"
			cmd.do(pdbformat)
			saveFormat = "save "+pdb_id+".pdb"
			cmd.do(saveFormat)
			os.rename(pdb_id+".pdb","PyProtif/"+pdb_id+"/"+pdb_id+".pdb")
			loadedPDBs.append(pdb_id)
			return 1
	#the files do not exist localy and we have to fetch them online, 
	# if we have the pdb format of the file online
	if returnCode==1:
		pdbformat  = "fetch " + pdb_id + " , type= pdb"
		cmd.do(pdbformat)
		value=0
		while pdb_id+".pdb" not in os.listdir():
			value+=1
		os.rename(pdb_id+".pdb","PyProtif/"+pdb_id+"/"+pdb_id+".pdb")

	# if we have only cif format of the file online
	if returnCode==2:
		ciformat  = "fetch " + pdb_id + " , type= cif"
		cmd.do(ciformat)
		saveFormat = "save "+pdb_id+".pdb"
		cmd.do(saveFormat)
		value=0
		while pdb_id+".pdb" not in os.listdir():
			value+=1
		os.rename(pdb_id+".cif","PyProtif/"+pdb_id+"/"+pdb_id+".cif")
		os.rename(pdb_id+".pdb","PyProtif/"+pdb_id+"/"+pdb_id+".pdb")
		returnCode=fetchPDBOffline(pdb_id)
	#if at any stage the return code is different then 1 means that the pdb file was not fetched and the
	#pdbID given is invalid
	if returnCode==1:
		validPDB = True
		loadedPDBs.append(pdb_id)
	else:
		validPDB = False
	return returnCode

def testExistence(pdb_id):
	#method to test if we are able to find the pdb id online or not :
	# return code 1 if it is found as pdb, return code 2 if found as .cif, return code 3 if we were not able to fetch
	pdb_id=pdb_id.lower()
	url = "https://files.rcsb.org/view/{:s}.pdb".format(pdb_id)
	try:
		response = urllib.request.urlopen(url)
		pdb_cont = response.read().decode('utf-8')
		return 1
	except:
		url = "https://files.rcsb.org/view/{:s}.cif".format(pdb_id)
		try:
			response = urllib.request.urlopen(url)
			cif_cont = response.read().decode('utf-8')
			with open(pdb_id+".cif",'w') as f:
				f.write(cif_cont)
				f.close()
			return 2
		except:
			return 3

def fetchPDBOffline(pdb_id):
	#check localy for the existence of .pdb files
	#return 1 if a .pdb file was found 3 in case not
	global validPDB
	pdb_id=pdb_id.lower()
	filename  = pdb_id + ".pdb"
	if pdb_id.casefold() in os.listdir("PyProtif"):
		if filename.casefold() in os.listdir("PyProtif"+fileSeparator+pdb_id):
			pdbformat  = "load PyProtif"+fileSeparator+pdb_id+fileSeparator+pdb_id + ".pdb"
			cmd.do(pdbformat)
			validPDB = True
			return 1
	else:
		validPDB = False
		return 3
def getFastaChain(pdb_id,chain,initial = 0):
	if initial == 0:
		filename  = pdb_id + "_adjusted_temp.fasta"
	else:
		filename  = pdb_id + "_adjusted.fasta"
	#pdb found localy
	if pdb_id.casefold() in os.listdir("PyProtif"):
		if filename.casefold() in os.listdir("PyProtif"+fileSeparator+pdb_id):
			openfile  =  open("PyProtif"+fileSeparator+pdb_id+fileSeparator+filename, "r")
			pdb_fasta = openfile.read()
			openfile.close()
	else:
		pdb_fasta=getFastaTemp(pdb_id)
	headsAndChains = pdb_fasta.split("\n>")
	for headAndChain in headsAndChains:
		if "_"+str(chain) in headAndChain:
			headNchain=headAndChain.split("\n")
			return headNchain[1].strip()
	return "not found"
def findResiduesToRemoveFromFasta(pdb_id,chainAdjustDict):
	oneLetterCodeDictionary={"ala":"A","arg":"R","asn":"N","asp":"D","asx":"B","cys":"C","glu":"E","gln":"Q","glx":"Z","gly":"G","his":"H","ile":"I","leu":"L","lys":"K","met":"M","phe":"F","pro":"P","ser":"S","thr":"T","trp":"W","tyr":"Y","val":"V"}
	filename  = pdb_id + ".pdb"
	#pdb found localy
	if pdb_id.casefold() in os.listdir("PyProtif"):
		if filename.casefold() in os.listdir("PyProtif"+fileSeparator+pdb_id):
			openfile  =  open("PyProtif"+fileSeparator+pdb_id+fileSeparator+filename, "r")
			pdb_cont = openfile.read()
			openfile.close()
	else:
		url = "https://files.rcsb.org/view/{:s}.pdb".format(pdb_id)
		response = urllib.request.urlopen(url)
		pdb_cont = response.read().decode('utf-8')
	#finding the DBREF / DBREF1 lines
	pattern = re.compile(r"\nREMARK 465.*")
	text = pattern.findall(pdb_cont)
	seq=""
	dictio = {}
	for i in text:
		split = i.split()
		if split[0] == "REMARK" and split[1] == "465" and len(split) ==5:
			split=split[2:]
			chainID = str(split[1].strip())
			residueName = str(split[0].lower().strip())
			if residueName not in oneLetterCodeDictionary.keys():
				continue
			resiudeSingleLetter = str(oneLetterCodeDictionary[residueName]).strip()
			residueIndexToRemove = int(split[2].strip())
			chainStartIndex = int(chainAdjustDict[chainID])
			originalChain = str(getFastaChain(pdb_id,chainID))
			chainSeqOrigin = originalChain
			residueFirstOccurence = chainSeqOrigin.index(resiudeSingleLetter)
			newChainSeq = chainSeqOrigin
			if chainID not in dictio.keys():
				dictio[chainID]=[{},{}]
				dictio[chainID][0][chainID]=0 # residues before the start of the chain
				dictio[chainID][1][chainID]=0 #residues removed after the start of the chain
			if residueIndexToRemove < chainStartIndex:
				newChainSeq = newChainSeq[residueFirstOccurence+1:]
				dictio[chainID][0][chainID]+=1
			else:
				#to split into two substrings and remove the element in between
				if dictio[chainID][1][chainID] ==0 :
					newChainSeq = chainSeqOrigin[dictio[chainID][0][chainID]:]  # remove before start residues
					newChainSeq2 = newChainSeq # create a copy of that string to split it  				
					residuesRemovedBefore = dictio[chainID][0][chainID]
					residuesRemovedAfter = dictio[chainID][1][chainID]
					indexToRemove= (residueIndexToRemove - chainStartIndex) + (residuesRemovedBefore - residuesRemovedAfter)
					#residue index - start + residues before the start - residues after the chain
					newChainSeq = newChainSeq2[:indexToRemove] + newChainSeq[indexToRemove+1: ]
					dictio[chainID][1][chainID]+=1
					#strObj = strObj[0 : index : ] + strObj[index + 1 : :]
				else:
					newChainSeq = dictio[chainID][-1]["AdjustedChain"]
					newChainSeq2 = newChainSeq # create a copy of that string to split it  				
					residuesRemovedBefore = dictio[chainID][0][chainID]
					residuesRemovedAfter = dictio[chainID][1][chainID]
					indexToRemove= (residueIndexToRemove - chainStartIndex) + (residuesRemovedBefore - residuesRemovedAfter)
					#residue index - start + residues before the start - residues after the chain
					newChainSeq = newChainSeq2[:indexToRemove] + newChainSeq[indexToRemove+1: ]
					dictio[chainID][1][chainID]+=1
			dictio[chainID].append({"residueToRemove":resiudeSingleLetter,"residueIndexToRemove":residueIndexToRemove,"chainStart":chainStartIndex,"fastaChain":originalChain,"AdjustedChain":str(newChainSeq.strip())})
	return dictio
def adjustForResidueDifference(pdb_id):
	#check adjustment to pdb files for different chain/end start position
	#returns a dictionary of chainsID as keys and start position as values
	filename  = pdb_id + ".pdb"
	#pdb found localy
	if pdb_id.casefold() in os.listdir("PyProtif"):
		if filename.casefold() in os.listdir("PyProtif"+fileSeparator+pdb_id):
			openfile  =  open("PyProtif"+fileSeparator+pdb_id+fileSeparator+filename, "r")
			pdb_cont = openfile.read()
			openfile.close()
	else:
		url = "https://files.rcsb.org/view/{:s}.pdb".format(pdb_id)
		response = urllib.request.urlopen(url)
		pdb_cont = response.read().decode('utf-8')
	#finding the DBREF / DBREF1 lines
	pattern = re.compile(r"\nDBREF.*")
	text = pattern.findall(pdb_cont)
	dictio = {}
	for i in text:
		split = i.split()
		if split[0] == "DBREF1" or split[0] == "DBREF":
			dictio[split[2]] = split[3]
	if dictio == {}:
		driver.close()
		cmd.delete("unselect")
		cmd.delete(pdb_id)
		unselect = "select unselect, none"
		cmd.do(unselect)
		raise Exception
	return dictio

def createDirectory(pathi):
	#create a directory
	if sys.version_info[0]<3:
		try: 
			os.makedirs(pathi)
		except OSError:
			if not os.path.isdir(pathi):
				raise
	else:
		os.makedirs(pathi, exist_ok=True)

def isConnected():
	#check if client is connected or not to internet
	#return a True if yes False if not
	try:
		url  = "https://www.rcsb.org/"
		urllib.request.urlopen(url)
		return True
	except:
		return False

def getSimilarChains(pdb_id):
	fastaFile  = getFastaTemp(pdb_id)
	lines = fastaFile.split("\n")
	SimilarChains={}
	for line in lines:
		if ">" in line:
			chains=line.split("_")[1:]
			for chain in chains:
				SimilarChains[chain]=chains
	return SimilarChains
########################################################################

################## Pymol managment methods  ############################

#Chain selection and grouping
def selectChain(name, elemA, elemB):
	# 6vww.A,6vww,A
	selection= "select "+str(name)+ " , "+str(elemA)+ " and chain "+str(elemB)
	cmd.do(selection)
def groupChainSelection(elemA, name,elemB,elemC,count=1):
	#6vww.chains, 6vww.A,6vww,A
	if count>0:
		selectChain(name,elemB,elemC)
	group(elemA,name)
# beta_hairpins, bab_motifs ,psi_loops, disulphides  selection
def selectAcrossMultiResiduesMotifs(name,residues,chain,pdb):
	residues=[int(elem) for elem in residues]
	miniRes=min(residues)
	maxiRes=max(residues)
	selectTwoResiduesMotifs(name,miniRes,maxiRes,chain,pdb)
def selectMultiResiduesMotifs(name,residues,chain,pdb):
	selection="select "+ name   +" , " 
	for residue in residues:
		if residue != residues[-1]:
			selection+= "(resi " +  str(residue[0]) + "-" +str(residue[1]) + " and chain " + str(chain) + " and "  + pdb+") or "
		else:
			selection+="(  resi " +str(residue[0]) +"-" + str(residue[1]) + " and chain " + str(chain) + " and "  + pdb+")"
	cmd.do(selection)
def selectFourResiduesMotifs(name,residue1,residue2,residue3,residue4,chain,pdb):
	selection="select "+ name   +" , " + "(resi " +  str(residue1) + "-" +str(residue2) + " and chain " + chain + " and "+pdb+" ) or ( resi " +str(residue3) +"-" + str(residue4) + " and chain " + chain + " and " +pdb+" )"
	cmd.do(selection)
# disulphideS:
def selectDisulphidesMotifs(name,residue1,chain1,residue2,chain2,pdb):
	selection="select "+ name   +" , " + "(resi " +  str(residue1) + " and chain " +str(chain1) + " and " +str(pdb) + "  ) or  (  resi " +str(residue2) +" and chain " + str(chain2) + " and "  + str(pdb)+")"
	cmd.do(selection)
# gamma_turns, strands ,helices, beta_turns  selection
def selectBySequence(name,sequences,chain,pdb):
	if chain=="":
		selection= "select "+ name   +" , "
		for sequence in sequences:
			if sequence != sequences[-1]:
				selection+= "( ( pepseq "+str(sequence.upper()) + ") and "+str(pdb)+" ) or"
			else:
				selection+= "( ( pepseq "+str(sequence.upper()) + ") and "+str(pdb)+" )"
	else:
		selection= "select "+ name   +" , "
		for sequence in sequences:
			if sequence != sequences[-1]:
				selection+= "( ( pepseq "+str(sequence.upper()) + ") and chain "+str(chain)+" and "+str(pdb)+" ) or"
			else:
				selection+= "( ( pepseq "+str(sequence.upper()) + ") and chain "+str(chain)+" and "+str(pdb)+" )"
	cmd.do(selection)
def selectTwoResiduesMotifs(name,residue1,residue2,chain,pdb):
	if chain=="":
		selection= "select "+ name   +" , (resi " +  str(residue1) + "-" +str(residue2) + " and " + pdb+" )"
	else:
		selection= "select "+ name   +" , (resi " +  str(residue1) + "-" +str(residue2) + " and chain " + str(chain) +" and " + pdb+" )"
	cmd.do(selection)
def selectStrandsMotifs(name,strands,pdb,chain):
	selection= "select "+ name   +" , "
	for strand in strands:
		resi1 = strand["Residues"]["0"][0]
		resi2 = strand["Residues"]["0"][1]
		selection+=	"(resi " +  str(resi1) + "-" +str(resi2) + " and chain " + str(chain) + " and "+pdb +")"
	cmd.do(selection)
def selectHelixInteractionMotifs(name,helicesInter,pdb,chain):
	selection= "select "+ name   +" , "
	for helInter in helicesInter:
		resi1 = helInter["Residues"]["0"][0]
		resi2 = helInter["Residues"]["0"][1]
		selection+=	"(resi " +  str(resi1) + "-" +str(resi2) + " and chain " + str(chain) + " and "+pdb +")"
	cmd.do(selection)
# groupping any two elements
def group(elemA, elemB):
	groupping="group "+str(elemA)+" , "+str(elemB)+" , add"
	cmd.do(groupping)
###########################################################################################################
def structuralMotifs(pdb_id):
	global similarityChain
	global visitedChains
	global motifDicti
	pdb_id=pdb_id.lower()
	returnCodeIni=fetchPDB(pdb_id)
	visitedChains=[]
	if returnCodeIni==3:
		return 1
	motifDicti=make_dictio()
	promotifText = getPromotifData(pdb_id)
	if promotifText=="":
		return 1
	result =  AnalyzePromotifData(promotifText,pdbId=pdb_id) 
	chains = findChains(pdb_id)
	for chain in chains[1]:
		group(pdb_id+"."+chain,pdb_id+"."+chain+".structural_motifs")
		selectChain(pdb_id+"."+chain,pdb_id,chain)
		group(pdb_id+".chains", pdb_id+"."+chain)
		for key in motifDicti.keys():
			motifDicti[key][chain]=[]
	group("protein_"+ pdb_id, pdb_id+".chains")
	for motifs in result:
		if motifs == "beta_hairpins":
			representBeta_Hairpins(result,pdb_id)
		if motifs == "gamma_turns":
			representGamma_Turns(result,pdb_id)
		if motifs == "strands":
			representStrands(result,pdb_id)
		if motifs == "helices":
			representHelices(result,pdb_id)
		if motifs == "beta_turns":
			representBeta_Turns(result,pdb_id)
		if motifs == "bab_motifs":
			representBab_Motifs(result,pdb_id)
		if motifs == "psi_loops":
			representPsi_Loops(result,pdb_id)
		if motifs == "disulphides":
			representDisulphides(result,pdb_id)
		if motifs == "sheets":
			representSheets(result,pdb_id)
		if motifs == "beta_bulges":
			representBeta_Bulges(result,pdb_id)
		if motifs == "helix_interactions":
			representHelix_Interactions(result,pdb_id)
	visitedChains=list(set(visitedChains))
	visitedChains.sort()
	chainsToDel=[]
	for chain in chains[1]:
		if chain not in visitedChains:
			chainsToDel.append(chain)
			cmd.delete(pdb_id+"."+chain+".structural_motifs")
	#cmd.do("order *, yes")
	return 0

def structuralMotifsOffline(pdb_id):
	global similarityChain
	global countGroupping
	global visitedChains
	global motifDicti
	pdb_id=pdb_id.lower()
	returnCodeIni = fetchPDB(pdb_id)
	if returnCodeIni==3:
		return 1
	result={}
	chains=[]
	visitedChains=[]
	motifDicti=make_dictio()
	with open("PyProtif/"+pdb_id.lower()+"/Structural_motifs/"+pdb_id.lower()+"_promotif_Analyzed.json") as json_file:
		result = json.load(json_file)
		json_file.close()
	findChains(pdb_id)
	chains=list(set(readChains(pdb_id)[1]))
	chains.sort()
	for chain in chains:
		countGroupping-=1
		group(pdb_id+"."+chain,pdb_id+"."+chain+".structural_motifs")
		groupChainSelection(pdb_id+".chains", pdb_id+"."+chain,pdb_id,chain,countGroupping)
		for key in motifDicti.keys():
			motifDicti[key][chain]=[]
	group("protein_"+ pdb_id, pdb_id+".chains")
	for motifs in result:
		if motifs == "beta_hairpins":
			representBeta_Hairpins(result,pdb_id)
		if motifs == "gamma_turns":
			representGamma_Turns(result,pdb_id)
		if motifs == "strands":
			representStrands(result,pdb_id)
		if motifs == "helices":
			representHelices(result,pdb_id)
		if motifs == "beta_turns":
			representBeta_Turns(result,pdb_id)
		if motifs == "bab_motifs":
			representBab_Motifs(result,pdb_id)
		if motifs == "psi_loops":
			representPsi_Loops(result,pdb_id)
		if motifs == "disulphides":
			representDisulphides(result,pdb_id)
		if motifs == "sheets":
			representSheets(result,pdb_id)
		if motifs == "beta_bulges":
			representBeta_Bulges(result,pdb_id)
	#cmd.do("order *, yes")
	visitedChains=list(set(visitedChains))
	visitedChains.sort()
	chainsToDel=[]
	for chain in chains:
		if chain not in visitedChains:
			chainsToDel.append(chain)
			cmd.delete(pdb_id+"."+chain+".structural_motifs")
	return 0

def summarizeStruct(pdb_id):
	global similarityChain
	returnCodeIni=fetchPDB(pdb_id)
	if returnCodeIni==3:
		return 1
	chains=readChains(pdb_id)[1]
	chains.sort()
	processedChains=[]
	stringOutput=""
	result={}
	with open("PyProtif/"+pdb_id.lower()+"/Structural_motifs/"+pdb_id.lower()+"_promotif_Analyzed.json") as json_file:
		result = json.load(json_file)
		json_file.close()
	for chain in chains:
		for SimilarChains in similarityChain[chain]:
			if SimilarChains not in processedChains:
				processedChains.append(SimilarChains)
				countTotal=0
				string="\nChain "+str(SimilarChains)+":\t\t'count_total' motifs found\n"
				for motif in result.keys():
					count=0
					string+="\t"+motif+":\t\t"
					for motifElement in result[motif]:
						if SimilarChains in motifElement["Chain"]:
							count+=1
					countTotal+=count
					string+=str(count)+"\n"
				string=string.replace("'count_total'",str(countTotal))
				stringOutput+=string
	return stringOutput



###########################################################################################################
def functionalMotifs(pdb_id,eValPfam=1.0,eValCdd=1.0,database="pfam"):
	global driver
	# we have never ran this search before then we call this method ( if the file is preexisting we directly skip to process the result file instead)
	##7RN1_pfam_cdd_1_0.2.txt example
	pdb_id = pdb_id.lower()
	eValPfam=float(eValPfam)
	eValCdd=float(eValCdd)
	createDirectory("PyProtif/"+pdb_id+"/")
	createDirectory("PyProtif/"+pdb_id+"/Functional_motifs")
	options = webdriver.ChromeOptions()
	prefs = {"download.default_directory": os.getcwd()+fileSeparator+"PyProtif"+fileSeparator+pdb_id+fileSeparator+"Functional_motifs","download.directory_upgrade": True}
	options.add_experimental_option("prefs",prefs)
	if "window" in system.lower():
		service=Service()
		service.creationflags = CREATE_NO_WINDOW
		driver = webdriver.Chrome(options=options,service=service)
	else:
		service=Service()
		driver = webdriver.Chrome(options=options,service=service)
	driver.maximize_window()
	driver.set_window_position(-10000,0)
	driver.get("https://www.genome.jp/tools/motif/")
	fastaTextField= driver.find_element(by=By.XPATH,value="/html/body/table[2]/tbody/tr[2]/td/form/table/tbody/tr[3]/td/table/tbody/tr[3]/td[2]/textarea")
	fasta_seq = getFasta(pdb_id)
	fastaTextField.clear()
	fastaTextField.send_keys(fasta_seq)
	
	pfamCheck= driver.find_element(by=By.XPATH,value="/html/body/table[2]/tbody/tr[2]/td/form/table/tbody/tr[4]/td/table/tbody/tr[2]/td[1]/input")
	cddCheck= driver.find_element(by=By.XPATH,value="/html/body/table[2]/tbody/tr[2]/td/form/table/tbody/tr[4]/td/table/tbody/tr[3]/td[1]/input")
	pfamEval=driver.find_element(by=By.XPATH,value="/html/body/table[2]/tbody/tr[2]/td/form/table/tbody/tr[4]/td/table/tbody/tr[2]/td[3]/input")
	cddEval=driver.find_element(by=By.XPATH,value="/html/body/table[2]/tbody/tr[2]/td/form/table/tbody/tr[4]/td/table/tbody/tr[3]/td[3]/input")
	submitButton=driver.find_element(by=By.XPATH,value="/html/body/table[2]/tbody/tr[2]/td/form/table/tbody/tr[2]/td/input[1]")
	if database != "pfam":
		if "pfam" not in database:
			#ncbi-cdd
			pfamCheck.click()
			cddCheck.click()
			cddEval.clear()
			cddEval.send_keys(eValCdd)
		else:
			cddCheck.click()
			cddEval.clear()
			cddEval.send_keys(eValCdd)
			pfamEval.clear()
			pfamEval.send_keys(eValPfam)
	else:
		pfamEval.clear()
		pfamEval.send_keys(eValPfam)
	submitButton.click()
	time.sleep(2)
	try: 
		motifCounter= driver.find_element(by=By.XPATH,value="/html/body/b[1]").text.split(" ")[-1].strip()
		if "0" == motifCounter:
			return {}
		results=driver.find_element(by=By.XPATH,value="/html/body/a")
		results.click()
		while not path.exists(os.getcwd()+"/PyProtif/"+pdb_id+"/Functional_motifs/result.txt"):
			time.sleep(1)
			#7RN1_pfam_cdd_1_0.2.txt
		os.rename("PyProtif/"+pdb_id+"/Functional_motifs/result.txt", "PyProtif/"+pdb_id+"/Functional_motifs/"+pdb_id+'_'+database.replace(",","_")+'_'+str(eValPfam)+'_'+str(eValCdd)+'.txt')
		#7RN1_pfam_cdd_1_0.2.png
		#with open("PyProtif/"+pdb_id+"/Functional_motifs/"+pdb_id+'_'+database.replace(",","_")+'_'+str(eValPfam)+'_'+str(eValCdd)+'.png', 'wb') as file:
		#	l = driver.find_element(by=By.XPATH,value="/html/body/img")
		#	file.write(l.screenshot_as_png)
		return processResultsTxt(pdb_id,eValPfam,eValCdd,database)
	except:
		try:
			invalidFastaSequence = driver.find_element(by=By.XPATH,value="/html/body/div/table/tbody/tr/td/span/h3")
			header = invalidFastaSequence.text.split("(")[1].split(")")[0]
			if "must be a protein sequence." in invalidFastaSequence.text:
				messagebox.showerror("Error","Error processing the pdb: "+pdb_id+" for functional motifs due to a none protein fasta sequence found with the following header "+header)
		except:
			messagebox.showerror("Error","Error processing the pdb: "+pdb_id+" for functional motifs")

def processResultsTxt(pdb_id,eValPfam=1.0,eValCdd=1.0,database="pfam"):
	# this method should work with local files too if they exist (we skip the method functionalMotifs then)
	global countGroupping
	pdb_id=pdb_id.lower()
	returnCodeIni = fetchPDB(pdb_id)
	if returnCodeIni == 3:
		return 1
	Adjustment=adjustForResidueDifference(pdb_id)
	with open("PyProtif/"+pdb_id+"/Functional_motifs/"+pdb_id+'_'+database.replace(",","_")+'_'+str(eValPfam)+'_'+str(eValCdd)+'.txt') as f:
		functionalMotifs=OrderedDict()
		functionalMotifs["pfam"]=[]
		functionalMotifs["ncbi-cdd"]=[]
		for line in f:
			if "#" in line:
				pass
			else:
				info=line.split("	")
				data={}
				chains = info[0].split("_")[1:]
				chains.sort()
				chainsInString=""
				for elem in chains:
					if elem==chains[-1]:
						chainsInString+=str(elem)
					else:
						chainsInString+=str(elem)+"_"
				data["chain"]=chainsInString
				data["ID"]=pdb_id+"_"+chainsInString+"/"+info[2].strip()
				data["database"]=info[1].strip()
				data["name"]=info[2].strip()
				data["sequence"]=info[4].replace("-","").upper().strip()
				data["score"]=info[6].strip()
				data["Eval"]=info[7].strip()
				data["chainSequence"] = getFastaChain(pdb_id,chains[0],initial=1)
				#adjust positions if necessary

				newposition = info[3].strip()
				arrayOfPositions=[]
				if "," in info[3]:#in case a motif is in multiple positions	
					twoPositions=newposition.split(",")
					newPositions=""
					for position in twoPositions:
						for key in Adjustment.keys():
							if key in chains:
								start=int(position.strip().split("..")[0])+int(Adjustment[key])-1
								end=int(position.strip().split("..")[1])+int(Adjustment[key])-1
						arrayOfPositions.append([str(start),str(end)])
					for position in arrayOfPositions:
						if position==arrayOfPositions[-1]:
							newPositions+=str(position[0])+".."+str(position[1])
						else:
							newPositions+=str(position[0])+".."+str(position[1])+","

					data["position"]=newPositions
				else:
					for key in Adjustment.keys():
						if key in chains:
							start=int(info[3].strip().split("..")[0])+int(Adjustment[key])-1
							end=int(info[3].strip().split("..")[1])+int(Adjustment[key])-1
							newposition=str(start)+".."+str(end)
					data["position"]=newposition
				sequences=[]
				splitSequencesbycomma = info[4].split(",")
				splitSequencesbycomma = [x for x in splitSequencesbycomma if x]
				splitSequencesByDash=[]
				for elem in splitSequencesbycomma:
					intermed = elem.split("-")
					for elem2 in intermed:
						sequences.append(elem2.strip())
				sequences = [x.upper() for x in sequences if x]
				functionalMotifs[data["database"]].append(data)
				if len(arrayOfPositions)==0:
					for chain in chains:
						selector  = pdb_id+"."+str(chain)+ "."+str(info[2].strip())+ "."+ data["position"].split("..")[0] +"_" +  data["position"].split("..")[1]
						selectBySequence(selector,sequences,chain,pdb_id)
						if data["database"]=="pfam":
							group(pdb_id+ "."+chain+".pfam.functional_motifs",selector)
						else:
							group(pdb_id+ "."+chain+".cdd.functional_motifs",selector)
						group(pdb_id+ "."+chain+".functional_motifs",pdb_id+ "."+chain+".pfam.functional_motifs")
						group(pdb_id+ "."+chain+".functional_motifs",pdb_id+ "."+chain+".cdd.functional_motifs")
						group(pdb_id+ "."+chain,pdb_id+ "."+chain+".functional_motifs")
						groupChainSelection(pdb_id+".chains", pdb_id+"."+chain,pdb_id,chain,countGroupping)
						countGroupping-=1
				else:
					#selectMultiResiduesMotifs
					for chain in chains:
						selector  = pdb_id+"."+str(chain)+ "."+str(info[2].strip())+ "."
						for position in arrayOfPositions:
							selector+="."+str(position[0])+"_"+str(position[1])
						selectBySequence(selector,sequences,chain,pdb_id)
						if data["database"]=="pfam":
							group(pdb_id+ "."+chain+".pfam.functional_motifs",selector)
						else:
							group(pdb_id+ "."+chain+".cdd.functional_motifs",selector)
						group(pdb_id+ "."+chain+".functional_motifs",pdb_id+ "."+chain+".pfam.functional_motifs")
						group(pdb_id+ "."+chain+".functional_motifs",pdb_id+ "."+chain+".cdd.functional_motifs")
						group(pdb_id+ "."+chain,pdb_id+ "."+chain+".functional_motifs")
						groupChainSelection(pdb_id+".chains", pdb_id+"."+chain,pdb_id,chain,countGroupping)	
						countGroupping-=1
				group("protein_"+ pdb_id, pdb_id+".chains")
				json_string = json.dumps(functionalMotifs)
				with open("PyProtif/"+pdb_id+"/Functional_motifs/"+pdb_id+'_'+database.replace(",","_")+'_'+str(eValPfam)+'_'+str(eValCdd)+'_Analyzed_Motifs.json', 'w') as outfile:
					outfile.write(json_string)
					outfile.close()
	return 0

def searchForFunctionalMotifs(keyword):
	global driver
	#https://pfam.xfam.org/family/PF00104#tabview=tab8
	#http://pfam.xfam.org/search#tabview=tab2
	#//*[@id="fstrucph"]
	keyword=keyword.lower()
	if "window" in system.lower():
		service=Service()
		service.creationflags = CREATE_NO_WINDOW
		driver = webdriver.Chrome(service=service)
	else:
		service=Service()
		driver = webdriver.Chrome(service=service)
	driver.maximize_window()
	driver.set_window_position(-10000,0)
	driver.get("http://pfam-legacy.xfam.org/search")
	keywordSearch=WebDriverWait(driver, 20).until(EC.presence_of_element_located((By.XPATH, "//*[@id=\"searchKeywordBlockSelector\"]/a")))
	keywordSearch.click()
	keywordField=driver.find_element(by=By.XPATH,value="//*[@id=\"kwQuery\"]")
	keywordField.clear()
	keywordField.send_keys(keyword)
	submitButton= WebDriverWait(driver, 20).until(EC.presence_of_element_located((By.XPATH, "//*[@id=\"keywordSearchForm\"]/div[2]/input[1]")))
	submitButton.click()
	try:
		noResult=driver.find_element(by=By.XPATH, value="/html/body/div[5]/h1")
		if noResult.text.lower()=="no results":
			json_string = json.dumps({keyword:["no pfam found"]})
			with open("PyProtif/"+keyword+"_pfamscan.json", 'w') as outfile:
				outfile.write(json_string)
				outfile.close()
			driver.close()
			return ""
	except:
		try:
			#//*[@id="pdbBlockSelector"]/a
			structuresTab=  WebDriverWait(driver, 20).until(EC.presence_of_element_located((By.XPATH, "//*[@id=\"pdbBlockSelector\"]/a")))
			try:
				structuresTab.click()
			except:
				fmSearch[pfamId]=["no stucture was found"]
				pfamId=driver.current_url.split("/")[-1].split("#")[0]
				#"https://pfam.xfam.org/family/"+str(pfamId)+"#tabview=tab8)" 
				#fmSearch["Link"]= "https://www.genome.jp/dbget-bin/www_bget?pf:"+pfamId
				fmSearch["Link"]="https://pfam-legacy.xfam.org/family/"+str(pfamId)+"#tabview=tab8"
				json_string = json.dumps(fmSearch)
				with open("PyProtif/"+keyword+"_pfamscan.json", 'w') as outfile:
					outfile.write(json_string)
					outfile.close()
				driver.close()
				return "https://pfam-legacy.xfam.org/family/"+str(pfamId)+"#tabview=tab8"
			fmSearch={}
			pfamId=driver.current_url.split("/")[-1].split("#")[0]
			fmSearch[pfamId]=[]
			structures=[]
			#WebDriverWait(driver, 20).until(EC.presence_of_element_located((By.XPATH, "//*[@id=\"structuresTable\"]/tbody/tr"))
			loading = driver.find_element(by=By.XPATH, value="//*[@id=\"fstrucph\"]")
			while loading:
				try:
					loading = driver.find_element(by=By.XPATH, value="//*[@id=\"fstrucph\"]")
				except:
					break
			structuresRow = driver.find_elements(by=By.XPATH, value="//*[@id=\"structuresTable\"]/tbody/tr")
			for structureRow in range(0,len(structuresRow)):
				if len(driver.find_elements(by=By.XPATH, value="//*[@id=\"structuresTable\"]/tbody/tr["+str(structureRow+1)+"]/td")) < 6:
					pass
				else:
					stucture=driver.find_element(by=By.XPATH, value="//*[@id=\"structuresTable\"]/tbody/tr["+str(structureRow+1)+"]/td[1]/a")
					stucture2=driver.find_element(by=By.XPATH, value="//*[@id=\"structuresTable\"]/tbody/tr["+str(structureRow+1)+"]/td[1]")
					actions=ActionChains(driver)
					actions.move_to_element(stucture2).perform()
					structures.append(stucture.get_attribute("text"))
			fmSearch[pfamId]=structures
			fmSearch["Link"]="https://pfam-legacy.xfam.org/family/"+str(pfamId)+"#tabview=tab8"
			json_string = json.dumps(fmSearch)
			with open("PyProtif/"+keyword+"_pfamscan.json", 'w') as outfile:
				outfile.write(json_string)
				outfile.close()
			driver.close()
			return "https://pfam-legacy.xfam.org/family/"+str(pfamId)+"#tabview=tab8"
			#https://www.genome.jp/dbget-bin/www_bget?pf:PF17990
		except:
			try:
				intermediate=driver.find_element(by=By.XPATH, value="//*[@id=\"resultTable\"]/tbody/tr[1]/td[2]/a")
				intermediate.click()
				structuresTab=  WebDriverWait(driver, 20).until(EC.presence_of_element_located((By.XPATH, "//*[@id=\"pdbBlockSelector\"]/a")))
				structuresTab.click()
			except:
				fmSearch[pfamId]=["no stucture was found"]
				pfamId=driver.current_url.split("/")[-1].split("#")[0]
				fmSearch["Link"]= "https://pfam-legacy.xfam.org/family/"+str(pfamId)+"#tabview=tab8"
				json_string = json.dumps(fmSearch)
				with open("PyProtif/"+keyword+"_pfamscan.json", 'w') as outfile:
					outfile.write(json_string)
					outfile.close()
				driver.close()
				return "https://pfam-legacy.xfam.org/family/"+str(pfamId)+"#tabview=tab8"
			fmSearch={}
			pfamId=driver.current_url.split("/")[-1].split("#")[0]
			fmSearch[pfamId]=[]
			structures=[]
			loading = driver.find_element(by=By.XPATH, value="//*[@id=\"fstrucph\"]")
			while loading:
				try:
					loading = driver.find_element(by=By.XPATH, value="//*[@id=\"fstrucph\"]")
				except:
					break
			structuresRow = driver.find_elements(by=By.XPATH, value="//*[@id=\"structuresTable\"]/tbody/tr")
			for structureRow in range(0,len(structuresRow)):
				if len(driver.find_elements(by=By.XPATH, value="//*[@id=\"structuresTable\"]/tbody/tr["+str(structureRow+1)+"]/td")) < 6:
					pass
				else:
					stucture=driver.find_element(by=By.XPATH, value="//*[@id=\"structuresTable\"]/tbody/tr["+str(structureRow+1)+"]/td[1]/a")
					stucture2=driver.find_element(by=By.XPATH, value="//*[@id=\"structuresTable\"]/tbody/tr["+str(structureRow+1)+"]/td[1]")
					actions=ActionChains(driver)
					actions.move_to_element(stucture2).perform()
					structures.append(stucture.get_attribute("text"))
			fmSearch[pfamId]=structures
			fmSearch["Link"]= "https://pfam-legacy.xfam.org/family/"+str(pfamId)+"#tabview=tab8"
			json_string = json.dumps(fmSearch)
			with open("PyProtif/"+keyword+"_pfamscan.json", 'w') as outfile:
				outfile.write(json_string)
				outfile.close()
			driver.close()
			return "https://pfam-legacy.xfam.org/family/"+str(pfamId)+"#tabview=tab8"
	
def loadweb(url):
	new = 1
	webbrowser.open(url,new=new)
def openweb():
	global url
	new = 1
	webbrowser.open(url,new=new)
def summarizeFMSearch(keyword):
	keyword=keyword.lower()
	previous={}
	with open("PyProtif/"+keyword+"_pfamscan.json") as json_file:
			previous = json.load(json_file)
			json_file.close()
	url=""
	string =""
	for key in previous.keys():
		if key == "Link":
			url=previous[key]
		else:
			summary="\nPFAM:\t\t"+key
			for elem in previous[key]:
				summary+="\n\t"+str(elem)
			string+=summary
	return [string,url]
def summarizeFunc(pdb_id,eValPfam=1.0,eValCdd=1.0,database="pfam"):
	#ncbi-cdd
	result={}
	with open("PyProtif/"+pdb_id+"/Functional_motifs/"+pdb_id+'_'+database.replace(",","_")+'_'+str(eValPfam)+'_'+str(eValCdd)+'_Analyzed_Motifs.json') as json_file:
		result = json.load(json_file)
		json_file.close()
	summary=""
	for database in result.keys():
		if len(result[database])>0:
			summaryMin=database+":"
			for funct in result[database]:
				summaryMin+="\n\tID: "+str(funct["ID"])+"\n\t\tPosition: "+str(funct["position"])+"\n\t\tScore: "+str(funct["score"])+"\n\t\tEval: "+str(funct["Eval"])
			summary+=summaryMin+"\n"
	return summary
def getFunctionalMotifOrdering(pdb_id,eValPfam=1.0,eValCdd=1.0,database="pfam"):
	result={}
	with open("PyProtif/"+pdb_id+"/Functional_motifs/"+pdb_id+'_'+database.replace(",","_")+'_'+str(eValPfam)+'_'+str(eValCdd)+'_Analyzed_Motifs.json') as json_file:
		result = json.load(json_file)
		json_file.close()
	summary={}
	for database in result.keys():
		summary[database]=""
		if len(result[database])>0:
			for funct in result[database]:
				summary[database]+="*"+str(pdb_id)+".*."+str(funct["name"])+".* "
	return summary
def getData(pdb_id,database,eValPfam,eValCdd):
	result=OrderedDict()
	with open("PyProtif/"+pdb_id+"/Functional_motifs/"+pdb_id+'_'+database.replace(",","_")+'_'+str(eValPfam)+'_'+str(eValCdd)+'_Analyzed_Motifs.json') as json_file:
		result = json.load(json_file)
		json_file.close()
	return result
def getChainsInfo(pdb_id):
	result={}
	finalResult={}
	f =  open("PyProtif/"+pdb_id+"/"+pdb_id+'_adjusted.fasta')
	for line in f:
		if ">" in line:
			chain = line[6:]
			chains = sorted(chain.split("_"))
			strKey=""
			for chaini in chains:
				if chaini ==chains[-1]:
					strKey+=chaini.strip()
				else:
					strKey+=chaini.strip()+"_"
			result[strKey]=0
		else:
			result[strKey]=[len(line.strip()),line.strip()]
	return result
def minimizeString(string):
	adjustedStr=""
	for index in range(0,len(string)):
		if (index+1)%40==0:
			adjustedStr+=string[index]+"<br>"
		else:
			adjustedStr+=string[index]
	return adjustedStr


		 
########################################### End of BackEnd ##############################################

########################################### Start of FrontEnd ###########################################
def checkValidPDB(pdb_id):
	global validPDB
	url = "https://www.rcsb.org/fasta/entry/{:s}/display".format(pdb_id)
	try:
		response = urllib.request.urlopen(url)
		fastaFile = response.read().decode('utf-8')
		if "No valid PDB IDs were submitted." == fastaFile.strip():
			validPDB=False
		else:
			validPDB=True
	except:
		validPDB=False
def searchFileLocaly(String,path=""):
	global decision
	allDatas=[]
	foundLocally = False
	for allInfo in os.scandir("PyProtif"+path):
		if not allInfo.is_dir():
			allDatas.append(allInfo.path)
	for allData in allDatas:
		if String in allData:
			if decision!=0:
				foundLocally = True
	return foundLocally	
def searchLocaly(String,path=""):
	global decision
	allDatas=[]
	foundLocally = False
	for allInfo in os.scandir("PyProtif"+path):
		if allInfo.is_dir():
			allDatas.append(allInfo.path)
	for allData in allDatas:
		if String in allData:
			if decision!=0:
				foundLocally = True
	return foundLocally
def cascade():
	global foundPDBLocally
	global foundPDBStructLocally
	global foundPDBFunctLocally
	global summary
	global funcSummary
	global runCounts
	global driver
	global countGroupping
	global similarityChain
	foundPDBLocally=False
	foundPDBStructLocally=False
	foundPDBFunctLocally=False
	if pdbIDField.get()=="":
		messagebox.showerror("Error", "Please enter at least one PDB ID to start")
	else:
		cmd.delete("unselect")
		ordering={}
		ordering["ncbi-cdd"]=""
		ordering["pfam"]=""
		pdb_ids = pdbIDField.get().split(" ")
		pdb_ids=list(set([pdb_i.lower() for pdb_i in pdb_ids]))
		pdb_ids.sort()
		invalidPDBs=[]
		if  searchFileLocaly("invalidPdbList.txt",path=""):
			f= open("PyProtif"+fileSeparator+"invalidPdbList.txt",'r')
			for line in f:
				if line.strip() in pdb_ids:
					invalidPDBs.append(line.strip())
					pdb_ids.remove(line.strip())
			f.close()
		for pdb_id in pdb_ids:
			countGroupping=1
			pdb_id=pdb_id.strip().lower()
			createDirectory("PyProtif/"+pdb_id)
			createDirectory("PyProtif/"+pdb_id+"/Structural_motifs")
			createDirectory("PyProtif/"+pdb_id+"/Functional_motifs")
			fetchPDB(pdb_id)
			if not validPDB:
				invalidPDBs.append(pdb_id)
				continue
			similarityChain=getSimilarChains(pdb_id)
			for key in similarityChain.keys():
				similarityChain[key].sort()
			foundPDBLocally=searchLocaly(pdb_id)
			if foundPDBLocally:
				foundPDBStructLocally = searchLocaly("Structural_motifs",path="/"+pdb_id)
				foundPDBFunctLocally = searchLocaly("Functional_motifs",path="/"+pdb_id)
				if decision==1 or decision==3:
					if foundPDBStructLocally:
						foundPDBStructLocalActualy = searchFileLocaly(pdb_id+"_promotif_Analyzed.json",path="/"+pdb_id+"/Structural_motifs")
						if foundPDBStructLocalActualy:
							code2=structuralMotifsOffline(pdb_id)
						else:
							code2=structuralMotifs(pdb_id)
					else:
						code2=structuralMotifs(pdb_id)
					if code2==0:
						summary=summarizeStruct(pdb_id)
						summaryPopUp("Structural Motifs: "+pdb_id,summary,showButt=False)
					else:
						runCounts-=1
						messagebox.showerror("Error", "Error loading the pdb: "+pdb_id+" file")
				if decision==2 or decision==3:
					if foundPDBFunctLocally:
						imageFile=""
						if decisionDB==1:
							foundPDBFunctLocalActualy = searchFileLocaly(pdb_id+'_pfam_'+str(float(pfamEvalField.get()))+'_1.0.txt',path="/"+pdb_id+"/Functional_motifs")
							if foundPDBFunctLocalActualy:
								code = processResultsTxt(pdb_id,eValPfam=float(pfamEvalField.get()))
								ordering=getFunctionalMotifOrdering(pdb_id,eValPfam=float(pfamEvalField.get()))
							else:
								code = functionalMotifs(pdb_id,eValPfam=float(pfamEvalField.get()))
								ordering=getFunctionalMotifOrdering(pdb_id,eValPfam=float(pfamEvalField.get()))
							if code==0:
								funcSummary=summarizeFunc(pdb_id,eValPfam=float(pfamEvalField.get()))
								ordering=getFunctionalMotifOrdering(pdb_id,eValPfam=float(pfamEvalField.get()))
								imageFile="PyProtif/"+pdb_id+"/Functional_motifs/"+pdb_id+'_pfam_'+str(float(pfamEvalField.get()))+'_1.0.png'
						if decisionDB==2:
							foundPDBFunctLocalActualy = searchFileLocaly(pdb_id+'_cdd_1.0_'+str(float(cddEvalField.get()))+'.txt',path="/"+pdb_id+"/Functional_motifs")
							if foundPDBFunctLocalActualy:
								code = processResultsTxt(pdb_id,eValCdd=float(cddEvalField.get()),database="cdd")
								ordering=getFunctionalMotifOrdering(pdb_id,eValCdd=float(cddEvalField.get()),database="cdd")
							else:
								code =functionalMotifs(pdb_id,eValCdd=float(cddEvalField.get()),database="cdd")
								ordering=getFunctionalMotifOrdering(pdb_id,eValCdd=float(cddEvalField.get()),database="cdd")
							if code==0:
								funcSummary=summarizeFunc(pdb_id,eValCdd=float(cddEvalField.get()),database="cdd")
								ordering=getFunctionalMotifOrdering(pdb_id,eValCdd=float(cddEvalField.get()),database="cdd")
								imageFile="PyProtif/"+pdb_id+"/Functional_motifs/"+pdb_id+'_cdd_1.0_'+str(float(cddEvalField.get()))+'.png'
						if decisionDB==3:
							foundPDBFunctLocalActualy = searchFileLocaly(pdb_id+'_pfam_cdd_'+str(float(pfamEvalField.get()))+'_'+str(float(cddEvalField.get()))+'.txt',path="/"+pdb_id+"/Functional_motifs")
							if foundPDBFunctLocalActualy:
								code = processResultsTxt(pdb_id,eValPfam=float(pfamEvalField.get()),eValCdd=float(cddEvalField.get()),database="pfam,cdd")
								ordering=getFunctionalMotifOrdering(pdb_id,eValPfam=float(pfamEvalField.get()),eValCdd=float(cddEvalField.get()),database="pfam,cdd")
							else:
								code = functionalMotifs(pdb_id,eValPfam=float(pfamEvalField.get()),eValCdd=float(cddEvalField.get()),database="pfam,cdd")
								ordering=getFunctionalMotifOrdering(pdb_id,eValPfam=float(pfamEvalField.get()),eValCdd=float(cddEvalField.get()),database="pfam,cdd")
							if code ==0:
								funcSummary=summarizeFunc(pdb_id,eValPfam=float(pfamEvalField.get()),eValCdd=float(cddEvalField.get()),database="pfam,cdd")
								ordering=getFunctionalMotifOrdering(pdb_id,eValPfam=float(pfamEvalField.get()),eValCdd=float(cddEvalField.get()),database="pfam,cdd")
								imageFile="PyProtif/"+pdb_id+"/Functional_motifs/"+pdb_id+'_pfam_cdd_'+str(float(pfamEvalField.get()))+'_'+str(float(cddEvalField.get()))+'.png'
					else:
						imageFile=""
						if decisionDB==1:
							code = functionalMotifs(pdb_id,eValPfam=float(pfamEvalField.get()))
							ordering=getFunctionalMotifOrdering(pdb_id,eValPfam=float(pfamEvalField.get()))
							if code ==0:
								funcSummary=summarizeFunc(pdb_id,eValPfam=float(pfamEvalField.get()))
								ordering=getFunctionalMotifOrdering(pdb_id,eValPfam=float(pfamEvalField.get()))
								imageFile="PyProtif/"+pdb_id+"/Functional_motifs/"+pdb_id+'_pfam_'+pfamEvalField.get()+'_1.0.png'
						if decisionDB==2:
							code = functionalMotifs(pdb_id,eValCdd=float(cddEvalField.get()),database="cdd")
							ordering=getFunctionalMotifOrdering(pdb_id,eValCdd=float(cddEvalField.get()),database="cdd")
							if code ==0:
								funcSummary=summarizeFunc(pdb_id,eValCdd=float(cddEvalField.get()),database="cdd")
								ordering=getFunctionalMotifOrdering(pdb_id,eValCdd=float(cddEvalField.get()),database="cdd")
								imageFile="PyProtif/"+pdb_id+"/Functional_motifs/"+pdb_id+'_cdd_1.0_'+cddEvalField.get()+'.png'
						if decisionDB==3:
							code = functionalMotifs(pdb_id,eValPfam=float(pfamEvalField.get()),eValCdd=float(cddEvalField.get()),database="pfam,cdd")
							ordering=getFunctionalMotifOrdering(pdb_id,eValPfam=float(pfamEvalField.get()),eValCdd=float(cddEvalField.get()),database="pfam,cdd")
							if code ==0:
								funcSummary=summarizeFunc(pdb_id,eValPfam=float(pfamEvalField.get()),eValCdd=float(cddEvalField.get()),database="pfam,cdd")
								ordering=getFunctionalMotifOrdering(pdb_id,eValPfam=float(pfamEvalField.get()),eValCdd=float(cddEvalField.get()),database="pfam,cdd")
								imageFile="PyProtif/"+pdb_id+"/Functional_motifs/"+pdb_id+'_pfam_cdd_'+pfamEvalField.get()+'_'+cddEvalField.get()+'.png'
					if code ==0:
						summaryPopUp("Functional Motifs: "+pdb_id,funcSummary,notShowImage=False,imgFile=imageFile,showButt=False)
					else:
						runCounts-=1
						messagebox.showerror("Error", "Error loading the pdb:  "+pdb_id+" file")
				if decision==0:
					messagebox.showwarning("Warning", "Please select at least one of the options")
			else:
				checkValidPDB(pdb_id)
				if not validPDB:
					invalidPDBs.append(pdb_id)
					continue
				else:
					#validpdb
					if decision==1 or decision==3:
						code = structuralMotifs(pdb_id)
						if code ==0:
							summary=summarizeStruct(pdb_id)
							summaryPopUp("Structural Motifs: "+pdb_id,summary,showButt=False)
						else:
							runCounts-=1
							messagebox.showerror("Error", "Error loading the pdb:  "+pdb_id+" file")
					if decision==2 or decision==3:
						imageFile=""
						if decisionDB==1:
							code2=functionalMotifs(pdb_id,eValPfam=float(pfamEvalField.get()))
							ordering=getFunctionalMotifOrdering(pdb_id,eValPfam=float(pfamEvalField.get()))
							imageFile="PyProtif/"+pdb_id+"/Functional_motifs/"+pdb_id+'_pfam_'+pfamEvalField.get()+'_1.0.png'
							if code2==0:
								funcSummary=summarizeFunc(pdb_id,eValPfam=float(pfamEvalField.get()))
						if decisionDB==2:
							code2=functionalMotifs(pdb_id,eValCdd=float(cddEvalField.get()),database="cdd")
							ordering=getFunctionalMotifOrdering(pdb_id,eValCdd=float(cddEvalField.get()),database="cdd")
							imageFile="PyProtif/"+pdb_id+"/Functional_motifs/"+pdb_id+'_cdd_1.0_'+cddEvalField.get()+'.png'
							if code2==0:
								funcSummary=summarizeFunc(pdb_id,eValCdd=float(cddEvalField.get()),database="cdd")
						if decisionDB==3:
							code2=functionalMotifs(pdb_id,eValPfam=float(pfamEvalField.get()),eValCdd=float(cddEvalField.get()),database="pfam,cdd")
							ordering=getFunctionalMotifOrdering(pdb_id,eValPfam=float(pfamEvalField.get()),eValCdd=float(cddEvalField.get()),database="pfam,cdd")
							imageFile="PyProtif/"+pdb_id+"/Functional_motifs/"+pdb_id+'_pfam_cdd_'+pfamEvalField.get()+'_'+cddEvalField.get()+'.png'
							if code2==0:
								funcSummary=summarizeFunc(pdb_id,eValPfam=float(pfamEvalField.get()),eValCdd=float(cddEvalField.get()),database="pfam,cdd")
						if code2==0:
							summaryPopUp("Functional Motifs: "+pdb_id,funcSummary,notShowImage=False,imgFile=imageFile,showButt=False)
						else:
							runCounts-=1
							messagebox.showerror("Error", "Error loading the pdb:  "+pdb_id+" file")
					if decision==0:
						messagebox.showwarning("Warning", "Please select at least one of the options")
		#cmd.do("order *, yes")
			orderCddMotifs = ordering["ncbi-cdd"]
			orderPfamMotifs = ordering["pfam"]
			orderFunctionalMotifsDatabases = str(pdb_id)+"*.cdd.functional_* "+str(pdb_id)+"*.pfam.functional_*"
			orderStructMotifs= str(pdb_id)+"*.bab_motifs* "+str(pdb_id)+"*.beta_bulges* "+str(pdb_id)+"*.beta_hairpins* "+str(pdb_id)+"*.beta_turns* "+str(pdb_id)+"*.disulphides* "+str(pdb_id)+"*.gamma_turns* "+str(pdb_id)+"*.helices* "+str(pdb_id)+"*.helix_interactions* "+str(pdb_id)+"*.psi_loops* "+str(pdb_id)+"*.sheets* "+str(pdb_id)+"*.strands*"
			orderInnerStructMotifs=str(pdb_id)+"*.bab_motif.* "+str(pdb_id)+"*.beta_bulge.* "+str(pdb_id)+"*.beta_hairpin.* "+str(pdb_id)+"*.beta_turn.* "+str(pdb_id)+"*.disulphide.* "+str(pdb_id)+"*.gamma_turn.* "+str(pdb_id)+"*.helice.* "+str(pdb_id)+"*.helix_interaction.* "+str(pdb_id)+"*.psi_loop.* "+str(pdb_id)+"*.sheet.* "+str(pdb_id)+"*.strand.*"
			orderMotifCategory="*.structural_* *.functional_*"
			cmd.order("*, yes")
			cmd.order(orderMotifCategory)
			cmd.order(orderStructMotifs)
			cmd.order(orderInnerStructMotifs,"yes")
			cmd.order(orderFunctionalMotifsDatabases)
			cmd.order(orderPfamMotifs)
			cmd.order(orderCddMotifs)
		

		
		loadedPDBs=cmd.get_object_list('solvent')
		if len(loadedPDBs)>0:
			unselect = "select unselect, none"
			cmd.do(unselect)
		try:
			driver.close()
		except:
			print("")
		if len(invalidPDBs)>0:
			#we have some invalid pdbs
			stringWarning=" The following PDB IDs were not processed because they are invalid:\n"
			f = open("PyProtif"+fileSeparator+"invalidPdbList.txt","w")
			for invalidPDB in list(set(invalidPDBs)):
				stringWarning+=invalidPDB+" "
				f.write(invalidPDB+"\n")
				if searchLocaly(invalidPDB,path=fileSeparator):
					os.removedirs("PyProtif/"+invalidPDB+"/Functional_motifs")
					os.removedirs("PyProtif/"+invalidPDB+"/Structural_motifs")
			f.close()
			invalidPDB=[]
			messagebox.showwarning("Warning", stringWarning)

def get_selection():
	global decision
	if (structM.get() == 1) & (funcM.get() == 0):
		decision=1
		choicePFAM.config(state= "disabled")
		choiceCDD.config(state= "disabled")
		pfamEvalField.config(state= "disabled")
		cddEvalField.config(state= "disabled")
	elif (structM.get() == 0) & (funcM.get() == 1):
		decision=2
		choicePFAM.config(state="normal")
		choiceCDD.config(state="normal")
		pfamEvalField.config(state="normal")
		cddEvalField.config(state="normal")
	elif (structM.get() == 0) & (funcM.get() == 0):
		messagebox.showwarning("Warning", "Please select at least one of the options")
		decision=0
		structM.set(1)
		choicePFAM.config(state= "disabled")
		choiceCDD.config(state= "disabled")
		pfamEvalField.config(state= "disabled")
		cddEvalField.config(state= "disabled")
	else:
		decision=3
		choicePFAM.config(state="normal")
		choiceCDD.config(state="normal")
		pfamEvalField.config(state="normal")
		cddEvalField.config(state="normal")
def get_DBselection():
	global decisionDB
	if (pfamD.get() == 1) & (cddD.get() == 0):
		decisionDB=1
	elif (pfamD.get() == 0) & (cddD.get() == 1):
		decisionDB=2
	elif (pfamD.get() == 0) & (cddD.get() == 0):
		messagebox.showwarning("Warning", "Invalid configuration the default settings will be used: pfam with eval of 1.0\n")
		choicePFAM.toggle()
		pfamDEval.set("1.0")
		decisionDB=1
	else:
		decisionDB=3
def summaryPopUp(string,data,notShowImage=True,imgFile="",showButt=False):
	global runCounts
	global decision
	global photo
	global fileSeparator

	top= tk.Toplevel(root)
	top.wm_iconphoto(False, photo)
	top.geometry('550x700')
	top.title(string)
	Btn = Button(top, text = "More Info",command=openweb)
	Btn.pack()
	FrameBIG = Frame(top)
	Main = Canvas(FrameBIG,height = 550,width =600)
	Main.configure(scrollregion=Main.bbox("all"))

	text = Text(Main, state='disabled', width=50, height=500)
	scrollbar = ttk.Scrollbar(
		FrameBIG,
		orient='vertical',
		command=Main.yview
	)
	
	scrollbar.pack(side="right", fill="y")
	Main.pack()
	FrameBIG.pack(anchor = W, fill = "x")
	text.configure(state='normal')
	text.insert('end', data)
	text.configure(state='disabled')
	text.pack()
	text['yscrollcommand'] = scrollbar.set
	Main.configure(yscrollcommand=scrollbar.set)
	if not notShowImage and imgFile != "":
		dataImg = imgFile.split(fileSeparator)[-1].replace(".png","").strip().split("_")
		dataImg=dataImg[1:]
		dataImg[0]=dataImg[0].split("/")[-1]
		# if len(dataImg)>4:
		# 	runVisualizer(dataImg[0],eValPfam=dataImg[3],eValCdd=dataImg[4],database=dataImg[1]+","+dataImg[2])
		# else:
		# 	runVisualizer(dataImg[0],eValPfam=dataImg[2],eValCdd=dataImg[3],database=dataImg[1])
		if decision==3:
			runCounts+=1
	if not showButt:
		runCounts+=1
		Btn.pack_forget()
	Main.configure(scrollregion=Main.bbox("all"))

def runFMSearch():
	global url
	url=""
	if FMField.get()=="":
		messagebox.showerror("Error", "Please enter at least one Keyword to start")
	else:
		FMKeys = FMField.get().split(" ")
		for FMKey in FMKeys:
			FMKey=FMKey.lower().strip()
			if FMKey!="":
				if not searchFileLocaly(FMKey+"_pfamscan.json"):
					url = searchForFunctionalMotifs(FMKey)
					summary = summarizeFMSearch(FMKey)[0]
					summaryPopUp("Functional Motifs for keyword: "+FMKey,summary,showButt=True)
				else:
					summary = summarizeFMSearch(FMKey)[0]
					url = summarizeFMSearch(FMKey)[1]
					summaryPopUp("Functional Motifs for keyword: "+FMKey,summary,showButt=True)


def run():
	global pdbIDField
	global FMField
	global errorHandlingLabel1
	global structM
	global funcM
	global decision
	global decisionDB
	global pfamD
	global cddD
	global choiceCDD
	global choicePFAM
	global pfamEvalField
	global cddEvalField
	global lbl2
	global lbl3
	global pfamDEval
	global cddDEval
	global choiceSM
	global choiceFM
	global root
	global url
	global loadedPDBs
	global photo

	url=""
	decision=1
	decisionDB=1
	count=0
	loadedPDBs=cmd.get_object_list('solvent')
	# create root window
	root = tk.Tk()
	photo = tk.PhotoImage(file = os.getcwd()+fileSeparator+"PyProtif"+fileSeparator+"PyProtifIcon.png")
	root.wm_iconphoto(False, photo)
	# root window title and dimension
	root.title("PyProtif")
	# Set geometry(widthxheight)
	root.geometry('600x400')
	# adding menu bar in root window
	menu = Menu(root)
	root.grid_rowconfigure(0, weight=1)
	root.grid_columnconfigure(0, weight=3)
	root.grid_rowconfigure(1, weight=1)
	root.grid_columnconfigure(1, weight=3)
	root.grid_rowconfigure(2, weight=1)
	root.grid_columnconfigure(2, weight=3)
	root.grid_rowconfigure(3, weight=1)
	root.grid_rowconfigure(4, weight=1)
	root.grid_rowconfigure(5, weight=1)
	root.grid_rowconfigure(6, weight=1)
	root.grid_rowconfigure(7, weight=1)
	root.grid_rowconfigure(8, weight=1)
	root.grid_rowconfigure(9, weight=1)

	root.config(menu=menu)
	# adding a label to the root window
	lbl = Label(root, text = "Enter Your PDB ID")
	lbl.grid(column =1, row =0)
	# adding Entry Field
	pdbIDField = Entry(root)
	pdbIDField.grid(column =1, row =1)

	# lbl4 = Label(root, text = "Search for Related uniprotID on pfam database\n using functional motifs")
	# lbl4.grid(column =1, row =7)

	# FMField = Entry(root)
	# FMField.grid(column =1, row =8)

	structM = tk.IntVar(value=1)
	funcM = tk.IntVar()
	pfamD = tk.IntVar()
	cddD = tk.IntVar()
	pfamDEval = tk.StringVar()
	cddDEval = tk.StringVar()
	choiceSM = tk.Checkbutton(root, text='Structural Motifs',variable=structM, onvalue=1, offvalue=0, command=get_selection)
	choiceSM.grid(column =1, row =2)
	choiceFM = tk.Checkbutton(root, text='Functional Motifs',variable=funcM, onvalue=1, offvalue=0, command=get_selection)
	choiceFM.grid(column =1, row =3)

	choicePFAM = tk.Checkbutton(root, text='Pfam', state= "disabled",variable=pfamD, onvalue=1, offvalue=0, command=get_DBselection)
	choicePFAM.select()
	choicePFAM.grid(column =0, row =4)
	choiceCDD = tk.Checkbutton(root, text='cdd', state= "disabled",variable=cddD, onvalue=1, offvalue=0, command=get_DBselection)
	choiceCDD.grid(column =0, row =5)
	lbl2 = Label(root, text = "pfam Eval:")
	lbl2.grid(column =1 , row =4)
	pfamEvalField = Entry(root,textvariable=pfamDEval,state= "disabled")
	pfamEvalField.grid(column =2, row =4)
	pfamDEval.set("1.0")
	lbl3 = Label(root, text = "cdd Eval:")
	lbl3.grid(column =1, row =5)
	cddEvalField = Entry(root,textvariable=cddDEval,state= "disabled")
	cddEvalField.grid(column =2, row =5)
	cddDEval.set("1.0")
	# button widget with blue color text insid3
	btn = Button(root, text = "Run" ,fg = "blue", command=cascade)
	# btn2 = Button(root, text = "Go" ,fg = "blue", command=runFMSearch)
	# Set Button Grid
	btn.grid(column=1, row=6)
	# btn2.grid(column=1, row=9)
	if not isConnected() and count==0:
		messagebox.showwarning("Warning", "You are not connected to the internet restrict your usage to local files\n")
		count+=1
	root.mainloop()
########################################################################################################
def logoGen():
	global loadedPDBs
	loadedPDBs=[]
	cmd.fetch("1MH1",type='pdb')
	cmd.select("myprot","1MH1")
	cmd.color("cyan", "myprot")
	cmd.bg_color("white")
	time.sleep(1)
	cmd.png("PyProtif"+fileSeparator+'PyProtifIcon.png')
	cmd.delete("1MH1")
	cmd.delete("myprot")
	cmd.bg_color("black")
	os.remove("1MH1.pdb")
def loadIcon():
	try:
		f=open("PyProtif"+fileSeparator+'PyProtifIcon.png')
		f.close()
	except:
		logoGen()


#########################################################################################################

def initialize():
	global fileSeparator
	global runCounts
	global driver
	global system
	#check which system we are using and update fileSeparator variable accordingly
	system=platform.platform(terse=True)
	if "windows" in system.lower():
		fileSeparator = "\\"
	else:
		if "darwin" in system.lower():
			#:
			fileSeparator = ":"
		else:
			fileSeparator = "/"
	#initialize a PyProtif directory if it does not exist
	createDirectory("PyProtif")
	try:
		loadIcon()
		today = date.today()
		cmd.log_open("PyProtif"+fileSeparator+"log_"+str(today)+".txt","a")
		#cmd.feedback("disable","all","everything")
		#cmd.feedback("enable","all","errors")
		run()
		cmd.log_close()
	except:
		messagebox.showerror("Error launching PyProtif","Another instance of the app is running please close it in order to run another instance of PyProtif")
	
	
	

def __init__(self):
	global runCounts
	runCounts=0
	self.menuBar.addmenuitem('Plugin', 'command',
                             "PyProtif",
                             label='PyProtif',
                             command=lambda s=self : initialize())