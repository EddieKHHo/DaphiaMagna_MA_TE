
##############################------------------------------Run TEFLoN for each genotype
##############################	Shown is code for only the FASC assembly as reference, but this was repeated for all nine reference assemblies


####################--------------------TEFLoN Step 1

#####declare GT
GT="FA"
SC=$GT"SC"

#####Declare sample ID that will be processed (using A_ID)
for ID in FASC FBSC FCSC GASC GBSC GCSC IASC IBSC ICSC; do
python /mnt/lfs2/schaack/eddieho/TEFLoN/teflon.v0.4.py \
-wd /mnt/lfs2/schaack/eddieho/results/trElementAnalysis/teflon/results_SC/$GT/ \
-d /mnt/lfs2/schaack/eddieho/results/trElementAnalysis/teflon/prep_SC/$SC/$SC.prep_TF/ \
-s /mnt/lfs2/schaack/eddieho/results/trElementAnalysis/teflon/prep_SC/fileLoc/$GT.txt \
-i $ID \
-eb /opt/modules/biology/bwa/0.7.17/bin/bwa \
-es /opt/modules/biology/samtools/1.9/bin/samtools \
-l1 class \
-l2 class \
-q 30 \
-t 10 \
-sd 150 \
> /mnt/lfs2/schaack/eddieho/results/trElementAnalysis/teflon/results_SC/LOG/$ID.log \
2> /mnt/lfs2/schaack/eddieho/results/trElementAnalysis/teflon/results_SC/LOG/$ID.log
done

####################--------------------TEFLoN Step 2
GT="FA"
SC=$GT"SC"
python /mnt/lfs2/schaack/eddieho/TEFLoN/teflon_collapse.py \
-wd /mnt/lfs2/schaack/eddieho/results/trElementAnalysis/teflon/results_SC/$GT/ \
-d /mnt/lfs2/schaack/eddieho/results/trElementAnalysis/teflon/prep_SC/$SC/$SC.prep_TF/ \
-s /mnt/lfs2/schaack/eddieho/results/trElementAnalysis/teflon/prep_SC/fileLoc/$GT.txt \
-es /opt/modules/biology/samtools/1.9/bin/samtools \
-n1 10 -n2 10 -q 30 -t 10

####################--------------------TEFLoN Step 3

#####declare GT
GT="FA"
SC=$GT"SC"

#####Declare sample ID that will be processed (using A_ID)
for ID in FASC FBSC FCSC GASC GBSC GCSC IASC IBSC ICSC; do
python /mnt/lfs2/schaack/eddieho/TEFLoN/teflon_count.py \
-wd /mnt/lfs2/schaack/eddieho/results/trElementAnalysis/teflon/results_SC/$GT/ \
-d /mnt/lfs2/schaack/eddieho/results/trElementAnalysis/teflon/prep_SC/$SC/$SC.prep_TF/ \
-s /mnt/lfs2/schaack/eddieho/results/trElementAnalysis/teflon/prep_SC/fileLoc/$GT.EC.txt \
-i $ID \
-eb /opt/modules/biology/bwa/0.7.17/bin/bwa \
-es /opt/modules/biology/samtools/1.9/bin/samtools \
-l2 class -q 30 -t 10 \
> /mnt/lfs2/schaack/eddieho/results/trElementAnalysis/teflon/results_SC/LOG/$ID.STEP3.log \
2> /mnt/lfs2/schaack/eddieho/results/trElementAnalysis/teflon/results_SC/LOG/$ID.STEP3.log
done

####################--------------------TEFLoN Step 4
GT="FA"
SC=$GT"SC"	
python /mnt/lfs2/schaack/eddieho/TEFLoN/teflon_genotype.py \
-wd /mnt/lfs2/schaack/eddieho/results/trElementAnalysis/teflon/results_SC/$GT/ \
-d /mnt/lfs2/schaack/eddieho/results/trElementAnalysis/teflon/prep_SC/$SC/$SC.prep_TF/ \
-s /mnt/lfs2/schaack/eddieho/results/trElementAnalysis/teflon/prep_SC/fileLoc/$GT.EC.txt \
-dt pooled




##############################------------------------------Analyze output of TEFLoN to obtain TE insertion site polymorphisms

##########----------Import packages
from __future__ import division
import sys
import pandas as pd
import numpy as np
######Add to sys path and import denovoFunctions
sys.path.append('/mnt/lfs2/schaack/eddieho/scripts/python')
from DmagnaDict import *


####################----------Create a vcf-like file of read counts for each GT
####################	all genotypes.txt file of the same GT have the SAME teID column (confirmed above)
####################	so can use teID as the unique identifier for a TE region
####################	also adds contig length and pN

#####-----define column names of genotypes.txt files
colNames = ['CONTIG','BP5','BP3','sID','cID','STRAND','refID','sc5','sc3','PRES','ABS','AMB','FREQ','teID']

#####-----define list of GT
scList = ['FASC','FBSC','FCSC','GASC','GBSC','GCSC','IASC','IBSC','ICSC']
gtList = ['FA','FB','FC','GA','GB','GC','IA','IB','IC']

#####-----create an vcf-like file for read counts for each GT
#####	order of columns is always FASC, FBSC, ... IBSC, ICSC regardless of focal GT
for GT in gtList:
	print('\nProcessing {}'.format(GT))
	SC = ancLines[GT]
	#####-----read FASC data and add TOTAL column
	DIR='/mnt/lfs2/schaack/eddieho/results/trElementAnalysis/teflon/results_SC/'+GT+'/genotypes/'
	fascDATA = pd.read_table(DIR+'FASC'+'.genotypes.txt', header=None, names=colNames)
	fascDATA['TOTAL'] = fascDATA['PRES']+fascDATA['ABS']
	#####-----create FASC column as PRES:ABS:TOTAL
	fascDATA['FASC'] = fascDATA[['PRES','ABS','TOTAL']].apply(lambda x: ':'.join(x.astype(str).values), axis=1)
	#####-----dfData only keeps teID and FASC data
	dfData = fascDATA[['CONTIG','teID','FASC']]
	#####-----add MA anc EC data to dfData
	for LINE in ['FBSC','FCSC','GASC','GBSC','GCSC','IASC','IBSC','ICSC']:
		#####-----read LINE data
		scDATA = pd.read_table(DIR+LINE+'.genotypes.txt', header=None, names=colNames)
		scDATA['TOTAL'] = scDATA['PRES']+scDATA['ABS']
		#####-----create SC column as PRES:ABS:TOTAL
		scDATA[LINE] = scDATA[['PRES','ABS','TOTAL']].apply(lambda x: ':'.join(x.astype(str).values), axis=1)
		#####-----only keep teID and LINE
		temp0 = scDATA[['teID',LINE]]
		#####-----merge to dfData
		dfData = pd.merge(dfData, temp0, how = 'inner', on=['teID'], left_index=True)
	#####-----add ref file stats of focal SC to dfData
	refStats = pd.read_csv('/mnt/lfs2/schaack/eddieho/results/NewReference/final/111/bioPython/'+SC+'.counts.csv')
	dfData = pd.merge(dfData, refStats[['CONTIG','LEN','pN']], on=['CONTIG'], left_index=True)
	#####-----rearrange columns
	newCOL = ['teID','CONTIG','LEN','pN']+scList
	dfData = dfData[newCOL].reset_index(drop=True)
	print('# all TE, TE on contigs >= 5kb: {}, {}'.format(len(dfData),len(dfData[dfData['LEN']>=5000])))
	#####-----output
	OUTFILE='/mnt/lfs2/schaack/eddieho/results/trElementAnalysis/teflon/results_SC/analysis/'+GT+'.readInfo.csv'
	dfData.to_csv(path_or_buf = OUTFILE, index = False)
	



####################----------FILTER 1: CONTIG LEN and TOTAL reads 
####################	Only keep TE on contigs that are longer than 5kb
####################	Only keep TE were TOTAL reads (PRES+ABS) is >= minTotal for all SC and MA lines

#####-----function to filter by CONTIG length and TOTAL reads
def filterByTotal(GT, minTotal):
	#####-----define all lines
	ALL = ['FASC','FBSC','FCSC','GASC','GBSC','GCSC','IASC','IBSC','ICSC']
	#####-----read data file
	INFILE='/mnt/lfs2/schaack/eddieho/results/trElementAnalysis/teflon/results_SC/analysis/'+GT+'.readInfo.csv'
	DATA = pd.read_csv(INFILE)
	#####-----only keep TE on contigs >= 5kb
	DATA = DATA[DATA['LEN']>=5000].reset_index(drop=True)
	#####-----loop through each TE and determine if it passes
	newDATA = pd.DataFrame()
	listPass = []
	for i in range(len(DATA)):
		TE = DATA[ALL].loc[i].values
		nPassTotal = 0
		for j in range(len(TE)):
			#####get num of pres, abs and total reads
			pres,abs,total = map(int,TE[j].split(':'))
			#####if total pass, then increment
			if total >= minTotal:
				nPassTotal += 1
		#####determine is TE passes
		if nPassTotal == len(ALL):
			newDATA = newDATA.append(DATA.loc[i])
	#####-----output
	print('# TE before and after: {}, {}'.format(len(DATA),len(newDATA)))
	newDATA = newDATA.reset_index(drop=True)
	newDATA = newDATA[['teID','CONTIG','LEN','pN']+ALL]
	return(newDATA)

#####-----define list of GT
gtList = ['FA','FB','FC','GA','GB','GC','IA','IB','IC']

#####-----set parameter
minTotal = 20

#####-----loop through each GT
for GT in gtList:
	print('Processing {}, minTotal = {}'.format(GT,minTotal))
	filterDATA = filterByTotal(GT, minTotal)
	OUTFILE='/mnt/lfs2/schaack/eddieho/results/trElementAnalysis/teflon/results_SC/analysis/filterTotal/'+GT+'.t'+str(minTotal)+'.csv'
	filterDATA.to_csv(path_or_buf = OUTFILE, index = False)




####################--------------------FILTER 2: based on minF and maxN

##########----------Import packages
from __future__ import division
import sys
import pandas as pd
import numpy as np
######Add to sys path and import denovoFunctions
sys.path.append('/mnt/lfs2/schaack/eddieho/scripts/python')
from DmagnaDict import *

#####-----function to return state given num of PRES, ABS reads and FREQ of pres reads
#####	0 = abs/abs, 1 = pres/abs, 2 = pres,pres, -1 = cannot determine with confidence
#####	het require PRES>0, ABS>0 and FREQ of pres >=minF
def getState(PRES,ABS,FREQ,minF):
	if PRES==0 and ABS>0:
		state = 0
	elif PRES>0 and ABS==0:
		state = 2
	elif PRES>0 and ABS>0 and FREQ >= minF:
		state = 1
	else:
		state = -1
	return(state)


#####-----define parameters
minTotal = 20
maxN = 3

#####-----define list of GT
gtList = ['FA','FB','FC','GA','GB','GC','IA','IB','IC']

for minP in [5]:
	for GT in gtList:
		print('Processing {}, minP: {}'.format(GT, minP))
		#####-----define all lines
		ALL = ['FASC','FBSC','FCSC','GASC','GBSC','GCSC','IASC','IBSC','ICSC']
		nSCMA = 1+len(maLines[GT])
		#####-----read data file
		INFILE='/mnt/lfs2/schaack/eddieho/results/trElementAnalysis/teflon/results_SC/analysis/filterTotal/'+GT+'.t'+str(minTotal)+'.csv'
		DATA = pd.read_csv(INFILE)
		#####-----loop through each TE and get info
		dfNewData = pd.DataFrame(columns=['teID','nPres','nAbs','nTotal','n0','n1','n2','nN','n12','FA','FB','FC','GA','GB','GC','IA','IB','IC'])
		for i in range(len(DATA)):
			ID = DATA['teID'].iloc[i]
			TE = DATA[ALL].loc[i].values
			nPres, nAbs, nTotal, n0, n1, n2, nN, n12, listStates = 0,0,0,0,0,0,0,0,[]
			for j in range(len(TE)):
				#####get num of pres, abs and total reads
				pres,abs,total = map(int,TE[j].split(':'))
				nPres += pres
				nAbs += abs
				nTotal += total
				#####get freq if total >0, else -1
				freq = pres/total if total>0 else -1
				#####get states of lines
				minF =  5/100
				scState = getState(pres,abs,freq,minF)
				if scState == 0:
					listStates.append(0)
					n0 += 1
				elif scState == 1:
					listStates.append(1)
					n1 += 1
				elif scState ==2:
					listStates.append(2)
					n2 += 1
				else:
					listStates.append(-1)
					nN += 1
			#####add data to list
			n12 = n1+n2
			newData = [ID,nPres,nAbs,nTotal,n0,n1,n2,nN,n12]+listStates
			dfNewData.loc[len(dfNewData)] = newData
		#####-----check that total SC counts is equal to 9
		print(pd.unique(dfNewData[['n0','n1','n2','nN']].sum(axis=1)))
		#####-----add stats to DATA
		newDATA = pd.merge(DATA, dfNewData, how = 'inner', on = ['teID'], left_index =True)
		#####-----FILTER by maxN
		newDATA_maxN = newDATA[newDATA['nN'] <=maxN]
		#####-----count number of line with state 0 in each pop
		popF,popG,popI = ['FA','FB','FC'],['GA','GB','GC'],['IA','IB','IC']
		newDATA_maxN['F0'] = newDATA_maxN[popF].eq(0).sum(1)
		newDATA_maxN['G0'] = newDATA_maxN[popG].eq(0).sum(1)
		newDATA_maxN['I0'] = newDATA_maxN[popI].eq(0).sum(1)
		#####-----add TEFLON data from SC genotype file
		colNames = ['CONTIG','BP5','BP3','sID','cID','STRAND','refID','sc5','sc3','PRES','ABS','AMB','FREQ','teID']
		INFILE='/mnt/lfs2/schaack/eddieho/results/trElementAnalysis/teflon/results_SC/'+GT+'/genotypes/'+GT+'SC.genotypes.txt'
		scDATA = pd.read_table(INFILE, header=None, names=colNames)
		#####add TE CLASS column
		scDATA['FAMILY'] = scDATA.sID.copy()
		scDATA['CLASS'] = scDATA.sID.str.split('_').str[0]
		#####merge
		mergeDATA = pd.merge(newDATA_maxN,scDATA[['teID','CLASS','FAMILY','refID','BP5','BP3','sc5','sc3']], on='teID')
		#####-----output
		OUTFILE='/mnt/lfs2/schaack/eddieho/results/trElementAnalysis/teflon/results_SC/analysis/filterTotal/gtInfo/'+GT+'.t'+str(minTotal)+'.'+str(minP)+'.'+str(maxN)+'.gtInfo.csv'
		mergeDATA.to_csv(path_or_buf = OUTFILE, index = False)
		


####################--------------------Get TIPs info across all TE families
####################	additional filter of pN <= 0.5

##########----------Import packages
from __future__ import division
import sys
import pandas as pd
import numpy as np
######Add to sys path and import denovoFunctions
sys.path.append('/mnt/lfs2/schaack/eddieho/scripts/python')
from DmagnaDict import *


##########----------Get TE that are unique to one SC line (absent in all but one focal SC line)
def getUniqTE_SC(refGT, DATA):
	#####define list of SC
	scList = ['FA','FB','FC','GA','GB','GC','IA','IB','IC']
	#####defein df
	dfUniqTE = pd.DataFrame()
	#####loop through each SC and get uniq TE
	for SC in scList:
		#####define columns to keep
		uniqCol = ['teID', 'CONTIG', 'LEN', 'pN','CLASS','FAMILY','F0','G0','I0',SC]
		#####subset uniqe TE for SC
		#temp0 = DATA[(DATA['nN']==0) & (DATA['n12']==1) & ((DATA[SC]==1 | (DATA[SC]==2)))][uniqCol]
		temp0 = DATA[(DATA['n12']==1) & ((DATA[SC]==1) | (DATA[SC]==2))][uniqCol]
		#####rename state of focal SC (col 'SC') to scState
		temp0 = temp0.rename(columns={SC:'scState'})
		#####make new SC to indicate focal SC
		temp0['SC'] = SC
		temp0['refGT'] = refGT
		#####reorder columns
		temp0 = temp0[['refGT','SC','teID', 'CONTIG', 'LEN', 'pN','CLASS','FAMILY','F0','G0','I0','scState']]
		#####append to df
		dfUniqTE = dfUniqTE.append(temp0).reset_index(drop=True)
	#####return
	return(dfUniqTE)

##########----------Get TE that are unique to one POP (absent in all but one focal POP)
def getUniqTE_POP(refGT, DATA):
	#####define list of SC
	scList = ['FA','FB','FC','GA','GB','GC','IA','IB','IC']
	#####define pops that need to have absence of TEs for each focal pop
	#####	e.g. for TE to be specigfic for F, I need all SC in G and I to lack TE (so G0 = 3and I0 = 3)
	popDict = {'F':['G0','I0'],'G':['F0','I0'],'I':['F0','G0']}
	#####define columns to keep
	uniqCol = ['n12','teID', 'CONTIG', 'LEN', 'pN','CLASS','FAMILY','F0','G0','I0']
	#####defein df
	dfUniqTE = pd.DataFrame()
	#####loop through each POP and get uniq TE
	for POP in ['F','G','I']:
		for nPres in [2,3]:
			#####define populations that need to have no TE (state 0)
			pop0 = popDict[POP]
			#####subset uniqe TE for POP
			temp0 = DATA[(DATA['n12']==nPres) & (DATA[pop0[0]]==3) & (DATA[pop0[1]]==3)][uniqCol]
			#####make new SC to indicate focal SC
			temp0['POP'] = POP
			temp0['refGT'] = refGT
			#####reorder columns
			temp0 = temp0[['refGT','POP','n12','teID', 'CONTIG', 'LEN', 'pN','CLASS','FAMILY','F0','G0','I0']]
			#####append to df
			dfUniqTE = dfUniqTE.append(temp0).reset_index(drop=True)
	#####return
	return(dfUniqTE)
	

#####-----define parameters
gtList = ['FA','FB','FC','GA','GB','GC','IA','IB','IC']

minTotal = 20
minP = 5
maxN = 0
maxPN = 0.5

#####-----create df for smy data
smyCol = [
'minTotal','minP','maxN','GT','nTE','nPoly',
'P0','P1','P2','P3','P4','P5','P6','P7','P8','P9',
'uFA','uFB','uFC','uGA','uGB','uGC','uIA','uIB','uIC',
'uF_2','uG_2','uI_2','uF_3','uG_3','uI_3','uF','uG','uI']
dfSmy = pd.DataFrame(columns=smyCol)

allUniqSC, allUniqPOP = pd.DataFrame(), pd.DataFrame()

#####-----loop through at GT
for GT in gtList:
	#####read file given GT (GT represents the ref used)
	INFILE='/mnt/lfs2/schaack/eddieho/results/trElementAnalysis/teflon/results_SC/analysis/filterTotal/gtInfo/'+GT+'.t'+str(minTotal)+'.'+str(minP)+'.'+str(maxN)+'.gtInfo.csv'
	focDATA = pd.read_csv(INFILE)
	
	#####FILTER for pN<=0.5
	focDATA = focDATA[focDATA.pN<=maxPN]
	#####FILTER for good FAMILIES (not include unclassified and all SINEs)
	goodFAM = ['DNA_Academ-1', 'DNA_CMC-EnSpm', 'DNA_MULE-MuDR', 'DNA_Merlin', 'DNA_P', 'DNA_P-Fungi', 'DNA_PIF-Harbinger', 'DNA_PIF-ISL2EU', 'DNA_Sola-2', 'DNA_TcMar-Fot1', 'DNA_TcMar-Tc1', 'DNA_Zisupton', 'DNA_hAT-Ac', 'DNA_hAT-Tip100', 'DNA_hAT-hATm', 'LINE_I', 'LINE_I-Jockey', 'LINE_L1-Tx1', 'LINE_Penelope', 'LINE_R1', 'LINE_R2-NeSL', 'LTR_Copia', 'LTR_DIRS', 'LTR_Gypsy', 'LTR_Ngaro', 'LTR_Pao', 'RC_Helitron']
	focDATA = focDATA[focDATA['FAMILY'].isin(goodFAM)]
	
	#####get total num TE and num polymorphic
	nTE,nPoly = len(focDATA),len(focDATA[(focDATA['n12']>0) & (focDATA['n12']<9)])
	#####count TE with x num of lines possessing the TE (range from 0 to 9)
	#####	remember that x=0, x=9 are actually monomorphic for TE presence
	listN12 = [len(focDATA[(focDATA['n12']==x)]) for x in range(10)]
	
	#####get TE that are unique to an SC
	dfUniqSC = getUniqTE_SC(GT, focDATA)
	listUniqSC = [len(dfUniqSC[(dfUniqSC['SC']==x)]) for x in ['FA','FB','FC','GA','GB','GC','IA','IB','IC']]
	
	#####get TE that are unique to a POP
	#####	can be unique with n12=2 or n12=3, I sepearate their counts and also total them
	dfUniqPOP = getUniqTE_POP(GT, focDATA)
	listUniqPOP_2 = [len(dfUniqPOP[(dfUniqPOP['n12']==2)&(dfUniqPOP['POP']==x)]) for x in ['F','G','I']]
	listUniqPOP_3 = [len(dfUniqPOP[(dfUniqPOP['n12']==3)&(dfUniqPOP['POP']==x)]) for x in ['F','G','I']]
	listUniqPOP_Both = [x + y for x, y in zip(listUniqPOP_2, listUniqPOP_3)]
	
	#####output uniq files
	OUTFILE='/mnt/lfs2/schaack/eddieho/results/trElementAnalysis/teflon/results_SC/analysis/filterTotal/tipInfo/'+GT+'.t'+str(minTotal)+'.'+str(minP)+'.'+str(maxN)+'.tip.uniqSC.csv'
	dfUniqSC.to_csv(path_or_buf = OUTFILE, index = False)
	OUTFILE='/mnt/lfs2/schaack/eddieho/results/trElementAnalysis/teflon/results_SC/analysis/filterTotal/tipInfo/'+GT+'.t'+str(minTotal)+'.'+str(minP)+'.'+str(maxN)+'.tip.uniqPOP.csv'
	dfUniqPOP.to_csv(path_or_buf = OUTFILE, index = False)
	
	#####append uniq files
	allUniqSC = allUniqSC.append(dfUniqSC)
	allUniqPOP = allUniqPOP.append(dfUniqPOP)
	
	#####create summary data vector
	smyData = [minTotal,minP,maxN,GT,nTE,nPoly]+listN12+listUniqSC+listUniqPOP_2+listUniqPOP_3+listUniqPOP_Both
	#####add to dfSmy
	dfSmy.loc[len(dfSmy)] = smyData

#####-----output 
OUTFILE='/mnt/lfs2/schaack/eddieho/results/trElementAnalysis/teflon/results_SC/analysis/filterTotal/tipInfo/'+'ALL.t'+str(minTotal)+'.'+str(minP)+'.'+str(maxN)+'.tip.smy.csv'
dfSmy.to_csv(path_or_buf = OUTFILE, index = False)

OUTFILE='/mnt/lfs2/schaack/eddieho/results/trElementAnalysis/teflon/results_SC/analysis/filterTotal/tipInfo/'+'ALL.t'+str(minTotal)+'.'+str(minP)+'.'+str(maxN)+'.tip.uniqSC.csv'
allUniqSC.to_csv(path_or_buf = OUTFILE, index = False)

OUTFILE='/mnt/lfs2/schaack/eddieho/results/trElementAnalysis/teflon/results_SC/analysis/filterTotal/tipInfo/'+'ALL.t'+str(minTotal)+'.'+str(minP)+'.'+str(maxN)+'.tip.uniqPOP.csv'
allUniqPOP.to_csv(path_or_buf = OUTFILE, index = False)





####################--------------------Get TIPs info for a specific TE family
####################	additional filter of pN <= 0.5

##########----------Import packages
from __future__ import division
import sys
import pandas as pd
import numpy as np
######Add to sys path and import denovoFunctions
sys.path.append('/mnt/lfs2/schaack/eddieho/scripts/python')
from DmagnaDict import *


##########----------Get TE that are unique to one SC line (absent in all but one focal SC line)
def getUniqTE_SC(refGT, DATA):
	#####define list of SC
	scList = ['FA','FB','FC','GA','GB','GC','IA','IB','IC']
	#####defein df
	dfUniqTE = pd.DataFrame()
	#####loop through each SC and get uniq TE
	for SC in scList:
		#####define columns to keep
		uniqCol = ['teID', 'CONTIG', 'LEN', 'pN','CLASS','FAMILY','F0','G0','I0',SC]
		#####subset uniqe TE for SC
		#temp0 = DATA[(DATA['nN']==0) & (DATA['n12']==1) & ((DATA[SC]==1 | (DATA[SC]==2)))][uniqCol]
		temp0 = DATA[(DATA['n12']==1) & ((DATA[SC]==1) | (DATA[SC]==2))][uniqCol]
		#####rename state of focal SC (col 'SC') to scState
		temp0 = temp0.rename(columns={SC:'scState'})
		#####make new SC to indicate focal SC
		temp0['SC'] = SC
		temp0['refGT'] = refGT
		#####reorder columns
		temp0 = temp0[['refGT','SC','teID', 'CONTIG', 'LEN', 'pN','CLASS','FAMILY','F0','G0','I0','scState']]
		#####append to df
		dfUniqTE = dfUniqTE.append(temp0).reset_index(drop=True)
	#####return
	return(dfUniqTE)

##########----------Get TE that are unique to one POP (absent in all but one focal POP)
def getUniqTE_POP(refGT, DATA):
	#####define list of SC
	scList = ['FA','FB','FC','GA','GB','GC','IA','IB','IC']
	#####define pops that need to have absence of TEs for each focal pop
	#####	e.g. for TE to be specigfic for F, I need all SC in G and I to lack TE (so G0 = 3and I0 = 3)
	popDict = {'F':['G0','I0'],'G':['F0','I0'],'I':['F0','G0']}
	#####define columns to keep
	uniqCol = ['n12','teID', 'CONTIG', 'LEN', 'pN','CLASS','FAMILY','F0','G0','I0']
	#####defein df
	dfUniqTE = pd.DataFrame()
	#####loop through each POP and get uniq TE
	for POP in ['F','G','I']:
		for nPres in [2,3]:
			#####define populations that need to have no TE (state 0)
			pop0 = popDict[POP]
			#####subset uniqe TE for POP
			temp0 = DATA[(DATA['n12']==nPres) & (DATA[pop0[0]]==3) & (DATA[pop0[1]]==3)][uniqCol]
			#####make new SC to indicate focal SC
			temp0['POP'] = POP
			temp0['refGT'] = refGT
			#####reorder columns
			temp0 = temp0[['refGT','POP','n12','teID', 'CONTIG', 'LEN', 'pN','CLASS','FAMILY','F0','G0','I0']]
			#####append to df
			dfUniqTE = dfUniqTE.append(temp0).reset_index(drop=True)
	#####return
	return(dfUniqTE)
	

#####-----define parameters
gtList = ['FA','FB','FC','GA','GB','GC','IA','IB','IC']

minTotal = 20
minP = 5
maxN = 0
maxPN = 0.5


#####-----define focal TE family
#####	these are families that have > 90 TE clusters in the gtInfo file when avg across the nine refGT
listFocFam = ['DNA_hAT-Ac','LINE_I','LTR_Copia','LTR_DIRS','LTR_Gypsy','LTR_Pao']

#####-----loop through at focFam and GT
for focFam in listFocFam:
	#####-----create df for smy data
	smyCol = [
	'minTotal','minP','maxN','GT','nTE','nPoly',
	'P0','P1','P2','P3','P4','P5','P6','P7','P8','P9',
	'uFA','uFB','uFC','uGA','uGB','uGC','uIA','uIB','uIC',
	'uF_2','uG_2','uI_2','uF_3','uG_3','uI_3','uF','uG','uI']
	dfSmy = pd.DataFrame(columns=smyCol)
	allUniqSC, allUniqPOP = pd.DataFrame(), pd.DataFrame()
	
	for GT in gtList:
		#####read file given GT (GT represents the ref used)
		INFILE='/mnt/lfs2/schaack/eddieho/results/trElementAnalysis/teflon/results_SC/analysis/filterTotal/gtInfo/'+GT+'.t'+str(minTotal)+'.'+str(minP)+'.'+str(maxN)+'.gtInfo.csv'
		focDATA = pd.read_csv(INFILE)
		
		#####FILTER for pN<=0.5
		focDATA = focDATA[focDATA.pN<=maxPN]
		
		#####FILTER for focal family
		focDATA = focDATA[focDATA.FAMILY==focFam]
		
		#####get total num TE and num polymorphic
		nTE,nPoly = len(focDATA),len(focDATA[(focDATA['n12']>0) & (focDATA['n12']<9)])
		#####count TE with x num of lines possessing the TE (range from 0 to 9)
		#####	remember that x=0, x=9 are actually monomorphic for TE presence
		listN12 = [len(focDATA[(focDATA['n12']==x)]) for x in range(10)]
		
		#####get TE that are unique to an SC
		dfUniqSC = getUniqTE_SC(GT, focDATA)
		listUniqSC = [len(dfUniqSC[(dfUniqSC['SC']==x)]) for x in ['FA','FB','FC','GA','GB','GC','IA','IB','IC']]
		
		#####get TE that are unique to a POP
		#####	can be unique with n12=2 or n12=3, I sepearate their counts and also total them
		dfUniqPOP = getUniqTE_POP(GT, focDATA)
		listUniqPOP_2 = [len(dfUniqPOP[(dfUniqPOP['n12']==2)&(dfUniqPOP['POP']==x)]) for x in ['F','G','I']]
		listUniqPOP_3 = [len(dfUniqPOP[(dfUniqPOP['n12']==3)&(dfUniqPOP['POP']==x)]) for x in ['F','G','I']]
		listUniqPOP_Both = [x + y for x, y in zip(listUniqPOP_2, listUniqPOP_3)]
		
		#####append uniq files
		allUniqSC = allUniqSC.append(dfUniqSC)
		allUniqPOP = allUniqPOP.append(dfUniqPOP)
		
		#####create summary data vector
		smyData = [minTotal,minP,maxN,GT,nTE,nPoly]+listN12+listUniqSC+listUniqPOP_2+listUniqPOP_3+listUniqPOP_Both
		#####add to dfSmy
		dfSmy.loc[len(dfSmy)] = smyData
	
	#####-----output 
	OUTFILE='/mnt/lfs2/schaack/eddieho/results/trElementAnalysis/teflon/results_SC/analysis/filterTotal/tipInfo/family/'+focFam+'/'+'ALL.t'+str(minTotal)+'.'+str(minP)+'.'+str(maxN)+'.tip.smy.csv'
	dfSmy.to_csv(path_or_buf = OUTFILE, index = False)
	OUTFILE='/mnt/lfs2/schaack/eddieho/results/trElementAnalysis/teflon/results_SC/analysis/filterTotal/tipInfo/family/'+focFam+'/'+'ALL.t'+str(minTotal)+'.'+str(minP)+'.'+str(maxN)+'.tip.uniqSC.csv'
	allUniqSC.to_csv(path_or_buf = OUTFILE, index = False)
	OUTFILE='/mnt/lfs2/schaack/eddieho/results/trElementAnalysis/teflon/results_SC/analysis/filterTotal/tipInfo/family/'+focFam+'/'+'ALL.t'+str(minTotal)+'.'+str(minP)+'.'+str(maxN)+'.tip.uniqPOP.csv'
	allUniqPOP.to_csv(path_or_buf = OUTFILE, index = False)
	