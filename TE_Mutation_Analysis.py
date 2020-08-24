
##############################------------------------------Run TEFLoN for each genotype
##############################	Shown is code for only the FASC genotype, but this was repeated for all nine genotypes

####################--------------------Prepare TE library for TEFLON
####################	Remove: Satellite, Simple_repeat, snRNA, Low_complexity, Unknown
####################	Sequence name change format to TEFLON format: ID#family/order

##########Load packages
from __future__ import division
import pandas as pd
import numpy as np

##########----------Create new 
#####designate labels at G2 that are NOT TE
BADG2 = ['Satellite', 'Simple_repeat', 'snRNA', 'Low_complexity', 'Unknown']

#####define in and out file
INFILE='/mnt/lfs2/schaack/eddieho/results/trElementAnalysis/mapTELib/teLib/ALL.repeats.nonRedundant.long.fasta'
OUTFILE='/mnt/lfs2/schaack/eddieho/results/trElementAnalysis/teflon/teLib/ALL.repeats.nonRedundant.long.onlyTE.fasta'

#####initialize nLine and tePASS as 0
nLine = 0
tePASS = 0
ID = 1
with open('%s' %INFILE, 'r') as file, open('%s' %OUTFILE, 'w') as out:
	#####Loop through all lines
	for line in file:
		#####if line is sequence name then check if it is TE
		if line.startswith('>'):
			#####split line by white space
			field = line.split(':')
			#####G2 contains family/order info
			G2 = field[1]
			#####if G2 not exist in BADG2, then output
			#####	replace '/' with '_' because TEFLON will have erros with '/'
			#####	replace 'SINE?' with just 'SINE'
			if G2 not in BADG2:
				G2New = G2.replace('/','_')
				if G2=='SINE?':
					G2New = 'SINE'
				name = '>DM'+str(ID)+'#'+G2New+'\n'
				out.write(name)
				ID += 1
				tePASS = 1
		#####if line is sequence, print if tePASS = 1
		elif tePASS==1:
			out.write(line)
			tePASS = 0



####################--------------------Prepare for TEFLON
####################	This was performed for each SC line

##########----------TEFLoN teflon_prep_custom.py
module load repeatmasker
SC='FASC'
python /mnt/lfs2/schaack/eddieho/TEFLoN/teflon_prep_custom.py \
-wd /mnt/lfs2/schaack/eddieho/results/trElementAnalysis/teflon/prep_PO/$SC/ \
-e /opt/modules/biology/repeatmasker/4.0.6/bin/RepeatMasker \
-g /mnt/lfs2/schaack/eddieho/results/NewReference/final/111/$SC.fasta \
-l /mnt/lfs2/schaack/eddieho/results/trElementAnalysis/teflon/teLib/ALL.repeats.nonRedundant.long.onlyTE.fasta \
-p $SC -s 150 -t 10


##########----------BWA map SC, MA and EC lines to SC.mappingRef.fa

for ID in FASC FA1 FA5 FA7 FA10 FA12 FA13 FA14 FAEC1 FAEC2; do
#####declare SC
SC=${ID:0:2}SC

#####Declare TE LIB as reference
refFile="/mnt/lfs2/schaack/eddieho/results/trElementAnalysis/teflon/prep_PO/$SC/$SC.prep_MP/$SC.mappingRef.fa"

#####Declare prefix for output
PREFIX="/mnt/lfs2/schaack/eddieho/results/trElementAnalysis/teflon/prep_PO/$SC/bwa"

#####------Bwa alignment-----#####
readFile_R1="/mnt/lfs2/schaack/megadaph/fmacrae/cleaned_reads/$ID.R1.fastq.gz"
readFile_R2="/mnt/lfs2/schaack/megadaph/fmacrae/cleaned_reads/$ID.R2.fastq.gz"
outputFile=$PREFIX"/sam/$ID.paired.sam"
bwa mem -t 8 -Y \
$refFile \
$readFile_R1 \
$readFile_R2 \
> $outputFile

#####-----Picard: Convert .sam to .bam and sort-----#####
input=$PREFIX"/sam/$ID.paired.sam"
output=$PREFIX"/bam/$ID.paired.bam"
picard SortSam INPUT="$input" OUTPUT="$output" SORT_ORDER=coordinate

#####-----Picard: Mark duplicates and also merge with pair reads bam-----#####
input_paired=$PREFIX"/bam/$ID.paired.bam"
output=$PREFIX"/bam/markDup/$ID.bam"
metrics=$PREFIX"/bam/markDup/$ID.metrics.txt"
picard MarkDuplicates MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=100 INPUT="$input_paired" OUTPUT="$output" METRICS_FILE="$metrics"

#####-----Samtools: index bam file-----#####
input=$PREFIX"/bam/markDup/$ID.bam"
samtools index $input
done


####################--------------------TEFLoN Step 1

#####declare GT
GT="FA"
SC=$GT"SC"

#####Declare sample ID that will be processed (using A_ID)
for ID in FASC FA1 FA5 FA7 FA10 FA12 FA13 FA14 FAEC1 FAEC2; do
python /mnt/lfs2/schaack/eddieho/TEFLoN/teflon.v0.4.py \
-wd /mnt/lfs2/schaack/eddieho/results/trElementAnalysis/teflon/results_EC/$GT/ \
-d /mnt/lfs2/schaack/eddieho/results/trElementAnalysis/teflon/prep_PO/$SC/$SC.prep_TF/ \
-s /mnt/lfs2/schaack/eddieho/results/trElementAnalysis/teflon/prep_PO/fileLoc/$GT.EC.txt \
-i $ID \
-eb /opt/modules/biology/bwa/0.7.17/bin/bwa \
-es /opt/modules/biology/samtools/1.9/bin/samtools \
-l1 class \
-l2 class \
-q 30 \
-t 10 \
-sd 150 \
> /mnt/lfs2/schaack/eddieho/results/trElementAnalysis/teflon/results_EC/LOG/$ID.log \
2> /mnt/lfs2/schaack/eddieho/results/trElementAnalysis/teflon/results_EC/LOG/$ID.log
done

####################--------------------TEFLoN Step 2
GT="FA"
SC=$GT"SC"
python /mnt/lfs2/schaack/eddieho/TEFLoN/teflon_collapse.py \
-wd /mnt/lfs2/schaack/eddieho/results/trElementAnalysis/teflon/results_EC/$GT/ \
-d /mnt/lfs2/schaack/eddieho/results/trElementAnalysis/teflon/prep_PO/$SC/$SC.prep_TF/ \
-s /mnt/lfs2/schaack/eddieho/results/trElementAnalysis/teflon/prep_PO/fileLoc/$GT.EC.txt \
-es /opt/modules/biology/samtools/1.9/bin/samtools \
-n1 10 -n2 10 -q 30 -t 10

####################--------------------TEFLoN Step 3

#####declare GT
GT="FA"
SC=$GT"SC"

#####Declare sample ID that will be processed (using A_ID)
for ID in FASC FA1 FA5 FA7 FA10 FA12 FA13 FA14 FAEC1 FAEC2; do
python /mnt/lfs2/schaack/eddieho/TEFLoN/teflon_count.py \
-wd /mnt/lfs2/schaack/eddieho/results/trElementAnalysis/teflon/results_EC/$GT/ \
-d /mnt/lfs2/schaack/eddieho/results/trElementAnalysis/teflon/prep_PO/$SC/$SC.prep_TF/ \
-s /mnt/lfs2/schaack/eddieho/results/trElementAnalysis/teflon/prep_PO/fileLoc/$GT.EC.txt \
-i $ID \
-eb /opt/modules/biology/bwa/0.7.17/bin/bwa \
-es /opt/modules/biology/samtools/1.9/bin/samtools \
-l2 class -q 30 -t 10 \
> /mnt/lfs2/schaack/eddieho/results/trElementAnalysis/teflon/results_EC/LOG/$ID.STEP3.log \
2> /mnt/lfs2/schaack/eddieho/results/trElementAnalysis/teflon/results_EC/LOG/$ID.STEP3.log
done

####################--------------------TEFLoN Step 4
GT="FA"
SC=$GT"SC"	
python /mnt/lfs2/schaack/eddieho/TEFLoN/teflon_genotype.py \
-wd /mnt/lfs2/schaack/eddieho/results/trElementAnalysis/teflon/results_EC/$GT/ \
-d /mnt/lfs2/schaack/eddieho/results/trElementAnalysis/teflon/prep_PO/$SC/$SC.prep_TF/ \
-s /mnt/lfs2/schaack/eddieho/results/trElementAnalysis/teflon/prep_PO/fileLoc/$GT.EC.txt \
-dt pooled





##############################------------------------------Analyze output of TEFLoN to obtain MA and EC mutations
##############################	In this code, the names of different types of mutation differ from the manuscript
##############################	0 -> 1 gain = homoIns
##############################	2 -> 1 loss = homoDel
##############################	1 -> 2 gain = hetIns
##############################	1 -> 0 loss = hetDel


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
gtList = ['FA','FB','FC','GA','GB','GC','IA','IB','IC']

#####-----create an vcf-like file for read counts for each GT
for GT in gtList:
	print('\nProcessing {}'.format(GT))
	SC = ancLines[GT]
	#####-----read SC data and add TOTAL column
	DIR='/mnt/lfs2/schaack/eddieho/results/trElementAnalysis/teflon/results_EC/'+GT+'/genotypes/'
	scDATA = pd.read_table(DIR+SC+'.genotypes.txt', header=None, names=colNames)
	scDATA['TOTAL'] = scDATA['PRES']+scDATA['ABS']
	#####-----create SC column as PRES:ABS:TOTAL
	scDATA[SC] = scDATA[['PRES','ABS','TOTAL']].apply(lambda x: ':'.join(x.astype(str).values), axis=1)
	#####-----dfData only keeps teID and SC data
	dfData = scDATA[['CONTIG','teID',SC]]
	#####-----add MA anc EC data to dfData
	for LINE in maLines[GT]+extLines[GT]:
		#####-----read LINE data
		maDATA = pd.read_table(DIR+LINE+'.genotypes.txt', header=None, names=colNames)
		maDATA['TOTAL'] = maDATA['PRES']+maDATA['ABS']
		#####-----create SC column as PRES:ABS:TOTAL
		maDATA[LINE] = maDATA[['PRES','ABS','TOTAL']].apply(lambda x: ':'.join(x.astype(str).values), axis=1)
		#####-----only keep teID and LINE
		temp0 = maDATA[['teID',LINE]]
		#####-----merge to dfData
		dfData = pd.merge(dfData, temp0, how = 'inner', on=['teID'], left_index=True)
	#####-----add ref file stats to dfData
	refStats = pd.read_csv('/mnt/lfs2/schaack/eddieho/results/NewReference/final/111/bioPython/'+SC+'.counts.csv')
	dfData = pd.merge(dfData, refStats[['CONTIG','LEN','pN']], on=['CONTIG'], left_index=True)
	#####-----rearrange columns
	newCOL = ['teID','CONTIG','LEN','pN']+[SC]+maLines[GT]+extLines[GT]
	dfData = dfData[newCOL].reset_index(drop=True)
	print('# all TE, TE on contigs >= 5kb: {}, {}'.format(len(dfData),len(dfData[dfData['LEN']>=5000])))
	#####-----output
	OUTFILE='/mnt/lfs2/schaack/eddieho/results/trElementAnalysis/teflon/results_EC/analysis/'+GT+'.readInfo.csv'
	dfData.to_csv(path_or_buf = OUTFILE, index = False)




####################----------FILTER 1: CONTIG LEN and TOTAL reads 
####################	Only keep TE on contigs that are longer than 5kb
####################	Only keep TE were TOTAL reads (PRES+ABS) is >= minTotal for all SC and MA lines

#####-----function to filter by CONTIG length and TOTAL reads
def filterByTotal(GT, minTotal):
	#####-----define all lines
	ALL = [ancLines[GT]]+maLines[GT]+extLines[GT]
	#####-----read data file
	INFILE='/mnt/lfs2/schaack/eddieho/results/trElementAnalysis/teflon/results_EC/analysis/'+GT+'.readInfo.csv'
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
	OUTFILE='/mnt/lfs2/schaack/eddieho/results/trElementAnalysis/teflon/results_EC/analysis/filterTotal/'+GT+'.t'+str(minTotal)+'.csv'
	filterDATA.to_csv(path_or_buf = OUTFILE, index = False)




####################--------------------Examine GOOD heterozygous sites with EC lines
####################	loosley define heterozygotes as lines where P>0 and A>0
####################	a good het must be het in ALL lines  ( >= minHet)
####################		so although individual het is loosely defined, we have incrased confidence when ALL lines are het

##########----------Import packages
from __future__ import division
import pandas as pd
import numpy as np
import sys
######Add to sys path and import denovoFunctions
sys.path.append('/mnt/lfs2/schaack/eddieho/scripts/python')
from DmagnaDict import *


def getGoodHet(GT, minTotal, minHet):
	#####-----define all lines
	ALL = [ancLines[GT]]+maLines[GT]+extLines[GT]
	#####-----read data file
	INFILE='/mnt/lfs2/schaack/eddieho/results/trElementAnalysis/teflon/results_EC/analysis/filterTotal/'+GT+'.t'+str(minTotal)+'.csv'
	DATA = pd.read_csv(INFILE)
	#####-----loop through each TE, only save good heterozygotes
	minHet = len(ALL)
	goodHet = pd.DataFrame()
	goodHetFreq = []
	for i in range(len(DATA)):
		TE = DATA[ALL].loc[i].values
		nHet,listFreq = 0,[]
		for j in range(len(TE)):
			#####get num of pres, abs and total reads
			pres,abs,total = map(int,TE[j].split(':'))
			#####if het then increment nHet and calc freqP
			if pres>0 and abs>0:
				nHet += 1
				listFreq.append(pres/total)
		#####if all lines are het, then add to goodHet and add to goodHetFreq
		if nHet >= minHet:
			goodHet = goodHet.append(DATA.loc[i])
			goodHetFreq.extend(listFreq)
	#####-----return
	return (goodHet.reset_index(drop=True), goodHetFreq)

#####-----define list of GT
gtList = ['FA','FB','FC','GA','GB','GC','IA','IB','IC']
minTotal = 20

#####-----loop through each GT
dfHetFreqStat = pd.DataFrame(columns=['GT','p0','p5','p10','p15','p20','p25','p50','p75','p80','p85','p90','p95','p100'])
for GT in gtList:
	print('Processing {}'.format(GT))
	minHet = len([ancLines[GT]]+maLines[GT]+extLines[GT])
	goodHet, goodHetFreq = getGoodHet(GT, minTotal, minHet)
	print('Num good het = {}'.format(len(goodHet)))
	print('Mean {}'.format(np.mean(goodHetFreq)))
	dfHetFreqStat.loc[len(dfHetFreqStat)] = [GT]+np.percentile(goodHetFreq,[0,5,10,15,20,25,50,75,80,85,90,95,100]).tolist()
	OUTFILE='/mnt/lfs2/schaack/eddieho/results/trElementAnalysis/teflon/results_EC/analysis/filterTotal/goodHet/'+GT+'.t'+str(minTotal)+'.goodHet.csv'
	goodHet.to_csv(path_or_buf = OUTFILE, index = False)

#####output stat file
OUTFILE='/mnt/lfs2/schaack/eddieho/results/trElementAnalysis/teflon/results_EC/analysis/filterTotal/goodHet/ALL.t'+str(minTotal)+'.goodHet.stats.csv'
dfHetFreqStat.to_csv(path_or_buf = OUTFILE, index = False)






####################--------------------FILTER 2: based on dfHetFreqStat 
####################	using dfHetFreqStat from on TE_Analysis_PO.py (does NOT include EC lines het data)
####################		original code is above

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
def getState(PRES,ABS,FREQ,minF,maxF):
	if PRES==0 and ABS>0:
		state = 0
	elif PRES>0 and ABS==0:
		state = 2
	elif FREQ >= minF and FREQ <= maxF:
		state = 1
	else:
		state = -1
	return(state)


#####-----define parameters
minTotal = 20

#####-----read goodHet stats (require minTotal)
INFILE='/mnt/lfs2/schaack/eddieho/results/trElementAnalysis/teflon/results_EC/analysis/filterTotal/goodHet/ALL.t'+str(minTotal)+'.goodHet.stats.csv'
goodHetStats = pd.read_csv(INFILE)

#####-----define list of GT
gtList = ['FA','FB','FC','GA','GB','GC','IA','IB','IC']
for minPercent in [10]:
	#####-----define minP,maxP
	minP = 'p'+str(minPercent)
	maxP = 'p'+str(100-minPercent)
	for GT in gtList:
		print('Processing {}, minP, maxP: {}, {}'.format(GT, minP,maxP))
		#####-----define all lines
		ALL = [ancLines[GT]]+maLines[GT]+extLines[GT]
		nSCMA = 1+len(maLines[GT])
		#####-----read data file
		INFILE='/mnt/lfs2/schaack/eddieho/results/trElementAnalysis/teflon/results_EC/analysis/filterTotal/'+GT+'.t'+str(minTotal)+'.csv'
		DATA = pd.read_csv(INFILE)
		#####-----loop through each TE and get info
		minFreq = float(goodHetStats[goodHetStats['GT']==GT][minP])
		maxFreq = float(goodHetStats[goodHetStats['GT']==GT][maxP])
		dfNewData = pd.DataFrame(columns=['teID','nPres','nAbs','nTotal','SC','nMA0','nMA1','nMA2','nMAn','nEC0','nEC1','nEC2','nECn'])
		for i in range(len(DATA)):
			ID = DATA['teID'].iloc[i]
			TE = DATA[ALL].loc[i].values
			nPassTotal, nPres, nAbs, nTotal, nMA0, nMA1, nMA2, nMAn, nEC0, nEC1, nEC2, nECn = 0,0,0,0,0,0,0,0,0,0,0,0
			for j in range(len(TE)):
				#####get num of pres, abs and total reads
				pres,abs,total = map(int,TE[j].split(':'))
				nPres += pres
				nAbs += abs
				nTotal += total
				#####get freq if total >0, else -1
				freq = pres/total if total>0 else -1
				#####different steps for SC vs MA
				if j==0:
					scState = getState(pres,abs,freq,minFreq,maxFreq)
				elif j < nSCMA:
					maState = getState(pres,abs,freq,minFreq,maxFreq)
					if maState == 0:
						nMA0 += 1
					elif maState == 1:
						nMA1 += 1
					elif maState ==2:
						nMA2 += 1
					else:
						nMAn += 1
				else:
					ecState = getState(pres,abs,freq,minFreq,maxFreq)
					if ecState == 0:
						nEC0 += 1
					elif ecState == 1:
						nEC1 += 1
					elif ecState ==2:
						nEC2 += 1
					else:
						nECn += 1
			#####add data to list
			newData = [ID,nPres,nAbs,nTotal,scState,nMA0,nMA1,nMA2,nMAn,nEC0,nEC1,nEC2,nECn]
			dfNewData.loc[len(dfNewData)] = newData
		#####-----check that total EC counts is equal to 2
		print(pd.unique(dfNewData[['nEC0','nEC1','nEC2','nECn']].sum(axis=1)))
		#####-----add stats to DATA
		newDATA = pd.merge(DATA, dfNewData, how = 'inner', on = ['teID'], left_index =True)
		#####-----output
		OUTFILE='/mnt/lfs2/schaack/eddieho/results/trElementAnalysis/teflon/results_EC/analysis/filterTotal/teInfo/'+GT+'.t'+str(minTotal)+'.'+minP+'.teInfo.csv'
		newDATA.to_csv(path_or_buf = OUTFILE, index = False)



###################--------------------FILTER 3: Probability
################### 

##########----------Import packages
from __future__ import division
import sys
import pandas as pd
import numpy as np
from scipy.stats import binom

######Add to sys path and import denovoFunctions
sys.path.append('/mnt/lfs2/schaack/eddieho/scripts/python')
from DmagnaDict import *


##########----------get prob that focal hom is actually a het
def getProb_het(DATA, ALL, TYPE):
	listProb, listFocal = [], []
	#####loop through each TE
	for i in range(len(DATA)):
		tempTE = DATA[ALL].loc[i].tolist()
		listP, listQ, testTotal = [],[],0
		#####loop through each line to get listP and listQ for non-focal hets and total for focal hom
		FOCAL=''
		for j,LINE in enumerate(tempTE):
			#####get num of pres, abs and total reads
			pres,abs,total = map(int,LINE.split(':'))
			#####only get freqP for heterozygotes (i.e. those that have PRES and ABS)
			#####	else, it is homozygote and we need to keep 'total'
			if pres > 0 and abs > 0:
				listP.append(pres/total)
				listQ.append(abs/total)
			else:
				testTotal = total
				FOCAL = ALL[j]
		#####if TYPE is DEL or INS then THETA is mean of listP or listQ
		if TYPE=='DEL':
			THETA = np.mean(listP)
		elif TYPE=='INS':
			THETA = np.mean(listQ)
		else:
			print('Unrecognized TYPE: {}'.format(TYPE))
			break			
		#####binom prob of sampling 0 PRES or ABS given testTotal samples with prob THETA
		pHomo = binom.pmf(0,testTotal, THETA)
		listProb.append(pHomo)
		listFocal.append(FOCAL)
		#print(TYPE,testTotal,THETA,FOCAL,pHomo)
	#####-----add columns
	DATA['PROB'] = listProb
	DATA['FOCAL'] = listFocal
	return(DATA)

##########----------get prob that that none of the reads sampled in non-focal contain the alternate that exists in the focal het
def getProb_homo(DATA, ALL, TYPE):
	listProb, listFocal = [], []
	#####loop through each TE
	for i in range(len(DATA)):
		tempTE = DATA[ALL].loc[i].tolist()
		focalP, focalQ, testTotal = 0,0,0
		#####lopo through lines to get focalP and focalQ of focal het and total reads of all non-focal
		FOCAL=''
		for j,LINE in enumerate(tempTE):
			#####get num of pres, abs and total reads
			pres,abs,total = map(int,LINE.split(':'))
			#####only get P and Q for the single focal heterozygous line
			#####	else, it is non-focal homozygote and we need to keep 'total'
			if pres > 0 and abs > 0:
				focalP = pres/total
				focalQ = abs/total
				FOCAL = ALL[j]
			else:
				testTotal += total
		#####if TYPE is DEL or INS then THETA focalQ or focalP
		if TYPE=='DEL':
			THETA = focalQ
		elif TYPE=='INS':
			THETA = focalP
		else:
			print('Unrecognized TYPE: {}'.format(TYPE))
			break			
		#####binom prob of sampling 0 PRES or ABS given prob THETA across ALL reads of non-focal lines
		pHet = binom.pmf(0, testTotal, THETA)
		listProb.append(pHet)
		listFocal.append(FOCAL)
	#####-----add columns
	DATA['PROB'] = listProb
	DATA['FOCAL'] = listFocal
	return(DATA)




##########---------function to get good TE indel mutations
def getGoodIndel_MA(GT, minTotal, minPercent, maxProb):
	#####-----get passing INDELs
	#####define minP
	minP = 'p'+str(minPercent)
	#####get lines
	SC,MA,EC = ancLines[GT],maLines[GT],extLines[GT]
	ALL=[SC]+MA+EC
	nMA, nNF, nLine = len(MA),len(MA)-1,len(ALL)
	#####-----input
	INFILE='/mnt/lfs2/schaack/eddieho/results/trElementAnalysis/teflon/results_EC/analysis/filterTotal/teInfo/'+GT+'.t'+str(minTotal)+'.'+minP+'.teInfo.csv'
	DATA = pd.read_csv(INFILE)
	#####get new INS and new DEL when SC was homo
	homoINS = DATA[(DATA['SC']==0) & (DATA['nMA0']==nNF) & (DATA['nMA1']==1) & (DATA['nEC0']==2)].reset_index(drop=True)
	homoDEL = DATA[(DATA['SC']==2) & (DATA['nMA2']==nNF) & (DATA['nMA1']==1) & (DATA['nEC2']==2)].reset_index(drop=True)
	#####get new INDEL when SC was het
	hetINS = DATA[(DATA['SC']==1) & (DATA['nMA1']==nNF) & (DATA['nMA2']==1) & (DATA['nEC1']==2)].reset_index(drop=True)
	hetDEL = DATA[(DATA['SC']==1) & (DATA['nMA1']==nNF) & (DATA['nMA0']==1) & (DATA['nEC1']==2)].reset_index(drop=True)
	#####get PROB for each mut type
	homoINS = getProb_homo(homoINS,ALL,'INS')
	homoDEL = getProb_homo(homoDEL,ALL,'DEL')
	hetINS = getProb_het(hetINS,ALL,'INS')
	hetDEL = getProb_het(hetDEL,ALL,'DEL')
	#####concat all indels and label by TYPE
	allINDEL = pd.concat([homoINS,homoDEL,hetINS,hetDEL]).reset_index(drop=True)
	allINDEL['TYPE'] = ['homoINS']*len(homoINS)+['homoDEL']*len(homoDEL)+['hetINS']*len(hetINS)+['hetDEL']*len(hetDEL)
	#####filter by PROB
	allINDEL_Pass = allINDEL[allINDEL['PROB']<=maxProb].reset_index(drop=True)
	#####-----add info from genotype file of TEFLon
	#####	using SC becasue the info I need (sID, etc) are the same across all lines
	colNames = ['CONTIG','BP5','BP3','sID','cID','STRAND','refID','sc5','sc3','PRES','ABS','AMB','FREQ','teID']
	INFILE='/mnt/lfs2/schaack/eddieho/results/trElementAnalysis/teflon/results_EC/'+GT+'/genotypes/'+SC+'.genotypes.txt'
	scDATA = pd.read_table(INFILE, header=None, names=colNames)
	scDATA_less = scDATA[['teID','sID','refID','sc5','sc3']]
	scDATA_less.columns = ['teID','FAMILY','refID','sc5','sc3']
	#####merge with allINDEL_Pass
	allINDEL_Pass_Merge = pd.merge(allINDEL_Pass,scDATA_less,on='teID')
	#####-----return
	return(allINDEL_Pass_Merge)
	
	

def getGoodIndel_EC(GT, minTotal, minPercent, maxProb):
	#####-----get passing INDELs
	#####define minP
	minP = 'p'+str(minPercent)
	#####get lines
	SC,MA,EC = ancLines[GT],maLines[GT],extLines[GT]
	ALL=[SC]+MA+EC
	nMA, nNF, nLine = len(MA),len(MA)-1,len(ALL)
	#####-----input
	INFILE='/mnt/lfs2/schaack/eddieho/results/trElementAnalysis/teflon/results_EC/analysis/filterTotal/teInfo/'+GT+'.t'+str(minTotal)+'.'+minP+'.teInfo.csv'
	DATA = pd.read_csv(INFILE)
	#####get new INS and new DEL when SC was homo
	homoINS = DATA[(DATA['SC']==0) & (DATA['nMA0']==nMA) & (DATA['nEC0']==1) & (DATA['nEC1']==1)].reset_index(drop=True)
	homoDEL = DATA[(DATA['SC']==2) & (DATA['nMA2']==nMA) & (DATA['nEC2']==1) & (DATA['nEC1']==1)].reset_index(drop=True)
	#####get new INDEL when SC was het
	hetINS = DATA[(DATA['SC']==1) & (DATA['nMA1']==nMA) & (DATA['nEC2']==1) & (DATA['nEC1']==1)].reset_index(drop=True)
	hetDEL = DATA[(DATA['SC']==1) & (DATA['nMA1']==nMA) & (DATA['nEC0']==1) & (DATA['nEC1']==1)].reset_index(drop=True)
	#####get PROB for each mut type
	homoINS = getProb_homo(homoINS,ALL,'INS')
	homoDEL = getProb_homo(homoDEL,ALL,'DEL')
	hetINS = getProb_het(hetINS,ALL,'INS')
	hetDEL = getProb_het(hetDEL,ALL,'DEL')
	#####concat all indels and label by TYPE
	allINDEL = pd.concat([homoINS,homoDEL,hetINS,hetDEL]).reset_index(drop=True)
	allINDEL['TYPE'] = ['homoINS']*len(homoINS)+['homoDEL']*len(homoDEL)+['hetINS']*len(hetINS)+['hetDEL']*len(hetDEL)
	#####filter by PROB
	allINDEL_Pass = allINDEL[allINDEL['PROB']<=maxProb].reset_index(drop=True)
	#####-----add info from genotype file of TEFLon
	#####	using SC becasue the info I need (sID, etc) are the same across all lines
	colNames = ['CONTIG','BP5','BP3','sID','cID','STRAND','refID','sc5','sc3','PRES','ABS','AMB','FREQ','teID']
	INFILE='/mnt/lfs2/schaack/eddieho/results/trElementAnalysis/teflon/results_EC/'+GT+'/genotypes/'+SC+'.genotypes.txt'
	scDATA = pd.read_table(INFILE, header=None, names=colNames)
	scDATA_less = scDATA[['teID','sID','refID','sc5','sc3']]
	scDATA_less.columns = ['teID','FAMILY','refID','sc5','sc3']
	#####merge with allINDEL_Pass
	allINDEL_Pass_Merge = pd.merge(allINDEL_Pass,scDATA_less,on='teID')
	#####-----return
	return(allINDEL_Pass_Merge)


##########---------main
#####-----define parameters
maxProb = 0.005
minTotal = 20
gtList = ['FA','FB','FC','GA','GB','GC','IA','IB','IC']

#####-----loop through each combo of minPercent and GT
for GT in gtList:
	for minPercent in [10]:
		allINDEL = getGoodIndel_MA(GT,minTotal,minPercent,maxProb)
		allINDEL_EC = getGoodIndel_EC(GT,minTotal,minPercent,maxProb)
		#####-----output; MANUALLY change the maxProb name
		OUTFILE = '/mnt/lfs2/schaack/eddieho/results/trElementAnalysis/teflon/results_EC/analysis/filterTotal/filterProb/'+GT+'.MA.t'+str(minTotal)+'.p'+str(minPercent)+'.'+'p005'+'.csv'
		allINDEL.to_csv(path_or_buf = OUTFILE, index = False)
		OUTFILE_EC = '/mnt/lfs2/schaack/eddieho/results/trElementAnalysis/teflon/results_EC/analysis/filterTotal/filterProb/'+GT+'.EC.t'+str(minTotal)+'.p'+str(minPercent)+'.'+'p005'+'.csv'
		allINDEL_EC.to_csv(path_or_buf = OUTFILE_EC, index = False)
		




###################--------------------Get number of callable sites and copy number of TEs from SC lines
###################	1) All sites >= minTotal number of reads (PRES+ABS)
###################	2) SC determines if ancestral state is homozygous or heterozygous

##########----------Import packages
from __future__ import division
import sys
import pandas as pd
import numpy as np
from scipy.stats import binom

######Add to sys path and import denovoFunctions
sys.path.append('/mnt/lfs2/schaack/eddieho/scripts/python')
from DmagnaDict import *


##########----------get count of callable sites and copy number in SC for each TE class
def getCN_Class(SC,minMA,minEC,DATA):
	#####-----loop through each class
	listClass = ['DNA', 'LINE', 'LTR', 'RC', 'SINE']
	listCN = pd.DataFrame(columns=['SC','CLASS','nHomoPres','nHetPres','CN'])
	for CLA in listClass:
		#####subset DATA to specific TE class
		temp0 = DATA[DATA['CLASS']==CLA].reset_index(drop=True)
		#####get num of sites where SC is homo or het for PRESENCE allele
		nHomoPres = len(temp0[(temp0['SC']==2) & (temp0['nMAp']>=minMA) & (temp0['nECp']>=minEC)])
		nHetPres = len(temp0[(temp0['SC']==1) & (temp0['nMAp']>=minMA) & (temp0['nECp']>=minEC)])
		#####get copy number of SC
		scCN = nHomoPres*2 + nHetPres
		#####add to listCN
		listCN.loc[len(listCN)] = [SC,CLA,nHomoPres,nHetPres,scCN]
	#####-----add total row
	listCN.loc[len(listCN)] = [SC,'TOTAL',np.sum(listCN.nHomoPres),np.sum(listCN.nHetPres),np.sum(listCN.CN)]
	return(listCN)

##########----------get count of callable sites and copy number in SC for each TE family
def getCN_Family(SC,minMA,minEC,DATA):
	#####-----loop through each family
	listFamily = ['DNA', 'DNA_Academ-1', 'DNA_CMC-EnSpm', 'DNA_IS3EU', 'DNA_MULE-MuDR', 'DNA_Merlin', 'DNA_P', 'DNA_P-Fungi', 'DNA_PIF-Harbinger', 'DNA_PIF-ISL2EU', 'DNA_Sola-2', 'DNA_TcMar-Fot1', 'DNA_TcMar-Tc1', 'DNA_Zisupton', 'DNA_hAT', 'DNA_hAT-Ac', 'DNA_hAT-Charlie', 'DNA_hAT-Tip100', 'DNA_hAT-hATm', 'LINE_I', 'LINE_I-Jockey', 'LINE_L1-Tx1', 'LINE_Penelope', 'LINE_R1', 'LINE_R2-NeSL', 'LTR', 'LTR_Copia', 'LTR_DIRS', 'LTR_ERV1', 'LTR_Gypsy', 'LTR_Ngaro', 'LTR_Pao', 'RC_Helitron', 'SINE', 'SINE_ID', 'SINE_tRNA-V-CR1']
	listCN = pd.DataFrame(columns=['SC','FAMILY','nHomoPres','nHetPres','CN'])
	for FAM in listFamily:
		#####subset DATA to specific TE class
		temp0 = DATA[DATA['sID']==FAM].reset_index(drop=True)
		#####get num of sites where SC is homo or het for PRESENCE allele
		nHomoPres = len(temp0[(temp0['SC']==2) & (temp0['nMAp']>=minMA) & (temp0['nECp']>=minEC)])
		nHetPres = len(temp0[(temp0['SC']==1) & (temp0['nMAp']>=minMA) & (temp0['nECp']>=minEC)])
		#####get copy number of SC
		scCN = nHomoPres*2 + nHetPres
		#####add to listCN
		listCN.loc[len(listCN)] = [SC,FAM,nHomoPres,nHetPres,scCN]
	#####-----add total row
	listCN.loc[len(listCN)] = [SC,'TOTAL',np.sum(listCN.nHomoPres),np.sum(listCN.nHetPres),np.sum(listCN.CN)]
	return(listCN)

##########----------create mergeDATA which is a combo of teInfo and TEFLoN genitype file
def getMergeData(GT, minTotal, minPercent):
	#####-----read teInfo file
	INFILE='/mnt/lfs2/schaack/eddieho/results/trElementAnalysis/teflon/results_EC/analysis/filterTotal/teInfo/'+GT+'.t'+str(minTotal)+'.p'+str(minPercent)+'.teInfo.csv'
	DATA = pd.read_csv(INFILE)
	#####-----get number of genotypes (i.e. passed) MA lines and EC lines
	DATA['nMAp'] =  DATA[['nMA0','nMA1','nMA2']].sum(axis=1)
	DATA['nECp'] =  DATA[['nEC0','nEC1','nEC2']].sum(axis=1)
	#####-----merge with TE id info
	#####read file
	SC=ancLines[GT]
	colNames = ['CONTIG','BP5','BP3','sID','cID','STRAND','refID','sc5','sc3','PRES','ABS','AMB','FREQ','teID']
	INFILE='/mnt/lfs2/schaack/eddieho/results/trElementAnalysis/teflon/results_EC/'+GT+'/genotypes/'+SC+'.genotypes.txt'
	scDATA = pd.read_table(INFILE, header=None, names=colNames)
	#####add TE CLASS column
	scDATA['CLASS'] = scDATA.sID.str.split('_').str[0]
	#####merge
	mergeDATA = pd.merge(DATA[['teID','SC','nMAp','nECp']],scDATA[['teID','sID','CLASS']], on='teID')
	#####return
	return(mergeDATA)


##########---------main
#####-----define parameters
minTotal = 20
gtList = ['FA','FB','FC','GA','GB','GC','IA','IB','IC']

#####-----loop through each combo of minPercent and GT
for minPercent in [10]:
	dfAllCN_Class = pd.DataFrame()
	dfAllCN_Family = pd.DataFrame()
	for GT in gtList:
		#####-----get merged data
		mergeDATA = getMergeData(GT, minTotal, minPercent)
		#####-----get CN by class and family
		#####	specify minMA = 3 to be fair to FB
		SC, nMA, nEC = ancLines[GT],len(maLines[GT]),len(extLines[GT])
		minMA, minEC = nMA, nEC
		cnClass = getCN_Class(SC, minMA, minEC, mergeDATA)
		cnFamily = getCN_Family(SC, minMA, minEC, mergeDATA)
		#####-----add to df
		dfAllCN_Class = pd.concat([dfAllCN_Class, cnClass])
		dfAllCN_Family = pd.concat([dfAllCN_Family, cnFamily])
	#####-----output
	OUTFILE = '/mnt/lfs2/schaack/eddieho/results/trElementAnalysis/teflon/results_EC/analysis/filterTotal/scCN/ALL.t'+str(minTotal)+'.p'+str(minPercent)+'.MA'+str(minMA)+'.EC'+str(minEC)+'.scCN.class.csv'
	dfAllCN_Class.to_csv(path_or_buf = OUTFILE, index = False)
	OUTFILE = '/mnt/lfs2/schaack/eddieho/results/trElementAnalysis/teflon/results_EC/analysis/filterTotal/scCN/ALL.t'+str(minTotal)+'.p'+str(minPercent)+'.MA'+str(minMA)+'.EC'+str(minEC)+'.scCN.family.csv'
	dfAllCN_Family.to_csv(path_or_buf = OUTFILE, index = False)







###################--------------------Get TE mutation rate for Mutation Accumulation (MA) lines
###################	There are multiple mutation rates
###################		novel insertion, novel deletion, all novel indel
###################		existing insertion, existing deletion, all existing indel (these are liklly recomb or gene conversion rather than TE machinery)
###################		all movement (total of all novel and all existing)
###################	Rates can be calculated by 
###################		TE family, TE class, all TE

##########----------Import packages
from __future__ import division
import sys
import pandas as pd
import numpy as np
from scipy.stats import binom

######Add to sys path and import denovoFunctions
sys.path.append('/mnt/lfs2/schaack/eddieho/scripts/python')
from DmagnaDict import *



#########----------count the occurence of each mut type
def ctMutType(dfMA):
	#####defin mut types
	listType = ['homoINS','homoDEL','hetINS','hetDEL']
	#####get counts of each type
	listCount = []
	for x in listType:
		tempType = dfMA[dfMA['TYPE']==x]
		listCount.append(len(tempType))
	#####calculate counts
	nNV_GAIN,nNV_LOSS,nLOH_GAIN,nLOH_LOSS = listCount[0],listCount[1],listCount[2],listCount[3]
	nNV = listCount[0]+listCount[1]
	nLOH = listCount[2]+listCount[3]
	nGAIN = listCount[0]+listCount[2]
	nLOSS = listCount[1]+listCount[3]
	nMU = np.sum(listCount)
	#####return
	return([nNV_GAIN,nNV_LOSS,nLOH_GAIN,nLOH_LOSS,nNV,nLOH,nGAIN,nLOSS,nMU])
	

##########----------function to get count of mutations of each CLASS for GT
##########	assumes several list and dataframes defined OUTSIDE function:
##########		listClass, genDict, cnClass,
def getMutCount_Class(minTotal, minPercent, GT):
	#####-----define variables
	POP,SC,MA = GT[0],ancLines[GT],maLines[GT]
	#####-----read TE mutation file for GT
	INFILE = '/mnt/lfs2/schaack/eddieho/results/trElementAnalysis/teflon/results_EC/analysis/filterTotal/filterProb/'+GT+'.MA.t'+str(minTotal)+'.p'+str(minPercent)+'.'+'p005'+'.csv'
	DATA = pd.read_csv(INFILE)
	DATA['CLASS'] = DATA.FAMILY.str.split('_').str[0]
	#####-----get CN of each class (exclude TOTAL)
	cnClass_SC = cnClass[(cnClass['SC']==SC) & (cnClass['CLASS']!='TOTAL')]
	nClass = len(cnClass_SC)
	#####-----loop through each MA line and get mut count for each class
	dfMutCount = pd.DataFrame()
	for ma in MA:
		#####subset to MA line and get num gen MA
		tempMA = DATA[DATA['FOCAL']==ma]
		GEN = genDict[ma]
		#####if line has mut then count otherwise give zero
		if len(tempMA)>0:
			#####count number of mut types for each class
			mutCount = pd.DataFrame(columns=['POP','GT','LINE','CLASS','scCN','maGEN','NV_GAIN','NV_LOSS','LOH_GAIN','LOH_LOSS','NV','LOH','GAIN','LOSS','MU'])
			for CLA in listClass:
				#####subset MA info for class
				tempCla = tempMA[tempMA['CLASS']==CLA]
				#####get SC CN for class
				scCN = int(cnClass_SC[cnClass_SC['CLASS']==CLA]['CN'])
				#####get mut counts
				if len(tempCla)>0:
					mutCount.loc[len(mutCount)] = [POP,GT,ma,CLA,scCN,GEN]+ctMutType(tempCla)
				else:
					mutCount.loc[len(mutCount)] = [POP,GT,ma,CLA,scCN,GEN,0,0,0,0,0,0,0,0,0]
		else:
			mutCount = pd.DataFrame({'POP':[POP]*nClass,'GT':[GT]*nClass,'LINE':[ma]*nClass,'CLASS':cnClass_SC.CLASS,'scCN':cnClass_SC.CN,'maGEN':[GEN]*nClass,
			'NV_GAIN':[0]*nClass,'NV_LOSS':[0]*nClass,'LOH_GAIN':[0]*nClass,'LOH_LOSS':[0]*nClass,
			'NV':[0]*nClass,'LOH':[0]*nClass,'GAIN':[0]*nClass,'LOSS':[0]*nClass,'MU':[0]*nClass})
		#####append to dfMutCount
		dfMutCount = dfMutCount.append(mutCount).reset_index(drop=True)
	#####-----reorder columns and return
	dfMutCount = dfMutCount[['POP','GT','LINE','CLASS','scCN','maGEN','NV_GAIN','NV_LOSS','LOH_GAIN','LOH_LOSS','NV','LOH','GAIN','LOSS','MU']]
	dfMutCount['DENOM'] = dfMutCount['scCN']*dfMutCount['maGEN']
	return(dfMutCount)

##########----------function to get count of mutations of each FAMILY for GT
##########	assumes several list and dataframes defined OUTSIDE function:
##########		listClass, genDict, cnClass,
def getMutCount_Family(minTotal, minPercent, GT):
	#####-----define variables
	POP,SC,MA = GT[0],ancLines[GT],maLines[GT]
	#####-----read TE mutation file for GT
	INFILE = '/mnt/lfs2/schaack/eddieho/results/trElementAnalysis/teflon/results_EC/analysis/filterTotal/filterProb/'+GT+'.MA.t'+str(minTotal)+'.p'+str(minPercent)+'.'+'p005'+'.csv'
	DATA = pd.read_csv(INFILE)
	DATA['CLASS'] = DATA.FAMILY.str.split('_').str[0]
	#####-----get CN of each family (exclude TOTAL)
	cnFamily_SC = cnFamily[(cnFamily['SC']==SC) & (cnFamily['FAMILY']!='TOTAL')]
	nFamily = len(cnFamily_SC)
	#####-----loop through each MA line and get mut count for each family
	dfMutCount = pd.DataFrame()
	for ma in MA:
		#####subset to MA line and get num gen MA
		tempMA = DATA[DATA['FOCAL']==ma]
		GEN = genDict[ma]
		#####if line has mut then count otherwise give zero
		if len(tempMA)>0:
			#####count number of mut types for each family
			mutCount = pd.DataFrame(columns=['POP','GT','LINE','FAMILY','scCN','maGEN','NV_GAIN','NV_LOSS','LOH_GAIN','LOH_LOSS','NV','LOH','GAIN','LOSS','MU'])
			for FAM in listFamily:
				#####subset MA info for family
				tempFam = tempMA[tempMA['FAMILY']==FAM]
				#####get SC CN for family
				scCN = int(cnFamily_SC[cnFamily_SC['FAMILY']==FAM]['CN'])
				#####get mut counts
				if len(tempFam)>0:
					mutCount.loc[len(mutCount)] = [POP,GT,ma,FAM,scCN,GEN]+ctMutType(tempFam)
				else:
					mutCount.loc[len(mutCount)] = [POP,GT,ma,FAM,scCN,GEN,0,0,0,0,0,0,0,0,0]
		else:
			mutCount = pd.DataFrame({'POP':[POP]*nFamily,'GT':[GT]*nFamily,'LINE':[ma]*nFamily,'FAMILY':cnFamily_SC.FAMILY,'scCN':cnFamily_SC.CN,'maGEN':[GEN]*nFamily,
			'NV_GAIN':[0]*nFamily,'NV_LOSS':[0]*nFamily,'LOH_GAIN':[0]*nFamily,'LOH_LOSS':[0]*nFamily,
			'NV':[0]*nFamily,'LOH':[0]*nFamily,'GAIN':[0]*nFamily,'LOSS':[0]*nFamily,'MU':[0]*nFamily})
		#####append to dfMutCount
		dfMutCount = dfMutCount.append(mutCount).reset_index(drop=True)
	#####-----reorder columns and return
	dfMutCount = dfMutCount[['POP','GT','LINE','FAMILY','scCN','maGEN','NV_GAIN','NV_LOSS','LOH_GAIN','LOH_LOSS','NV','LOH','GAIN','LOSS','MU']]
	dfMutCount['DENOM'] = dfMutCount['scCN']*dfMutCount['maGEN']
	return(dfMutCount)
	


##########-----Define number MA gen based on LINE
genDATA = pd.read_csv('/mnt/lfs2/schaack/eddieho/config/genMA.csv')
genDict = dict(zip(genDATA['LINE'], genDATA['GEN']))

##########----------define parameters
listClass = ['DNA', 'LINE', 'LTR', 'RC', 'SINE']
listFamily = ['DNA', 'DNA_Academ-1', 'DNA_CMC-EnSpm', 'DNA_IS3EU', 'DNA_MULE-MuDR', 'DNA_Merlin', 'DNA_P', 'DNA_P-Fungi', 'DNA_PIF-Harbinger', 'DNA_PIF-ISL2EU', 'DNA_Sola-2', 'DNA_TcMar-Fot1', 'DNA_TcMar-Tc1', 'DNA_Zisupton', 'DNA_hAT', 'DNA_hAT-Ac', 'DNA_hAT-Charlie', 'DNA_hAT-Tip100', 'DNA_hAT-hATm', 'LINE_I', 'LINE_I-Jockey', 'LINE_L1-Tx1', 'LINE_Penelope', 'LINE_R1', 'LINE_R2-NeSL', 'LTR', 'LTR_Copia', 'LTR_DIRS', 'LTR_ERV1', 'LTR_Gypsy', 'LTR_Ngaro', 'LTR_Pao', 'RC_Helitron', 'SINE', 'SINE_ID', 'SINE_tRNA-V-CR1']
minTotal = 20
minPercent = 10
minMA = 8
minEC = 2

##########----------read SC CN file based on minTotal, minPercent, minMA, minEC
INFILE='/mnt/lfs2/schaack/eddieho/results/trElementAnalysis/teflon/results_EC/analysis/filterTotal/scCN/ALL.t'+str(minTotal)+'.p'+str(minPercent)+'.MA'+str(minMA)+'.EC'+str(minEC)+'.scCN.class.csv'
cnClass = pd.read_csv(INFILE)
INFILE='/mnt/lfs2/schaack/eddieho/results/trElementAnalysis/teflon/results_EC/analysis/filterTotal/scCN/ALL.t'+str(minTotal)+'.p'+str(minPercent)+'.MA'+str(minMA)+'.EC'+str(minEC)+'.scCN.family.csv'
cnFamily = pd.read_csv(INFILE)


##########----------Get mut counts and mut rates for each GT by CLASS and by FAMILY
gtList = ['FA','FB','FC','GA','GB','GC','IA','IB','IC']
allCount_Class,allCount_Family,allRate_Class,allRate_Family = pd.DataFrame(),pd.DataFrame(),pd.DataFrame(),pd.DataFrame()
for GT in gtList:
	print('Processing: {}'.format(GT))
	#####get mut counts
	mutCount_Class = getMutCount_Class(minTotal,minPercent,GT)
	mutCount_Family = getMutCount_Family(minTotal,minPercent,GT)
	#####get mut rates (exclude cases where DENOM==0, since SC CN = 0)
	mutRate_Class = mutCount_Class[mutCount_Class['DENOM']>0]
	mutRate_Family = mutCount_Family[mutCount_Family['DENOM']>0]
	mutRate_Class[['nvgRate','nvlRate','lohgRate','lohlRate','nvRate','lohRate','gRate','lRate','muRate']] = mutRate_Class[['NV_GAIN','NV_LOSS','LOH_GAIN','LOH_LOSS','NV','LOH','GAIN','LOSS','MU']].div(mutRate_Class['DENOM'],axis=0)
	mutRate_Family[['nvgRate','nvlRate','lohgRate','lohlRate','nvRate','lohRate','gRate','lRate','muRate']] = mutRate_Family[['NV_GAIN','NV_LOSS','LOH_GAIN','LOH_LOSS','NV','LOH','GAIN','LOSS','MU']].divide(mutRate_Family['DENOM'],axis=0)
	#####append
	allCount_Class = allCount_Class.append(mutCount_Class).reset_index(drop=True)
	allCount_Family = allCount_Family.append(mutCount_Family).reset_index(drop=True)
	allRate_Class = allRate_Class.append(mutRate_Class).reset_index(drop=True)
	allRate_Family = allRate_Family.append(mutRate_Family).reset_index(drop=True)

#####-----output
DIR='/mnt/lfs2/schaack/eddieho/results/trElementAnalysis/teflon/results_EC/analysis/filterTotal/mutRate/'
allCount_Class.to_csv(path_or_buf = DIR+'ALL.MA.t'+str(minTotal)+'.p'+str(minPercent)+'.p005.count.class.csv', index = False)
allCount_Family.to_csv(path_or_buf = DIR+'ALL.MA.t'+str(minTotal)+'.p'+str(minPercent)+'.p005.count.family.csv', index = False)
allRate_Class.to_csv(path_or_buf = DIR+'ALL.MA.t'+str(minTotal)+'.p'+str(minPercent)+'.p005.rate.class.csv', index = False)
allRate_Family.to_csv(path_or_buf = DIR+'ALL.MA.t'+str(minTotal)+'.p'+str(minPercent)+'.p005.rate.family.csv', index = False)






###################--------------------get TE mutation rate for Extant Control (EC) LINES
###################	There are multiple mutation rates
###################		novel insertion, novel deletion, all novel indel
###################		existing insertion, existing deletion, all existing indel (these are liklly recomb or gene conversion rather than TE machinery)
###################		all movement (total of all novel and all existing)
###################	Rates can be calculated by 
###################		TE family, TE class, all TE



##########----------Import packages
from __future__ import division
import sys
import pandas as pd
import numpy as np
from scipy.stats import binom

######Add to sys path and import denovoFunctions
sys.path.append('/mnt/lfs2/schaack/eddieho/scripts/python')
from DmagnaDict import *



#########----------count the occurence of each mut type
def ctMutType(dfMA):
	#####defin mut types
	listType = ['homoINS','homoDEL','hetINS','hetDEL']
	#####get counts of each type
	listCount = []
	for x in listType:
		tempType = dfMA[dfMA['TYPE']==x]
		listCount.append(len(tempType))
	#####calculate counts
	nNV_GAIN,nNV_LOSS,nLOH_GAIN,nLOH_LOSS = listCount[0],listCount[1],listCount[2],listCount[3]
	nNV = listCount[0]+listCount[1]
	nLOH = listCount[2]+listCount[3]
	nGAIN = listCount[0]+listCount[2]
	nLOSS = listCount[1]+listCount[3]
	nMU = np.sum(listCount)
	#####return
	return([nNV_GAIN,nNV_LOSS,nLOH_GAIN,nLOH_LOSS,nNV,nLOH,nGAIN,nLOSS,nMU])
	

##########----------function to get count of mutations of each CLASS for GT
##########	assumes several list and dataframes defined OUTSIDE function:
##########		listClass, genDict, cnClass,
def getMutCount_Class(minTotal, minPercent, GT):
	#####-----define variables
	POP,SC,EC = GT[0],ancLines[GT],extLines[GT]
	#####-----read TE mutation file for GT
	INFILE = '/mnt/lfs2/schaack/eddieho/results/trElementAnalysis/teflon/results_EC/analysis/filterTotal/filterProb/'+GT+'.EC.t'+str(minTotal)+'.p'+str(minPercent)+'.'+'p005'+'.csv'
	DATA = pd.read_csv(INFILE)
	DATA['CLASS'] = DATA.FAMILY.str.split('_').str[0]
	#####-----get CN of each class (exclude TOTAL)
	cnClass_SC = cnClass[(cnClass['SC']==SC) & (cnClass['CLASS']!='TOTAL')]
	nClass = len(cnClass_SC)
	#####-----loop through each EC line and get mut count for each class
	dfMutCount = pd.DataFrame()
	for ec in EC:
		#####subset to MA line and get num gen EC
		tempMA = DATA[DATA['FOCAL']==ec]
		GEN = genDict[ec]
		#####if line has mut then count otherwise give zero
		if len(tempMA)>0:
			#####count number of mut types for each class
			mutCount = pd.DataFrame(columns=['POP','GT','LINE','CLASS','scCN','maGEN','NV_GAIN','NV_LOSS','LOH_GAIN','LOH_LOSS','NV','LOH','GAIN','LOSS','MU'])
			for CLA in listClass:
				#####subset EC info for class
				tempCla = tempMA[tempMA['CLASS']==CLA]
				#####get SC CN for class
				scCN = int(cnClass_SC[cnClass_SC['CLASS']==CLA]['CN'])
				#####get mut counts
				if len(tempCla)>0:
					mutCount.loc[len(mutCount)] = [POP,GT,ec,CLA,scCN,GEN]+ctMutType(tempCla)
				else:
					mutCount.loc[len(mutCount)] = [POP,GT,ec,CLA,scCN,GEN,0,0,0,0,0,0,0,0,0]
		else:
			mutCount = pd.DataFrame({'POP':[POP]*nClass,'GT':[GT]*nClass,'LINE':[ec]*nClass,'CLASS':cnClass_SC.CLASS,'scCN':cnClass_SC.CN,'maGEN':[GEN]*nClass,
			'NV_GAIN':[0]*nClass,'NV_LOSS':[0]*nClass,'LOH_GAIN':[0]*nClass,'LOH_LOSS':[0]*nClass,
			'NV':[0]*nClass,'LOH':[0]*nClass,'GAIN':[0]*nClass,'LOSS':[0]*nClass,'MU':[0]*nClass})
		#####append to dfMutCount
		dfMutCount = dfMutCount.append(mutCount).reset_index(drop=True)
	#####-----reorder columns and return
	dfMutCount = dfMutCount[['POP','GT','LINE','CLASS','scCN','maGEN','NV_GAIN','NV_LOSS','LOH_GAIN','LOH_LOSS','NV','LOH','GAIN','LOSS','MU']]
	dfMutCount['DENOM'] = dfMutCount['scCN']*dfMutCount['maGEN']
	return(dfMutCount)

##########----------function to get count of mutations of each FAMILY for GT
##########	assumes several list and dataframes defined OUTSIDE function:
##########		listClass, genDict, cnClass,
def getMutCount_Family(minTotal, minPercent, GT):
	#####-----define variables
	POP,SC,EC = GT[0],ancLines[GT],extLines[GT]
	#####-----read TE mutation file for GT
	INFILE = '/mnt/lfs2/schaack/eddieho/results/trElementAnalysis/teflon/results_EC/analysis/filterTotal/filterProb/'+GT+'.EC.t'+str(minTotal)+'.p'+str(minPercent)+'.'+'p005'+'.csv'
	DATA = pd.read_csv(INFILE)
	DATA['CLASS'] = DATA.FAMILY.str.split('_').str[0]
	#####-----get CN of each family (exclude TOTAL)
	cnFamily_SC = cnFamily[(cnFamily['SC']==SC) & (cnFamily['FAMILY']!='TOTAL')]
	nFamily = len(cnFamily_SC)
	#####-----loop through each EC line and get mut count for each family
	dfMutCount = pd.DataFrame()
	for ec in EC:
		#####subset to EC line and get num gen EC
		tempMA = DATA[DATA['FOCAL']==ec]
		GEN = genDict[ec]
		#####if line has mut then count otherwise give zero
		if len(tempMA)>0:
			#####count number of mut types for each family
			mutCount = pd.DataFrame(columns=['POP','GT','LINE','FAMILY','scCN','maGEN','NV_GAIN','NV_LOSS','LOH_GAIN','LOH_LOSS','NV','LOH','GAIN','LOSS','MU'])
			for FAM in listFamily:
				#####subset EC info for family
				tempFam = tempMA[tempMA['FAMILY']==FAM]
				#####get SC CN for family
				scCN = int(cnFamily_SC[cnFamily_SC['FAMILY']==FAM]['CN'])
				#####get mut counts
				if len(tempFam)>0:
					mutCount.loc[len(mutCount)] = [POP,GT,ec,FAM,scCN,GEN]+ctMutType(tempFam)
				else:
					mutCount.loc[len(mutCount)] = [POP,GT,ec,FAM,scCN,GEN,0,0,0,0,0,0,0,0,0]
		else:
			mutCount = pd.DataFrame({'POP':[POP]*nFamily,'GT':[GT]*nFamily,'LINE':[ec]*nFamily,'FAMILY':cnFamily_SC.FAMILY,'scCN':cnFamily_SC.CN,'maGEN':[GEN]*nFamily,
			'NV_GAIN':[0]*nFamily,'NV_LOSS':[0]*nFamily,'LOH_GAIN':[0]*nFamily,'LOH_LOSS':[0]*nFamily,
			'NV':[0]*nFamily,'LOH':[0]*nFamily,'GAIN':[0]*nFamily,'LOSS':[0]*nFamily,'MU':[0]*nFamily})
		#####append to dfMutCount
		dfMutCount = dfMutCount.append(mutCount).reset_index(drop=True)
	#####-----reorder columns and return
	dfMutCount = dfMutCount[['POP','GT','LINE','FAMILY','scCN','maGEN','NV_GAIN','NV_LOSS','LOH_GAIN','LOH_LOSS','NV','LOH','GAIN','LOSS','MU']]
	dfMutCount['DENOM'] = dfMutCount['scCN']*dfMutCount['maGEN']
	return(dfMutCount)




##########-----Define number MA gen based on LINE
genDATA = pd.read_csv('/mnt/lfs2/schaack/eddieho/config/genEC.csv')
genDict = dict(zip(genDATA['LINE'], genDATA['GEN']))

##########----------define parameters
listClass = ['DNA', 'LINE', 'LTR', 'RC', 'SINE']
listFamily = ['DNA', 'DNA_Academ-1', 'DNA_CMC-EnSpm', 'DNA_IS3EU', 'DNA_MULE-MuDR', 'DNA_Merlin', 'DNA_P', 'DNA_P-Fungi', 'DNA_PIF-Harbinger', 'DNA_PIF-ISL2EU', 'DNA_Sola-2', 'DNA_TcMar-Fot1', 'DNA_TcMar-Tc1', 'DNA_Zisupton', 'DNA_hAT', 'DNA_hAT-Ac', 'DNA_hAT-Charlie', 'DNA_hAT-Tip100', 'DNA_hAT-hATm', 'LINE_I', 'LINE_I-Jockey', 'LINE_L1-Tx1', 'LINE_Penelope', 'LINE_R1', 'LINE_R2-NeSL', 'LTR', 'LTR_Copia', 'LTR_DIRS', 'LTR_ERV1', 'LTR_Gypsy', 'LTR_Ngaro', 'LTR_Pao', 'RC_Helitron', 'SINE', 'SINE_ID', 'SINE_tRNA-V-CR1']
minTotal = 20
minPercent = 10
minMA = 8
minEC = 2

##########----------read SC CN file based on minTotal, minPercent, minMA, minEC
INFILE='/mnt/lfs2/schaack/eddieho/results/trElementAnalysis/teflon/results_EC/analysis/filterTotal/scCN/ALL.t'+str(minTotal)+'.p'+str(minPercent)+'.MA'+str(minMA)+'.EC'+str(minEC)+'.scCN.class.csv'
cnClass = pd.read_csv(INFILE)
INFILE='/mnt/lfs2/schaack/eddieho/results/trElementAnalysis/teflon/results_EC/analysis/filterTotal/scCN/ALL.t'+str(minTotal)+'.p'+str(minPercent)+'.MA'+str(minMA)+'.EC'+str(minEC)+'.scCN.family.csv'
cnFamily = pd.read_csv(INFILE)


##########----------Get mut counts and mut rates for each GT by CLASS and by FAMILY
gtList = ['FA','FB','FC','GA','GB','GC','IA','IB','IC']
allCount_Class,allCount_Family,allRate_Class,allRate_Family = pd.DataFrame(),pd.DataFrame(),pd.DataFrame(),pd.DataFrame()
for GT in gtList:
	print('Processing: {}'.format(GT))
	#####get mut counts
	mutCount_Class = getMutCount_Class(minTotal,minPercent,GT)
	mutCount_Family = getMutCount_Family(minTotal,minPercent,GT)
	#####get mut rates (exclude cases where DENOM==0, since SC CN = 0)
	mutRate_Class = mutCount_Class[mutCount_Class['DENOM']>0]
	mutRate_Family = mutCount_Family[mutCount_Family['DENOM']>0]
	mutRate_Class[['nvgRate','nvlRate','lohgRate','lohlRate','nvRate','lohRate','gRate','lRate','muRate']] = mutRate_Class[['NV_GAIN','NV_LOSS','LOH_GAIN','LOH_LOSS','NV','LOH','GAIN','LOSS','MU']].div(mutRate_Class['DENOM'],axis=0)
	mutRate_Family[['nvgRate','nvlRate','lohgRate','lohlRate','nvRate','lohRate','gRate','lRate','muRate']] = mutRate_Family[['NV_GAIN','NV_LOSS','LOH_GAIN','LOH_LOSS','NV','LOH','GAIN','LOSS','MU']].divide(mutRate_Family['DENOM'],axis=0)
	#####append
	allCount_Class = allCount_Class.append(mutCount_Class).reset_index(drop=True)
	allCount_Family = allCount_Family.append(mutCount_Family).reset_index(drop=True)
	allRate_Class = allRate_Class.append(mutRate_Class).reset_index(drop=True)
	allRate_Family = allRate_Family.append(mutRate_Family).reset_index(drop=True)

#####-----output
DIR='/mnt/lfs2/schaack/eddieho/results/trElementAnalysis/teflon/results_EC/analysis/filterTotal/mutRate/'
allCount_Class.to_csv(path_or_buf = DIR+'ALL.EC.t'+str(minTotal)+'.p'+str(minPercent)+'.p005.count.class.csv', index = False)
allCount_Family.to_csv(path_or_buf = DIR+'ALL.EC.t'+str(minTotal)+'.p'+str(minPercent)+'.p005.count.family.csv', index = False)
allRate_Class.to_csv(path_or_buf = DIR+'ALL.EC.t'+str(minTotal)+'.p'+str(minPercent)+'.p005.rate.class.csv', index = False)
allRate_Family.to_csv(path_or_buf = DIR+'ALL.EC.t'+str(minTotal)+'.p'+str(minPercent)+'.p005.rate.family.csv', index = False)


