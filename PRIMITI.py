
#import and load
import os
import csv
import sys
import time
import math

import joblib
import pandas as pd
import numpy as np
import argparse
import pickle
import subprocess
import scipy as sp
import datetime

from PIL import Image, ImageDraw, ImageFont

import warnings
warnings.filterwarnings("ignore")


from subprocess import Popen, PIPE
from Bio import SeqIO
from Bio.Seq import transcribe
from Bio.Seq import Seq

import re
import multiprocessing as mp

from joblib import load

from config import MainConfig


#PRIMITI-TS 1. Prepare the list for miRNA-transcript pairs ------------------------------------------------------------------------------------##


#PRIMITI-TS 2. Run a code here to get a first set of features as a CSV file -------------------------------------------------------------------##


#PRIMITI-TS 3. Generate csv files for a list of miRNA and a list of transcripts ---------------------------------------------------------------##


#PRIMITI-TS 4. Call a python 3 iLearn script that use csv files to generate a file contain miRNA & a file contain mRNA representations --------##


#PRIMITI-TS 5. Run a code to retreive informations from 2 files and integrate as a part of features -> get a final set of features ------------##


#PRIMITI-TS 6. Standardize the data using pre-trained stardardizer ----------------------------------------------------------------------------##


#PRIMITI-TS 7. Predict the probability for each miRNA-target site pairs using a pre-trained model and store it in the memory ------------------##


#PRIMITI-TM 8. Calculate the PRIMITI-TM feature for each miRNA-transcript pairs based on probability of each target site ----------------------##


#PRIMITI-TM 9. Predict the probability for each miRNA-transcript pairs using a pre-trained model. ---------------------------------------------##


#PRIMITI-TM 10. Process them into an output format that contains a overall score & score for each site ----------------------------------------##



#####################FOR STEP1######################################################################################################################


def write_error_file(error_inputs, error_type, mode, job_id):
    errors = list(set(error_inputs))
    file_errors_path = "result" + "/errors/" + str(job_id) + "_errors_miRNAs_and_transcripts.txt"
    file_errors = open(file_errors_path, mode)    
    if(len(errors) > 0):
    	message = "# The following " + error_type + "(s) was (were) not supported"
    	if(error_type == "input_1"):
    		message += " --- The types of inputs are not recognized. Please check if the inputs are in correct forms"
        if(error_type == "input_2"):
            message += " --- The types of inputs are recognized; however, no matched supported ENSTs are found. Please check if the inputs specify protein-coding transcripts"
        if(error_type == "input_3"):
            message += " --- The provided ENST are not supported. Please check if the inputs are in correct forms and specify protein-coding transcripts"
    	message += ": \n"
        file_errors.write(message)
        file_errors.write(str(errors[0]))
        i = 1

        while (i < len(errors)):
            file_errors.write(", " + str(errors[i]))
            i += 1
        file_errors.write("\n\n")    
    else:
         file_errors.write("# No " + error_type + " errors were found. \n")

    file_errors.close()

    return file_errors_path


def check_empty_list(list_mirna_or_transcript):

    if not type(list_mirna_or_transcript) is np.ndarray:
        if list_mirna_or_transcript == 'none':
            return "False"
    
    
    if type(list_mirna_or_transcript) is np.ndarray:
        if list_mirna_or_transcript.shape[0] == 0:
            return "False"

        else:
            return "True"


def check_input_type_transcript(input_list, map_ENSG_to_ENST, map_RefSeq_to_ENST, map_GeneSymbol_to_ENSG, map_HGNCGene_to_ENSG, map_MIMGene_to_ENSG, map_GeneID_to_ENSG):

    input_type_list = []

    for input_transcript in input_list:

        if 'ENST' in input_transcript:
            input_type_list.append('ENST')

        elif input_transcript in map_ENSG_to_ENST.keys():
            input_type_list.append('ENSG') 

        elif input_transcript in map_RefSeq_to_ENST.keys():
            input_type_list.append('RefSeq') 

        elif input_transcript in map_GeneSymbol_to_ENSG.keys():
            input_type_list.append('GeneSymbol') 

        elif 'NCBIGeneID:' in input_transcript:
            input_type_list.append('NCBIGeneID')

        elif 'HGNC:' in input_transcript:
            input_type_list.append('HGNC') 

        elif 'MIM:' in input_transcript:
            input_type_list.append('MIMGene') 

        else:
            input_type_list.append('error') 

    return input_type_list


def check_input_type_miRNA(miRNA_list, map_miRNA_ID_to_name):
    input_type_list = []

    for input_miRNA in miRNA_list:
        if 'hsa' in input_miRNA:
            input_type_list.append('miRNA_name')
        elif input_miRNA in map_miRNA_ID_to_name.keys():
            input_type_list.append('miRNA_ID')
        else:
            input_type_list.append('error')

    return input_type_list

def get_list_miRNAs_transcripts(miRNA_list, transcript_list, supported_miRNA_list, supported_transcript_list):
    
    input_list = transcript_list

    new_miRNA_list = np.array([], dtype=str)
    transcript_list = np.array([],dtype=str)

    error_miRNAs = np.array([],dtype=str)
    error_transcripts_type1 = np.array([], dtype=str) #For bad input
    error_transcripts_type2 = np.array([], dtype=str) #For good input, but no supported ENST
    error_transcripts_type3 = np.array([], dtype=str) 

    #Check input type for miRNAs & Change them all to miRNA name
    map_miRNA_ID_to_name_path = "code/Data_file/Mapping/ID_map/miRNA_ID_to_name_dict.pkl"
    
    with open(map_miRNA_ID_to_name_path, 'rb') as file:
        map_miRNA_ID_to_name = pickle.load(file)

    input_type_list = check_input_type_miRNA(miRNA_list, map_miRNA_ID_to_name)

    i = 0
    for input_type in input_type_list:
        if input_type == 'miRNA_name':
            miRNA_name = miRNA_list[i]
            new_miRNA_list = np.append(new_miRNA_list, miRNA_name)

        elif input_type == 'miRNA_ID':
            miRNA_ID = miRNA_list[i]
            miRNA_name = map_miRNA_ID_to_name.get(miRNA_ID, 'error')
            new_miRNA_list = np.append(new_miRNA_list, miRNA_name)


        elif input_type == 'error':
            error_miRNAs = np.append(error_miRNAs, miRNA_list[i])

        else:
            error_miRNAs = np.append(error_miRNAs, miRNA_list[i])

        i += 1  

    miRNA_list = new_miRNA_list

    #Check input type for transcript

    map_ENSG_to_ENST_path = "code/Data_file/Mapping/ID_map/ENSG_to_ENST_dict.pkl"
    map_RefSeq_to_ENST_path = "code/Data_file/Mapping/ID_map/NCBI_RefSeq_to_ENST_dict.pkl"
    map_GeneSymbol_to_ENSG_path = "code/Data_file/Mapping/ID_map/GeneSymbol_to_ENSG_dict.pkl"
    map_HGNCGene_to_ENSG_path = "code/Data_file/Mapping/ID_map/HGNCGene_to_ENSG_dict.pkl"
    map_MIMGene_to_ENSG_path = "code/Data_file/Mapping/ID_map/MIMGene_to_ENSG_dict.pkl"
    map_GeneID_to_ENSG_path = "code/Data_file/Mapping/ID_map/NCBIGeneID_to_ENSG_dict.pkl"

    with open(map_ENSG_to_ENST_path, 'rb') as file:
        map_ENSG_to_ENST = pickle.load(file)

    with open(map_RefSeq_to_ENST_path, 'rb') as file:
        map_RefSeq_to_ENST = pickle.load(file)

    with open(map_GeneSymbol_to_ENSG_path, 'rb') as file:
        map_GeneSymbol_to_ENSG = pickle.load(file)

    with open(map_HGNCGene_to_ENSG_path, 'rb') as file:
        map_HGNCGene_to_ENSG = pickle.load(file)

    with open(map_MIMGene_to_ENSG_path, 'rb') as file:
        map_MIMGene_to_ENSG = pickle.load(file)

    with open(map_GeneID_to_ENSG_path, 'rb') as file:
        map_GeneID_to_ENSG = pickle.load(file)

    input_type_list = check_input_type_transcript(input_list, map_ENSG_to_ENST, map_RefSeq_to_ENST, map_GeneSymbol_to_ENSG, map_HGNCGene_to_ENSG, map_MIMGene_to_ENSG, map_GeneID_to_ENSG)

    map_ENST_to_input = {}

    i = 0
    for input_type in input_type_list:
        if input_type == 'ENST':
            transcript = input_list[i]
            if transcript in supported_transcript_list:
                transcript_list = np.append(transcript_list, transcript)
                map_ENST_to_input[transcript] = transcript
            else:
                error_transcripts_type3 = np.append(error_transcripts_type3, input_list[i])


        elif input_type == 'ENSG':
            ENSG = input_list[i]
            transcript_match = map_ENSG_to_ENST.get(ENSG, 'error')
            found = 0
            for transcript in transcript_match:
                if transcript in supported_transcript_list:
                    found += 1
                    transcript_list = np.append(transcript_list, transcript)
                    map_ENST_to_input[transcript] = input_list[i]
            if found == 0:
                error_transcripts_type2 = np.append(error_transcripts_type2, input_list[i])

        elif input_type == 'RefSeq':
            transcript_match = map_RefSeq_to_ENST.get(input_list[i], 'error')
            found = 0
            for transcript in transcript_match:
                if transcript in supported_transcript_list:
                    found += 1
                    transcript_list = np.append(transcript_list, transcript)
                    map_ENST_to_input[transcript] = input_list[i]
            if found == 0:
                error_transcripts_type2 = np.append(error_transcripts_type2, input_list[i])

        elif input_type == 'GeneSymbol':
            ENSG = map_GeneSymbol_to_ENSG.get(input_list[i], 'error')
            found = 0
            for member in ENSG:
                transcript_match = map_ENSG_to_ENST.get(member, 'error')
                for transcript in transcript_match:
                    if transcript in supported_transcript_list:
                        found += 1
                        transcript_list = np.append(transcript_list, transcript)
                        map_ENST_to_input[transcript] = input_list[i]
            if found == 0: 
                error_transcripts_type2 = np.append(error_transcripts_type2, input_list[i])

        elif input_type == 'NCBIGeneID':
            query = input_list[i].split('NCBIGeneID:')[1]
            ENSG = map_GeneID_to_ENSG.get(int(query), 'error')
            if ENSG == 'error':
                error_transcripts_type1 = np.append(error_transcripts_type1, input_list[i])
                i += 1
                continue
            else:
                found = 0
                for member in ENSG:
                    transcript_match = map_ENSG_to_ENST.get(member, 'error')
                    for transcript in transcript_match:
                        if transcript in supported_transcript_list:
                            found += 1
                            transcript_list = np.append(transcript_list, transcript)
                            map_ENST_to_input[transcript] = input_list[i]
            if found == 0:
                error_transcripts_type2 = np.append(error_transcripts_type2, input_list[i])


        elif input_type == 'HGNC':
            query = input_list[i].split('HGNC:')[1]
            ENSG = map_HGNCGene_to_ENSG.get(query, 'error')
            if ENSG == 'error':
                error_transcripts_type1 = np.append(error_transcripts_type1, input_list[i])
                i += 1
                continue
            else:
                found = 0
                for member in ENSG:
                    transcript_match = map_ENSG_to_ENST.get(member, 'error')
                    for transcript in transcript_match:
                        if transcript in supported_transcript_list:
                            found += 1
                            transcript_list = np.append(transcript_list, transcript)
                            map_ENST_to_input[transcript] = input_list[i]
            if found == 0:
                error_transcripts_type2 = np.append(error_transcripts_type2, input_list[i])

        elif input_type == 'MIMGene':
            query = input_list[i].split('MIM:')[1]
            ENSG = map_MIMGene_to_ENSG.get(query, 'error')
            if ENSG == 'error':
                error_transcripts_type1 = np.append(error_transcripts_type1, input_list[i])
                i += 1
                continue
            else:
                found = 0
                for member in ENSG:
                    transcript_match = map_ENSG_to_ENST.get(member, 'error')
                    for transcript in transcript_match:
                        if transcript in supported_transcript_list:
                            found += 1
                            transcript_list = np.append(transcript_list, transcript)
                            map_ENST_to_input[transcript] = input_list[i]
            if found == 0:
                error_transcripts_type2 = np.append(error_transcripts_type2, input_list[i])

        elif input_type == 'error':
            error_transcripts_type1 = np.append(error_transcripts_type1, input_list[i])
        
        else:
            error_transcripts_type1 = np.append(error_transcripts_type1, input_list[i])
        
        i += 1

    i = 0
    for member in miRNA_list:
        if not member in supported_miRNA_list:
            miRNA_list = np.delete(miRNA_list, i)
            error_miRNAs = np.append(error_miRNAs, member)
            i -= 1
        i += 1


    miRNA_list = miRNA_list.astype('str')
    transcript_list = transcript_list.astype('str')

    miRNA_list = np.unique(miRNA_list)
    transcript_list = np.unique(transcript_list)

    return miRNA_list, transcript_list, map_ENST_to_input, error_miRNAs, error_transcripts_type1, error_transcripts_type2, error_transcripts_type3


#####################FOR STEP2######################################################################################################################

#set variables for step 2
PC30_path_pt1 = "code/Data_file/Mapping/conserv_map/conserv_dataset/phastcons/30way/"
PC100_path_pt1 = "code/Data_file/Mapping/conserv_map/conserv_dataset/phastcons/100way/"
PP30_path_pt1 = "code/Data_file/Mapping/conserv_map/conserv_dataset/phyloP/30way/"
PP100_path_pt1 = "code/Data_file/Mapping/conserv_map/conserv_dataset/phyloP/100way/"
mRNA_seq_path_pt1 = "code/Data_file/Mapping/mRNA_seq_map/seq_dataset_str/"

list_df = []

def load_all_pickle():

    map_ENST_Pickle_conserv_path = "code/Data_file/Mapping/conserv_map/conserv_dataset/map_ENSTtoPickle_Conserv.pkl"


    miRNA_seq_dict_path = "code/Data_file/Mapping/miRNA_seq_map/miRNA_seq_dict.pkl"

    map_ENST_Pickle_seq_path = "code/Data_file/Mapping/mRNA_seq_map/map_ENSTtoPickle_seq.pkl"

    DRSNP_path = "code/Data_file/Mapping/SNP_map/DRSNP_UTR3.pkl"


    with open(map_ENST_Pickle_conserv_path, 'rb') as file:
        map_ENST_Pickle_conserv = pickle.load(file)

    with open(miRNA_seq_dict_path, 'rb') as file:
        miRNA_seq_dict = pickle.load(file)

    with open(map_ENST_Pickle_seq_path, 'rb') as file:
        map_ENST_Pickle_seq = pickle.load(file)


    with open(DRSNP_path, 'rb') as file:
        DRSNP_dict = pickle.load(file)
        
    return map_ENST_Pickle_conserv, miRNA_seq_dict, map_ENST_Pickle_seq, DRSNP_dict


def find_6seed_index(miRNA_6seed_recom, UTR_3):
    index = []
    if miRNA_6seed_recom in UTR_3:
        for m in re.finditer(str(miRNA_6seed_recom), str(UTR_3)):
            index.append(m.start())
    if len(index) > 0:
        return index
    else:
        return 'not found'

def predict_duplex(mRNA_cdna, miRNA_seq, num_seed, site_start_cdna_temp, site_stop_cdna_temp, miRNA_length):
        #Calculate miRNA-mRNA hybrid structure
        try:
            if len(mRNA_cdna) > 30000: #Because the limitation for the argument send to subprocess is 32,768 character in windows
                limit_range = len(mRNA_cdna) - 30000 
                option_1 = '--target=' + str(mRNA_cdna[limit_range:])
            else:    
                option_1 = '--target=' + str(mRNA_cdna)
            
            option_2 = '--query=' + str(miRNA_seq)
            option_3 = '--outmode=C'

            option_5 = "--outCsvCols=start1,end1,start2,end2,hybridDPfull,hybridDBfull,E,ED1,ED2,Pu1,Pu2,E_hybrid,seedE,seedStart2,seedEnd2,seedPu1,seedPu2,subseqDP"

            #seed restriction
            option_6 = '--seedBP=' + str(num_seed)
            option_7 = '--seedQMaxUP=0' # the maximal number of unpaired bases within the query's seed region
            option_8 = '--seedTMaxUP=0' # the maximal number of unpaired bases within the target's seed region
            if num_seed == 6:
                option_9 = '--seedQRange=2-7'
            else:
                option_9 = '--seedQRange=2-8'
            option_10 = '--seedTRange=' + str(site_start_cdna_temp) +'-'+ str(site_stop_cdna_temp)

            #Enable accessibility computation using selected energy model
            #using a window technique to reduce computation workload
            #1) for query (miRNA)
                #0 sets it to the whole sequence length
            option_13 = '--qAccW=0' #sliding window length
            option_14 = '--qAccL=0' #maximal length of considered intramolecular base pairs (the maximal number of position enclosed by a base pair)
            #2) for target (mRNA)
            option_15 = '--tAccW=100'
            option_16 = '--tAccL=50'
        
            #Further restriction to improve quality of predicted structure
            option_17 = '--tIntLenMax=60' #Restrict overall length of an interaction in target mRNA
            option_18 = '--qIntLenMax=' + str(miRNA_length)
            
            acc_prof_miRNA = ""

            for i in range(0, miRNA_length, 1):
                if i == [8, 9, 10]:
                    acc_prof_miRNA = acc_prof_miRNA + "b" #nt position 9 is blocked, not use in ED computation
                else:
                    acc_prof_miRNA = acc_prof_miRNA + "x" 
        
            option_22 = '--qAccConstr='+ acc_prof_miRNA #Exluce nt 9 from interaction (being blocked in Ago2)
            option_23 = '--mode=M' #Exact mode (not heuristic)
        
            option_25 = '--windowWidth=6000' #use window-based prediction to reduce memory requirement
            option_26 = '--windowOverlap=2000'
            option_27 = '--threads=0'

            p = Popen([os.path.join(MAIN_ENV_PATH, 'bin', 'IntaRNA'), option_1, option_2, option_3,option_5, option_6, option_7, option_8, option_9, option_10,
                    option_13, option_14, option_15, option_16, option_17, option_18, option_22, option_23, option_25, option_26, option_27
                    ], stdin=PIPE, stdout=PIPE, stderr=PIPE)
            
            ans = p.communicate()   
            out = str(ans[0]).split("\n")[1].split(";")
            return out
        except:
            return "'"    

def get_candidates(row, map_ENST_Pickle_conserv, miRNA_seq_dict, map_ENST_Pickle_seq, DRSNP_dict):

    result_df = pd.DataFrame(columns=['miRNA',
                    "miRNA sequence 5'-3'",
                    'Provided ID', 
                    'Transcript ENST', 
                    'Gene ENSG',
                    'binded mRNA sequence', 
                    'miRNA_length',
                    'site_start_UTR_3', 
                    'site_stop_UTR_3',
                    'site_start_cdna', 
                    'site_stop_cdna', 
                    'site_type',
                    'Rep_dict',

                    '6-mer', 
                    '8-mer',

                    'Overall_interaction_energy',

                    'PP30_nt2_7',        
                    'PP100_nt8', 
                    'PC30_nt8', 
                    'PC100_nt2_7',

                    'Num_com_bases_nt10_12',

                    'dist_ter_3UTR', 
                    'ratio_dist_ter_3UTR', 

                    'global_3UTR_AU_content',

                    "drsnp_nt5",
                    "drsnp_nt8",

                    'miRNA_iFeature_Kmer_1', 
                    'miRNA_iFeature_DAC_4', 
                    'miRNA_iFeature_PseDNC_10', 
                    'miRNA_iFeature_PseDNC_11', 

                    'mRNA_iFeature_Kmer_3', 
                    'mRNA_iFeature_Kmer_4', 
                    'mRNA_iFeature_Kmer_6', 
                    'mRNA_iFeature_Kmer_7', 
                    'mRNA_iFeature_DAC_10'])
    miRNA = row[3]
    Provided_ID = row[1]
    mRNA_ENST = row[2]
    mRNA_ENSG = row[0]  
    
    miRNA_seq = Seq(str(Seq(str(miRNA_seq_dict.get(miRNA))).transcribe()))
    
    miRNA_seq_recom = miRNA_seq.reverse_complement()
    
    miRNA_6seed_recom = miRNA_seq_recom[-7:-1] #these ones not equal 6-mer, 7-mer, 8-mer!
    miRNA_7seed_recom = miRNA_seq_recom[-8:-1]
    
    miRNA_length = len(miRNA_seq)
    

    #2) mRNA 
        
    #Get mRNA sequence
    
    mRNA_seq_path = mRNA_seq_path_pt1 + str(map_ENST_Pickle_seq.get(mRNA_ENST)) + '.pkl'
    with open(mRNA_seq_path, 'rb') as file:
        mRNA_seq_file = pickle.load(file)
        
    mRNA_seq = mRNA_seq_file.get(mRNA_ENST)
    
    mRNA_cdna = Seq(str(Seq(str(mRNA_seq.get('cdna'))).transcribe()))
    mRNA_cds = Seq(str(Seq(str(mRNA_seq.get('cds'))).transcribe()))
    
    UTR_3 = Seq(str(mRNA_cdna)[str(mRNA_cdna).find(str(mRNA_cds)):].replace(str(mRNA_cds), '')) #the first [] was used to trim 5'UTR
    

    #find a position of binding in selected sequence of 3'UTR 
    if find_6seed_index(miRNA_6seed_recom, UTR_3) == 'not found':
        return 'no_perfect_complementary_sequence_found'
            
    for position in find_6seed_index(miRNA_6seed_recom, UTR_3):
            
        if miRNA_6seed_recom == UTR_3[position:position+6]:
            if miRNA_7seed_recom == UTR_3[position-1:position+6]:
                num_seed = 7 #2-8
            else:
                num_seed = 6 #2-7  
                    
        #create temporary mRNA_seq_binded
        mRNA_seq_binded_temp = UTR_3[position+7-miRNA_length:position+7]
            
        site_start_cdna_temp = mRNA_cdna.transcribe().find(mRNA_seq_binded_temp) + 1 #plus 1 because position should start at 1
        site_stop_cdna_temp = site_start_cdna_temp + miRNA_length - 1
    
        out = predict_duplex(mRNA_cdna, miRNA_seq, num_seed, site_start_cdna_temp, site_stop_cdna_temp, miRNA_length)


        #Data from IntaRNA (Energy in kcal/mol)

        if out == ["'"]: # if there is no predicted structure from IntaRNA in this miRNA-mRNA hybrid
            continue
                
        try: #if it can not tranform to int, there is no predicted structure
            int(out[0])
        except:
            continue
        
        start1 = int(out[0]) #start index of hybrid in seq1
        end1 = int(out[1]) #end index of hybrid in seq1
        #start2 = int(out[2]) #start index of hybrid in seq2
        #end2 = int(out[3]) #end index of hybrid in seq2
        #hybridDP = out[4] #hybrid in VRNA dot-bracket notation (interaction sites only)
        hybridDB = out[5] #hybrid in dot-bar notation (interactin sites only)
        E_overall = out[6] #overall interaction energy (ED1 + ED2 + E_hybrid)
        #ED1 = out[7] #ED value of seq1
        #ED2 = out[8] #ED value of seq2
        #Pu1 = out[9] #probability to be accessible for seq1
        #Pu2 = out[10] #probability to be accessible for seq2
        #E_hybrid = out[11] #energy of hybridization
        #E_seed = out[12] #only overall energy of the seed (including seedED)
        #seedStart2 = out[13] #start index of seed in seq2
        #seedEnd2 = out[14] #stop index of seed in seq2
        #seedPu1 = out[15] #probability of seed region to be accessible for seq1
        #seedPu2 = out[16] #probability of seed region to be accessible for seq2 
        subsetDP = out[17]
            
    
        #Get hybrid interaction (this sequence is trimmed)
        #int_mRNA = hybridDB.split("&")[0][end1+1-miRNA_length: end1+1] #5'-3'
        int_miRNA_5 = hybridDB.split("&")[1][1:] #5'-3'
        #int_miRNA_3 = int_miRNA_5[::-1] #3'-5'

        #mRNA_hyb_seq = subsetDP.split("&")[0]
        miRNA_hyb_seq_5 = subsetDP.split("&")[1]
        #miRNA_hyb_seq_3 = miRNA_hyb_seq_5[::-1]

        site_start_cdna = start1 #start at 1
        site_stop_cdna = end1 #start at 1
        #site_start_cds = mRNA_cdna.transcribe().find(UTR_3) + site_start_cdna
        #site_stop_cds = mRNA_cdna.transcribe().find(UTR_3) + site_stop_cdna
    
        #Get interaction profile (Match 1 , Mismatch 0, GU-wobble 2, Bulge 3)

        int_nt_mRNA = [] 
        int_nt_miRNA = []

        index = 0
        for i in hybridDB.split('&')[0]:
            if i == '|':
                int_nt_mRNA.append(index)
            index += 1

        int_nt_mRNA = int_nt_mRNA[::-1]

        index = 0
        for i in hybridDB.split('&')[1]:
            if i == '|':
                int_nt_miRNA.append(index)
            index += 1
        
    
        index_bind = []
        if num_seed == 6:
            wanted_list = [2, 3, 4, 5, 6, 7]
        if num_seed == 7:
            wanted_list = [2, 3, 4, 5, 6, 7, 8]
    
        for i in range(0, len(int_nt_miRNA), 1):
            if int_nt_miRNA[i] in wanted_list:
                index_bind.append(i)
            
        int_nt_mRNA_sel = [int_nt_mRNA[i] for i in index_bind]
            
        mRNA_seq_binded = mRNA_cdna[max(int_nt_mRNA_sel)+1-miRNA_length: max(int_nt_mRNA_sel)+1].transcribe()
 
        #Get true site start/stop for cdna/cds
        if not UTR_3.find(mRNA_seq_binded) == -1:  #sometimes, the temp part are located in 3'UTR, but new part are located in cds
            site_start_cdna = mRNA_cdna.transcribe().find(mRNA_seq_binded) + 1 #plus 1 because position should start at 1
            site_stop_cdna = site_start_cdna + len(mRNA_seq_binded) - 1
        else:
            continue
    
        site_start_UTR_3 = UTR_3.find(mRNA_seq_binded) + 1 #plus 1 because position should start at 1
        site_stop_UTR_3 = site_start_UTR_3 + len(mRNA_seq_binded) - 1
            
            
        #Check point if sequence exist
        if mRNA_seq_binded == Seq(''):
            continue
            
        #Feature generation
        #1) canonical site types (8-mer, 7-mer-m8, 7-mer-A1, 6-mer) & The presence of A in first position in mRNA
        if mRNA_seq_binded[-1] == 'A':
            if num_seed == 6:
                site_type = '7-mer-A1'
            elif num_seed == 7:
                site_type = '8-mer'
        else:
            if num_seed == 6:
                site_type = '6-mer'
            elif num_seed == 7:
                site_type = '7-mer-m8'


        if site_type == '6-mer':
            mer6 = 1
            mer8 = 0
        elif site_type == '8-mer':
            mer8 = 1
            mer6 = 0
        else:
            mer6 = 0
            mer8 = 0
            

        #4) interaction profile for mRNA & miRNA
        temp_dict = {}


        for pos_miRNA in range(1, miRNA_length+1, 1):
            if pos_miRNA in int_nt_miRNA: #If bind with mRNA nucleotide
                pos_mRNA_binded = int_nt_mRNA[int_nt_miRNA.index(pos_miRNA)]
                pos_mRNA_notbinded = '-'
            elif pos_miRNA == 1: #If it's the first position -> easy calculation
                pos_mRNA_binded = '-'
                pos_mRNA_notbinded = int_nt_mRNA[0]+1
            elif pos_miRNA == miRNA_length: #If it's the last one
                pos_mRNA_binded = '-'
                if temp_dict[pos_miRNA-1][1] != '-':
                    pos_mRNA_notbinded = temp_dict[pos_miRNA - 1][1] - 1
                else:
                    pos_mRNA_notbinded = temp_dict[pos_miRNA - 1][2] - 1
            
            else:
                pos_mRNA_binded = '-'
                if (pos_miRNA+1) in int_nt_miRNA:
                    if not int_nt_mRNA[int_nt_miRNA.index(pos_miRNA+1)] == '-':
                        pos_mRNA_notbinded = int_nt_mRNA[int_nt_miRNA.index(pos_miRNA+1)]+1
                else:
                    if temp_dict[pos_miRNA-1][1] != '-': #Previous have position
                        pos_mRNA_notbinded = temp_dict[pos_miRNA - 1][1] - 1
                    else:
                        pos_mRNA_notbinded = temp_dict[pos_miRNA - 1][2] - 1
            temp_dict[pos_miRNA] = (pos_miRNA, pos_mRNA_binded, pos_mRNA_notbinded)
    

        #Update with the buldge in mRNA
        temp_dict_ver2 = {}
        plus_value = 0

        for key in temp_dict.keys():
            if key == 1:
                temp_dict_ver2[key] = (temp_dict.get(key)[0], temp_dict.get(key)[1], temp_dict.get(key)[2])
                continue
        
            current_position = temp_dict.get(key)[1]
            previous_position = temp_dict.get(key-1)[1]
    
            if current_position == '-':
                current_position = temp_dict.get(key)[2]

            if previous_position == '-':
                previous_position = temp_dict.get(key-1)[2]
    
            if current_position == previous_position:
                temp_dict_ver2[key+plus_value] = (temp_dict.get(key)[0], temp_dict.get(key)[1], '-')
            elif (previous_position - current_position) > 1:
                for i in range(0, previous_position - current_position - 1, 1):
                    temp_dict_ver2[key+plus_value] = ('-', '-', previous_position - 1 - i)
                    plus_value += 1
                temp_dict_ver2[key+plus_value] = (temp_dict.get(key)[0], temp_dict.get(key)[1], temp_dict.get(key)[2])
                plus_value += 1
            else:
                temp_dict_ver2[key+plus_value] = (temp_dict.get(key)[0], temp_dict.get(key)[1], temp_dict.get(key)[2])

        Rep_dict = {}
    
        for key in temp_dict_ver2.keys():
            ntm_A, ntm_C, ntm_U, ntm_G, ntm_X = 0, 0, 0, 0, 0
            ntmi_A, ntmi_C, ntmi_U, ntmi_G, ntmi_X = 0, 0, 0, 0, 0
            ntp_0, ntp_1, ntp_2 = 0, 0, 0
            p_status = 0
    
            ntmRNA_notbind_index = temp_dict_ver2[key][2]
            ntmRNA_bind_index = temp_dict_ver2[key][1]
            ntmiRNA_index = temp_dict_ver2[key][0]
        
            if ntmRNA_bind_index != '-':
                ntmRNA = mRNA_cdna[ntmRNA_bind_index-1]
                if ntmRNA == 'A':
                    ntm_A = 1
                    p_status = 1
                elif ntmRNA == 'C':
                    ntm_C = 1
                    p_status = 1
                elif ntmRNA == 'G':
                    ntm_G = 1
                    p_status = 1
                elif ntmRNA == 'U':
                    ntm_U = 1
                    p_status = 1
            elif ntmRNA_notbind_index != '-':
                ntmRNA = mRNA_cdna[ntmRNA_notbind_index-1]
                if ntmRNA == 'A':
                    ntm_A = 1
                elif ntmRNA == 'C':
                    ntm_C = 1
                elif ntmRNA == 'G':
                    ntm_G = 1
                elif ntmRNA == 'U':
                    ntm_U = 1
            else:
                ntm_X = 1 
            
        
            if not ntmiRNA_index == '-':
                ntmiRNA = miRNA_seq[ntmiRNA_index-1]
                if ntmiRNA == 'A':
                    ntmi_A = 1
                elif ntmiRNA == 'C':
                    ntmi_C = 1
                elif ntmiRNA == 'G':
                    ntmi_G = 1
                elif ntmiRNA == 'U':
                    ntmi_U = 1
            else:
                ntmi_X = 1        
    
    
            if (((ntm_A == 1 and ntmi_U) or (ntm_U == 1 and ntmi_A == 1)) and p_status == 1):
                ntp_1 = 1
            elif (((ntm_C == 1 and ntmi_G) or (ntm_G == 1 and ntmi_C == 1)) and p_status == 1):
                ntp_1 = 1
            elif (((ntm_G == 1 and ntmi_U) or (ntm_U == 1 and ntmi_G == 1)) and p_status == 1):
                ntp_2 = 1
            else:
                ntp_0 = 1
        
            if ntmiRNA_index == 1 or ntmiRNA_index ==9:
                ntp_0 = 1
                ntp_1 = 0
                ntp_2 = 0
    
            Rep_dict[key] = [ntm_A, ntm_C, ntm_U, ntm_G, ntm_X, ntmi_A, ntmi_C, ntmi_U, ntmi_G, ntmi_X, ntp_0, ntp_1, ntp_2]
    
        #5) 3'supplementary pairing 

        num_sup_center = 0

        for i in range(10, 13, 1):
            if i in int_nt_miRNA:
                num_sup_center += 1 
        
        #6) distance from start and terminal of 3'UTR
        length_UTR_3 = len(UTR_3)
        dist_ter_UTR_3 = length_UTR_3 - site_stop_UTR_3    

        ratio_dist_ter_3UTR = float(float(dist_ter_UTR_3)/float(length_UTR_3))
    
        #Global AU content
        global_3UTR_AU_content = float((UTR_3.count("A") + UTR_3.count("U")))/len(UTR_3)
    
        float(1060/1375)

        #8)Conservation score
        chromosome = map_ENST_Pickle_conserv.get(mRNA_ENST)['chromosome']
        file_name = map_ENST_Pickle_conserv.get(mRNA_ENST)['file']

        PC30_path = PC30_path_pt1 + file_name + '.pkl'
        PC100_path = PC100_path_pt1 + file_name + '.pkl'
        PP30_path = PP30_path_pt1 + file_name + '.pkl'
        PP100_path = PP100_path_pt1 + file_name + '.pkl'
            
        with open(PC30_path, 'rb') as file:
            PC30_file = pickle.load(file).get(mRNA_ENST)

        with open(PC100_path, 'rb') as file:
            PC100_file = pickle.load(file).get(mRNA_ENST)

        with open(PP30_path, 'rb') as file:
            PP30_file = pickle.load(file).get(mRNA_ENST)

        with open(PP100_path, 'rb') as file:
            PP100_file = pickle.load(file).get(mRNA_ENST)

            
        PC30_nt8 = PC30_file.get(site_stop_UTR_3-7).get('PC30')

        PP100_nt8 = PP100_file.get(site_stop_UTR_3-7).get('PP100')

        PC100_nt2 = PC100_file.get(site_stop_UTR_3-1).get('PC100')
        PC100_nt3 = PC100_file.get(site_stop_UTR_3-2).get('PC100')
        PC100_nt4 = PC100_file.get(site_stop_UTR_3-3).get('PC100')
        PC100_nt5 = PC100_file.get(site_stop_UTR_3-4).get('PC100')
        PC100_nt6 = PC100_file.get(site_stop_UTR_3-5).get('PC100')
        PC100_nt7 = PC100_file.get(site_stop_UTR_3-6).get('PC100')

        PP30_nt2 = PP30_file.get(site_stop_UTR_3-1).get('PP30')
        PP30_nt3 = PP30_file.get(site_stop_UTR_3-2).get('PP30')
        PP30_nt4 = PP30_file.get(site_stop_UTR_3-3).get('PP30')
        PP30_nt5 = PP30_file.get(site_stop_UTR_3-4).get('PP30')
        PP30_nt6 = PP30_file.get(site_stop_UTR_3-5).get('PP30')
        PP30_nt7 = PP30_file.get(site_stop_UTR_3-6).get('PP30')

        gencoor_nt5 = PP30_file.get(site_stop_UTR_3-4).get('gencoor')
        gencoor_nt8 = PP30_file.get(site_stop_UTR_3-7).get('gencoor')

        PC100_nt2_7 = np.median([PC100_nt2, PC100_nt3, PC100_nt4, PC100_nt5, PC100_nt6, PC100_nt7])

        PP30_nt2_7 = np.median([PP30_nt2, PP30_nt3, PP30_nt4, PP30_nt5, PP30_nt6, PP30_nt7])


        #Get DRSNP (need genomic coordinate of each site (1-8)) 
        drsnp_nt5, drsnp_nt8= 0, 0

        if gencoor_nt5 in DRSNP_dict.get(chromosome):
            drsnp_nt5 = 1

        if gencoor_nt8 in DRSNP_dict.get(chromosome):
            drsnp_nt8 = 1


        result_df = result_df.append({'miRNA': miRNA,
                    "miRNA sequence 5'-3'" : miRNA_seq,
                    'Provided ID' : Provided_ID, 
                    'Transcript ENST' : mRNA_ENST, 
                    'Gene ENSG' : mRNA_ENSG,
                    'binded mRNA sequence' : mRNA_seq_binded, 
                    'miRNA_length' : miRNA_length,
                    'site_start_UTR_3': site_start_UTR_3, 
                    'site_stop_UTR_3' : site_stop_UTR_3,
                    'site_start_cdna' : site_start_cdna, 
                    'site_stop_cdna' : site_stop_cdna, 
                    'site_type' : site_type,
                    'Rep_dict' : Rep_dict,

                    '6-mer' : mer6, 
                    '8-mer' : mer8,

                    'Overall_interaction_energy' : E_overall,

                    'PP30_nt2_7' : PP30_nt2_7,        
                    'PP100_nt8': PP100_nt8, 
                    'PC30_nt8' : PC30_nt8, 
                    'PC100_nt2_7' : PC100_nt2_7,

                    'Num_com_bases_nt10_12': num_sup_center,

                    'dist_ter_3UTR' : dist_ter_UTR_3, 
                    'ratio_dist_ter_3UTR' : ratio_dist_ter_3UTR, 

                    'global_3UTR_AU_content' : global_3UTR_AU_content,

                    "drsnp_nt5": drsnp_nt5,
                    "drsnp_nt8": drsnp_nt8,

                    'miRNA_iFeature_Kmer_1': np.nan, 
                    'miRNA_iFeature_DAC_4' : np.nan, 
                    'miRNA_iFeature_PseDNC_10': np.nan, 
                    'miRNA_iFeature_PseDNC_11' : np.nan, 

                    'mRNA_iFeature_Kmer_3' : np.nan, 
                    'mRNA_iFeature_Kmer_4': np.nan, 
                    'mRNA_iFeature_Kmer_6': np.nan, 
                    'mRNA_iFeature_Kmer_7' : np.nan, 
                    'mRNA_iFeature_DAC_10' : np.nan
                    }, ignore_index=True)

    if result_df.shape[0] > 0:
        if type(result_df) == pd.core.frame.DataFrame:
            return result_df

def run_process_by_row(i, row, map_ENST_Pickle_conserv, miRNA_seq_dict, map_ENST_Pickle_seq, DRSNP_dict):
    result = get_candidates(row, map_ENST_Pickle_conserv, miRNA_seq_dict, map_ENST_Pickle_seq, DRSNP_dict)
    return(i, result)

def collect_result(result):
    global list_df
    list_df.append(result)
    

#####################FOR STEP8######################################################################################################################

def rankdata(a, method='ordinal'):
    if method not in ('average', 'min', 'max', 'dense', 'ordinal'):
        raise ValueError('unknown method "{0}"'.format(method))

    arr = np.ravel(np.asarray(a))
    algo = 'mergesort' if method == 'ordinal' else 'quicksort'
    sorter = np.argsort(arr, kind=algo)

    inv = np.empty(sorter.size, dtype=np.intp)
    inv[sorter] = np.arange(sorter.size, dtype=np.intp)

    if method == 'ordinal':
        return inv + 1

    arr = arr[sorter]
    obs = np.r_[True, arr[1:] != arr[:-1]]
    dense = obs.cumsum()[inv]

    if method == 'dense':
        return dense

    # cumulative counts of each unique value
    count = np.r_[np.nonzero(obs)[0], len(obs)]

    if method == 'max':
        return count[dense]

    if method == 'min':
        return count[dense - 1] + 1

    # average method
    return .5 * (count[dense] + count[dense - 1] + 1)



def create_TM_df(target_list, miRNA_TS_df, primiti_ts_proba, miRNA_TS_dist_df):
    #Set a default value (0 or np.nan) for columns
    target_list.loc[:,'num_site'] = 0
    target_list.loc[:,'miRNA length'] = 0
    target_list.loc[:, 'Site_1_prob'] = np.nan
    target_list.loc[:,'Site_1_pos'] = np.nan
    target_list.loc[:,'Site_1_index'] = np.nan
    target_list.loc[:,'Site_2_prob'] = np.nan
    target_list.loc[:,'Site_2_pos'] = np.nan
    target_list.loc[:,'Site_2_index'] = np.nan
    target_list.loc[:,'Site_3_prob'] = np.nan
    target_list.loc[:,'Site_3_pos'] = np.nan
    target_list.loc[:,'Site_3_index'] = np.nan
    target_list.loc[:,'Site_4_prob'] = np.nan
    target_list.loc[:,'Site_4_pos'] = np.nan
    target_list.loc[:,'Site_4_index'] = np.nan
    target_list.loc[:,'Site_5_prob'] = np.nan
    target_list.loc[:,'Site_5_pos'] = np.nan
    target_list.loc[:,'Site_5_index'] = np.nan
    target_list.loc[:,'Site_6_prob'] = np.nan
    target_list.loc[:,'Site_6_pos'] = np.nan
    target_list.loc[:,'Site_6_index'] = np.nan
    target_list.loc[:,'Site_7_prob'] = np.nan
    target_list.loc[:,'Site_7_pos'] = np.nan
    target_list.loc[:,'Site_7_index'] = np.nan
    target_list.loc[:,'Site_8_prob'] = np.nan
    target_list.loc[:,'Site_8_pos'] = np.nan
    target_list.loc[:,'Site_8_index'] = np.nan
    target_list.loc[:,'Site_9_prob'] = np.nan
    target_list.loc[:,'Site_9_pos'] = np.nan
    target_list.loc[:,'Site_9_index'] = np.nan
    target_list.loc[:,'Site_10_prob'] = np.nan
    target_list.loc[:,'Site_10_pos'] = np.nan
    target_list.loc[:,'Site_10_index'] = np.nan
    target_list.loc[:,'Site_11_prob'] = np.nan
    target_list.loc[:,'Site_11_pos'] = np.nan
    target_list.loc[:,'Site_11_index'] = np.nan
    target_list.loc[:,'Site_12_prob'] = np.nan
    target_list.loc[:,'Site_12_pos'] = np.nan
    target_list.loc[:,'Site_12_index'] = np.nan
    target_list.loc[:,'Site_13_prob'] = np.nan
    target_list.loc[:,'Site_13_pos'] = np.nan
    target_list.loc[:,'Site_13_index'] = np.nan
    target_list.loc[:,'Site_14_prob'] = np.nan
    target_list.loc[:,'Site_14_pos'] = np.nan
    target_list.loc[:,'Site_14_index'] = np.nan
    target_list.loc[:,'Site_15_prob'] = np.nan
    target_list.loc[:,'Site_15_pos'] = np.nan
    target_list.loc[:,'Site_15_index'] = np.nan


    for index, row in target_list.iterrows():
        miRNA = row[0]
        mRNA_ENST = row[2]
        site_indexes = miRNA_TS_df[(miRNA_TS_df['miRNA']==miRNA) & (miRNA_TS_df['Transcript ENST']== mRNA_ENST)].index
        
        num_site = 1
        for site_index in site_indexes:
            site_prob = primiti_ts_proba[site_index]
            site_pos = miRNA_TS_dist_df.iloc[site_index]["site_start_UTR_3"]
            col_name_prob = 'Site_' + str(num_site) + '_prob'
            col_name_pos = 'Site_' + str(num_site) + '_pos'
            col_name_index = 'Site_' + str(num_site) + '_index'
            target_list.loc[index, col_name_prob] = site_prob[1]
            target_list.loc[index, col_name_pos] = int(site_pos)
            target_list.loc[index, col_name_index] = site_index
            num_site += 1
        
        miRNA_length = miRNA_TS_df[miRNA_TS_df['miRNA']==miRNA].iloc[0]['miRNA_length']
        
        target_list.loc[index, 'miRNA length'] = miRNA_length
            
        target_list.loc[index, 'num_site'] = num_site-1

    #Get rank based on position
    cs_df = pd.DataFrame({})
    overall_num_site_target = int((target_list.shape[1]-6)/3)

    for index, row in target_list.iterrows():
        num_site = row['num_site']
        miRNA_length = row['miRNA length']
        
        default = 0
        prob_list = []
        for ind_site in range(1, overall_num_site_target+1, 1):
            site_prob = row['Site_'+str(ind_site)+'_pos'] #change _prob & _pos to change what to rank
            prob_list.append(site_prob)
            
        rank = rankdata(prob_list)
        
        rank_prob = list([])
        rank_pos = list([])
        rank_index = list([])
        
        if overall_num_site_target >= 15:
            for i in range(1, 16, 1): #get 15 site with most score
                rank_prob.append(row['Site_'+str(np.where(rank == i)[0][0]+1)+'_prob'])
                rank_pos.append(row['Site_'+str(np.where(rank == i)[0][0]+1)+'_pos'])
                rank_index.append(row['Site_'+str(np.where(rank == i)[0][0]+1)+'_index'])
        else:
            for i in range(1, 16, 1):
                try:
                    rank_prob.append(row['Site_'+str(np.where(rank == i)[0][0]+1)+'_prob'])
                    rank_pos.append(row['Site_'+str(np.where(rank == i)[0][0]+1)+'_pos'])
                    rank_index.append(row['Site_'+str(np.where(rank == i)[0][0]+1)+'_index'])
                except:
                    rank_prob.append(np.nan)
                    rank_pos.append(np.nan)
                    rank_index.append(np.nan)                
            
        #If 2 sites are overlapped, selected those with higher probability

        deleted_pos = []

        for j in range(1, 15, 1):
            if rank_pos[j] < rank_pos[j-1]+miRNA_length:
                if rank_prob[j] > rank_prob[j-1]:
                    deleted_pos.append(j-1)
                elif rank_prob[j] <= rank_prob[j-1]:
                    deleted_pos.append(j)
                
        for k in deleted_pos:
            if k == 0:
                rank_prob = rank_prob[1:] + [rank_prob[0]]
                rank_pos = rank_pos[1:] + [rank_pos[0]]
                rank_prob[-1] = np.nan
                rank_pos[-1] = np.nan
                rank_index[-1] = np.nan
            elif k != 0:
                rank_prob = rank_prob[0:k] + rank_prob[k+1:] + [rank_prob[k]]
                rank_pos = rank_pos[0:k] + rank_pos[k+1:] + [rank_pos[k]]
                rank_prob[-1] = np.nan
                rank_pos[-1] = np.nan
                rank_index[-1] = np.nan
            
        cs_df = cs_df.append({'num_site': num_site,
                                    'p1': rank_prob[0], 'ind1': rank_index[0],
                                    'p2': rank_prob[1], 'ind2': rank_index[1],
                                    'p3': rank_prob[2], 'ind3': rank_index[2],
                                    'p4': rank_prob[3], 'ind4': rank_index[3],
                                    'p5': rank_prob[4], 'ind5': rank_index[4],
                                    'p6': rank_prob[5], 'ind6': rank_index[5],
                                    'p7': rank_prob[6], 'ind7': rank_index[6],
                                    'p8': rank_prob[7], 'ind8': rank_index[7],
                                    'p9': rank_prob[8], 'ind9': rank_index[8],
                                    'p10': rank_prob[9], 'ind10': rank_index[9],
                                    'p11': rank_prob[10], 'ind11': rank_index[10],
                                    'p12': rank_prob[11],  'ind12': rank_index[11],
                                    'p13': rank_prob[12],  'ind13': rank_index[12],
                                    'p14': rank_prob[13], 'ind14': rank_index[13],
                                    'p15': rank_prob[14], 'ind15': rank_index[14]}, ignore_index=True)
        
        if (rank_pos[1]-(rank_pos[0]+miRNA_length-1)-1) < 0:
            continue
        
    cs_df = cs_df.fillna(0)
    cs_df['index'] = range(0, cs_df.shape[0], 1)
    cs_df = cs_df.set_index('index')
    cs_df['p_sum'] = 0
    cs_df['p_m1'] = 0
    cs_df['p_m1_index'] = 0
    cs_df['p_m2'] = 0
    cs_df['p_m2_index'] = 0
    cs_df['p_m3'] = 0
    cs_df['p_m3_index'] = 0
    
    #Create p_max, p_max2, p_max3
    list_p = ['p1', 'p2', 'p3', 'p4', 'p5', 'p6', 'p7', 'p8', 'p9', 'p10', 'p11', 'p12', 'p13', 'p14', 'p15']

    for index, row in cs_df.iterrows():
        p_sum = row[list_p].sum()
        list_rank_p = row[list_p].values.argsort()
        
        p_m1_index = list_rank_p[-1]
        p_m2_index = list_rank_p[-2]
        p_m3_index = list_rank_p[-3]
        
        p_m1 = row[list_p][p_m1_index]
        p_m2 = row[list_p][p_m2_index]
        p_m3 = row[list_p][p_m3_index]
            
        #modify data in the row
        cs_df.loc[index,'p_sum']= p_sum
        cs_df.loc[index,'p_m1'] = p_m1
        cs_df.loc[index,'p_m1_index'] = int(row['ind' + str(p_m1_index+1)]) 
        cs_df.loc[index,'p_m2'] = p_m2
        cs_df.loc[index,'p_m2_index'] = int(row['ind' + str(p_m2_index+1)]) 
        cs_df.loc[index,'p_m3'] = p_m3
        cs_df.loc[index,'p_m3_index'] = int(row['ind' + str(p_m3_index+1)]) 


    cs_df = cs_df.drop(list_p, axis=1)

    list_drop = ['ind1', 'ind2', 'ind3', 'ind4', 'ind5', 'ind6', 'ind7', 'ind8', 'ind9', 'ind10', 'ind11', 'ind12', 'ind13', 'ind14', 'ind15', 'p_m1_index', 'p_m2_index', 'p_m3_index']
    cs_df = cs_df.drop(list_drop, axis=1)

    #Standardization
    list_of_standardized_col = ['p_sum', 'p_m1', 'p_m2', 'p_m3']

    pre_trained_scaler_path =  "code/Data_file/Scaler/pre_trained_scaler_PRIMITI-TM.pickle"
    pre_trained_scaler = load(pre_trained_scaler_path)

    cs_df.loc[:,list_of_standardized_col] = pre_trained_scaler.transform(cs_df.loc[:,list_of_standardized_col])

    return cs_df

#####################FOR STEP10#####################################################################################################################

def make_result_table1(miRNA_TM_df, primiti_tm_pred, primiti_tm_proba, target_list, map_ENST_to_input, job_id):
    
    map_ENST_to_ENSG_path = "code/Data_file/Mapping/ID_map/ENST_to_ENSG_dict.pkl"
    
    with open(map_ENST_to_ENSG_path, 'rb') as file:
        map_ENST_to_ENSG = pickle.load(file)
    
    def map_ENST_to_Provided_ID_func(ENST):
        Provided_ID = map_ENST_to_input.get(ENST)
        return Provided_ID
    

    def map_ENST_to_ENSG_func(ENST):
        ENSG = map_ENST_to_ENSG.get(ENST)
        return ENSG
    
    def assign_YorN(x):
        if x >= 0.9:
            return 'Yes'
        else:
            return 'No'
    
    Result_table1 = pd.DataFrame(columns=["miRNA", 
        "Provided ID", 
        "Transcript ENST", 
        "Gene ENSG", 
        "Interaction",
        "Interaction Prob. Score [0-1]", 
        "No. of target sites",
        ])
    
    
    Result_table1['miRNA'] = target_list['miRNA']
    Result_table1['Transcript ENST'] = target_list['Transcript ENST']
    Result_table1['No. of target sites'] = target_list['num_site']
    Result_table1['No. of target sites'] = Result_table1['No. of target sites'].astype(int)
    Result_table1['Interaction'] = primiti_tm_proba[:,1]
    Result_table1['Interaction'] = Result_table1['Interaction'].astype(float)
    Result_table1['Interaction'] = Result_table1['Interaction'].apply(assign_YorN)

    Result_table1['Interaction Prob. Score [0-1]'] = primiti_tm_proba[:,1]
    Result_table1['Interaction Prob. Score [0-1]'] = Result_table1['Interaction Prob. Score [0-1]'].round(3)

    Result_table1['Provided ID'] = Result_table1['Transcript ENST'].apply(map_ENST_to_Provided_ID_func)
    
    Result_table1['Gene ENSG'] = Result_table1['Transcript ENST'].apply(map_ENST_to_ENSG_func)
    
    Result_table1_path = "result/" + str(job_id) + "_target_mRNA_result_table.csv"
    
    Result_table1.to_csv(Result_table1_path)

    ##For table 1 save and return table & path
    return Result_table1, Result_table1_path

###########Drawing################

A_color = '#FF1F5B'
U_color = '#009ADE'
C_color = '#AF58BA'
G_color = '#FFC61E'
Bond_color = '#2A363B'
GU_Bond_color = '#999DA0'

A1_color = '#A6761D'

font_path = "code/Data_file/Font/" + "CourierPrime-Regular.ttf"
font_bold_path = "code/Data_file/Font/" + "CourierPrime-Bold.ttf"


def draw_interaction(index, job_id, miRNA_seq, seqstruct, transcript_seq, path):

    max_len = max(len(miRNA_seq), len(seqstruct), len(transcript_seq))
    
    #Specify the filename

    filename = path + "/" + str(index) + ".png"
    

    #Start drawing a picture
    img = Image.new('RGB', (2200, 400), color=(0,0,0))

    
    fnt = ImageFont.truetype(font_path, 75)
    fnt_bold = ImageFont.truetype(font_bold_path, 75)

    d = ImageDraw.Draw(img)

    #Add miRNA and mRNA label

    start_point_x = 50
    start_point_y_miRNA = 50
    start_point_y_mRNA = 250

    color = (1, 0, 0)
    d.text((start_point_x, start_point_y_miRNA), 'miRNA', font=fnt, fill=color)

    d.text((start_point_x, start_point_y_mRNA), 'mRNA', font=fnt, fill=color)

    start_point_x = 300

    d.text((start_point_x, start_point_y_miRNA), "5'", font=fnt, fill=color)

    d.text((start_point_x, start_point_y_mRNA), "3'", font=fnt, fill=color)

    start_point_x = 400
    start_point_y = 50

    for nucleotide in miRNA_seq:
        if nucleotide == 'A':
            color = A_color
        elif nucleotide == 'U':
            color = U_color
        elif nucleotide == 'C':
            color = C_color
        elif nucleotide == 'G':
            color = G_color
        d.text((start_point_x, start_point_y), nucleotide, font=fnt, fill=color)
        start_point_x += 50

    start_point_y += 100
    start_point_x = 400
    for interaction in seqstruct:
        if interaction == ' ':
            color = '#FFFFFF'
        elif interaction == '|':
            color = Bond_color
        elif interaction == ':':
            color = GU_Bond_color

        d.text((start_point_x, start_point_y), interaction, font=fnt, fill=color)
        start_point_x += 50   

    start_point_y += 100
    start_point_x = 400
    pos = 0
    for nucleotide in transcript_seq:
        if nucleotide == 'A':
            color = A_color
            if pos == 0:
                color = A1_color
                d.text((start_point_x, start_point_y), nucleotide, font=fnt_bold, fill=color)
                start_point_x += 50
                pos += 1
                continue
        elif nucleotide == 'U':
            color = U_color
        elif nucleotide == 'C':
            color = C_color
        elif nucleotide == 'G':
            color = G_color
        d.text((start_point_x, start_point_y), nucleotide, font=fnt, fill=color)
        start_point_x += 50    
        pos += 1


    rgba = img.convert("RGBA")
    datas = rgba.getdata()

    newData = []
    for item in datas:
        if item[0] == 0 and item[1] == 0 and item[2] == 0:  # finding black colour by its RGB value
            # storing a transparent value when we find a black colour
            newData.append((255, 255, 255, 0))
        else:
            newData.append(item)  # other colours remain unchanged

    rgba.putdata(newData)
    rgba.save(filename, "PNG")



#################################
def make_result_table2(Result_table1, miRNA_TS_df, primiti_ts_pred, primiti_ts_proba, job_id):
    Result_table2 = pd.DataFrame(columns=["miRNA", 
        "Provided ID", "Transcript ENST", "Gene ENSG", "Interaction",
        "Interaction Prob. Score [0-1]", "No. of target sites", "Target site no.", 
        "Position in 3'-UTR", "Position in cDNA", "Binding", "Binding Prob. Score [0-1]", "Binding type",
        "miRNA sequence", "Bound transcript sequence", "Predicted structure by IntaRNA",
        "Col1", "Col2","Col3","Col4","Col5","Col6","Col7","Col8",  "Col9","Col10","Col11","Col12","Col13", "Col14", "Col15","Col16", 
        "Col17","Col18", "Col19","Col20", "Col21", "Col22", "Col23","Col24",  "Col25","Col26","Col27","Col28","Col29",
        "Col30","Col31", "Col32","Col33", "Col34","Col35","Col36","Col37","Col38", "Col39","Col40", "Col41",
        "Col42","Col43","Col44","Col45","Col46","Col47","Col48","Col49","Col50" ])

    path = "result/" + str(job_id) + "_graphics"
    os.mkdir(path)

    map_index_TM_to_index_TS = {}

    index_TM = 0
    for row in Result_table1.values:
        miRNA = row[0]
        Provided_ID = row[1]
        Transcript_ENST = row[2]
        Gene_ENSG = row[3]
        Interaction = row[4]
        Interaction_prob = row[5]
        Num_target_sites = row[6]
        
         
        sel_df = miRNA_TS_df[(miRNA_TS_df['miRNA']==miRNA)&(miRNA_TS_df['Transcript ENST']==Transcript_ENST)]

        target_site_num = 1

        list_index_TS = []
        for index, row in sel_df.iterrows():

            list_index_TS.append(index)

            Binding = primiti_ts_pred[index]
            if Binding == 1:
                Binding = 'Yes'
            else:
                Binding = 'No'
            Binding_prob = round(primiti_ts_proba[index, 1], 3)
            site_start_UTR_3 = row['site_start_UTR_3']
            site_stop_UTR_3 = row['site_stop_UTR_3']
            site_start_cdna = row['site_start_cdna']
            site_stop_cdna = row['site_stop_cdna']
            position_UTR_3 = str(site_start_UTR_3) + '-' + str(site_stop_UTR_3)
            position_cdna = str(site_start_cdna) + '-' + str(site_stop_cdna)
            Binding_type = row['site_type']

            #Generating a sequence
            Rep_dict = row['Rep_dict']
            miRNA_seq = row["miRNA sequence 5'-3'"]
            binded_mRNA_seq = row['binded mRNA sequence']

            Row1_dict = {}
            Row2_dict = {}
            Row3_dict = {}

            Row1_string = ''
            Row2_string = ''
            Row3_string = ''

            for i in range(0, 50, 1):
                Row1_dict[i] = ''
                Row2_dict[i] = ''
                Row3_dict[i] = ''
            
            num = 0
            for key in Rep_dict.keys():
                [ntm_A, ntm_C, ntm_U, ntm_G, ntm_X, ntmi_A, ntmi_C, ntmi_U, ntmi_G, ntmi_X, ntp_0, ntp_1, ntp_2] = Rep_dict[key]

                #miRNA
                if ntmi_A == 1:
                    Row1_dict[num] = 'A'
                    Row1_string += 'A'
                elif ntmi_C == 1:
                    Row1_dict[num] = 'C'
                    Row1_string += 'C'
                elif ntmi_U == 1:
                    Row1_dict[num] = 'U'
                    Row1_string += 'U'
                elif ntmi_G == 1:
                    Row1_dict[num] = 'G'
                    Row1_string += 'G'
                elif ntmi_X == 1:
                    Row1_dict[num] = ' '
                    Row1_string += ' '

                #Interaction
                if key == 0: #Force that the binding in the first position will not happen :) (Biologically make sense)
                    Row2_dict[num] = ' '
                    Row2_string += ' '
#             
                else:
                    if ntp_0 == 1:
                        Row2_dict[num] = ' '
                        Row2_string += ' '
                    elif ntp_1 == 1:
                        Row2_dict[num] = '|'
                        Row2_string += '|'
                    elif ntp_2 == 1:
                        Row2_dict[num] = ':'
                        Row2_string += ':'
                    else:
                        Row2_dict[num] = ' '
                        Row2_string += ' '

                #mRNA
                if ntm_A == 1:
                    Row3_dict[num] = 'A'
                    Row3_string += 'A'
                elif ntm_C == 1:
                    Row3_dict[num] = 'C'
                    Row3_string += 'C'
                elif ntm_U == 1:
                    Row3_dict[num] = 'U'
                    Row3_string += 'U'
                elif ntm_G == 1:
                    Row3_dict[num] = 'G'
                    Row3_string += 'G'
                elif ntm_X == 1:
                    Row3_dict[num] = ' '
                    Row3_string += ' '

                num += 1


            int_table = pd.DataFrame()

            #Draw an interaction :D
            draw_interaction(index, job_id, Row1_string, Row2_string, Row3_string, path)

            #Row1
            Result_table2 = Result_table2.append({ "miRNA" : miRNA, 
                                                    "Provided ID" : Provided_ID, 
                                                    "Transcript ENST" : Transcript_ENST, 
                                                    "Gene ENSG": Gene_ENSG, 
                                                    "Interaction" : Interaction,
                                                    "Interaction Prob. Score [0-1]" : Interaction_prob, 
                                                    "No. of target sites" : Num_target_sites, 
                                                    "Target site no." : target_site_num, 
                                                    "Position in 3'-UTR" : position_UTR_3,
                                                    "Position in cDNA" : position_cdna,
                                                    "Binding" : Binding,
                                                    "Binding Prob. Score [0-1]" : Binding_prob, 
                                                    "Binding type" : Binding_type, 
                                                    "miRNA sequence" : miRNA_seq,
                                                    "Bound transcript sequence" : binded_mRNA_seq,
                                                    "Col1": Row1_dict[0],
                                                    "Col2": Row1_dict[1],
                                                    "Col3": Row1_dict[2],
                                                    "Col4": Row1_dict[3],
                                                    "Col5": Row1_dict[4],
                                                    "Col6": Row1_dict[5],
                                                    "Col7": Row1_dict[6],
                                                    "Col8": Row1_dict[7],                                                   
                                                    "Col9": Row1_dict[8],
                                                    "Col10": Row1_dict[9],
                                                    "Col11": Row1_dict[10],
                                                    "Col12": Row1_dict[11],
                                                    "Col13": Row1_dict[12],
                                                    "Col14": Row1_dict[13],
                                                    "Col15": Row1_dict[14],
                                                    "Col16": Row1_dict[15], 
                                                    "Col17": Row1_dict[16],
                                                    "Col18": Row1_dict[17],
                                                    "Col19": Row1_dict[18],
                                                    "Col20": Row1_dict[19],
                                                    "Col21": Row1_dict[20],
                                                    "Col22": Row1_dict[21],
                                                    "Col23": Row1_dict[22],
                                                    "Col24": Row1_dict[23],                                                    
                                                    "Col25": Row1_dict[24],
                                                    "Col26": Row1_dict[25],
                                                    "Col27": Row1_dict[26],
                                                    "Col28": Row1_dict[27],
                                                    "Col29": Row1_dict[28],
                                                    "Col30": Row1_dict[29],
                                                    "Col31": Row1_dict[30],
                                                    "Col32": Row1_dict[31],                                                  
                                                    "Col33": Row1_dict[32],
                                                    "Col34": Row1_dict[33],
                                                    "Col35": Row1_dict[34],
                                                    "Col36": Row1_dict[35],
                                                    "Col37": Row1_dict[36],
                                                    "Col38": Row1_dict[37],
                                                    "Col39": Row1_dict[38],
                                                    "Col40": Row1_dict[39], 
                                                    "Col41": Row1_dict[40],
                                                    "Col42": Row1_dict[41],
                                                    "Col43": Row1_dict[42],
                                                    "Col44": Row1_dict[43],
                                                    "Col45": Row1_dict[44],
                                                    "Col46": Row1_dict[45],
                                                    "Col47": Row1_dict[46],
                                                    "Col48": Row1_dict[47],
                                                    "Col49": Row1_dict[48],
                                                    "Col50": Row1_dict[49]                                              
                                                    },
                                                    ignore_index=True)
            #Row2
            Result_table2 = Result_table2.append({ "miRNA" : '', 
                                                    "Provided ID" : '', 
                                                    "Transcript ENST" : '', 
                                                    "Gene ENSG": '', 
                                                    "Interaction" : '',
                                                    "Interaction Prob. Score [0-1]" : '', 
                                                    "No. of target sites" : '', 
                                                    "Target site no." : '', 
                                                    "Position in 3'-UTR" : '',
                                                    "Position in cDNA" : '',
                                                    "Binding" : '',
                                                    "Binding Prob. Score [0-1]" : '', 
                                                    "Binding type" : '', 
                                                    "miRNA sequence" : '',
                                                    "Bound transcript sequence" : '',
                                                    "Col1": Row2_dict[0],
                                                    "Col2": Row2_dict[1],
                                                    "Col3": Row2_dict[2],
                                                    "Col4": Row2_dict[3],
                                                    "Col5": Row2_dict[4],
                                                    "Col6": Row2_dict[5],
                                                    "Col7": Row2_dict[6],
                                                    "Col8": Row2_dict[7],                                                   
                                                    "Col9": Row2_dict[8],
                                                    "Col10": Row2_dict[9],
                                                    "Col11": Row2_dict[10],
                                                    "Col12": Row2_dict[11],
                                                    "Col13": Row2_dict[12],
                                                    "Col14": Row2_dict[13],
                                                    "Col15": Row2_dict[14],
                                                    "Col16": Row2_dict[15], 
                                                    "Col17": Row2_dict[16],
                                                    "Col18": Row2_dict[17],
                                                    "Col19": Row2_dict[18],
                                                    "Col20": Row2_dict[19],
                                                    "Col21": Row2_dict[20],
                                                    "Col22": Row2_dict[21],
                                                    "Col23": Row2_dict[22],
                                                    "Col24": Row2_dict[23],                                                    
                                                    "Col25": Row2_dict[24],
                                                    "Col26": Row2_dict[25],
                                                    "Col27": Row2_dict[26],
                                                    "Col28": Row2_dict[27],
                                                    "Col29": Row2_dict[28],
                                                    "Col30": Row2_dict[29],
                                                    "Col31": Row2_dict[30],
                                                    "Col32": Row2_dict[31],                                                  
                                                    "Col33": Row2_dict[32],
                                                    "Col34": Row2_dict[33],
                                                    "Col35": Row2_dict[34],
                                                    "Col36": Row2_dict[35],
                                                    "Col37": Row2_dict[36],
                                                    "Col38": Row2_dict[37],
                                                    "Col39": Row2_dict[38],
                                                    "Col40": Row2_dict[39], 
                                                    "Col41": Row2_dict[40],
                                                    "Col42": Row2_dict[41],
                                                    "Col43": Row2_dict[42],
                                                    "Col44": Row2_dict[43],
                                                    "Col45": Row2_dict[44],
                                                    "Col46": Row2_dict[45],
                                                    "Col47": Row2_dict[46],
                                                    "Col48": Row2_dict[47],
                                                    "Col49": Row2_dict[48],
                                                    "Col50": Row2_dict[49]  
                                                    },
                                                    ignore_index=True)
            #Row3
            Result_table2 = Result_table2.append({ "miRNA" : '', 
                                                    "Provided ID" : '', 
                                                    "Transcript ENST" : '', 
                                                    "Gene ENSG": '', 
                                                    "Interaction" : '',
                                                    "Interaction Prob. Score [0-1]" : '', 
                                                    "No. of target sites" : '', 
                                                    "Target site no." : '', 
                                                    "Position in 3'-UTR" : '',
                                                    "Position in cDNA" : '',
                                                    "Binding" : '',
                                                    "Binding Prob. Score [0-1]" : '', 
                                                    "Binding type" : '', 
                                                    "miRNA sequence" : '',
                                                    "Bound transcript sequence" : '',
                                                    "Col1": Row3_dict[0],
                                                    "Col2": Row3_dict[1],
                                                    "Col3": Row3_dict[2],
                                                    "Col4": Row3_dict[3],
                                                    "Col5": Row3_dict[4],
                                                    "Col6": Row3_dict[5],
                                                    "Col7": Row3_dict[6],
                                                    "Col8": Row3_dict[7],                                                   
                                                    "Col9": Row3_dict[8],
                                                    "Col10": Row3_dict[9],
                                                    "Col11": Row3_dict[10],
                                                    "Col12": Row3_dict[11],
                                                    "Col13": Row3_dict[12],
                                                    "Col14": Row3_dict[13],
                                                    "Col15": Row3_dict[14],
                                                    "Col16": Row3_dict[15], 
                                                    "Col17": Row3_dict[16],
                                                    "Col18": Row3_dict[17],
                                                    "Col19": Row3_dict[18],
                                                    "Col20": Row3_dict[19],
                                                    "Col21": Row3_dict[20],
                                                    "Col22": Row3_dict[21],
                                                    "Col23": Row3_dict[22],
                                                    "Col24": Row3_dict[23],                                                    
                                                    "Col25": Row3_dict[24],
                                                    "Col26": Row3_dict[25],
                                                    "Col27": Row3_dict[26],
                                                    "Col28": Row3_dict[27],
                                                    "Col29": Row3_dict[28],
                                                    "Col30": Row3_dict[29],
                                                    "Col31": Row3_dict[30],
                                                    "Col32": Row3_dict[31],                                                  
                                                    "Col33": Row3_dict[32],
                                                    "Col34": Row3_dict[33],
                                                    "Col35": Row3_dict[34],
                                                    "Col36": Row3_dict[35],
                                                    "Col37": Row3_dict[36],
                                                    "Col38": Row3_dict[37],
                                                    "Col39": Row3_dict[38],
                                                    "Col40": Row3_dict[39], 
                                                    "Col41": Row3_dict[40],
                                                    "Col42": Row3_dict[41],
                                                    "Col43": Row3_dict[42],
                                                    "Col44": Row3_dict[43],
                                                    "Col45": Row3_dict[44],
                                                    "Col46": Row3_dict[45],
                                                    "Col47": Row3_dict[46],
                                                    "Col48": Row3_dict[47],
                                                    "Col49": Row3_dict[48],
                                                    "Col50": Row3_dict[49]   
                                                    },
                                                    ignore_index=True)
            target_site_num += 1
        
        map_index_TM_to_index_TS[index_TM] = list_index_TS

        index_TM += 1

    Result_table2_path = "result/" + str(job_id) + "_target_site_result_table.csv"
    ##For table 2 only return a save path :)
    Result_table2, Result_table2.to_csv(Result_table2_path, index=False)


    return Result_table2, Result_table2_path, map_index_TM_to_index_TS





































###########################Step: 1 - 10#############################################################################################################################


def step1_primiti(miRNA_list, transcript_list, supported_miRNA_list, supported_transcript_list):
    #Get miRNAs and transcripts from an input file
    miRNA_list, transcript_list, map_ENST_to_input, error_miRNAs, error_transcripts_type1, error_transcripts_type2, error_transcripts_type3 = get_list_miRNAs_transcripts(miRNA_list, transcript_list, supported_miRNA_list, supported_transcript_list)

    error_miRNAs = np.unique(error_miRNAs)
    error_transcripts_type1 = np.unique(error_transcripts_type1)
    error_transcripts_type2 = np.unique(error_transcripts_type2)    
    error_transcripts_type3 = np.unique(error_transcripts_type3)



    if miRNA_list.shape[0] == 0:
        miRNA_list = 'none'

    if transcript_list.shape[0] == 0:
        transcript_list = 'none'

    return miRNA_list, transcript_list, map_ENST_to_input, error_miRNAs, error_transcripts_type1, error_transcripts_type2, error_transcripts_type3



def step2_primiti(miRNA_list, transcript_list, map_ENST_Pickle_conserv, miRNA_seq_dict, map_ENST_Pickle_seq, DRSNP_dict):
    #only generate the features we need in the final model :)

    list_columns = ['miRNA','Provided ID', 'Transcript ENST', 'Gene ENSG',]

    #create a pair list from two list
    pair_list_df = pd.DataFrame(columns=list_columns)

    for miRNA in miRNA_list:
        i = 0
        for transcript in transcript_list:
            original_input = np.nan
            d = {'miRNA' : [miRNA], 'Provided ID': [original_input], 'Transcript ENST' : [transcript]}
            df = pd.DataFrame(data=d)
            pair_list_df = pair_list_df.append(df, ignore_index=True)
            i += 1 

    #Get candidate canonical target sites

    #cores = mp.cpu_count()-7 ######JUST FOR TEST IN OUR PC :) can increase more !########
    cores = 4
    pool = mp.Pool(cores)
    
    global list_df
    list_df = []

    for i, row in enumerate(pair_list_df.values):
        pool.apply_async(run_process_by_row, args=(i, row, map_ENST_Pickle_conserv, miRNA_seq_dict, map_ENST_Pickle_seq, DRSNP_dict), callback=collect_result)
    
    pool.close()
    pool.join()

    list_df.sort(key=lambda x: x[0])
    list_df = [r for i, r in list_df]

    new_list_df = []
    for i in list_df:
        if type(i) == pd.core.frame.DataFrame:
            new_list_df.append(i)

    miRNA_TS_df = pd.concat(new_list_df)

    miRNA_TS_df = miRNA_TS_df.drop_duplicates(subset=['miRNA', 'Transcript ENST', 'site_start_UTR_3', 'site_stop_UTR_3'])

    miRNA_TS_df['index'] = range(0, miRNA_TS_df.shape[0], 1)
    miRNA_TS_df.set_index('index', inplace=True)
    
    return miRNA_TS_df



def step3_primiti(miRNA_TS_df, job_id):
    found_miRNA_list = np.unique(miRNA_TS_df['miRNA'].values) #Different from miRNA_list as some miRNAs may not found any interaction site with any transcript :(
    
    #Create miRNA name to sequence dict
    miRNA_name_to_seq_dict = {val1: val2["miRNA sequence 5'-3'"].values[0] for val1,val2 in miRNA_TS_df.groupby(["miRNA"])}
    miRNA_to_miRNA_index_dict = {}
    index = 0
    for miRNA in found_miRNA_list:
        miRNA_to_miRNA_index_dict[miRNA] = index
        index += 1

    #Create index to mRNA binded sequence dict
    index_to_seq_dict = miRNA_TS_df["binded mRNA sequence"].to_dict()

    #Create special fasta format for miRNA

    save_miRNA_ilearn_format_path = 'result/temp/iLearn/miRNA_ilearn_format_input_' + str(job_id) + '.txt'

    file_handle = open(save_miRNA_ilearn_format_path,"w") 

    for name in found_miRNA_list:
    #for name in miRNA_name_to_seq_dict.keys(): #####################Check for this one !!! It's a serious issue :(
        seq = miRNA_name_to_seq_dict.get(name)
        first_line = '>'+ name + '|1|training'
        second_line = seq

        file_handle.write(str(first_line+'\n')) 
        file_handle.write(str(second_line+'\n')) 
        
    file_handle.close()

    #Create special fasta format for mRNA binded sequence
    save_mRNA_ilearn_format_path = 'result/temp/iLearn/mRNA_ilearn_format_input' + str(job_id) + '.txt'

    file_handle = open(save_mRNA_ilearn_format_path,"w") 

    for index in index_to_seq_dict.keys():
        seq = index_to_seq_dict.get(int(index))
        first_line = '>'+ str(index) + '|1|training'
        second_line = seq
        file_handle.write(str(first_line+'\n')) 
        file_handle.write(str(second_line+'\n')) 
        
    file_handle.close()

    return save_miRNA_ilearn_format_path, save_mRNA_ilearn_format_path, miRNA_to_miRNA_index_dict


def step4_primiti(miRNA_ilearn_input, mRNA_ilearn_input, job_id):

    mRNA_output_basic = 'result/temp/iLearn/mRNA_ilearn_format_output_basic_' + str(job_id) + '.txt'
    mRNA_output_DACC = 'result/temp/iLearn/mRNA_ilearn_format_output_ACC_' + str(job_id) + '.txt'
    miRNA_output_basic = 'result/temp/iLearn/miRNA_ilearn_format_output_basic_' + str(job_id) + '.txt'
    miRNA_output_DACC = 'result/temp/iLearn/miRNA_ilearn_format_output_ACC_' + str(job_id) + '.txt'
    miRNA_output_Pse = 'result/temp/iLearn/miRNA_ilearn_format_output_Pse_' + str(job_id) + '.txt'

    cmd = SUB_ENV_PATH + "/bin/python3 " + 'code/iLearn/iLearn-nucleotide-basic.py' + " --file " + mRNA_ilearn_input + ' --method Kmer ' + '--format svm ' + '--out ' + mRNA_output_basic
    os.system(cmd)

    cmd = SUB_ENV_PATH + "/bin/python3 " + 'code/iLearn/iLearn-nucleotide-Pse.py' + " --file " + mRNA_ilearn_input + ' --method PseDNC ' + ' --type RNA ' + ' --format svm ' + '--out ' + mRNA_output_DACC
    os.system(cmd)

    cmd = SUB_ENV_PATH + "/bin/python3 " + 'code/iLearn/iLearn-nucleotide-basic.py' + " --file " + miRNA_ilearn_input + ' --method Kmer ' + '--format svm ' + '--out ' + miRNA_output_basic
    os.system(cmd)

    cmd = SUB_ENV_PATH + "/bin/python3 " + 'code/iLearn/iLearn-nucleotide-acc.py' + " --file " + miRNA_ilearn_input + ' --method DAC ' + ' --type RNA ' + '--format svm ' + '--out ' + miRNA_output_DACC
    os.system(cmd)

    cmd = SUB_ENV_PATH + "/bin/python3 " + 'code/iLearn/iLearn-nucleotide-Pse.py' + " --file " + miRNA_ilearn_input + ' --method PseDNC ' + ' --type RNA ' + '--format svm ' + '--out ' + miRNA_output_Pse
    os.system(cmd)

    return mRNA_output_basic, mRNA_output_DACC, miRNA_output_basic, miRNA_output_DACC, miRNA_output_Pse


def step5_primiti(miRNA_TS_df, mRNA_output_basic, mRNA_output_ACC, miRNA_output_basic, miRNA_output_ACC, miRNA_output_Pse, miRNA_to_miRNA_index_dict):

    def get_val(x):
        if ':' in str(x):
            val = x.split(':')[1]
            return val
        else:
            return x

    basic_miRNA_ilearn = pd.read_csv(miRNA_output_basic, sep='  ', engine='python', header=None)
    basic_mRNA_ilearn = pd.read_csv(mRNA_output_basic, sep='  ', engine='python', header=None)

    for key in basic_miRNA_ilearn.keys():
        basic_miRNA_ilearn[key] = basic_miRNA_ilearn[key].apply(get_val)
        
    for key in basic_mRNA_ilearn.keys():
        basic_mRNA_ilearn[key] = basic_mRNA_ilearn[key].apply(get_val)

    acc_miRNA_ilearn = pd.read_csv(miRNA_output_ACC, sep='  ', engine='python', header=None)
    acc_mRNA_ilearn = pd.read_csv(mRNA_output_ACC, sep=' ', engine='python', header=None)

    for key in acc_miRNA_ilearn.keys():
        acc_miRNA_ilearn[key] = acc_miRNA_ilearn[key].apply(get_val)

    for key in acc_mRNA_ilearn.keys():
        acc_mRNA_ilearn[key] = acc_mRNA_ilearn[key].apply(get_val)
        
    pse_miRNA_ilearn = pd.read_csv(miRNA_output_Pse, sep='  ', engine='python', header=None)

    for key in pse_miRNA_ilearn.keys():
        pse_miRNA_ilearn[key] = pse_miRNA_ilearn[key].apply(get_val)


    basic_miRNA_ilearn_dict = basic_miRNA_ilearn.to_dict() #16 feat
    basic_mRNA_ilearn_dict = basic_mRNA_ilearn.to_dict()

    acc_miRNA_ilearn_dict = acc_miRNA_ilearn.to_dict() #12 feat
    acc_mRNA_ilearn_dict = acc_mRNA_ilearn.to_dict()

    pse_miRNA_ilearn_dict = pse_miRNA_ilearn.to_dict() #18 feat

    def get_ilearn_feat_miRNA(x):
        miRNA_index = miRNA_to_miRNA_index_dict.get(x)
        feat = data_dict.get(feat_num).get(miRNA_index)
        return feat

    def get_ilearn_feat_mRNA(x):
        feat = data_dict.get(feat_num).get(x)
        return feat

    #Basic Kmer iLearn miRNA
    feature_name_list = ['miRNA_iFeature_Kmer_1', 'miRNA_iFeature_Kmer_2', 'miRNA_iFeature_Kmer_3', 'miRNA_iFeature_Kmer_4', 'miRNA_iFeature_Kmer_5', 'miRNA_iFeature_Kmer_6', 
                'miRNA_iFeature_Kmer_7', 'miRNA_iFeature_Kmer_8', 'miRNA_iFeature_Kmer_9', 'miRNA_iFeature_Kmer_10', 'miRNA_iFeature_Kmer_11', 'miRNA_iFeature_Kmer_12', 
                'miRNA_iFeature_Kmer_13', 'miRNA_iFeature_Kmer_14', 'miRNA_iFeature_Kmer_15', 'miRNA_iFeature_Kmer_16']
    data_dict = basic_miRNA_ilearn_dict

    feat_num = 1
    for feature_name in feature_name_list:
        if feat_num in [1]:
            miRNA_TS_df[feature_name] = miRNA_TS_df['miRNA'].apply(get_ilearn_feat_miRNA)
        feat_num += 1

    #ACC DAC iLearn miRNA
    feature_name_list = ['miRNA_iFeature_DAC_1', 'miRNA_iFeature_DAC_2', 'miRNA_iFeature_DAC_3', 'miRNA_iFeature_DAC_4', 'miRNA_iFeature_DAC_5', 'miRNA_iFeature_DAC_6', 
                'miRNA_iFeature_DAC_7', 'miRNA_iFeature_DAC_8', 'miRNA_iFeature_DAC_9', 'miRNA_iFeature_DAC_10', 'miRNA_iFeature_DAC_11', 'miRNA_iFeature_DAC_12']
    data_dict = acc_miRNA_ilearn_dict
    feat_num = 1
    for feature_name in feature_name_list:
        if feat_num in [4]:
            miRNA_TS_df[feature_name] = miRNA_TS_df['miRNA'].apply(get_ilearn_feat_miRNA)
        feat_num += 1

    #Pse Pse iLearn miRNA
    feature_name_list = ['miRNA_iFeature_PseDNC_1', 'miRNA_iFeature_PseDNC_2', 'miRNA_iFeature_PseDNC_3', 'miRNA_iFeature_PseDNC_4', 'miRNA_iFeature_PseDNC_5', 'miRNA_iFeature_PseDNC_6', 
                'miRNA_iFeature_PseDNC_7', 'miRNA_iFeature_PseDNC_8', 'miRNA_iFeature_PseDNC_9', 'miRNA_iFeature_PseDNC_10', 'miRNA_iFeature_PseDNC_11', 'miRNA_iFeature_PseDNC_12', 
                'miRNA_iFeature_PseDNC_13', 'miRNA_iFeature_PseDNC_14', 'miRNA_iFeature_PseDNC_15', 'miRNA_iFeature_PseDNC_16', 'miRNA_iFeature_PseDNC_17', 'miRNA_iFeature_PseDNC_18']
    data_dict = pse_miRNA_ilearn_dict
    feat_num = 1
    for feature_name in feature_name_list:
        if feat_num in [10, 11]:
            miRNA_TS_df[feature_name] = miRNA_TS_df['miRNA'].apply(get_ilearn_feat_miRNA)
        feat_num += 1


    #######################################
    index_series = pd.Series(miRNA_TS_df.index.map(int))

    #Basic Kmer iLearn mRNA
    
    feature_name_list = ['mRNA_iFeature_Kmer_1', 'mRNA_iFeature_Kmer_2', 'mRNA_iFeature_Kmer_3', 'mRNA_iFeature_Kmer_4', 'mRNA_iFeature_Kmer_5', 'mRNA_iFeature_Kmer_6', 
                    'mRNA_iFeature_Kmer_7', 'mRNA_iFeature_Kmer_8', 'mRNA_iFeature_Kmer_9', 'mRNA_iFeature_Kmer_10', 'mRNA_iFeature_Kmer_11', 'mRNA_iFeature_Kmer_12', 
                    'mRNA_iFeature_Kmer_', 'mRNA_iFeature_Kmer_14', 'mRNA_iFeature_Kmer_15', 'mRNA_iFeature_Kmer_16']
    data_dict = basic_mRNA_ilearn_dict
    feat_num = 1
    for feature_name in feature_name_list:
        if feat_num in [3, 4, 6, 7]:
            miRNA_TS_df[feature_name] = index_series.apply(get_ilearn_feat_mRNA)
        feat_num += 1
        

    #ACC DAC iLearn mRNA
    feature_name_list = ['mRNA_iFeature_DAC_1', 'mRNA_iFeature_DAC_2', 'mRNA_iFeature_DAC_3', 'mRNA_iFeature_DAC_4', 'mRNA_iFeature_DAC_5', 'mRNA_iFeature_DAC_6', 
                    'mRNA_iFeature_DAC_7', 'mRNA_iFeature_DAC_8', 'mRNA_iFeature_DAC_9', 'mRNA_iFeature_DAC_10', 'mRNA_iFeature_DAC_11', 'mRNA_iFeature_DAC_12']
    data_dict = acc_mRNA_ilearn_dict
    feat_num = 1
    for feature_name in feature_name_list:
        if feat_num in [10]:
            miRNA_TS_df[feature_name] = index_series.apply(get_ilearn_feat_mRNA)
        feat_num += 1

    return miRNA_TS_df


def step6_primiti(miRNA_TS_df):

    list_of_standarized_col = ['6-mer' , 
                    'miRNA_iFeature_DAC_4', 
                    'Overall_interaction_energy',
                    'dist_ter_3UTR', 
                    'PC100_nt2_7', 
                    'miRNA_iFeature_PseDNC_11', 
                    'ratio_dist_ter_3UTR', 
                    'mRNA_iFeature_Kmer_3', 
                    'PP30_nt2_7', 
                    'miRNA_iFeature_PseDNC_10', 
                    'mRNA_iFeature_Kmer_6', 
                    '8-mer',
                    'mRNA_iFeature_Kmer_4', 
                    'global_3UTR_AU_content', 
                    'miRNA_iFeature_Kmer_1', 
                    'PP100_nt8', 
                    'PC30_nt8', 
                    "drsnp_nt5",
                    "drsnp_nt8",
                    'mRNA_iFeature_Kmer_7', 
                    'mRNA_iFeature_DAC_10',
                    'Num_com_bases_nt10_12']

    #Save location in 3'-UTR before trained for showing later 
    
    sel_col = ["site_start_UTR_3", "site_stop_UTR_3", "site_start_cdna", "site_stop_cdna", "miRNA_length"]
    miRNA_TS_dist_df = miRNA_TS_df.loc[:,sel_col]

    #Load pre-trained scaler
    pre_trained_scaler_path = "code/Data_file/Scaler/" + "pre_trained_scaler_PRIMITI-TS.pickle"
    pre_trained_scaler = load(pre_trained_scaler_path)

    miRNA_TS_df.loc[:,list_of_standarized_col] = pre_trained_scaler.transform(miRNA_TS_df.loc[:,list_of_standarized_col])

    return miRNA_TS_df, miRNA_TS_dist_df
    

def step7_primiti(miRNA_TS_df):

    list_feat_col = ['6-mer' , 
                    'miRNA_iFeature_DAC_4', 
                    'Overall_interaction_energy',
                    'dist_ter_3UTR', 
                    'PC100_nt2_7', 
                    'miRNA_iFeature_PseDNC_11', 
                    'ratio_dist_ter_3UTR', 
                    'mRNA_iFeature_Kmer_3', 
                    'PP30_nt2_7', 
                    'miRNA_iFeature_PseDNC_10', 
                    'mRNA_iFeature_Kmer_6', 
                    '8-mer',
                    'mRNA_iFeature_Kmer_4', 
                    'global_3UTR_AU_content', 
                    'miRNA_iFeature_Kmer_1', 
                    'PP100_nt8', 
                    'PC30_nt8', 
                    "drsnp_nt5",
                    "drsnp_nt8",
                    'mRNA_iFeature_Kmer_7', 
                    'mRNA_iFeature_DAC_10',
                    'Num_com_bases_nt10_12']

    #Load pre-trained model
    model_path = "code/models/" + "pre_trained_PRIMITI_TS_XGBoost_py2.sav"

    model = pickle.load(open(model_path, 'rb'))

    X = miRNA_TS_df.loc[:, list_feat_col]

    primiti_ts_pred = model.predict(X)
    primiti_ts_proba = model.predict_proba(X)

    return primiti_ts_pred, primiti_ts_proba


def step8_primiti(miRNA_TS_df, primiti_ts_proba, miRNA_TS_dist_df):
    sel_col = ['miRNA', 'Provided ID', 'Transcript ENST', 'Gene ENSG']
    target_list = miRNA_TS_df[sel_col].copy()
    target_list = target_list.drop_duplicates()

    miRNA_TM_df = create_TM_df(target_list, miRNA_TS_df, primiti_ts_proba, miRNA_TS_dist_df)

    return miRNA_TM_df, target_list


def step9_primiti(miRNA_TM_df):
    list_feat_col = ['p_sum', 'p_m1', 'p_m2', 'p_m3']

    #Load pre-trained model
    model_path = "code/models/pre_trained_PRIMITI_TM_XGBoost_py2.sav"
    model = pickle.load(open(model_path, 'rb'))

    X = miRNA_TM_df.loc[:, list_feat_col]

    primiti_tm_pred = model.predict(X)
    primiti_tm_proba = model.predict_proba(X)

    return primiti_tm_pred, primiti_tm_proba 


def step10_primiti(miRNA_TS_df, primiti_ts_pred, primiti_ts_proba, miRNA_TM_df, primiti_tm_pred, primiti_tm_proba, target_list, map_ENST_to_input, job_id):
    
    Result_table1, Result_table1_path = make_result_table1(miRNA_TM_df, primiti_tm_pred, primiti_tm_proba, target_list, map_ENST_to_input, job_id)

    Result_table2, Result_table2_path, map_index_TM_to_index_TS = make_result_table2(Result_table1, miRNA_TS_df, primiti_ts_pred, primiti_ts_proba, job_id)

    return Result_table1, Result_table1_path, Result_table2, Result_table2_path, map_index_TM_to_index_TS


##############################################################################################
def get_list(input_csv):
    input_list = pd.read_csv(input_csv, quotechar='\"', delimiter=',', header=0)
    input_list = input_list.values

    new_input_list = []
    for input in input_list:
        new_input_list.append(input[0])
    return new_input_list


if __name__ == '__main__':

    print("**********Started running JOB***********")
    print(datetime.datetime.now())
    print("****************************************")
    parser = argparse.ArgumentParser(description='ex) python PRIMITI.py input_miRNA.csv input_transcript.csv output.csv')
    parser.add_argument("input_miRNA_csv", help="Choose csv file contain a list of miRNAs")
    parser.add_argument("input_transcript_csv", help="Choose csv file contain a list of transcripts")
    parser.add_argument("job_id", help="Select job_id for outputs - the results with job id can be found in the result folder. They include 1) miRNA-target mRNA interactions table 2) miRNA-target site binding table 3) graphical representation for binding interaction as a subfolder 4) if there is an error(s), the report will be found in errors folder")
    
    args = parser.parse_args()
    input_miRNA_csv = args.input_miRNA_csv
    input_transcript_csv = args.input_transcript_csv
    job_id = args.job_id

    MAIN_ENV_PATH = MainConfig.MAIN_ENV_PATH
    SUB_ENV_PATH = MainConfig.SUB_ENV_PATH

    print('job_id = ', job_id)

    #PRIMITI-TS 1. Prepare the list for miRNA-transcript pairs ------------------------------------------------------------------------------------##
    try:
        print("Retrieving miRNA and disease lists from files\n")

        miRNA_list = get_list(input_miRNA_csv)
        transcript_list = get_list(input_transcript_csv)

        print('a list of input miRNA', miRNA_list)
        print('a list of input transcript/gene', transcript_list)

        #Supported miRNA and transcripts in the analysis
        supported_miRNA_list_pkl_path = "code/Data_file/Mapping/Supported_list/List_of_supported_miRNA.pkl"
        with open(supported_miRNA_list_pkl_path, 'rb') as file:
            supported_miRNA_list = pickle.load(file)

        supported_transcript_list_pkl_path = "code/Data_file/Mapping/Supported_list/List_of_supported_ENST.pkl"
        with open(supported_transcript_list_pkl_path, 'rb') as file:
            supported_transcript_list = pickle.load(file)

        miRNA_list, transcript_list, map_ENST_to_input, error_miRNAs, error_transcripts_type1, error_transcripts_type2, error_transcripts_type3 = step1_primiti(miRNA_list, transcript_list, supported_miRNA_list, supported_transcript_list)

        write_error_file(error_miRNAs, "miRNA", "w", job_id)
        write_error_file(error_transcripts_type1, "input_1", "a", job_id)
        write_error_file(error_transcripts_type2, "input_2", "a", job_id)
        write_error_file(error_transcripts_type3, "input_3", "a", job_id)

        check1 = check_empty_list(miRNA_list)
        check2 = check_empty_list(transcript_list)
        list_of_checks = [check1, check2]

        if ("False" in list_of_checks):
            sys.exit("SystemExit: error found")


        #PRIMITI-TS 2. Run a code here to get a first set of features as a CSV file -------------------------------------------------------------------##

        map_ENST_Pickle_conserv, miRNA_seq_dict, map_ENST_Pickle_seq, DRSNP_dict = load_all_pickle()

        miRNA_TS_df = step2_primiti(miRNA_list, transcript_list, map_ENST_Pickle_conserv, miRNA_seq_dict, map_ENST_Pickle_seq, DRSNP_dict)

        #PRIMITI-TS 3. Generate files for a list of miRNA and a list of transcripts ---------------------------------------------------------------##

        miRNA_ilearn_input, mRNA_ilearn_input, miRNA_to_miRNA_index_dict = step3_primiti(miRNA_TS_df, job_id)

        #PRIMITI-TS 4. Call a python 3 iLearn script that use csv files to generate a file contain miRNA & a file contain mRNA representations --------##

        mRNA_output_basic, mRNA_output_ACC, miRNA_output_basic, miRNA_output_ACC, miRNA_output_Pse = step4_primiti(miRNA_ilearn_input, mRNA_ilearn_input, job_id)

        #PRIMITI-TS 5. Run a code to retreive informations from 2 files and integrate as a part of features -> get a final set of features ------------##

        miRNA_TS_df = step5_primiti(miRNA_TS_df, mRNA_output_basic, mRNA_output_ACC, miRNA_output_basic, miRNA_output_ACC, miRNA_output_Pse, miRNA_to_miRNA_index_dict)

        #PRIMITI-TS 6. Standardize the data using pre-trained stardardizer ----------------------------------------------------------------------------##

        miRNA_TS_df, miRNA_TS_dist_df = step6_primiti(miRNA_TS_df)

        #PRIMITI-TS 7. Predict the probability for each miRNA-target site pairs using a pre-trained model and store it in the memory ------------------##

        primiti_ts_pred, primiti_ts_proba = step7_primiti(miRNA_TS_df)

        #PRIMITI-TM 8. Calculate the PRIMITI-TM feature for each miRNA-transcript pairs based on probability of each target site ----------------------##

        miRNA_TM_df, target_list = step8_primiti(miRNA_TS_df, primiti_ts_proba, miRNA_TS_dist_df)

        #PRIMITI-TM 9. Predict the probability for each miRNA-transcript pairs using a pre-trained model. ---------------------------------------------##

        primiti_tm_pred, primiti_tm_proba = step9_primiti(miRNA_TM_df)

        #PRIMITI-TM 10. Process them into an output format that contains a overall score & score for each site ----------------------------------------##
 
        Result_table1, Result_table1_path, Result_table2, Result_table2_path, map_index_TM_to_index_TS = step10_primiti(miRNA_TS_df, primiti_ts_pred, primiti_ts_proba, miRNA_TM_df, primiti_tm_pred, primiti_tm_proba, target_list, map_ENST_to_input, job_id)

        #Try remove temp file if allowed
        try:
            path = "result/temp/iLearn/"
            for file_name in os.listdir(path):
                file = path + file_name
                if os.path.isfile(file):
                    print('Deleting file:', file)
                    os.remove(file)
            print('successfuly remove temp files in', path)
        except:
            print('failed to remove file in', path)
  
        print("****************Finished running JOB****************")
        print(datetime.datetime.now())
        print("****************************************************")    

        
    except Exception as error_log:
        print(error_log)

