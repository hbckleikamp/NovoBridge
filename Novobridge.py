#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  28 17:32:32 2020

@author: hugokleikamp
"""
#%% clear variables and console

try:
    from IPython import get_ipython
    get_ipython().magic('clear')
    get_ipython().magic('reset -f')
except:
    pass

#%% parameters Part 1: Unipept submission
starting_vars=dir()

ALC_cutoff=40      # miniumum required ALC score (Peaks score)
Score_cutoff=-0.1  # mininum required score cutoff (DeepNovo score)

ppm_cutoff=20      # maximum allowed ppm
length_cutoff=7    # mininum required peptide length
Area_cutoff=0      # mininum required peak area
Intensity_cutoff=0 # mininum required intensity

#get variables for writing to output
current_vars=set(dir()).difference(starting_vars)
parameters=[i for i in locals().items() if i[0] in current_vars]

#%% parameters Part 2: Composition 
starting_vars=dir()

#(defaults slightly more stringent filtering)
comp_ALC_cutoff=70      # miniumum required ALC score (Peaks score)
comp_Score_cutoff=-0.1 # mininum required score cutoff (DeepNovo score)

comp_ppm_cutoff=15      # maximum allowed ppm
comp_length_cutoff=7    # mininum required peptide length
comp_Area_cutoff=0      # mininum required peak area
comp_Intensity_cutoff=0 # mininum required intensity

cutbranch=3    # minimum amount of unique peptides per taxonomic branch in denoising

normalize=False # normalize quantification to total for that ranks

#get variables for writing to output
current_vars=set(dir()).difference(starting_vars)
comp_parameters=[i for i in locals().items() if i[0] in current_vars]

#%% Parameters Part 3: Function
starting_vars=dir()

#which pathways to include
Pathways=['09100 Metabolism', 
          '09120 Genetic Information Processing'
          '09130 Environmental Information Processing'
          '09140 Cellular Processes'] 

#which levels to include
cats=["cat1","cat2","cat3","cat4"]                                    

current_vars=set(dir()).difference(starting_vars)
fun_parameters=comp_parameters+[i for i in locals().items() if i[0] in current_vars]

#%% Modules
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import random, re, requests
import threading, time, string
import pickle


from Bio.SeqUtils import molecular_weight
from itertools import chain, groupby
from collections import Counter
from openpyxl import load_workbook 

#%% change directory to script directory (should work on windows and mac)
import os
from pathlib import Path
from inspect import getsourcefile
os.chdir(str(Path(os.path.abspath(getsourcefile(lambda:0))).parents[0]))
print(os.getcwd())

#remove chained assignment warnings
#pd.options.mode.chained_assignment = None  # default='warn'
#%% load peptide file 
#Important! the only required column for analysis is a list of peptides with the header: "Peptide". 

"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
                                                    "Part 1: Annotation with unipept" 
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#input folder
pathin="input_peaks"
for filename in os.listdir(pathin): 
    if filename[0].isalnum(): # check if not a hidden file, filename should start with alphanumeric
        for Randomize in [False, 'scramble']: 
        
            if filename.endswith('.csv'):                               xlsdf=pd.read_csv(str(Path(pathin,filename)),sep=",")
            if filename.endswith('.xlsx') or filename.endswith('.xls'): xlsdf=pd.read_excel(str(Path(pathin,filename)))   
            
            xlsdf=xlsdf.fillna("0") #replace nans with zero
            
            #%% if DeepNovo output, convert to PEAKS output
            if 'predicted_sequence' in xlsdf.columns:
                
                xlsdf["Peptide"]=xlsdf["predicted_sequence"]           
                xlsdf=xlsdf[xlsdf["Peptide"]!="0"]
            
                #calculate ppm shift from m/z
                mass=np.zeros((1,len(xlsdf)))
                mass+=xlsdf["Peptide"].str.count("(Carbamidomethylation)").fillna(0)*57.021463
                mass+=xlsdf["Peptide"].str.count("(Oxidation)").fillna(0)*15.994915

                xlsdf['Peptide']=xlsdf['Peptide'].apply(lambda x: re.sub("[\(\[].*?[\)\]]", "", x).replace(",","")) #remove ptms in peptides
                xlsdf["Tag Length"]=xlsdf['Peptide'].apply(len)

                
                #%% calculate peptide mass differently
                
                std_aa_mass = {'G': 57.02146, 'A': 71.03711, 'S': 87.03203, 'P': 97.05276, 'V': 99.06841,
                               'T': 101.04768,'C': 103.00919,'L': 113.08406,'I': 113.08406,'J': 113.08406,
                               'N': 114.04293,'D': 115.02694,'Q': 128.05858,'K': 128.09496,'E': 129.04259,
                               'M': 131.04049,'H': 137.05891,'F': 147.06841,'U': 150.95364,'R': 156.10111,
                               'Y': 163.06333,'W': 186.07931,'O': 237.14773}    

                def mass_calc(x,std_aa_mass=std_aa_mass):
                    return sum(std_aa_mass.get(aa) for aa in x)+18.01056
                
                xlsdf['calculated_mass']=mass[0]+xlsdf['Peptide'].apply(lambda x: mass_calc(x)).values
                xlsdf['precursor_mass']=xlsdf['precursor_mz']*xlsdf['precursor_charge']-xlsdf['precursor_charge']*1.007277                
                xlsdf["ppm"]=(1000000/xlsdf['calculated_mass'])*(xlsdf['calculated_mass']-xlsdf['precursor_mass'])

                #rename columns
                if "feature_id" in xlsdf.columns: xlsdf=xlsdf.rename(columns={"feature_id":"Scan"})  
                if "feature_area" in xlsdf.columns: xlsdf=xlsdf.rename(columns={"feature_area":"Area"})  
                if "feature_intensity" in xlsdf.columns: xlsdf=xlsdf.rename(columns={"feature_intensity":"Intensity"})  
            
            #set datatypes as float
            for i in ['Tag Length','ALC (%)','predicted_score','ppm','Area','Intensity']:
                if i in xlsdf.columns:
                    xlsdf[i]=xlsdf[i].astype(float) 
            
            #add Scan if not present, add Scan
            if 'Scan' not in xlsdf.columns: 
                xlsdf['Scan']=list(range(0,len(xlsdf)))
                        
            # recalibrate ppm   
            if 'ppm' in xlsdf.columns:
                hist, bin_edges = np.histogram(xlsdf['ppm'].tolist(),bins=100)
                xlsdf['ppm']=abs(xlsdf['ppm'] -(bin_edges[np.where(hist==np.max(hist))[0]]
                                               +bin_edges[np.where(hist==np.max(hist))[0]+1])/2)
            #remove ptms in peptides
            peplist=xlsdf['Peptide'].tolist()
            xlsdf['Peptide']=list(map(lambda x: re.sub("[\(\[].*?[\)\]]", "", x),peplist)) 
            
            #%% filter by ALC scores, ppm and peptide lengths 
            filt="on"
            if filt=="on":    
                if 'Tag Length' in xlsdf.columns: 
                    xlsdf=xlsdf[xlsdf['Tag Length']>=length_cutoff]
                
                if 'ALC (%)' in xlsdf.columns:
                    xlsdf=xlsdf[xlsdf['ALC (%)']>=ALC_cutoff] #scoring Peaks 
                
                if 'predicted_score' in xlsdf.columns:
                    xlsdf=xlsdf[xlsdf['predicted_score']<=Score_cutoff] #scoring DeepNovo 
                
                if 'ppm' in xlsdf.columns: 
                    xlsdf=xlsdf[xlsdf['ppm']<=ppm_cutoff]
                
                if 'Area' in xlsdf.columns: 
                    xlsdf=xlsdf[xlsdf['Area']>=Area_cutoff]
            
                if 'Intensity' in xlsdf.columns: 
                    xlsdf=xlsdf[xlsdf['Intensity']>=Intensity_cutoff]
                                
            #%% randomization (optional)
            if   Randomize=='scramble': #scramble all in front of cleavage site
                 xlsdf['Peptide']=[''.join(random.sample(i[:len(i)-1], len(i)-1)+[i[-1]]) for i in xlsdf['Peptide']]
                   
            #%% submit peptides to unipept
            
            #urls
            twl='http://api.unipept.ugent.be/api/v1/pept2lca.json?input[]=';
            twr='&equate_il=true&extra=true&names=true';
            fwl='http://api.unipept.ugent.be/api/v1/pept2funct.json?input[]=';
            fwr='&equate_il=true';
            
            fields=["peptide",
                    "taxon_name",
                    "superkingdom_name",
                    "phylum_name",
                    "class_name",
                    "order_name",
                    "family_name",
                    "genus_name"];
              
            unipeps=np.unique(xlsdf['Peptide'])
            batchsize=100
            steps=list(range(0,len(unipeps),batchsize))
            
            # Threading unipept
            def unipept_scrape(r,url):
                
                while True:
                    try:
                        r.extend(requests.get(url,stream=True).json())
                        break
                    except:
                        "sleeping"
                        time.sleep(2)
                        
                        
            # is used to divide a list into "chunks" with a generator
            def chunks(lst,n):
                for i in range(0,len(lst),n):
                    yield lst[i:i+n]
            

            taxalist=list()
            funlist=list()

            threads=[]
            counter=0
            for chunk in chunks(unipeps,batchsize):
                counter+=1
                print(counter)
                query="&input[]=".join(chunk)
                    
                #taxonomy
                turl=twl+query+twr 
                t=threading.Thread(target=unipept_scrape, args=[taxalist,turl])
                t.start()
                threads.append(t)
                
                #function
                furl=fwl+query+fwr 
                t=threading.Thread(target=unipept_scrape, args=[funlist,furl])
                t.start()
                threads.append(t)
                 
            for thread in threads:
                thread.join()
            
            #%% post processing
            #remove non-regular taxonomic ranks
            taxa=pd.DataFrame(taxalist)
            [taxa.pop(x) for x in taxa.columns if x not in fields]
            
            #parse function dataframe
            funs=pd.DataFrame(funlist)
            funs["ec"]= funs["ec"].apply( lambda x: " ".join(pd.json_normalize(x)["ec_number"]) if x else [])
            funs["go"]= funs["go"].apply( lambda x: " ".join(pd.json_normalize(x)["go_term"]) if x else [])
            funs["ipr"]=funs["ipr"].apply(lambda x: " ".join(pd.json_normalize(x)["code"]) if x else [])

            #merge annotations
            taxa.set_index('peptide',inplace=True, drop=True)
            funs.set_index('peptide',inplace=True, drop=True)
            xlsdf.set_index('Peptide',inplace=True, drop=False)
            xlsdf=xlsdf.join(taxa,how='left')
            xlsdf=xlsdf.join(funs,how='left')
            
            xlsdf=xlsdf.mask(xlsdf.applymap(str).eq('[]'))
            xlsdf=xlsdf.fillna("")
            #%% write result

            pathout="output_unipept"
            if not os.path.exists(pathout): os.makedirs(pathout)
            
            xlsfilename=str(Path(pathout,"unipept_"+filename.replace(os.path.splitext(filename)[1], '.xlsx')))
            writer = pd.ExcelWriter(xlsfilename, engine='xlsxwriter')
            xlsdf.to_excel(writer, sheet_name='Output')
            pd.DataFrame(parameters,columns=["Name","Value"]).to_excel(writer, sheet_name='Parameters')
            writer.save()
            writer.close()
    
    #%%
            """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
                                                                "Part 2: Compositional analysis" 
            """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
            
            
            #%% Pre_processing filtering & denoising
            
            filt="on"
            if filt=="on":    
                if 'Tag Length' in xlsdf.columns: 
                    xlsdf=xlsdf[xlsdf['Tag Length']>=comp_length_cutoff]
                
                if 'ALC (%)' in xlsdf.columns:
                    xlsdf=xlsdf[xlsdf['ALC (%)']>=comp_ALC_cutoff]
                    
                if 'predicted_score' in xlsdf.columns:
                    xlsdf=xlsdf[xlsdf['predicted_score']<=comp_Score_cutoff] #scoring DeepNovo 
                
                if 'ppm' in xlsdf.columns: 
                    xlsdf=xlsdf[xlsdf['ppm']<=comp_ppm_cutoff]
                
                if 'Area' in xlsdf.columns: 
                    xlsdf=xlsdf[xlsdf['Area']>=comp_Area_cutoff]
            
                if 'Intensity' in xlsdf.columns: 
                    xlsdf=xlsdf[xlsdf['Intensity']>=comp_Intensity_cutoff]
            
            taxa=["superkingdom_name",
                    "phylum_name",
                    "class_name",
                    "order_name",
                    "family_name",
                    "genus_name"]
              
            denoise="on"
            if denoise=="on":
                unirows=xlsdf[["Peptide"]+taxa].drop_duplicates()[taxa].astype(str).values.tolist()
                jbranch=["#*".join(i) for i in unirows]
                fbranch=[branch for branch, counts in Counter(jbranch).items() if counts >= cutbranch]
                allowed_taxids=set(chain(*[i.split("#*") for i in fbranch]))
                for i in taxa:
                    xlsdf.loc[~xlsdf[i].isin(allowed_taxids),i]="" 
                        
            #%% krona plot (only coded for spectral counting)
            
            pathout="output_composition"
            if not os.path.exists(pathout): os.makedirs(pathout)
            
            if "krona_template.xlsm" in os.listdir():
                
                grouped=xlsdf.groupby(taxa)
                vals=grouped.size() 
                vals=vals.reset_index(name="Count")  
                vals=vals[~(vals[taxa]=="").all(axis=1)] #remove empty rows
                
                #results
                branches=vals[taxa].values.tolist()  
                counts=  vals["Count"].tolist()
                
                #fill gaps
                for i in range(len(branches)):
                    j=[j for j in range(len(branches[i])) if branches[i][j]!=""]
                    for l in range(0,max(j)):
                        if branches[i][l]=="": branches[i][l]="annotation_gap"     
                vals[taxa]=branches
                
                branchdf=vals

                kronafilename=str(Path(pathout,filename.replace(os.path.splitext(filename)[1], '.xlsm')))
                letters=list(string.ascii_uppercase)
                wb = load_workbook("krona_template.xlsm",read_only=False, keep_vba=True)
                ws = wb.active
                for i in range(0,len(branches)):
                
                    ws['A{0}'.format(4+i)].value=filename
                    ws['B{0}'.format(4+i)].value=counts[i]
                    for j in range(0,len(taxa)):
                        ws['{0}{1}'.format(letters[j+3],4+i)].value=branches[i][j]
                wb.save(kronafilename)   
                
            else:
                print("No Krona template found, proceding without generating krona plots")
                
            #%% make stacked bars
            
            def stacked_bar(taxa,df,ylabel,pathout,filename): 
            
                labels=[i.split("_name")[0] for i in taxa];
                countcols=[i for i in df.columns if "count" in i]
                absdat=df[countcols] 
                absdat.columns=labels
                normdat= absdat/np.nansum(absdat.to_numpy(),axis=0)
                
                figure, axes = plt.subplots(1, 2)
                ax1=absdat.T.plot(ax=axes[0],kind='bar', stacked=True, figsize=(10, 6), legend=False)
                ax1.set_ylabel(ylabel)
                ax1.set_xlabel('taxonomic ranks')
                ax1.set_xticklabels(labels, rotation=30)
                ax1.set_title("Absolute")
                
                ax2=normdat.T.plot(ax=axes[1],kind='bar', stacked=True, figsize=(10, 6), legend=False)
                ax2.set_ylabel(ylabel)
                ax2.set_xlabel('taxonomic ranks')
                ax2.set_xticklabels(labels, rotation=30)
                ax2.set_title("Normalized")
                
                figname=str(Path(pathout,(ylabel+"_"+filename.replace(os.path.splitext(filename)[1], '.png'))))
                plt.savefig(figname)
                return figname      
            
            #%% quantification and visual outputs, write output to xlsx
            pathout="output_composition"
            if not os.path.exists(pathout): os.makedirs(pathout)
            
            xlsfilename=str(Path(pathout,"composition_"+filename.replace(os.path.splitext(filename)[1], '.xlsx')))
            if Randomize=="scramble":
                xlsfilename=str(Path(pathout,"Rand_composition_"+filename.replace(os.path.splitext(filename)[1], '.xlsx')))
            
            writer = pd.ExcelWriter(xlsfilename, engine='xlsxwriter')
            pd.DataFrame(comp_parameters,columns=["Name","Value"]).to_excel(writer, sheet_name='Parameters')
            letters=list(string.ascii_uppercase)
                        
            if Randomize==False:            
                #Spectral counts
                Spectral_Counts=pd.DataFrame()
                for i in taxa:
                    counts=Counter(xlsdf[i].astype(str))
                    if "" in counts.keys(): counts.pop("") #ignore unassigned peptides
                    freqs=pd.Series(counts).reset_index().sort_values(by="index", axis=0, ascending=True, inplace=False).reset_index(drop=True)
                    freqs.columns=[i,i+"_count"]
                    
                    if normalize: freqs[i+"_count"]=freqs[i+"_count"]/freqs[i+"_count"].sum()*100 #normalize to 100% for easier comparison to metadata
                    #alphabetic sorting
                    Spectral_Counts = pd.concat([Spectral_Counts, freqs], axis=1) 
    
                namecols=[i for i in Spectral_Counts.columns  if "count" not in i] 
                countcols=[i for i in Spectral_Counts.columns  if "count" in i] 
                Spectral_Counts=Spectral_Counts.fillna(0)
                
                #write spectral counts    
                Spectral_Counts[namecols].to_excel(writer, sheet_name='TAX_LINEAGES')
                Spectral_Counts[countcols].to_excel(writer, sheet_name='TAX_COUNTS')
                stacked_bar(taxa,Spectral_Counts,'spectral counts',pathout, filename)
                writer.save()   
                writer.close()
                
            if Randomize=="scramble":            
                #Spectral counts
                Spectral_Counts=pd.DataFrame()
                counts=[sum(xlsdf[i]!="") for i in taxa]
                ranks=[i.split("_")[0] for i in taxa]                          
                #write spectral counts    
                pd.Series(counts,index=ranks).to_excel(writer, sheet_name='RANDOM_Cts')
                writer.save()   
                writer.close()
            #%%
    
            """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
                                                                "Part 3: Functional analysis" 
            """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

            #%% pre processing
            
            #preselect xlsdf to only have rows that have ec numbers
            xlsdf=xlsdf[xlsdf["ec"].notnull()]
            #calculte total area/intensity/spectral counts for normalization
            totSpect=len(xlsdf)
            if "Area" in xlsdf.columns:      totArea=sum(xlsdf["Area"])
            if "Intensity" in xlsdf.columns: totInt =sum(xlsdf["Intensity"])
            
            #separate ec annotations for exact matching
            xlsdf["ec"]=xlsdf["ec"].apply(lambda x: str(x).split(" ")) 
            expdf=xlsdf.explode('ec')    
            
           
            #%% quantification
            
            #check if kegg database is present
            if not any("keg.pkl" in i for i in os.listdir()):
                print("no local kegg database found, please run download_utilities.py before functional annotation can be done")
                
            else: #Load local kegg database
                for pkl in os.listdir():
                    if pkl.endswith('keg.pkl'):  
                        infile = open(pkl,'rb')
                        keggdf  = pickle.load(infile)
                        infile.close()
                
                #only select pathways that are in the parameters
                keggdf=keggdf[keggdf.isin(Pathways).any(axis=1)] 
                #only select cats that are in the parameters
                keggdf=keggdf[cats+["ECs"]]
                
                #always quantify on spectral counts, but if present: quantify based on area/intensity
                columns=["Area","Intensity"] 
                methods=["total"] 
                
                for j in range(len(cats)):
                    df=pd.DataFrame()
                    print(j)
                    counter=0
                    for i in np.unique(keggdf[cats[j]].to_numpy()):
                        counter+=1
                        print(counter)
                        ECs=keggdf.loc[i==keggdf[cats[j]],"ECs"]
                        uniscans=np.unique(expdf.loc[expdf["ec"].isin(ECs),"Scan"].to_numpy())
                        if len(uniscans):
                            sr=pd.Series() #every new series contains a single pathwayand its quantification
                            sr[cats[j]]=i  #name of category (used for joining)
                            sr[cats[j]+"_"+"Spectral_counts"]=len(uniscans) #spectral counts
                            
                            # Area or Intensity quantification
                            for c in columns:               
                                if c in xlsdf.columns:
                                    values=xlsdf.loc[xlsdf["Scan"].isin(uniscans),c]
                                    for m in methods:        
                                        if  m=="total": sr[cats[j]+"_"+c+"_"+m]=sum(values)                                      
                            
                            #append series to dataframe 
                            df=df.append(sr,ignore_index=True)
                            df = df[sr.index] #reorder based on series index
                    if len(df):
                        keggdf=keggdf.merge(df,on=cats[j]) #join quantifiactions with keggdf        
            
                #write to xlsx
                keggdf=keggdf.sort_values(cats)
                
                #normalize to total spectra/area/intensity, only normalize this way when total counts are used
                for i in keggdf.columns:
                    if "Spectral"   in i: keggdf[i]=keggdf[i]/totSpect
                    if "Area_total" in i: keggdf[i]=keggdf[i]/totArea
                    if "Intensity"  in i: keggdf[i]=keggdf[i]/totInt

                pathout="output_function"
                if not os.path.exists(pathout): os.makedirs(pathout)
                
                xlsfilename=str(Path(pathout,"function_"+filename.replace(os.path.splitext(filename)[1], '.xlsx')))
                writer = pd.ExcelWriter(xlsfilename, engine='xlsxwriter')
                pd.DataFrame(fun_parameters,columns=["Name","Value"]).to_excel(writer, sheet_name='Parameters')
    
                for i in cats:
                    cols=[j for j in keggdf.columns if i in j]
                    #add previous cols
                    c=0; prev=list()
                    while i!=cats[c]:
                        prev.append(cats[c])
                        c+=1
                    cols=prev+cols
                    keggdf[cols].drop_duplicates().to_excel(writer, sheet_name=i) 
                 
                writer.save()   
                writer.close()
                        
#%% Cleanup (remove .png files)
cleanup=True
if cleanup:
    [os.remove(i) for i in os.listdir() if i.endswith(".png")]



