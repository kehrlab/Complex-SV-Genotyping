# written by Clementine Döhring
# 
# Takes VCF files from gnomAD-SV as input and outputs variant descriptions
# for GGTyper in a json file. Additonal files for variant frequencies and more are generated.

from sys import argv   
import pysam            
from pysam import VariantFile
import json
import csv


######## functions

def SplitInterval(CPXinterval):
    region_interval= CPXinterval.split("_")
    chromWindow = region_interval[1].split(":")
    chrm_intervall = chromWindow[0]
    BeginEnd = chromWindow[1].split("-")
    start_intervall = int(BeginEnd[0])
    end_intervall = int(BeginEnd[1])
    return chrm_intervall, start_intervall, end_intervall

def getName(rec):
    SV_Name= str(rec.chrom)+'_'+str(rec.pos)+ '_'+ str(rec.stop)+'_'+str(rec.info['CPX_TYPE'])
    k = 0
    while True:
        if SV_Name not in SV_Variants:
            return SV_Name
        SV_Name= str(rec.chrom)+'_'+str(rec.pos)+ '_'+ str(rec.stop)+'_'+str(rec.info['CPX_TYPE'])+'_'+ str(k)
        k = k +1

def getNameCanonical(rec):
    SV_Name= str(rec.chrom)+'_'+str(rec.pos)+ '_'+ str(rec.stop)+'_'+str(rec.info['SVTYPE'])
    return SV_Name

def getVariantType(rec):
    if rec.info['SVTYPE']== 'CPX':
        VariantType = str(rec.info['CPX_TYPE'])
    else:
        VariantType = str(rec.info['SVTYPE'])   #######DEL = DEL :ME
    return VariantType

def getMaximumAF(rec):
    if rec.alts[0] != '<CNV>':
        VariantType = getVariantType(rec)
        if rec.info['AF_afr'][0]> max_AF['afr'][VariantType][0]:
                    max_AF['afr'][VariantType][0] =rec.info['AF_afr'][0]
                    max_AF['afr'][VariantType][1] = rec.pos
                    #print("afr",rec.pos, VariantType, max_AF['afr'][VariantType])
        if rec.info['AF_eas'][0]> max_AF['eas'][VariantType][0]:
                    max_AF['eas'][VariantType][0] =rec.info['AF_eas'][0]
                    max_AF['eas'][VariantType][1] = rec.pos
                    #print("eas",rec.pos, rec.info['CPX_TYPE'], max_AF['eas'][VariantType])
        if (rec.info['AF_fin'][0] or rec.info['AF_nfe'][0]) > max_AF['fin/nfe'][VariantType][0]:
                    max_AF['fin/nfe'][VariantType][0] = max(rec.info['AF_fin'][0], rec.info['AF_nfe'][0])
                    max_AF['fin/nfe'][VariantType][1] = rec.pos
                    #print("fin/nfe",rec.pos, rec.info['CPX_TYPE'], max_AF['fin/nfe'][VariantType])

def checkChromosome(rec): #checks if all intervalls & start are on the same chromosome
    chr=['','','']
    for i in range(len(rec.info['CPX_INTERVALS'])):
        chr[i],startInt, endInt= SplitInterval(rec.info['CPX_INTERVALS'][i])
    if len(rec.info['CPX_INTERVALS'])==2:
        control = chr[0] ==chr[1]
    if len(rec.info['CPX_INTERVALS'])==3:
        control = chr[0] ==chr[1] == chr[2]
    if len(rec.info['CPX_INTERVALS'])==1:
        control = chr[0] == rec.chrom
    if control == False:
        #print(control, rec.info['CPX_TYPE'])
        if (rec.info['CPX_TYPE']!= 'dDUP'):
            if (rec.info['CPX_TYPE']!= 'dDUP_iDEL') :
                print('! Check Chromosomes in Variant description in: '+ rec.info['CPX_TYPE'])

def getAFfreq(rec):
    AF_entry = [getName(rec), rec.info['AF_fin'][0],rec.info['AN_fin'],rec.info['AF_nfe'][0],rec.info['AN_nfe'],rec.info['AF_eas'][0],rec.info['AN_eas'],rec.info['AF_afr'][0],rec.info['AN_afr']]
    AF_freqs.append(AF_entry)

def createJSON(rec):
    if rec.alts[0] == '<DEL>':
        SV_Name= getNameCanonical(rec)
        DEL= {
            SV_Name : {
                "VAR":{
                    rec.chrom:{
                        "1":{
                            "rNameLeft": rec.chrom,
                            "rNameRight": rec.chrom,
                            "xLeft": rec.pos,
                            "xRight": rec.stop,
                            "directionLeft": "left",
                            "directionRight": "right"
                        }
                    }
                }
                
            }
        }
        countVariants['countDEL'] = countVariants['countDEL'] +1
        SV_Variants.update(DEL) 
        #getMaximumAF(rec)

    elif rec.alts[0] == '<DUP>':
        SV_Name= getNameCanonical(rec)
        DUP= {
            SV_Name : {
                "VAR":{
                    rec.chrom:{
                        "1":{
                            "rNameLeft": rec.chrom,
                            "rNameRight": rec.chrom,
                            "xLeft": rec.stop,
                            "xRight": rec.pos,
                            "directionLeft": "left",
                            "directionRight": "right"
                        },
                        "2":{
                            "rNameLeft": rec.chrom,
                            "rNameRight": rec.chrom,
                            "xLeft": rec.stop,
                            "xRight": rec.stop,
                            "directionLeft": "left",
                            "directionRight": "right"
                        }
                    }
                }
                
            }
        }  
        countVariants['countDUP'] = countVariants['countDUP'] +1
        SV_Variants.update(DUP) 
        #getMaximumAF(rec)

    elif rec.alts[0] == '<INS>':
        SV_Name= getNameCanonical(rec)
        countVariants['countINS'] = countVariants['countINS'] +1
        #getMaximumAF(rec)
        ###### skip insertion
        # INS= {
        #     SV_Name : {
        #         "VAR":{
        #             rec.chrom:{
        #                 "1":{
        #                     "rNameLeft": rec.chrom,
        #                     "rNameRight": rec.chrom,
        #                     "xLeft": rec.stop,
        #                     "xRight": rec.pos,
        #                     "directionLeft": "left",
        #                     "directionRight": "right"
        #                 },
        #                 "2":{
        #                     "rNameLeft": rec.chrom,
        #                     "rNameRight": rec.chrom,
        #                     "xLeft": rec.stop,
        #                     "xRight": rec.stop,
        #                     "directionLeft": "left",
        #                     "directionRight": "right"
        #                 }
        #             }
        #         }
                
        #     }
        # }  
        # SV_Variants.update(INS) 

    elif rec.alts[0] == '<INV>':
        SV_Name= getNameCanonical(rec)
        INV= {
            SV_Name : {
                "VAR":{
                    rec.chrom:{
                        "1":{
                            "rNameLeft": rec.chrom,
                            "rNameRight": rec.chrom,
                            "xLeft": rec.pos-1,
                            "xRight": rec.stop-1,
                            "directionLeft": "left",
                            "directionRight": "left"
                        },
                        "2":{
                            "rNameLeft": rec.chrom,
                            "rNameRight": rec.chrom,
                            "xLeft": rec.pos,
                            "xRight": rec.stop,
                            "directionLeft": "right",
                            "directionRight": "right"
                        }
                    }
                }
                
            }
        }  
        countVariants['countINV'] = countVariants['countINV'] +1
        SV_Variants.update(INV)
        #getMaximumAF(rec)

    elif rec.alts[0] == '<CTX>':
        SV_Name= getNameCanonical(rec)
        CTX= {
            SV_Name : {
                "VAR":{
                    rec.chrom:{
                        "1":{
                            "rNameLeft": rec.chrom,
                            "rNameRight": rec.info['CHR2'],
                            "xLeft": rec.pos-1,
                            "xRight": rec.info['POS2'],
                            "directionLeft": "left",
                            "directionRight": "right"
                        },
                    },
                    rec.info['CHR2']:{
                        "1":{
                            "rNameLeft": rec.info['CHR2'],
                            "rNameRight": rec.chrom,
                            "xLeft": rec.info['POS2']-1,
                            "xRight": rec.pos,
                            "directionLeft": "left",
                            "directionRight": "right"
                        }
                    }
                }
                
            }
        }  
        countVariants['countCTX'] = countVariants['countCTX'] +1
        SV_Variants.update(CTX)
        #getMaximumAF(rec)

    elif rec.alts[0] == '<BND>':
        countVariants['countBND'] = countVariants['countBND'] +1
        #getMaximumAF(rec)


    elif rec.alts[0] == '<CPX>':
        #getMaximumAF(rec)
        #print(rec.alts, rec.info['CPX_TYPE'], rec.pos)
        #print(rec.chrom,rec.pos,rec.info['CPX_TYPE'],rec.info['CPX_INTERVALS'],rec.stop)
        checkChromosome(rec)
        
        ########## FOR CPX_TYPE: delINV
        if rec.info['CPX_TYPE']== 'delINV':
            SV_Name= str(rec.chrom)+'_'+str(rec.pos)+'_delINV'#getName(rec)     #
            k = 0
            while True:
                if SV_Name not in SV_Variants:
                    break
                SV_Name= str(rec.chrom)+'_'+str(rec.pos)+ '_'+str(rec.info['CPX_TYPE'])+'_'+ str(k)
                k = k +1
        #    # print(SV_Name)
        #     interval_1= rec.info['CPX_INTERVALS'][0].split("_")
        #     chromWindow1 = interval_1[1].split(":")
        #     chrm1 = chromWindow1[0]
        #     #print(chromWindow1)
        #     BeginEnd1 = chromWindow1[1].split("-")
        #     start1 = int(BeginEnd1[0])
        #     end1 = int(BeginEnd1[1])

        #     interval_2= rec.info['CPX_INTERVALS'][1].split("_")
        #     chromWindow2 = interval_2[1].split(":")
        #     chrm2 = chromWindow2[0]
        #     BeginEnd2 = chromWindow2[1].split("-")
        #     start2 = int(BeginEnd2[0])
        #     end2 = int(BeginEnd2[1])
            chrm1,start1, end1=SplitInterval(rec.info['CPX_INTERVALS'][0])
            chrm2,start2, end2=SplitInterval(rec.info['CPX_INTERVALS'][1])
            delINV= {
                SV_Name : {
                    "VAR":{
                        rec.chrom:{
                            "1":{
                                "rNameLeft": chrm1,
                                "rNameRight": chrm2,
                                "xLeft": start1,
                                "xRight": end2,
                                "directionLeft": "left",
                                "directionRight": "left"
                            },
                            "2" : {
                                "rNameLeft": chrm2,
                                "rNameRight": chrm2,
                                "xLeft": start2+1,
                                "xRight": rec.stop +1,
                                "directionLeft": "right",
                                "directionRight": "right",
                                "check END<INVend+10":end2<= rec.stop < end2 +10,
                                "End Inversion": end2,
                                "End Variant":rec.stop
                            }
                        }
                    }
                    
                }
            }
            #print(delINV)
            countVariants['countdelINV'] = countVariants['countdelINV'] +1
            svNames[rec.id] = SV_Name
            SV_Variants.update(delINV)

        ########## FOR CPX_TYPE: INVdel
        elif rec.info['CPX_TYPE']== 'INVdel':
            SV_Name= str(rec.chrom)+'_'+str(rec.pos)+'_INVdel' #getName(rec)     #
            k = 0
            while True:
                if SV_Name not in SV_Variants:
                    break
                SV_Name= str(rec.chrom)+'_'+str(rec.pos)+ '_'+str(rec.info['CPX_TYPE'])+'_'+ str(k)
                k = k +1
            
            chrm1,start1, end1=SplitInterval(rec.info['CPX_INTERVALS'][0])
            chrm2,start2, end2=SplitInterval(rec.info['CPX_INTERVALS'][1])
            INVdel= {
                SV_Name : {
                    "VAR":{
                        rec.chrom:{
                            "1":{
                                "rNameLeft": chrm1,
                                "rNameRight": chrm1,
                                "xLeft": rec.pos-1,
                                "xRight": end1-1,
                                "directionLeft": "left",
                                "directionRight": "left"
                            },
                            "2" : {
                                "rNameLeft": chrm1,
                                "rNameRight": chrm2,
                                "xLeft": start1,
                                "xRight": rec.stop,
                                "directionLeft": "right",
                                "directionRight": "right",
                                "check start-10<pos":start1 -10 < rec.pos <= start1
                            }
                        }
                    }
                    
                }
            }
            countVariants['countINVdel'] = countVariants['countINVdel'] +1
            SV_Variants.update(INVdel)
            svNames[rec.id] = SV_Name
            ##getMaximumAF(rec)

        ######## FOR CPX_TYPE: dDUP
        elif rec.info['CPX_TYPE']== 'dDUP':
            SV_Name= str(rec.chrom)+'_'+str(rec.pos)+'_dDUP'#getName(rec)     #
            k = 0
            while True:
                if SV_Name not in SV_Variants:
                    break
                SV_Name= str(rec.chrom)+'_'+str(rec.pos)+ '_'+str(rec.info['CPX_TYPE'])+'_'+ str(k)
                k = k +1
            
            
            ####### with Inversion
            if len(rec.info['CPX_INTERVALS']) ==2:
                chrm1,start1, end1=SplitInterval(rec.info['CPX_INTERVALS'][0])
                dDUP= {
                    SV_Name : {
                        "VAR":{
                            rec.chrom:{
                                "1":{
                                    "rNameLeft": rec.chrom,
                                    "rNameRight": chrm1,
                                    "xLeft": rec.pos,
                                    "xRight": end1,
                                    "directionLeft": "left",
                                    "directionRight": "left"
                                },
                                "2" : {
                                    "rNameLeft": chrm1,
                                    "rNameRight": rec.chrom,
                                    "xLeft": start1,
                                    "xRight": rec.stop+1,
                                    "directionLeft": "right",
                                    "directionRight": "right",
                                    "start > end -10":rec.pos > rec.stop -10
                                }
                            },
                            
                        }
                        
                    }
                }   

            ####### without Inversion
            elif len(rec.info['CPX_INTERVALS']) ==1:
                chrm1,start1, end1=SplitInterval(rec.info['CPX_INTERVALS'][0])
                dDUP= {
                    SV_Name : {
                        "VAR":{
                            rec.chrom:{
                                "1":{
                                    "rNameLeft": rec.chrom,
                                    "rNameRight": chrm1,
                                    "xLeft": rec.pos,
                                    "xRight": start1,
                                    "directionLeft": "left",
                                    "directionRight": "right"
                                },
                                "2" : {
                                    "rNameLeft": chrm1,
                                    "rNameRight": rec.chrom,
                                    "xLeft": end1+1,
                                    "xRight": rec.stop+1,
                                    "directionLeft": "left",
                                    "directionRight": "right",
                                    "start > end -10":rec.pos > rec.stop -10
                                }
                            },
                            
                        }
                        
                    }
                }
            else:
                print("OTHER CASE",rec.info['CPX_INTERVALS'])
            countVariants['countdDUP'] =countVariants['countdDUP'] +1
            SV_Variants.update(dDUP)
            svNames[rec.id] = SV_Name
        
        ######## FOR CPX_TYPE: dDUP_iDEL
        elif rec.info['CPX_TYPE']== 'dDUP_iDEL':
            SV_Name= str(rec.chrom)+'_'+str(rec.pos)+'_dDUP_iDEL'#getName(rec)     #
            k = 0
            while True:
                if SV_Name not in SV_Variants:
                    break
                SV_Name= str(rec.chrom)+'_'+str(rec.pos)+ '_'+str(rec.info['CPX_TYPE'])+'_'+ str(k)
                k = k +1
            
            
            ####### with Inversion
            if len(rec.info['CPX_INTERVALS']) ==3:
                chrm1,start1, end1=SplitInterval(rec.info['CPX_INTERVALS'][0])
                dDUP_iDEL= {
                    SV_Name : {
                        "VAR":{
                            rec.chrom:{
                                "1":{
                                    "rNameLeft": rec.chrom,
                                    "rNameRight": chrm1,
                                    "xLeft": rec.pos,
                                    "xRight": end1,
                                    "directionLeft": "left",
                                    "directionRight": "left"
                                },
                                "2" : {
                                    "rNameLeft": chrm1,
                                    "rNameRight": rec.chrom,
                                    "xLeft": start1,
                                    "xRight": rec.stop+1,
                                    "directionLeft": "right",
                                    "directionRight": "right",
                                }
                            }
                        }
                        
                    }
                }   

            ####### without Inversion
            elif len(rec.info['CPX_INTERVALS']) ==2:
                chrm1,start1, end1=SplitInterval(rec.info['CPX_INTERVALS'][0])
                dDUP_iDEL= {
                    SV_Name : {
                        "VAR":{
                            rec.chrom:{
                                "1":{
                                    "rNameLeft": rec.chrom,
                                    "rNameRight": chrm1,
                                    "xLeft": rec.pos,
                                    "xRight": start1,
                                    "directionLeft": "left",
                                    "directionRight": "right"
                                },
                                "2" : {
                                    "rNameLeft": chrm1,
                                    "rNameRight": rec.chrom,
                                    "xLeft": end1+1,
                                    "xRight": rec.stop+1,
                                    "directionLeft": "left",
                                    "directionRight": "right",
                                }
                            }
                        }
                        
                    }
                }
            countVariants['countdDUP_iDEL'] = countVariants['countdDUP_iDEL'] +1
            SV_Variants.update(dDUP_iDEL)
            svNames[rec.id] = SV_Name

        ########## FOR CPX_TYPE: delINVdel
        elif rec.info['CPX_TYPE']== 'delINVdel':
            SV_Name= str(rec.chrom)+'_'+str(rec.pos)+'_delINVdel'#getName(rec)     #
            k = 0
            while True:
                if SV_Name not in SV_Variants:
                    break
                SV_Name= str(rec.chrom)+'_'+str(rec.pos)+ '_'+str(rec.info['CPX_TYPE'])+'_'+ str(k)
                k = k +1
            
            chrm1,start1, end1=SplitInterval(rec.info['CPX_INTERVALS'][0])
            chrm2,start2, end2=SplitInterval(rec.info['CPX_INTERVALS'][1])
            chrm3,start3, end3=SplitInterval(rec.info['CPX_INTERVALS'][2])
            delINVdel= {
                SV_Name : {
                    "VAR":{
                        rec.chrom:{
                            "1":{
                                "rNameLeft": chrm1,
                                "rNameRight": chrm2,
                                "xLeft": start1-1,
                                "xRight": end2-1,
                                "directionLeft": "left",
                                "directionRight": "left"
                            },
                            "2" : {
                                "rNameLeft": chrm2,
                                "rNameRight": chrm1,
                                "xLeft": start2,
                                "xRight": rec.stop,
                                "directionLeft": "right",
                                "directionRight": "right",
                            }
                        }
                    }
                    
                }
            }
            countVariants['countdelINVdel'] = countVariants['countdelINVdel'] +1
            SV_Variants.update(delINVdel)
            svNames[rec.id] = SV_Name  

        ########## FOR CPX_TYPE: dupINVdup
        elif rec.info['CPX_TYPE']== 'dupINVdup':
            SV_Name=getName(rec)
            chrm1,start1, end1=SplitInterval(rec.info['CPX_INTERVALS'][0])
            chrm2,start2, end2=SplitInterval(rec.info['CPX_INTERVALS'][1])
            chrm3,start3, end3=SplitInterval(rec.info['CPX_INTERVALS'][2])
            dupINVdup= {
                SV_Name : {
                    "VAR":{
                        rec.chrom:{
                            "1":{
                                "rNameLeft": chrm1,
                                "rNameRight": chrm3,
                                "xLeft": end1-1,
                                "xRight": end3-1,
                                "directionLeft": "left",
                                "directionRight": "left"
                            },
                            "2" : {
                                "rNameLeft": chrm1,
                                "rNameRight": chrm3,
                                "xLeft": start1,
                                "xRight": start3,
                                "directionLeft": "right",
                                "directionRight": "right",
                            }
                        }
                    }
                    
                }
            }
            countVariants['countdupINVdup'] = countVariants['countdupINVdup'] +1
            SV_Variants.update(dupINVdup)
            svNames[rec.id] = SV_Name
            #if rec.info['afr_AF'][0]> max_AF_afr_dupINVdup:
            #   max_AF_afr_dupINVdup =rec.info['afr_AF'][0]
            #  print(rec.pos, rec.info['CPX_TYPE'], max_AF_afr_dupINVdup)
            #print(rec.info['afr_AF'])
            ##getMaximumAF(rec)
                
        ########## FOR CPX_TYPE: dupINVdel
        elif rec.info['CPX_TYPE']== 'dupINVdel':
            SV_Name=getName(rec)
            chrm1,start1, end1=SplitInterval(rec.info['CPX_INTERVALS'][0])
            chrm2,start2, end2=SplitInterval(rec.info['CPX_INTERVALS'][1])
            #chrm3,start3, end3=SplitInterval(rec.info['CPX_INTERVALS'][2])
            dupINVdel= {
                SV_Name : {
                    "VAR":{
                        rec.chrom:{
                            "1":{
                                "rNameLeft": chrm1,
                                "rNameRight": chrm2,
                                "xLeft": end1-1,
                                "xRight": end2-1,
                                "directionLeft": "left",
                                "directionRight": "left"
                            },
                            "2" : {
                                "rNameLeft": chrm2,
                                "rNameRight": chrm1,
                                "xLeft": start2,
                                "xRight":rec.stop,
                                "directionLeft": "right",
                                "directionRight": "right",
                            }
                        }
                    }
                    
                }
            }
            countVariants['countdupINVdel'] = countVariants['countdupINVdel'] +1
            SV_Variants.update(dupINVdel)
            svNames[rec.id] = SV_Name

        ########## FOR CPX_TYPE: delINVdup
        elif rec.info['CPX_TYPE']== 'delINVdup':
            SV_Name=getName(rec)
            chrm1,start1, end1=SplitInterval(rec.info['CPX_INTERVALS'][0])
            chrm2,start2, end2=SplitInterval(rec.info['CPX_INTERVALS'][1])
            chrm3,start3, end3=SplitInterval(rec.info['CPX_INTERVALS'][2])
            delINVdup= {
                SV_Name : {
                    "VAR":{
                        rec.chrom:{
                            "1":{
                                "rNameLeft": chrm1,
                                "rNameRight": chrm2,
                                "xLeft": rec.pos-1,
                                "xRight": end2-1,
                                "directionLeft": "left",
                                "directionRight": "left"
                            },
                            "2" : {
                                "rNameLeft": chrm2,
                                "rNameRight": chrm1,
                                "xLeft": start2,
                                "xRight":start3,
                                "directionLeft": "right",
                                "directionRight": "right",
                            }
                        }
                    }
                    
                }
            }
            countVariants['countdelINVdup'] = countVariants['countdelINVdup'] +1
            SV_Variants.update(delINVdup)
            svNames[rec.id] = SV_Name

        ########## FOR CPX_TYPE: INVdup
        elif rec.info['CPX_TYPE']== 'INVdup':
            SV_Name=getName(rec)
            chrm1,start1, end1=SplitInterval(rec.info['CPX_INTERVALS'][0])
            chrm2,start2, end2=SplitInterval(rec.info['CPX_INTERVALS'][1])
            INVdup= {
                SV_Name : {
                    "VAR":{
                        rec.chrom:{
                            "1":{
                                "rNameLeft": chrm1,
                                "rNameRight": chrm1,
                                "xLeft": rec.pos-1,
                                "xRight": end1-1,
                                "directionLeft": "left",
                                "directionRight": "left"
                            },
                            "2" : {
                                "rNameLeft": chrm1,
                                "rNameRight": chrm1,
                                "xLeft": start1,
                                "xRight":start2,
                                "directionLeft": "right",
                                "directionRight": "right",
                            }
                        }
                    }
                    
                }
            }
            countVariants['countINVdup'] = countVariants['countINVdup'] +1
            SV_Variants.update(INVdup)
            svNames[rec.id] = SV_Name

        ########## FOR CPX_TYPE: dupINV
        elif rec.info['CPX_TYPE']== 'dupINV':
            SV_Name=getName(rec)
            chrm1,start1, end1=SplitInterval(rec.info['CPX_INTERVALS'][0])
            chrm2,start2, end2=SplitInterval(rec.info['CPX_INTERVALS'][1])
            dupINV= {
                SV_Name : {
                    "VAR":{
                        rec.chrom:{
                            "1":{
                                "rNameLeft": chrm1,
                                "rNameRight": chrm2,
                                "xLeft": end1-1,
                                "xRight": end2-1,
                                "directionLeft": "left",
                                "directionRight": "left"
                            },
                            "2" : {
                                "rNameLeft": chrm2,
                                "rNameRight": chrm1,
                                "xLeft": start2,
                                "xRight":rec.stop,
                                "directionLeft": "right",
                                "directionRight": "right",
                            }
                        }
                    }
                    
                }
            }
            countVariants['countdupINV'] = countVariants['countdupINV'] +1
            SV_Variants.update(dupINV)
            svNames[rec.id] = SV_Name

        ########## FOR CPX_TYPE: INS_iDEL
        elif rec.info['CPX_TYPE']== 'INS_iDEL':
            countVariants['countINS_iDEL'] = countVariants['countINS_iDEL'] +1
            

    elif rec.alts[0] == '<INS:ME:ALU>':
        countVariants['countALUINS'] = countVariants['countALUINS'] +1
        #getMaximumAF(rec)

    elif rec.alts[0] == '<INS:ME:LINE1>':
        countVariants['countLINEINS'] = countVariants['countLINEINS']+1

    elif rec.alts[0] == '<DEL:ME:ALU>':
        countVariants['countALUDEL'] = countVariants['countALUDEL'] +1
        SV_Name= str(getNameCanonical(rec))+"_ME"
        DEL= {
            SV_Name : {
                "VAR":{
                    rec.chrom:{
                        "1":{
                            "rNameLeft": rec.chrom,
                            "rNameRight": rec.chrom,
                            "xLeft": rec.pos,
                            "xRight": rec.stop,
                            "directionLeft": "left",
                            "directionRight": "right"
                        }
                    }
                }
                
            }
        }
        SV_Variants.update(DEL) 
        #getMaximumAF(rec)


    elif rec.alts[0] == '<DEL:ME:LINE1>':
        countVariants['countLINEDEL'] = countVariants['countLINEDEL'] +1
        SV_Name= str(getNameCanonical(rec))+"_ME"
        DEL= {
            SV_Name : {
                "VAR":{
                    rec.chrom:{
                        "1":{
                            "rNameLeft": rec.chrom,
                            "rNameRight": rec.chrom,
                            "xLeft": rec.pos,
                            "xRight": rec.stop,
                            "directionLeft": "left",
                            "directionRight": "right"
                        }
                    }
                }
                
            }
        }
        SV_Variants.update(DEL) 
        #getMaximumAF(rec)

    elif rec.alts[0] == '<DEL:ME:SVA>':
        countVariants['countSVADEL'] = countVariants['countSVADEL'] +1
        SV_Name= str(getNameCanonical(rec))+"_ME"
        DEL= {
            SV_Name : {
                "VAR":{
                    rec.chrom:{
                        "1":{
                            "rNameLeft": rec.chrom,
                            "rNameRight": rec.chrom,
                            "xLeft": rec.pos,
                            "xRight": rec.stop,
                            "directionLeft": "left",
                            "directionRight": "right"
                        }
                    }
                }
                
            }
        }
        SV_Variants.update(DEL) 
        #getMaximumAF(rec)

    elif rec.alts[0] == '<DEL:ME:HERVK>':
        countVariants['countHERVKDEL'] = countVariants['countHERVKDEL'] +1
        SV_Name= str(getNameCanonical(rec))+"_ME"
        DEL= {
            SV_Name : {
                "VAR":{
                    rec.chrom:{
                        "1":{
                            "rNameLeft": rec.chrom,
                            "rNameRight": rec.chrom,
                            "xLeft": rec.pos,
                            "xRight": rec.stop,
                            "directionLeft": "left",
                            "directionRight": "right"
                        }
                    }
                }
                
            }
        }
        SV_Variants.update(DEL) 
        #getMaximumAF(rec)

    elif rec.alts[0] == '<CNV>':
        countVariants['countCNV'] =countVariants['countCNV'] +1
        ##getMaximumAF(rec)

    else:
        #print( rec.chrom, rec.pos, rec.alts[0])
        SV_Name= str(rec.chrom)+'_'+str(rec.pos)+ '_'+ str(rec.stop)
        var = {
            SV_Name:{
                "chrom":rec.chrom,
                "pos":rec.pos,
                "type":rec.alts[0]
            }
        }
        #restvariants.append( rec.chrom, rec.pos, rec.alts[0])
        restvariants.update(var)







#### MAIN

if len(argv) != 3:
    print('Usage: ' + argv[0] + ' <vcffile.vcf> <outfile.json>')
    exit(1)

## Get the command line arguments
filepath = argv[1]
outpath = argv[2]   #### im Moment nch nicht verwendet
SV_Variants ={}
svNames = {}
AF_freqs =[['Name','fin_AF','fin_AN','nfe_AF','nfe_AN','eas_AF','eas_AN','afr_AF','afr_AN']]

countVariants={
    'countBND' :0,
    'countALUINS':0,
    'countLINEINS' :0,
    'countSVAINS' : 0,
    'countALUDEL':0,
    'countLINEDEL' :0,
    'countSVADEL' : 0,
    'countHERVKDEL' : 0,

    'countINS' : 0,
    'countDEL' : 0,
    'countDUP' : 0,
    'countCNV' : 0,
    'countINV' : 0,
    'countCTX' : 0,

    'countdupINVdup' : 0,
    'countdelINVdel' : 0,
    'countdelINVdup' : 0,
    'countdupINVdel' : 0,
    'countdelINV' : 0,
    'countINVdel' : 0,
    'countdDUP_iDEL' : 0,
    'countdupINV' : 0,
    'countINVdup' : 0,
    'countdDUP' : 0,
    'countINS_iDEL' :0
}



max_AF_afr_dupINVdup = 0
max_AF = {
    'afr':{
        'Population':'afr',
        'chr': '',
        'dupINVdup': [0,0],
        'delINVdel': [0,0],
        'delINVdup' : [0,0],
        'dupINVdel' : [0,0],
        'delINV' : [0,0],
        'INVdel' : [0,0],
        'dDUP_iDEL' : [0,0],
        'INS_iDEL':[0,0],
        'dupINV' : [0,0],
        'INVdup' : [0,0],
        'dDUP': [0,0],
        'DEL': [0,0],
        'INS':[0,0],
        'INV':[0,0],
        'CTX':[0,0],
        'CNV':[0,0],
        'BND':[0,0],
        'DUP':[0,0]
        },
    'eas':{
        'Population':'afr',
        'chr': '',
        'dupINVdup': [0,0],
        'delINVdel': [0,0],
        'delINVdup' : [0,0],
        'dupINVdel' : [0,0],
        'delINV' : [0,0],
        'INVdel' : [0,0],
        'dDUP_iDEL' : [0,0],
        'INS_iDEL':[0,0],
        'dupINV' : [0,0],
        'INVdup' : [0,0],
        'dDUP': [0,0],
        'DEL': [0,0],
        'INS':[0,0],
        'INV':[0,0],
        'CTX':[0,0],
        'CNV':[0,0],
        'BND':[0,0],
        'DUP':[0,0]
        },
    'fin/nfe':{
        'Population':'afr',
        'chr': '',
        'dupINVdup': [0,0],
        'delINVdel': [0,0],
        'delINVdup' : [0,0],
        'dupINVdel' : [0,0],
        'delINV' : [0,0],
        'INVdel' : [0,0],
        'dDUP_iDEL' : [0,0],
        'INS_iDEL':[0,0],
        'dupINV' : [0,0],
        'INVdup' : [0,0],
        'dDUP': [0,0],
        'DEL': [0,0],
        'INS':[0,0],
        'INV':[0,0],
        'CTX':[0,0],
        'CNV':[0,0],
        'BND':[0,0],
        'DUP':[0,0]
        }
    }

restvariants = {}


window = 'chr5:28932400- 28932900'
window = 'chr5:148173475-148173476'  ###### SPINK14
window = 'chr1:187495690- 187497640' ###### sv 1
window = 'chr2:124294170-124294190'  ###### sv 2
window = 'chr3:80013710-80013730'    ###### sv 3
window = 'chr5:28932500-28935600'    ###### sv 4
window = 'chr9:62806329-62806340'    ###### sv 5
window = 'chr9:86539640-86539645'    ###### sv 6
window = 'chr10:57497180-57497195'   ###### sv 7
window = 'chr2:61473500-61473700'    ###### sv 8
window = 'chr9:105055050-105055067'  ###### sv 9
window = 'chr12:77990000-77990298'   ###### sv 10
window = 'chr2:10685910-10685940'    ###### NA 2
window = 'chr3:95746900-95746980'    ###### NA 3
window = 'chr5:79753700-79753770'    ###### NA 4
window = 'chr8:99145580-99145590'    ###### NA 5
window = 'chr10:125508653-125508680' ###### NA 6
#window = 'chr17:5691360-5691390'     ###### NA 7




samplelist=[]
if filepath.endswith('.vcf.gz') or filepath.endswith('.vcf'):
    samplelist = [filepath]
if filepath.endswith('.txt'):
    file = open(filepath, 'r')
    for line in file: 
        #print(line)
        samplelist.append(line.strip())

for name in samplelist:
    print("go through file:", name)
    vcf_in = VariantFile(name)  # auto-detect input format
    for rec in vcf_in.fetch():
        ##### optional Filter:
        #if rec.alts[0] == '<CPX>' and ((rec.info['fin_AF'][0] or rec.info['nfe_AF'][0] or rec.info['fin_AF'][0] or rec.info['nfe_AF'][0]) > 0.9):
    #for rec in vcf_in.fetch(region=window):
        chromosome= rec.chrom
        getMaximumAF(rec)
        
        if rec.alts[0] == '<CPX>':
            getAFfreq(rec)
            createJSON(rec)

#max_AF['fin/nfe']['chr'] = chromosome
#max_AF['eas']['chr'] = chromosome
#max_AF['afr']['chr'] = chromosome

chromosome= filepath
prefix=outpath[:-5] # assumes output file has '.json' filetype

with open(prefix + "_name_assignments.tsv", 'w') as outFile:
    for k in svNames:
        outFile.write(k + "\t" + svNames[k] + "\n")

json_object = json.dumps(SV_Variants, indent=4)
with open(outpath, "w") as outfile:
    outfile.write(json_object)

#### Variants that are not categorized (hopefully none)
json_object2 = json.dumps(restvariants, indent=4)
with open(prefix+"_restvariants.json", "w") as outfile:
    outfile.write(json_object2)

json_object3 = json.dumps(max_AF, indent=4)
with open(prefix+"_maximum_AF.json", "w") as outfile:
    outfile.write(json_object3)

########### FREQUENCY
fields =['Variant_Type', 'frequency']
rows =[
    ['Breakends',countVariants['countBND'] ],
    ['INS:ME:ALU',countVariants['countALUINS']],
    ['INS:ME:LINE1',countVariants['countLINEINS'] ],
    ['INS:ME:SVA',countVariants['countSVAINS'] ],
    ['DEL:ME:ALU',countVariants['countALUDEL']],
    ['DEL:ME:LINE1',countVariants['countLINEDEL'] ],
    ['DEL:ME:SVA',countVariants['countSVADEL']],
    ['DEL:ME:HERVK',countVariants['countHERVKDEL'] ],
    ['INS', countVariants['countINS']],#countINS = 0
    ['DEL',countVariants['countDEL'] ],
    ['DUP',countVariants['countDUP'] ],
    ['CNV', countVariants['countCNV']],#countCNV = 0
    ['INV',countVariants['countINV'] ],
    ['CTX',countVariants['countCTX'] ],

    ['dupINVdup',countVariants['countdupINVdup'] ],
    ['delINVdel',countVariants['countdelINVdel'] ],
    ['delINVdup',countVariants['countdelINVdup'] ],
    ['dupINVdel',countVariants['countdupINVdel'] ],
    ['delINV ',countVariants['countdelINV'] ],
    ['INVdel',countVariants['countINVdel'] ],
    ['dDUP_iDEL',countVariants['countdDUP_iDEL'] ],
    ['dupINV',countVariants['countdupINV'] ],
    ['INVdup',countVariants['countINVdup'] ],
    ['dDUP',countVariants['countdDUP'] ],
    ['INS_iDEL', countVariants['countINS_iDEL']]
]
#print("INS",countVariants['countINS'], "CNV", countVariants['countCNV'])
#print(max_AF)
print("INS_iDEL: ", countVariants['countINS_iDEL'])
        
with open(prefix+"_frequency.csv", 'w') as f:
     
    # using csv.writer method from CSV package
    write = csv.writer(f)
     
    write.writerow(fields)
    write.writerows(rows)


with open(prefix+'_max_AF.csv', 'w') as f:  # You will need 'wb' mode in Python 2.x
    w = csv.DictWriter(f, max_AF['afr'].keys())
    w.writeheader()
    w.writerow(max_AF['afr'])
    w = csv.DictWriter(f, max_AF['eas'].keys())
    w.writerow(max_AF['eas'])
    w = csv.DictWriter(f, max_AF['fin/nfe'].keys())
    w.writerow(max_AF['fin/nfe'])
      
with open(prefix+'_AF.csv','w') as f:
    write = csv.writer(f)
    write.writerows(AF_freqs)
