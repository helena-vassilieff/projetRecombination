#!/usr/bin/env python3
import os
from decimal import Decimal
import numpy as np
from statistics import mean, median, mode, stdev
from scipy import stats
from scipy.stats import norm
from datetime import datetime
import argparse
import sys
import TestClass as ov




class identity:
    """
    identity is compsed by an individual
    id its id after plink2linkphase.py its position in fam and fam.orig files
    Attributes:
    ----------
    - individual : str 'X'
              one individual
    - sex : {UNKNOWN, MALE, FEMALE}
              individual sex
    - mother : str 'X'
             
    - father : str 'X' 
             
    - meioseIntervals : list of : obj : 'RecombIntervals'
              list of intervals for individuals that are parents
    - intervals : list of : obj : 'RecombIntervals'
              list of intervals for individuals
    - offsprings : list of : offsprings : str 'X'
    """

    def __init__(self, id_fam, pos, id_fam_orig = None):
        self.id_fam = id_fam      
        self.position = pos 
        self.id_fam_orig = id_fam_orig



id_count_dico = {}
count_line = 0
with open("RID_genea_merged_12.fam", "r") as fin:
    for line in fin :
        fields = line.rstrip().split()
        id_count_dico[str(count_line)] = identity(fields[1],str(count_line))
        count_line += 1


link_dico = {}
duo_dico = {}

count_line = 0
with open("RID_genea_merged_12.fam.orig", "r") as fin:
    for line in fin :
        fields = line.rstrip().split()
        id_count_dico[str(count_line)].id_fam_orig = fields[1]
        link_dico[id_count_dico[str(count_line)].id_fam] = id_count_dico[str(count_line)]
        duo_dico[id_count_dico[str(count_line)].id_fam_orig] = id_count_dico[str(count_line)]

        count_line += 1




"""
def link2duo (link_indiv, link_list, duo_indiv,  duo_list):
    #print("link2duo")
    print("indiv sexe overlap CImax CImin intervals inter_shared")
    for indiv in link_indiv :
        if link_dico[link_indiv[indiv].indiv].id_fam_orig in duo_dico:
            if len(link_indiv[indiv].offsprings) != 0:  
                 max_res = 100000
                 max_duo_intervals = [inter for inter in duo_indiv[link_dico[link_indiv[indiv].indiv].id_fam_orig].meioseIntervals if len(duo_indiv[link_dico[link_indiv[indiv].indiv].id_fam_orig].offsprings) >1] #if inter.size < max_res]
                 max_resintervals = [ inter for inter in link_indiv[indiv].meioseIntervals if len(link_indiv[indiv].offsprings) >1 ]  #  if inter.size < max_res ]
                 if len(max_resintervals) != 0:
                     p_overlap = 0
                     p_overlap = ov.RecombInterval.percent_overlap(max_resintervals, max_duo_intervals)
                     n = len(max_resintervals)
                     n2 =len(max_duo_intervals)
                     p = p_overlap*0.01
                     ci = 1.96*np.sqrt((p*(1-p))/n)
                     CImax = p + ci
                     CImin = p - ci
                     #ligne_overlap += str(individuals[indiv].indiv)+"\t"+str(ov.sexs(individuals[indiv].sex))+"\t"+str(p_overlap)+"\t"+str(CImax*100)+"\t"+str(CImin*100)+"\t"+str(n)+"\n"
                     print(link_indiv[indiv].indiv," ",ov.sexs(link_indiv[indiv].sex)," ",p_overlap," ",CImax*100," ",CImin*100," ",n,n/100*p_overlap)




"""
def duo2link (duo_indiv, duo_list, link_indiv, link_list):
    #print("duo2link")
    print("indivD indivL sexe overlap CImax CImin intervals inter_shared")
    all_inter = [inter.size for inter in duo_list ]
    print(len(all_inter))
    for indiv in duo_indiv :
        if duo_indiv[indiv].indiv in duo_dico:
            if duo_dico[duo_indiv[indiv].indiv].id_fam in link_indiv:
                if len(duo_indiv[indiv].meioseIntervals) != 0:  
                     max_res = 100000
                     max_resintervals = [ inter for inter in duo_indiv[indiv].meioseIntervals if len(duo_indiv[indiv].offsprings) >1 ]#if inter.size < max_res ]
                     max_link_intervals = [inter for inter in link_indiv[duo_dico[duo_indiv[indiv].indiv].id_fam].meioseIntervals if len(link_indiv[duo_dico[duo_indiv[indiv].indiv].id_fam].offsprings) >1 and  len(link_indiv[duo_dico[duo_indiv[indiv].indiv].id_fam].meioseIntervals) != 0] #if inter.size < max_res]
                     if len(max_resintervals) != 0:
                         p_overlap = 0
                         p_overlap = ov.RecombInterval.percent_overlap(max_resintervals, max_link_intervals)
                         n = len(max_resintervals)
                         n2 =len(max_link_intervals)
                         p = p_overlap*0.01
                         ci = 1.96*np.sqrt((p*(1-p))/n)
                         CImax = p + ci
                         CImin = p - ci
                         #ligne_overlap += str(individuals[indiv].indiv)+"\t"+str(ov.sexs(individuals[indiv].sex))+"\t"+str(p_overlap)+"\t"+str(CImax*100)+"\t"+str(CImin*100)+"\t"+str(n)+"\n"
                         print(duo_indiv[indiv].indiv," ",duo_dico[duo_indiv[indiv].indiv].id_fam," ",ov.sexs(duo_indiv[indiv].sex)," ",p_overlap," ",CImax*100," ",CImin*100," ",n,n/100*p_overlap)












