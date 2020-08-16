#!/usr/bin/env python3
import os
from collections import defaultdict
from collections import OrderedDict
from decimal import Decimal
import numpy as np
import matplotlib.pyplot as plt
from statistics import mean, median, mode, stdev
import getopt
from scipy import stats
from scipy.stats import norm
from datetime import datetime
import pybedtools as pybed
from pybedtools import BedTool
import argparse
import sys
import copy 
import cProfile
import re
import random



UNKNOWN = "0"
MALE = "1"
FEMALE = "2"

# Utility functions
def eprint(*args, **kwargs):
    print(*args,  file=sys.stderr, **kwargs)


## 1/ Une classe pour représenter les intervalles
#######

class RecombInterval:
    """
    An interval in which is detected a recombination (crossover) event.
    Attributes
    ----------
    - chrom : str 'chrX'
              chromosome
    - left :  int
              left coordinate
    - right : int
              right coordinate
    - meiosis : {UNKNOWN, MALE, FEMALE}
              Sex of the individual that produced the recombination event
    """
    def __init__(self, chrom, left, right, meiosis = UNKNOWN):
        self.chrom = chrom
        self.left = left
        self.right = right
        self.meiosis = meiosis
        self._is_overlapping = None
        self._pc_overlap = None
        
        if self.left > self.right :
            raise Exception("left marker is bigger than right marker")
        if not isinstance(self.chrom,str ) :
            raise Exception("your chromosome should be chrx ")
	
    @property
    def size(self):
        """ 
        Size of the interval (in bp)
        """
        return self.right - self.left

    def to_bed(self):
        interval_list = list(map(str, [self.chrom, self.left, self.right]))
        return pybed.create_interval_from_list(interval_list)


    def random_shift(self, sd):
        """Creates a new interval with position randomly shifted
        The new interval is a copy of the object with coordinates shifted
        by a random deviation r ~ N(0, sd).
        Parameters
        ----------
        sd : float, positive
            Standard deviation of the Gaussian distribution
        
        Returns
        -------
        obj:`RecombInterval`
            The shifted interval 
        """
        assert sd >= 0
        r = norm.rvs(loc = 0, scale = sd)
        while self.left+r < 0 or self.right+r < 0:
            r = norm.rvs(loc = 0, scale = sd)
        new_interval = copy.deepcopy(self)
        new_interval.left = int(self.left+r)
        new_interval.right = int(self.right+r)
        return new_interval

    def is_overlapping(self, others=None):
        """Does the interval overlaps one in `others` ?
        
        Parameters
        ----------
        others : list :obj:`RecombIntervals`
                 List of intervals
        
        Returns
        -------
        bool
             True if the interval overlaps one in `others`
        """
        if others is None:
            assert self._is_overlapping is not None
        else:
            self._is_overlapping = RecombInterval.percent_overlap([self], others)==100
        return self._is_overlapping

    def pc_random_overlap(self, others=None, sd=2e5, nsim=1000):
        """Proportion (%) of overlap by chance with `others`
        Parameters
        ----------
        others : list of :obj:`RecombIntervals`
            List of intervals
        sd : float, optional
            Standard deviation for random shift
        nsim : int, optional
            Number of shifts to perform   
        Returns
        -------
        float
            Calculated proportion
        """
        ## TODO implementer renvoi _pc_overlap si others == None
        if others is None:
            assert self._is_overlapping is not None
        else:
            rnd_perturbations = [ self.random_shift(sd) for i in range(nsim)]
            self._pc_overlap = RecombInterval.percent_overlap(rnd_perturbations, others)
        return self._pc_overlap

    @staticmethod
    def intervals_to_pybed(intervals):
        """create a list of bed intervals from intervals objects
        Parameters
        ----------
        intervals : list of :obj:`RecombIntervals`
       
        Returns
        -------
        list of : bed : intervals
            
        """
        pybed_intervals = []
        for inter in intervals:
            pybed_intervals.append(inter.to_bed())
        return pybed.BedTool(pybed_intervals).sort()

    @staticmethod
    def percent_overlap(interv_list_a, interv_list_b):
        """
        This function computes the proportion of intervals in interv_list_a 
        that overlap one or more intervals in interv_list_b
        Parameters
        ----------
        interv_list_a : list of :bed :intervals
           example :  List of intervals
        interv_list_a : list of :bed :intervals
           example :  List of hotspots
        Returns
        -------
        float
            Calculated proportion of overlap
        """
        pybed_a = RecombInterval.intervals_to_pybed(interv_list_a)
        pybed_b = RecombInterval.intervals_to_pybed(interv_list_b)
        overlap = pybed_a.intersect(pybed_b,u=True)#, sort =True)

        return float((overlap.count()/pybed_a.count())*100)





# 2/ Une classe pour représentées les marqueurs
########
class Marker:
    """
    An marker used to detect a recombination (crossover) event.
    Attributes:
    ----------
    - chrom : str 'chrX'
              chromosome
    - identity :  str 'x'
              marker identity
    - name : str 'oar_xxx'
              marker name
    - position : int 
             true position in bp
    """   
    def __init__(self, chrom, identity, name, position):
        self.chrom = chrom
        self.identity = identity
        self.name = name
        self.position = position




# 4/ Une classe pour les individus
#######

class Individual:
    """
    individuals is compsed by an individual
    id its sex and a list of intervals objects
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

    def __init__(self, individual, sex, mom=None, dad=None, rexIntervals = 0 ):
        self.indiv = individual
        self.sex = sex
        self.mother = mom
        self.father = dad 
        self.meioseIntervals = []
        self.intervals = []
        self.offsprings = []
        self._num_max_resIntervals = rexIntervals



    def find_num_maxIntervals(self, max_res = 3e4):
        """
        calculate the length of set intervals for intervals
        under max_res
        ----------
        max_res: int: by default 30 000 pb
        """
        if self._num_max_resIntervals == 0 :
            max_resIntervals = [inter for inter in self.meioseIntervals if  inter.size < max_res]     
            self._num_max_resIntervals = len(max_resIntervals)
        return self._num_max_resIntervals


    def populate_intervals(self, interval):
        """
        full the intervals list
        ----------
        interval : obj: `RecombIntervals`
        """
        self.intervals.append(interval)
    
    def populate_MeioseIntervals(self, interval):
        """
        full the MeioseIntervals list
            ----------
        interval : obj: `RecombIntervals`
        """
        self.meioseIntervals.append(interval)


    def populate_offsprings(self, offspring):
        """
        full the offsprings list
        ----------
        offspring : str 'X'  
        """

        self.offsprings.append(offspring)



# 5/ Fonctions de lecture des données
########
def alpha_estimate(intervals, ldhot=None, npoints=101):
    """
    this function will try to find for one individual
    its best fit for alpha 
    Parameters
    ----------
    intervals : list of :obj:`RecombIntervals`
        List of intervals
    ldhot : list of :obj:`RecombIntervals` from a hotspot file
        Standard deviation for random shift
    npoints : int, optional
        density of the grid of tested alpha  
    size_max : int, optional
        size of intervals analysed
    Returns
    -------
    str 
        Calculated alpha, 95% confidence interval based on likelihood within two units of the maximum log likelihood
    """

    deltas = [ inter.is_overlapping(ldhot) for inter in intervals ]
    P_overlap_chance = [ inter.pc_random_overlap(ldhot) / 100 for inter in intervals ]

    alpha_vals = np.linspace(0,1,num=npoints)
    loglik = np.zeros(len(alpha_vals),dtype=np.float64)
    for i,alpha in enumerate(alpha_vals):
        for idx,inter in enumerate(intervals):
            P_over_hotspots = alpha + (1 - alpha) * P_overlap_chance[idx]
            likelihood = deltas[idx] * P_over_hotspots + (1 - deltas[idx]) * (1 - P_over_hotspots)
            with np.errstate(divide='ignore'): 
                loglik[i] += np.log(likelihood)
    idx_max = np.argmax(loglik)
    lmax = np.max(loglik)
    alpha_max = alpha_vals[idx_max]

    ci_alphas = np.argwhere( loglik >= lmax-2)
    ci_lo = min(alpha_vals[ci_alphas])
    ci_hi = max(alpha_vals[ci_alphas])
    
    return alpha_max, (ci_lo,ci_hi), lmax, alpha_vals, loglik



def read_ldhotspots( hotspot_file):
    """
    reads in a hotspot_file
    and returns a list of RecombInterval objects
    Parameters
        ----------
        hotspot_file
           
        -------
        list
            obj : `RecombIntervals`
        
    """
    hotspots = []
    with open (hotspot_file,"r") as fin:
        for line in fin:
            fields=line.rstrip().split()
            if fields[4] == 'TRUE':
                chrom = "chr%s" % fields[0]
                start, end = fields[1:3]
                start = int(start)
                end = int(end)
                hotspot = RecombInterval(chrom, start, end)
                hotspots.append(hotspot)

    return hotspots
    


def read_typMap(typMap_file, marker_dict, chrom):
    """
    reads typage.map (a tsv file with 
    the chromosome the marker name and the position
    and return a dictionnary with marker id
    as key and an object with marker informatioon
    as value
    Parameters
        ----------
        typage.map : tsv file
            contains the chromosome number, the marker name and the marker position
        marker_dict : dictionnary
            key : marker name
            value : marker identity
        chrom : str 'chrX'
              
        Returns
        -------
        dictionnary
            key : marker indentity
            value : obj : `Marker`
        
    """
    markers = defaultdict()
    with open(typMap_file,'r') as fin:
        for line in  fin:
            fields=line.rstrip().split()
            mrk_name = fields[1]
            if mrk_name in  marker_dict:
                mrk_id = marker_dict[mrk_name]
                pos = fields[3]
                markers[mrk_id] = Marker(chrom,mrk_id, mrk_name, pos)
                                       
    return markers
                


def read_map(map_file):
    """ 
   this function reads in map and 
    returns a dictionnary with
    marker name as key and marker id 
    as value
    Parameters
        ----------
        map_file : tsv file
            contains the marker identity and the marker name 
                      
        Returns
        -------
        dictionnary
            key : marker name
            value : marker identity
        
    """
    marker_name_id = defaultdict()
    with open (map_file,'r') as fin:
        for line in  fin:
            fields=line.rstrip().split()
            marker_name_id[fields[1]] = fields[0]

    return marker_name_id



def sexs(sex_indiv):
   if sex_indiv == FEMALE:
       return 'FEMELLE'
   if sex_indiv == MALE:
       return 'MALE'



def read_nrec(nrec_file):
    """
    This function reads in nrec_file a tsv file
    and populate the Individual class
    Parameters
    ---------
        nrec_file : tsv file
           contains offsprings parents and parents sex
        Returns
    -------
        dictionnary:
        key : individual identity (whether they are parents or offsprings) :str 'X'
        value : obj : `Individual`
        """
    individuals_id = defaultdict()
    parents = []
    with open (nrec_file,'r') as fin:
        for line in  fin:
            
            fields=line.rstrip().split()
            offspring = fields[0]
            parent = fields[1]
            sexParent = fields[2]
            
            if offspring not in individuals_id :
                one_individual = Individual(offspring, '0')
                if sexParent == FEMALE:
                    one_individual.mother = parent

                if  sexParent == MALE:
                    one_individual.father = parent 
                
                individuals_id[offspring] = one_individual
            else:
                if sexParent == FEMALE:
                     individuals_id[offspring].mother = parent
                if sexParent == MALE:
                     individuals_id[offspring].father = parent

            if parent not in  individuals_id and parent not in parents:
                the_parent = Individual(parent, sexParent)
                the_parent.populate_offsprings(offspring)
                parents.append(parent)
                individuals_id[parent]=the_parent
            else:
                individuals_id[parent].populate_offsprings(offspring)
                if individuals_id[parent].sex =='0':
                    individuals_id[parent].sex = sexParent
            

    return individuals_id




def read_fam(fam_file):
    """
    This function reads in fam_file a tsv file
    and populate the Individual class
    Parameters
    ---------
        fam_file : tsv file
           contains offsprings parents and parents sex
        Returns
    -------
        dictionnary:
        key : individual identity (whether they are parents or offsprings) :str 'X'
        value : obj : `Individual`
        """
    individuals_duo_id = defaultdict()
    parents = []
    with open (fam_file,'r') as fin:
        for line in  fin:
            
            fields=line.rstrip().split()
            offspring = fields[1]
            fath = fields[2]
            moth = fields[3]
            par=(fath, moth)
            if offspring != fields[0] :
                if offspring not in individuals_duo_id and offspring not in parents:
                    one_individual = Individual(offspring, '0')
                    if moth != "0":
                        one_individual.mother = moth
                    if fath != "0":
                        one_individual.father = fath 
                
                    individuals_duo_id[offspring] = one_individual
               
                for i in range(0,2):
                
                    if par[i] not in  individuals_duo_id and par[i] not in parents and par[i] != "0":
                        the_parent = Individual(par[i], str(int(i)+1))
                        the_parent.populate_offsprings(offspring)
                        parents.append(par[i])
                        individuals_duo_id[par[i]]=the_parent
                    elif par[i] != "0":
                        individuals_duo_id[par[i]].populate_offsprings(offspring)
                        if individuals_duo_id[par[i]].sex =='0':
                            individuals_duo_id[par[i]].sex = str(int(i)+1)
            

    return individuals_duo_id


def get_linkphase_filename(lkphase_dir, chrom, extension):
    return os.path.join(lkphase_dir, chrom, extension)

def read_intervals(linkphase_dir, chromosome):
    """
    This function reads in 'recombination_hmm" files 
    located in 'linkphase_dir' and subdirectories and returns a list 
    of RecombInterval objects and populate Individuals objects with intervals 
    
    Parameters
    ---------
         linkphase_dir: repertory
           contains chromosomes repertory that contains nrec_file, map_file , typage.map and recombinations_hmm
 
         chromosome : map
            chromosome number
        Returns
    -------
        dictionnary:
        key : individual identity (whether they are parents or offsprings) :str 'X'
        value : obj : `Individual`
        
        list: of obj : `RecombIntervals`
    """
    observed_intervals = []
    individuals_objects = defaultdict()
    parents_objects = defaultdict()
    chromosome2 = copy.deepcopy(chromosome)
    first_chromosome = min(chromosome2)
    nrec_file =  get_linkphase_filename(linkphase_dir, 'chr_'+first_chromosome+"/", "nrec_hmm.txt" )    
    individuals_information = read_nrec(nrec_file)

    for chrom in chromosome:
        chrom_directory='chr_'+chrom+"/"
        chrom_name='chr'+chrom
        map_file = get_linkphase_filename(linkphase_dir, chrom_directory, "map" )
        typageMap_file = get_linkphase_filename(linkphase_dir, chrom_directory, "typage.map" )
        recomb_file = get_linkphase_filename(linkphase_dir, chrom_directory, "recombinations_hmm" )
        marker_information1 = read_map(map_file)
        marker_information = read_typMap(typageMap_file, marker_information1, chrom_name)
        with open (recomb_file, 'r') as fin:
                for line in  fin:
                    fields=line.rstrip().split()
                    if marker_information[fields[2]].chrom == chrom_name and \
                    marker_information[fields[3]].chrom == chrom_name  :  
                        left_marker = int(marker_information[fields[2]].position)
                        right_marker = int(marker_information[fields[3]].position)
                        parent = fields[1]
                        offspring = fields[0]
                        one_interval =  RecombInterval(chrom_name, \
                        left_marker,right_marker, individuals_information[parent].sex)
                        if offspring in individuals_information :
                            individuals_information[offspring].populate_intervals(one_interval)   
                        if parent in individuals_information :
                            individuals_information[parent].populate_MeioseIntervals(one_interval)
                        
                        observed_intervals.append(one_interval)
    
    
    return observed_intervals, individuals_information



def get_duohmm_filename(duohmm_dir, chrom, extension):
    return os.path.join(duohmm_dir, chrom, extension)

def read_duo_intervals(duohmm_dir, chromosome, probability):
    """
    This function reads in 'chX-recombinations.txt" files 
    located in 'shapeit' and subdirectories and returns a list 
    of RecombInterval objects and populate Individuals objects with intervals 
    
    Parameters
    ---------
         shapeit: repertory
           contains chromosomes repertory that contains chrX-recombinations.txt, original chrX.fam file
 
         chromosome : map
            chromosome number
        Returns
    -------
        dictionnary:
        key : individual identity (whether they are parents or offsprings) :str 'X'
        value : obj : `Individual`
        
        list: of obj : `RecombIntervals`
    """
    observed_intervals = []
    individuals_objects = defaultdict()
    parents_objects = defaultdict()
    chromosome2 = copy.deepcopy(chromosome)
    first_chromosome = min(chromosome2)
    fam_file =  get_duohmm_filename(duohmm_dir, 'chr'+str(first_chromosome)+'/', "RID_genea_merged_12_chr"+str(first_chromosome)+".fam" )    
    individuals_information = read_fam(fam_file)
    for chrom in chromosome:
        
        chrom_directory='chr'+chrom+"/"
        chrom_name='chr'+chrom
        recomb_file = get_duohmm_filename(duohmm_dir, chrom_directory, str(chrom_name)+"-recombinations.txt")
        x_li = 0 
        with open (recomb_file, 'r') as fin:
                for line in  fin:
                    if x_li > 0:
                        fields=line.rstrip().split()
                        if float(fields[4]) >= probability :  
                            left_marker = int(fields[2])
                            right_marker = int(fields[3])
                            parent = fields[1]
                            offspring = fields[0]
                            if len(individuals_information[parent].offsprings)> 1:
                                one_interval =  RecombInterval(chrom_name, \
                            left_marker,right_marker, individuals_information[parent].sex)
                            
                                if offspring in individuals_information :
                                    individuals_information[offspring].populate_intervals(one_interval)   
                                if parent in individuals_information  :

                                    individuals_information[parent].populate_MeioseIntervals(one_interval)                    

                                observed_intervals.append(one_interval)
                    x_li += 1
                        
    return observed_intervals, individuals_information


# 5/
#######



def scramble(recomb_intervals, per):
    """ 
    This functions scrambles the position of each interval in
    the 'recomb_intervals' list using permutation 'per' and returns
    a list of new scrambled intervals
    Parameters
    ---------
         recomb_intervals: list : of obj : 'RecombIntervals'
 
         per : int 'X'
           a given perturbation in bp
        Returns
    -------
        list: of obj : deepcopy of : `RecombIntervals`
    """
    scrambled_intervals = []
    for inter in recomb_intervals:
        new_inter = inter.random_shift(per) 
        scrambled_intervals.append(new_inter)

    return scrambled_intervals


def catch_alpha(indiv,sexe, alpha,CI,alphas,loglik,sample_size,ligne_alpha,ligne_loglik):
    """ 
    This functions construct lines for tsv file with alpha information
    and confidendence informations
    Parameters
    ---------
         indiv : str (group of individual) or int (individual id)
 
         sexe : int '0' or '1'
           
         intervals : list of : intervals
 
         hot : list of : intervals
 
         ligne_alpha : str
         ligne_loglik : str
        Returns
    -------
        str : alpha with confidence intervals
              alpha with log likelihood
    """

    ligne1=''
    ligne2=''
    ligne1+=str(indiv)+"\t"+str(sexe)+"\t"+str(alpha*100)+"\t"+str(float(CI[0]*100))+"\t"+str(float(CI[1]*100))+"\t"+str(sample_size)+"\n"
    for i,lik in enumerate(loglik):
        ligne2+=(str(indiv)+'\t'+str(lik)+'\t'+str(alphas[i]*100)+'\n')

    ligne_alpha += ligne1
    ligne_loglik += ligne2  
    return ligne_alpha, ligne_loglik


def indiv_likelihood(individuals, sample, perm):
    """ 
    This functions calculate a log likelihood for a group of 
    individuals by constructing a false set of intervals
    for a given individual this set will have the length 
    of the true set of intervals under the maximum ofresolution



    Parameters
    ---------
         individuals : list of : obj : individuals
 
         sample : list of : obj : intervals
           
         
        Returns
    -------
        float : sum of maximum log likelihood of 
            individuals

    """
    start = 0
    cumulative_log = 0
    ligne = ''
    for i,indiv in enumerate(individuals):
        end = start + indiv.find_num_maxIntervals()
        new_meioses_intervals = sample[start:end]
        alpha,CI,indiv_loglik,alphas, loglik = alpha_estimate(new_meioses_intervals)           
        start = end 
        cumulative_log += indiv_loglik
        ligne += (str(indiv.indiv)+'\t'+str(indiv.sex)+'\t'+str(perm)+'\t'+str(alpha)+'\n')
    return cumulative_log, ligne

