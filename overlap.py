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
import Recombinations_analyses as ov
import duo2link as dl
import copy 

# Utility functions
def eprint(*args, **kwargs):
    print(*args,  file=sys.stderr, **kwargs)

UNKNOWN = "0"
MALE = "1"
FEMALE = "2"


def main(hotspot_file, link_phase_dir,duohmm_phase_dir, resolutions, perturbations, chromosomes, permutations, probability): 
    """this main will call function to construct intervals and
    analyses those intervals sampled by sex and resolution
    (and perturbations) to know their overlap percentage 
    with intervals from hotspots 
    Parameters
    ---------
    hotspot_file: tsv file
        contains intervals with positions 
 
    link_phase_dir : repertory
        contains repertory that contains files with intervals positions
    resolutions : list of : int, optional
        determine the size of analysed intervals
    perturbations :  list of : int, optional
        determine the shift of intervals 
    chromosomes : map of : str 'X' , optional
        chromosome id
    
        Returns
    -------
       float : 
          overlap percentages (according to intervals size, intervals perturbation)

       file : `repetOverlap.tsv`
           resolutions, perturbations , repetition until 100 , overlap percentage those parameters

       file : `overlap_caterpillar.tsv`
           individual, sex, overlap percentage, confidence interval, intervals number, offsprings number

    """
    ligne='resolution\tperturbation\trepetition\toverlapPercentage\n'
    ligne_overlap = 'indiv\tsexe\toverlap\tCImax\tCImin\tintervals\tdesc\n'   
    eprint("resolutions: %s" % ", ".join(map(str,resolutions)))

    eprint("perturbations: %s" % ", ".join(map(str,perturbations)))
    ldhotspots = ov.read_ldhotspots(hotspot_file)
    if link_phase_dir != None and duohmm_phase_dir != None :
        intervals_link , individuals_link= ov.read_intervals(link_phase_dir, chromosomes)
        intervals_duo, individuals_duo = ov.read_duo_intervals(duohmm_phase_dir, chromosomes, probability)  
        #dl.link2duo(individuals_link,intervals_link, individuals_duo,intervals_duo)
        dl.duo2link(individuals_duo, intervals_duo, individuals_link, intervals_link)
    else:
        if link_phase_dir != None and duohmm_phase_dir == None :    
            intervals , individuals = ov.read_intervals(link_phase_dir, chromosomes)
        if link_phase_dir == None and duohmm_phase_dir != None :
            intervals, individuals = ov.read_duo_intervals(duohmm_phase_dir, chromosomes, probability)  


        for indiv in individuals :

            if len(individuals[indiv].offsprings) != 1:  
                 max_res = int(max(resolutions))
                 max_resintervals = [ inter for inter in individuals[indiv].meioseIntervals if inter.size < max_res ]
                 p_overlap = 0
                 if len(max_resintervals) != 0 :
                     p_overlap = ov.RecombInterval.percent_overlap(max_resintervals,ldhotspots)
                     n = len(max_resintervals)
                     p = p_overlap*0.01
                     ci = 1.96*np.sqrt((p*(1-p))/n)
                     CImax = p + ci
                     CImin = p - ci

                     ligne_overlap += str(individuals[indiv].indiv)+"\t"+str(ov.sexs(individuals[indiv].sex))+"\t"+str(p_overlap)+"\t"+str(CImax*100)+"\t"+str(CImin*100)+"\t"+str(len(individuals[indiv].offsprings))+"\t"+str(n)+"\n"

       
        for res in  resolutions:
            res = int(res)
            resintervals = [ inter for inter in intervals if inter.size < res ]
            p_overlap = ov.RecombInterval.percent_overlap(resintervals, ldhotspots)
            print("%s %d %d %3.2f" % ( 0, res, 0, p_overlap))
            if res == max(resolutions):
                p = p_overlap*0.01
                n = len(resintervals)
                ci = 1.96*np.sqrt((p*(1-p))/n)
                CImax = p + ci
                CImin = p - ci
                ligne_overlap += ("overlap over\tall\t"+str(p_overlap)+"\t"+str(CImax*100)+"\t"+str(CImin*100)+"\t"+str(n)+"\n")

            if len(resintervals) != 0:
                for sex in ( MALE, FEMALE):
                    sexintervals = [ inter for inter in resintervals if inter.meiosis == sex ]
                    female_distribution = [ inter.size for inter in resintervals if inter.meiosis == FEMALE ]
                    male_distribution = [ inter.size for inter in resintervals if inter.meiosis == MALE ]
                    if len(sexintervals) > 0:
                        p_overlap = ov.RecombInterval.percent_overlap( sexintervals, ldhotspots)
                        print("%s %d %d %3.2f" % ( sex, res, 0, p_overlap))
                    if res == max(resolutions) : 
                        p = p_overlap*0.01
                        n = len(sexintervals)
                        ci = 1.96*np.sqrt((p*(1-p))/n)
                        CImax = p + ci
                        CImin = p - ci
                        ligne_overlap += ("overlap over"+str(sex)+"\t"+str(sex)+"\t"+str(p_overlap)+"\t"+str(CImax*100)+"\t"+str(CImin*100)+"\t"+str(n)+"\n")
                if len(female_distribution) > 0 and len(male_distribution) > 0:
                    print("%d %d %d %d" % (1, len(male_distribution), 2, len(female_distribution)))
                    kolmo=stats.ks_2samp(male_distribution,female_distribution)
                    print("%s %d %3.2f %3.2f" % ('k', res,kolmo[0],kolmo[1]))
                for perturb in perturbations:
                    p_overlap = 0
                    p_overlap2 = 0
                    for i in range(1,101):
                        scrambled_intervals = ov.scramble(resintervals, perturb)
                        p_overlap2 = ov.RecombInterval.percent_overlap(scrambled_intervals, ldhotspots) 
                        p_overlap += p_overlap2
                        ligne+=(str(res)+'\t'+str(perturb)+'\t'+str(i)+'\t'+str(p_overlap2)+'\n')
                    print("%s %d %d %3.2f" % ('s',res, perturb, p_overlap / 100))
        with open('repetOverlap.tsv','a') as f1 : 
            f1.write(ligne)
        with open('overlap_caterpillar.tsv','a') as f1 : 
            f1.write(ligne_overlap)




def parse_arguments():
    parser = argparse.ArgumentParser(
                        description='Detect overlaps between LD recombination '
                        'hotspots and familial recombination intervals.')
    parser.add_argument("-hot", "--hotspots", type=str,
                        required=True, help='The hotspot file')
    parser.add_argument("-lphase", "--link_phase_dir", type=str,
                        required=False, help='The linkphase dir')
    parser.add_argument("-duo", "--duohmm_phase_dir", type=str,
                        required=False, help='The duohmm dir')
    parser.add_argument("-res", "--resolutions", nargs='+',
                        required=False, help='Resolutions',
                        default=[5000, 10000, 20000, 30000])
    parser.add_argument("-per", "--perturbations", nargs='+',
                        required=False, help='Perturbations to perform',
                        default=[200000, 500000, 1000000])
    parser.add_argument("-c", "--chromosomes", nargs='+',
                        required=False, help='The chromosomes',
                        default=list(map(str, range(1, 27))))
    parser.add_argument("-perm", "--permutations", nargs='+',
                        required=False, help='Permutations to perform',
                        default = 10000)
    parser.add_argument("-prob", "--prob_recombination", nargs='+',
                        required=False, help='probability of recombinations',
                        default = 0.9)
    args = parser.parse_args( )
    return args


if __name__ == '__main__':
    args = parse_arguments()
    main(hotspot_file=args.hotspots,
         link_phase_dir=args.link_phase_dir,
         duohmm_phase_dir=args.duohmm_phase_dir,
         resolutions=args.resolutions,
         perturbations=args.perturbations,
         chromosomes=args.chromosomes,
         permutations=args.permutations,
         probability=args.prob_recombination
        )
