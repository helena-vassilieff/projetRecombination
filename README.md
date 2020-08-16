# projetRecombination
create a directory : data/family/Romane/
in this directory put the hotspot file and the linkphase.in repertory (in wich we can find chr_x repertory)
scripts overlap.py and alpha.py are int the code/ repertory with Recombinations_analyses.py (ESSENTIAL to run properly overlap.py and alpha.py)

to run the both scripts execute this command (hotspot file and linkphase.in repertory are needed):
python3 code/overlap.py (or alpha.py) -hot data/family/Romane/hotspot_600k.txt -lphase data/family/Romane/linkphase.in/ 

you can precise option normally by default

* the resolution :
python3 code/overlap.py -hot data/family/Romane/hotspot_600k.txt -lphase data/family/Romane/linkphase.in/ -res 40000

* the perturbation : 

python3 code/overlap.py -hot data/family/Romane/hotspot_600k.txt -lphase data/family/Romane/linkphase.in/ -per 80000

* the chromosome 

python3 code/Overlap.py -hot data/family/Romane/hotspot_600k.txt -lphase data/family/Romane/linkphase.in/ -c 26

* the permutations : 

python3 code/alpha.py -hot data/family/Romane/hotspot_600k.txt -lphase data/family/Romane/linkphase.in/ -perm 500
