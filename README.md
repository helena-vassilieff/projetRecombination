# projetRecombination
create a directory : data/family/Romane/
in this directory put the hotspot file and the linkphase.in (or duohmm) repertory (in wich we can find chr_x repertory with crossing over locations)
scripts overlap.py, alpha.py and duo2link.py are in the code/ repertory with Recombinations_analyses.py (ESSENTIAL to run properly overlap.py and alpha.py)

alpha.py estimate hotspots usage rates with maximum likelihood for population and individuals .
overlap.py estimate naives hotspots usage rates.
duo2link.py compare crossing over from linkphase and crossing over from duohmm (python3 code/oduo2link.py -hot data/family/Romane/hotspot_600k.txt -lphase data/family/Romane/linkphase.in/ -duo data/family/Romane/duohmm/)

to run the scripts execute this command (hotspot file and linkphase.in or duohmm repertories are needed):
python3 code/overlap.py (or alpha.py ) -hot data/family/Romane/hotspot_600k.txt -lphase data/family/Romane/linkphase.in/ 

you can precise option normally by default

* the resolution :
python3 code/overlap.py -hot data/family/Romane/hotspot_600k.txt -lphase data/family/Romane/linkphase.in/ -res 40000

* the perturbation : 

python3 code/overlap.py -hot data/family/Romane/hotspot_600k.txt -lphase data/family/Romane/linkphase.in/ -per 80000

* the chromosome :

python3 code/Overlap.py -hot data/family/Romane/hotspot_600k.txt -lphase data/family/Romane/linkphase.in/ -c 26

* the permutations : 

python3 code/alpha.py -hot data/family/Romane/hotspot_600k.txt -lphase data/family/Romane/linkphase.in/ -perm 500
