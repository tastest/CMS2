#################################################################
                          Code Organization
#################################################################

CORE/ - contains code common to everyone. Don't change anything 
      here unless you are absolutely sure what you do and what 
      effects it will have on all users

NtupleTools/ - tools to make new loopers from scratch

data/ - ROOT files with corrections, fake rates etc

Tools/ - common tools


#################################################################
                      Ntuple Location and Naming
#################################################################

At all sites a user needs to set CMS2_NTUPLE_LOCATION environment 
variable that points to data location. At this location ntuples
are organized like that:

cms2-V00-04-00/merge_tW.root  
cms2-V00-04-00/merge_WW.root
cms2-V00-04-01/merge_DY.root  
cms2-V00-04-01/merge_ttbar.root  
cms2-V00-04-01/merge_Wjet.root
cms2-V00-05-00/merge_WZ.root  
cms2-V00-05-00/merge_ZZ.root

it's what we have right now, in future it will be cleaned up.


#################################################################
                      Analysis Code and User Code
#################################################################

All other subdirectories represent analysis and user code. It's up
to users how to keep their code clean, but here is a general 
guidance:
  * at the top level only final code should be visible
  * all intermediate hacks, tests and developments have to be
    in user sub-directories (sandboxes)
  * brief HOWTO/README should be available for the top level 
    analysis code so that anyone can reproduce the results



#################################################################
                      HOWTO Run New Looper                      
#################################################################

NOTE: the make files need to be reorganized. Below is just an
      example how to use what is available.

Set CMS2_NTUPLE_LOCATION environment variable that points to data 
location. 

  setenv CMS2_NTUPLE_LOCATION /data/tmp/

Push the button to get baseline WW result:
  
  make WW_Results.tbl

10 mins later, print results. 

  cat WW_Results.tbl

After deciphering it for a bit, you should get something like that:

      |       ww      |      wz      |      zz     |      wjets   |     dyee    |      dymm    |      dytt    |       ttbar   |      tw      |       total    |
| all | 498.4 +/- 8.7 | 13.6 +/- 1.4 | 4.2 +/- 0.7 | 45.5 +/- 9.5 | 2.0 +/- 1.4 | 13.4 +/- 4.3 | 22.8 +/- 5.6 | 115.4 +/- 7.1 | 47.6 +/- 2.6 | 762.9 +/- 16.7 |
| mm  |  89.2 +/- 3.7 |  3.4 +/- 0.7 | 1.7 +/- 0.4 |  1.0 +/- 1.0 | 0.0 +/- 0.0 | 11.1 +/- 4.1 |  0.9 +/- 0.9 |  31.7 +/- 3.7 | 11.2 +/- 1.3 | 150.3 +/-  6.9 |
| em  | 347.4 +/- 7.3 |  7.9 +/- 1.0 | 0.3 +/- 0.2 | 40.4 +/- 9.2 | 0.0 +/- 0.0 |  2.2 +/- 1.6 | 21.9 +/- 5.6 |  67.1 +/- 5.4 | 30.1 +/- 2.1 | 517.4 +/- 14.4 |
| ee  |  61.7 +/- 3.1 |  2.3 +/- 0.6 | 2.1 +/- 0.5 |  4.1 +/- 2.0 | 2.0 +/- 1.4 |  0.0 +/- 0.0 |  0.0 +/- 0.0 |  16.6 +/- 2.7 |  6.4 +/- 0.9 |  95.3 +/-  4.9 |

Output as of 11/18/2008:
  * feb selection
  * nTrk>2
  * MET->tcMET
  * el with calo iso - 0.92
  * trkJet veto
  * extra muon veto

