#################################################################
                          Code Organization
#################################################################

CORE/ - contains code common to everyone. Don't change anything 
      here unless you are absolutely sure what you do and what 
      effects it will have on all users

NtupleTools/ - tools to make new loopers from scratch

Tools/ - common tools

Templates/ - a template directory for a simple dilepton
looper (to be extended by the user into a real analysis)

#################################################################

All other subdirectories represent analysis and user code. It's up
to users how to keep their code clean, but here is a general 
guidance:
  * at the top level only final code should be visible
  * all intermediate hacks, tests and developments have to be
    in user sub-directories (sandboxes)
  * brief HOWTO/README should be available for the top level 
    analysis code so that anyone can reproduce the results
    Generally, this is in form of another AAAREADME.txt file.

The following is a list if directories and their respective owners:

ConversionTrackStudies:  fkw/WA
- used for the met note on the 2009 data. Continues to be useful for a bit longer.

SSDilep:  fkw/SP
- used for same sign dilepton analysis since summer 2009





