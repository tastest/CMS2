import FWCore.ParameterSet.Config as cms

from CMS2.NtupleMaker.tcmetMaker_cfi import *

tcmetSequence = cms.Sequence( tcmetMaker )

