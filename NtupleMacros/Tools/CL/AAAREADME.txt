libCMS2CL is a library that allows you to use CMS2 selection functions from
TTree::Draw(), and eventually also from TTree::Scan().  (CMS2CL stands for CMS2
from the command line.)

Instructions:

- run make in this directory; this should give you libCMS2CL.so
- unfortunately, there are some dependencies on other ROOT libraries (such as
  libThread.so and libNet.so, because obviously those are sine qua nons for
  looping over trees...); you need to load all the libraries found in setup.C
- once you have all the libraries loaded, try out something like: 
    Events->Draw("hyp_p4.M()>>h3(150, 0, 150)", "hyp_p4.M2() > 0 && hyp_lt_id * hyp_ll_id < 0 && abs(hyp_lt_id) == 11 && electronSelection_cand01(hyp_lt_index) && abs(hyp_ll_id) == 11 && electronSelection_cand01(hyp_ll_index)")
  where Events could be the Z --> ee CMS2 ntuple, for example

