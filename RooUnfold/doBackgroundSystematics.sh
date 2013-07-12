#!/bin/bash

root -b -q macro_1DUnfolding_bkgsyst.C &> macro_1DUnfolding_bkgsyst.log &
root -b -q macro_2DUnfolding_bkgsyst.C &> macro_2DUnfolding_bkgsyst.log &
root -b -q macro_2DUnfolding_ttpt_bkgsyst.C &> macro_2DUnfolding_ttpt_bkgsyst.log &
root -b -q macro_2DUnfolding_ttrapidity2_bkgsyst.C &> macro_2DUnfolding_ttrapidity2_bkgsyst.log &
