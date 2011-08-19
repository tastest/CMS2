#!/usr/bin/env python

import sys,os,json

def convertRange(runs):
    result = {}
    for run in runs.keys():
        tmp = []
        for lumirange in runs[run]:
            tmp.extend(range(lumirange[0],lumirange[1]+1))
        result[run] = tmp
    return result
    
def convertArray(runs):
    result = {}
    for run in runs.keys():
        lumis = runs[run]
        lumis.sort()
        endRange = lumis[0]
        startRange = lumis[0]
        lumiRanges = []
        for lumi in range(1,len(lumis)) :
            if endRange+1 == lumis[lumi] :
                endRange = lumis[lumi]
            else :
                lumiRanges.append([startRange,endRange])
                endRange = lumis[lumi]
                startRange = lumis[lumi]
        lumiRanges.append([startRange,endRange])
        result[run] = lumiRanges
    return result

input_arrays = []
for argv in sys.argv[1:]:
    input_arrays.append(convertRange(json.load(open(argv))))
    
result_array = {}
for array in input_arrays:
    for run in array.keys():
        if run in result_array.keys():
            for lumi in array[run]:
                if lumi not in result_array[run]:
                    result_array[run].append(lumi)
        else :
            result_array[run] = array[run]            

result_range = convertArray(result_array)

output_handle = open("merged.json",'w')
json.dump(result_range,output_handle)
output_handle.close()


