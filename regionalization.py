import pysal
import numpy as np
import random
from constants import *
import csv

def regionlization(text):
    print_seperator(text,"Reginalization",FALSE)
    w=pysal.open("province_queen.gal").read()
    f = pysal.open("GDP3.csv")
    pci = np.array([f.by_col[str(y)]for y in range(1978, 2012)])
    pci.shape
    np.random.seed(100)
    random.seed(10)
    r=pysal.Maxp(w,pci,floor=5,floor_variable=np.ones((31,1)),initial=110)
    region_values = "r.regions" + str(r.regions)
    writeCalculations(text,region_values,False)
    #print ("r.regions",r.regions)
    header_values = "f.header" + str(f.header)
    writeCalculations(text,header_values,False)
    # print ("f.header",f.header)
    names=f.by_col('NAME')
    names=np.array(names)
    print (names)
    for region in r.regions:
        ids=map(int,region)
        ids=[i-1 for i in ids]
        print (names[ids])
    r.inference()
    pvalue = "r.pvalue" + str(r.pvalue)
    writeCalculations(text,pvalue,False)
    #print ("r.pvalue",r.pvalue)
    I=r.regions
    I2=[]
    for i in range(len(I)):
            I2.append([i+1]*len(I[i]))
    res=[]
    for i in range(len(I)):
            e1=I[i]
            e2=I2[i]
            r=zip(e1,e2)
            res=res+r
    res_values = "res" + str(res)
    writeCalculations(text, res_values, TRUE)
    #print "res", res
    # csvfile = open('C:/Python27/mydata/region2.csv', 'w')
    # csvwriter = csv.writer(csvfile)
    # csvwriter.writerow(['ID','Number'])
    # for e in res:
    #     csvwriter.writerow(list(e))
    # map(csvwriter.writerow, res)




