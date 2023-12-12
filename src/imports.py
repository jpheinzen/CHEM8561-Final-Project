import numpy as np
import matplotlib.pyplot as plt
import time



def importU():
    txtIE = "Importing Energy"
    print(txtIE, end="")
    file = open("C:/Users/jphei/Desktop/U.txt", "r")
    # txt = fo.read()

    # print(type(txt))
    # print(txt)
    # txt2 = txt.split()
    # print(len(txt))
    # print(txt2)
    # num = np.array(txt2,dtype=float)
    # print(num)
    i = 0

    sizeMat = file.readline().strip().split()

    sizeMat = [int(x) for x in sizeMat]         # Turn list of 3 strings into ints
    # arrSize = np.prod(sizeMat[:2])
    # print(arrSize)
    arrSize = sizeMat[0]//sizeMat[2]+sizeMat[1]*sizeMat[2]
    # print(arrSize)
    # print(sizeMat[0]//sizeMat[2])
    # print(sizeMat[0]//sizeMat[2]+sizeMat[1]*sizeMat[2])
    arr = np.zeros((arrSize,2))

    for line in file:
        txt = line.strip().split()
        # arr[i,0] = np.array(txt[0],dtype=int)
        arr[i,:] = np.array(txt,dtype=float)
        i += 1

    file.close()

    # Plotting
    resid = abs(arr[1:,1] - arr[:-1,1])

    time.sleep(1)
    print('\b'*len(txtIE), end="")
    print('Imported Energy ')
    return arr, resid


def importPts(fast = False):
    txtIP = "Importing Points: "
    print(txtIP, end="")
    file = open("C:/Users/jphei/Desktop/pts.txt", "r")
    
    i = -1
    j = 0

    sizeMat = file.readline().strip().split()   # Read the first line

    sizeMat = [int(x) for x in sizeMat]         # Turn list of 3 strings into ints
    nPts = sizeMat[1]
    arrSize = (sizeMat[0]//sizeMat[2] + 1)
    # print(arrSize)

    arr = np.zeros((arrSize,nPts,3))
    # iter, V, dp, dv
    stats = np.zeros((arrSize,4))
    
    
    fmt = ""
    if not fast:
        for line in file:
            txt = line.strip().split()
            # print(txt)
            txt = np.array(txt,dtype=float)
            if np.size(txt) == 4:
                i += 1
                stats[i,:] = txt
                j = 0
                # print(i)
            elif np.size(txt) == 3:
                arr[i,j,:] = txt
                j +=1
            
            if (i % (arrSize//10) == 0):
                erase = '\b'*len(fmt)
                fmt = "%i/%i" % (i+1,arrSize)
                print("%s%s" % (erase, fmt), end="")


    file.close()
    print('\b'*(len(fmt)+len(txtIP)),end="")
    print("Imported Points  ","%s"%(" "*(2*round(np.log10(arrSize))+2+1)))

    # arr *= 1e9  # to get to nanometers

    return arr, stats

if __name__ == "__main__":
    # importU()
    importPts()

    print((1005001//100 +1)* (108) + 1)

    print(round(np.log10(10051)))






