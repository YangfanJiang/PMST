import numpy as np

def loadFile(datafile):
    xEls = []
    ycs = []
    yrs = []
    ype = []
    ytotal = []
    with open(datafile) as f:
        lines = f.readlines()
        for l in lines[3:]:
            l = l.split('|')
            energy = float(l[0]) *1000 # convert Mev to Kev
            xEls.append(energy)
            csCo = float(l[2])
            ycs.append(csCo)
            rsCo = float(l[1])
            yrs.append(rsCo)
            peCo = float(l[3])
            ype.append(peCo)
            totalCoWithCo = float(l[4])
            ytotal.append(totalCoWithCo)
    return xEls, ycs, yrs, ype, ytotal
def interpolate(pEnergy, xEls, ycs, yrs, ype, ytotal):
    for i in range (0, len(xEls)-1):
        if pEnergy == xEls[i]:
            return ype[i], ycs[i], yrs[i], ytotal[i]
        if pEnergy > xEls[i] and  pEnergy< xEls[i+1]:
            ##for log interpolation
            # peCoLog = np.exp(np.log(pEnergy/xEls[i]) * (np.log(ype[i+1]/ype[i])/np.log(xEls[i+1]/xEls[i])) + np.log(ype[i]))
            # csCoLog = np.exp(
            #     np.log(pEnergy / xEls[i]) * (np.log(ycs[i + 1] / ycs[i]) / np.log(xEls[i + 1] / xEls[i])) + np.log(
            #         ycs[i]))
            # rsCoLog = np.exp(
            #     np.log(pEnergy / xEls[i]) * (np.log(yrs[i + 1] / yrs[i]) / np.log(xEls[i + 1] / xEls[i])) + np.log(
            #         yrs[i]))
            # totalCoLog = np.exp(
            #     np.log(pEnergy / xEls[i]) * (np.log(ytotal[i + 1] / ytotal[i]) / np.log(xEls[i + 1] / xEls[i])) + np.log(
            #         ytotal[i]))
            ##Linear interpolation
            peCo = ype[i] + (pEnergy - xEls[i]) * (ype[i+1] - ype[i]) / (xEls[i+1] - xEls[i])
            csCo = ycs[i] + (pEnergy - xEls[i]) * (ycs[i+1] - ycs[i]) / (xEls[i+1] - xEls[i])
            rsCo = yrs[i] + (pEnergy - xEls[i]) * (yrs[i+1] - yrs[i]) / (xEls[i+1] - xEls[i])
            totalCo = ytotal[i] + (pEnergy - xEls[i]) * (ytotal[i+1] - ytotal[i]) / (xEls[i+1] - xEls[i])
            return peCo, csCo, rsCo, totalCo
def interpolate1(pEnergy, datafile):
    xEls = []
    ycs = []
    yrs = []
    ype = []
    ytotal = []
    with open(datafile) as f:
        lines = f.readlines()
        for l in lines[3:]:
            l = l.split('|')
            energy = float(l[0]) *1000 # convert Mev to Kev
            xEls.append(energy)
            csCo = float(l[2])
            ycs.append(csCo)
            rsCo = float(l[1])
            yrs.append(rsCo)
            peCo = float(l[3])
            ype.append(peCo)
            totalCoWithCo = float(l[4])
            ytotal.append(totalCoWithCo)
    # print(ype, ycs, yrs, ytotal)

    # global xEls, ycs, yrs, ype, ytotal
    for i in range (0, len(xEls)-1):
        if pEnergy == xEls[i]:
            return ype[i], ycs[i], yrs[i], ytotal[i]
        if pEnergy > xEls[i] and  pEnergy< xEls[i+1]:
            ##for log interpolation
            # peCoLog = np.exp(np.log(pEnergy/xEls[i]) * (np.log(ype[i+1]/ype[i])/np.log(xEls[i+1]/xEls[i])) + np.log(ype[i]))
            # csCoLog = np.exp(
            #     np.log(pEnergy / xEls[i]) * (np.log(ycs[i + 1] / ycs[i]) / np.log(xEls[i + 1] / xEls[i])) + np.log(
            #         ycs[i]))
            # rsCoLog = np.exp(
            #     np.log(pEnergy / xEls[i]) * (np.log(yrs[i + 1] / yrs[i]) / np.log(xEls[i + 1] / xEls[i])) + np.log(
            #         yrs[i]))
            # totalCoLog = np.exp(
            #     np.log(pEnergy / xEls[i]) * (np.log(ytotal[i + 1] / ytotal[i]) / np.log(xEls[i + 1] / xEls[i])) + np.log(
            #         ytotal[i]))
            ##Linear interpolation
            peCo = ype[i] + (pEnergy - xEls[i]) * (ype[i+1] - ype[i]) / (xEls[i+1] - xEls[i])
            csCo = ycs[i] + (pEnergy - xEls[i]) * (ycs[i+1] - ycs[i]) / (xEls[i+1] - xEls[i])
            rsCo = yrs[i] + (pEnergy - xEls[i]) * (yrs[i+1] - yrs[i]) / (xEls[i+1] - xEls[i])
            totalCo = ytotal[i] + (pEnergy - xEls[i]) * (ytotal[i+1] - ytotal[i]) / (xEls[i+1] - xEls[i])
            return peCo, csCo, rsCo, totalCo
#
# re = interpolate(55)
# for co in re:
#     print(co)
