#! /usr/bin/env python

import numpy as np


def DipTestSig(xpdf, nboot):
    """Dip test with significance
    """

    (dip, xlow, xup, ifault, gcm, lcm, mn, mj) = DipTest(xpdf)

    bootDip = []
    for i in range(nboot):
        unifpdf = np.sort(np.random.rand(xpdf.shape[0]))
        bootDip.append(DipTest(unifpdf)[0])
    bootDip = np.array(bootDip)
    p = np.sum(np.less(dip, bootDip)) / float(nboot)

    return (dip, p, xlow, xup)


def DipTest(xpdf):
    """Hartigan's dip test. 

    This is a copy 
    """
    x = np.sort(xpdf)
    N = x.shape[0]
    mn = np.zeros(x.shape, dtype=x.dtype)
    mj = np.zeros(x.shape, dtype=x.dtype)
    lcm = np.zeros(x.shape, dtype=x.dtype)
    gcm = np.zeros(x.shape)

    # check if N is one
    if N == 1:
        xl = x[0]
        xu = x[0]
        dip = 0.0
        ifault = 2
        print('\nHartigansDipTest.    InputError :  ifault=%s' + str(ifault))
        return (dip, xl, xu, ifault, gcm, lcm, mn, mj)

    #check for case 1<N<4 or all identical values
    if N <= 4 or x[N - 1] == x[0]:
        xl = x[0]
        xu = x[0]
        dip = 0.0
        ifault = 4
        print('\nHartigansDipTest.    InputError :  ifault=%s' + str(ifault))
        return (dip, xl, xu, ifault, gcm, lcm, mn, mj)

    #check if x is perfectly unimodal
    xsign = -np.sign(np.diff(np.diff(xpdf)))
    posi = np.greater(xsign, 0.0)
    negi = np.less(xsign, 0.0)
    if np.sum(posi) == 0 or np.sum(negi) == 0 or np.sum(np.less(posi, np.min(negi))) == N:
        #A unimodal function is its own best unimodal approximation, 
        #with a zero corresponding dip
        xl = x[0]
        xu = x[N - 1]
        dip = 0.0
        ifault = 5
        print('\nHartigansDipTest.    InputError :  ifault=%s' + str(ifault))
        return (dip, xl, xu, ifault, gcm, lcm, mn, mj)

    # LOW  contains the index of the current estimate of the lower 
    # end of the modal interval
    # HIGH contains the index of the current estimate of the upper 
    # end of the modal interval
    fn = N
    low = 1
    high = N
    dip = 1. / fn
    xl = x[low - 1]
    xu = x[high - 1]

    # establish the indices over which combination is necessary 
    # for the convex minorant fit
    mn[0] = 1.
    for j in range(2, N + 1):
        mn[j - 1] = j - 1

        mnj = mn[j - 1]
        mnmnj = mn[mnj - 1]
        a = mnj - mnmnj
        b = j - mnj
        while not ((mnj == 1) or (x[j - 1] - x[mnj - 1]) * a < (x[mnj - 1] - x[mnmnj - 1]) * b):
            mn[j - 1] = mnmnj
            mnj = mn[j - 1]
            mnmnj = mn[mnj - 1]
            a = mnj - mnmnj
            b = j - mnj

    # establish the indices over which combination is necessary 
    # for the concave majorant fit
    mj[N - 1] = N
    na = N - 1
    for jk in range(1, na + 1):
        k = N - jk
        mj[k - 1] = k + 1

        mjk = mj[k - 1]
        mjmjk = mj[mjk - 1]
        a = mjk - mjmjk
        b = k - mjk
        while not ( (mjk == N) or (x[k - 1] - x[mjk - 1]) * a < (x[mjk - 1] - x[mjmjk - 1]) * b):
            mj[k - 1] = mjmjk
            mjk = mj[k - 1]
            mjmjk = mj[mjk - 1]
            a = mjk - mjmjk
            b = k - mjk

    itarate_flag = True

    while itarate_flag:
        ic = 1
        gcm[0] = high
        igcm1 = gcm[ic - 1]
        ic += 1
        gcm[ic - 1] = mn[igcm1 - 1]

        while gcm[ic - 1] > low:
            igcm1 = gcm[ic - 1]
            ic += 1
            gcm[ic - 1] = mn[igcm1 - 1]

        icx = ic

        # collect the change points for the LCM from LOW to HIGH
        ic = 1
        lcm[0] = low
        lcm1 = lcm[ic - 1]
        ic += 1
        lcm[ic - 1] = mj[lcm1 - 1]
        while lcm[ic - 1] < high:
            lcm1 = lcm[ic - 1]
            ic += 1
            lcm[ic - 1] = mj[lcm1 - 1]

        icv = ic

        # ICX, IX, IG are counters for the convex minorant
        # ICV, IV, IH are counters for the concave majorant
        ig = icx
        ih = icv

        # find the largest distance greater than 'DIP' 
        # between the GCM and the LCM from low to high

        ix = icx - 1
        iv = 2
        d = 0.0

        if not (icx != 2 or icv != 2):
            d = 1. / fn
        else:
            iterate_BP50 = True

            while iterate_BP50:
                igcmx = gcm[ix - 1]
                lcmiv = lcm[iv - 1]
                if not (igcmx > lcmiv):
                    # if the next point of either the GCM or 
                    # LCM is from the LCM then calculate distance 
                    #here OTHERWISE, GOTO BREAK POINT 55

                    lcmiv1 = lcm[iv - 1 - 1]
                    a = lcmiv - lcmiv1
                    b = igcmx - lcmiv1 - 1
                    dx = (x[igcmx - 1] - x[lcmiv1 - 1]) * a / (fn * (x[lcmiv - 1] - x[lcmiv1 - 1])) - b / fn
                    ix -= 1
                    if dx < d:
                        goto60 = True
                    else:
                        d = dx
                        ig = ix + 1
                        ih = iv
                        goto60 = True
                else:
                    # if the next point of either the GCM or 
                    # LCM is from the GCM then calculate distance 
                    # here CODE BREAK POINT 55
                    lcmiv = lcm[iv - 1]
                    igcm = gcm[ix - 1]
                    igcm1 = gcm[ix + 1 - 1]
                    a = lcmiv - igcm1 + 1
                    b = igcm - igcm1
                    dx = a / fn - ((x[lcmiv - 1] - x[igcm1 - 1]) * b) / (fn * (x[igcm - 1] - x[igcm1 - 1]))
                    iv += 1

                    if not dx < d:
                        d = dx
                        ig = ix + 1
                        ih = iv - 1

                    goto60 = True

                if goto60:
                    if ix < 1: ix = 1
                    if iv > icv: iv = icv
                    iterate_BP50 = gcm[ix - 1] != lcm[iv - 1]

        itarate_flag = not d < dip
        if itarate_flag:
            # if itarate_flag is true, then continue 
            # calculations and the great iteration cycle
            #if itarate_flag is NOT true, then stop 
            # calculations here, and break out of 
            #great iteration cycle to BREAK POINT 100

            # calculate the DIPs for the current LOW and HIGH

            #the DIP for the convex minorant
            dl = 0.
            if ig != icx:
                icxa = icx - 1
                for j in range(ig, icxa + 1):
                    temp = 1. / fn
                    jb = gcm[j + 1 - 1]
                    je = gcm[j - 1]
                    if not (je - jb <= 1):
                        if not (x[je - 1] == x[jb - 1]):
                            a = je - jb
                            const = a / (fn * (x[je - 1] - x[jb - 1]))
                            for jr in range(int(jb), int(je + 1)):
                                b = jr - jb + 1
                                t = b / fn - (x[jr - 1] - x[jb - 1]) * const
                                if t > temp:
                                    temp = t
                    if dl < temp:
                        dl = temp

            du = 0.
            if not (ih == icv):
                icva = icv - 1
                for k in range(ih, icva + 1):
                    temp = 1. / fn
                    kb = lcm[k - 1]
                    ke = lcm[k + 1 - 1]
                    if not (ke - kb <= 1):
                        if not (x[ke - 1] == x[kb - 1]):
                            a = ke - kb
                            const = a / (fn * (x[ke - 1] - x[kb - 1]))
                            for kr in range(int(kb), int(ke + 1)):
                                b = kr - kb - 1
                                t = (x[kr - 1] - x[kb - 1]) * const - b / fn
                                if t > temp: temp = t
                    if du < temp: du = temp

            dipnew = dl
            if du > dl: dipnew = du
            if dip < dipnew: dip = dipnew
            low = gcm[ig - 1]
            high = lcm[ih - 1]
            #end if itarate_flag

    dip *= 0.5

    return dip

