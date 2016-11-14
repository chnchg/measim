import numpy as np

# load spike data
def load_spk(fn):
    spkf = open(fn,'r')
    dat = np.loadtxt(spkf)
    spkf.close()
    return dat

# calculate sliding histogram
def slide_hist(dat,win):
    hst = [dat[0]-win]
    cnt = [0]
    i0 = 0
    i1 = 0
    ct = 0
    while i0<len(dat):
        if i1<len(dat) and dat[i1]<dat[i0]+win:
            i1 += 1
            t = dat[i1-1]-win/2
            hst.extend([t,t])
            cnt.extend([cnt[-1],i1-i0])
        else:
            i0 += 1
            t = dat[i0-1]+win/2
            hst.extend([t,t])
            cnt.extend([cnt[-1],i1-i0])
            if i0>=len(dat): break
    return hst,cnt

# calculate sliding average
def slide_average(ts,vs,win):
    hst = [ts[0]-win]
    ave = [0]
    i0 = 0
    i1 = 0
    ct = 0
    while i0<len(ts):
        if i1<len(ts) and ts[i1]<ts[i0]+win:
            i1 += 1
            t = ts[i1-1]-win/2
            hst.extend([t,t])
            ave.extend([ave[-1],np.mean(vs[i0:i1])])
        else:
            i0 += 1
            t = ts[i0-1]+win/2
            hst.extend([t,t])
            ave.extend([ave[-1],np.mean(vs[i0:i1])])
            if i0>=len(ts): break
    return hst,ave

# burst detection
def spk_bursts(dat):
    spkt = [x[0]/1000 for x in dat]
    twin = 0.1
    hst,cnt = slide_hist(spkt,twin)
    srt = [x/twin/1000 for x in cnt] # rate in spike per millisecond
    
    mx = max(srt)
    # some parameters for detection
    epsilon = mx * 0.04
    delta = mx * 0.2
    tau = 1.0
    b = False
    bs = []
    dr = []
    ton = hst[0]
    toff = hst[0]
    psrt = 0

    for i in range(len(srt)):
        if psrt < epsilon:
            if srt[i] >= epsilon: ton = hst[i]
        elif srt[i] < epsilon: toff = hst[i]
        if not b:
            if srt[i] > delta:
                # bursting starts
                bs.append(ton)
                b = True
        else:
            if srt[i] < epsilon and hst[i] - toff > tau:
                dr.append(toff)
                b = False
        psrt = srt[i]
    if b: dr.append(hst[-1]) # end last burst
    return bs,dr

# burst detection, assuming second for time unit
def burst_detect(spk,**kwargs):
    twin = kwargs.pop('twin',0.1)
    epsilon = kwargs.pop('epsilon',0.04)
    delta = kwargs.pop('delta',0.2)
    tau = kwargs.pop('tau',1.0)
    h,c = slide_hist(spk,twin)
    mx = max(c)
    epm = mx*epsilon
    dem = mx*delta
    b = False
    bs = []
    dr = []
    ton = h[0]
    toff = h[0]
    ps = 0
    for i in range(len(c)):
        if ps<epm:
            if c[i]>=epm: ton = h[i]
        elif c[i]<epm: toff = h[i]
        if not b:
            if c[i]>dem:
                bs.append(ton)
                b = True
        else:
            if c[i]<epm and h[i]-toff>tau:
                dr.append(toff)
                b = False
        ps = c[i]
    if b: dr.append(h[-1]) # last burst
    return bs,dr

# burst detection from histogram, assuming second for time unit
def hst_burst(hst,cnt,**kwargs):
    epsilon = kwargs.pop('epsilon',0.04)
    delta = kwargs.pop('delta',0.2)
    tau = kwargs.pop('tau',1.0)
    mx = max(cnt)
    epm = mx*epsilon
    dem = mx*delta
    b = False
    bs = []
    dr = []
    ton = hst[0]
    toff = hst[0]
    ps = 0
    for i in range(len(cnt)):
        if ps<epm:
            if cnt[i]>=epm: ton = hst[i]
        elif cnt[i]<epm: toff = hst[i]
        if not b:
            if cnt[i]>dem:
                bs.append(ton)
                b = True
        else:
            if cnt[i]<epm and hst[i]-toff>tau:
                dr.append(toff)
                b = False
        ps = cnt[i]
    if b: dr.append(hst[-1]) # last burst
    return bs,dr

# symmertric peak detection
def hst_peaks(hst,cnt,**kwargs):
    epsilon = kwargs.pop('epsilon',0.04)
    eps = max(cnt)*epsilon
    vmn = cnt[0]
    vmx = cnt[0]
    pni = 0
    pkt = []
    pkh = []
    for i in range(len(hst)):
        if vmn < eps: vmn = eps
        if cnt[i] < vmn:
            vmn = cnt[i]
            pni = i # where the minimum is
            vmx = vmn
        elif cnt[i] > vmx:
            vmx = cnt[i]
            pxi = i # where the maximum is
        elif vmx > 2 * eps and vmx > 2 * vmn and cnt[i] < vmx / 2:
            #if vmn <= eps:
            pkt.append(hst[pxi])
            pkh.append(cnt[pxi])
            vmn = cnt[i]
            vmx = cnt[i]
    return pkt,pkh

# symmertric peak detection
def hst_peaks1(hst,cnt,**kwargs):
    epsilon = kwargs.pop('epsilon',0.04)
    delta = kwargs.pop('delta',0.1)
    eps = max(cnt)*epsilon
    dlt = max(cnt)*delta
    vmn = cnt[0]
    vmx = cnt[0]
    pni = 0
    pkt = []
    pkh = []
    pks = []
    pke = []
    pcnt = 100000
    ppct = 100000
    lcmm = 0 # track previous local minimum
    for i in range(len(hst)):
        if pcnt < ppct and cnt[i] >= pcnt:
            lcmm = i-1
        ppct = pcnt
        pcnt = cnt[i]
        if vmn < eps:
            vmn = eps
            if len(pke) < len(pkt): pke.append(hst[i])
        if cnt[i] < vmn:
            vmn = cnt[i]
            if lcmm<0: pni = i # where the minimum is
            else: pni = lcmm
            vmx = vmn
        elif cnt[i] > vmx:
            vmx = cnt[i]
            pxi = i # where the maximum is
        elif vmx > 2 * eps and vmx > 2 * vmn and cnt[i] < vmx / 2:
            if len(pke) < len(pkt): pke.append(hst[pni])
            if cnt[pxi]>dlt:
                pkt.append(hst[pxi])
                pkh.append(cnt[pxi])
                pks.append(hst[pni])
            vmn = 100000
            vmx = cnt[i]
    if len(pke) < len(pkt): pke.append(hst[i])
    return pkt,pkh,pks,pke

# symmertric peak detection 2
def hst_peaks2(hst,cnt,**kwargs):
    epsilon = kwargs.pop('epsilon',0.04)
    delta = kwargs.pop('delta',0.1)
    eps = max(cnt)*epsilon
    dlt = max(cnt)*delta
    vmn = cnt[0]
    vmx = cnt[0]
    pni = 0
    pkt = []
    pkh = []
    pks = []
    pke = []
    pcnt = 100000
    ppct = 100000
    lcmm = 0 # track previous local minimum
    for i in range(len(hst)):
        if pcnt < ppct and cnt[i] >= pcnt:
            lcmm = i-1
        ppct = pcnt
        pcnt = cnt[i]
        if vmn < eps:
            vmn = eps
            if len(pke) < len(pkt): pke.append(hst[i])
            pni = i
        if cnt[i] < vmn:
            vmn = cnt[i]
            if lcmm<0: pni = i # where the minimum is
            else: pni = lcmm
            vmx = vmn
        elif cnt[i] > vmx:
            vmx = cnt[i]
            pxi = i # where the maximum is
        elif vmx > 2 * eps and vmx > 2 * vmn and cnt[i] < vmx / 2:
            if len(pke) < len(pkt): pke.append(hst[pni])
            if cnt[pxi]>dlt:
                pkt.append(hst[pxi])
                pkh.append(cnt[pxi])
                pks.append(hst[pni])
            vmn = 100000
            vmx = cnt[i]
    if len(pke) < len(pkt): pke.append(hst[i])
    return pkt,pkh,pks,pke

def calc_burst(fn):
    global dat
    global hst
    global srt
    global midx
    global bsm
    global drm
    global mhst
    global msrt
    global bs
    global dr
    dat = load_spk(fn)
    bs, dr = spk_bursts(dat)
    durs = np.array(dr)-bs
    midx = np.argsort(durs)[len(durs)/2]
    hst,cnt = slide_hist([x[0]/1000 for x in dat],0.01)
    srt = [x/0.01/1000 for x in cnt]
    bsm = bs[midx]
    drm = dr[midx]
    bsi = 0
    while bsi+1 < len(hst) and hst[bsi] < bsm: bsi = bsi + 1
    dri = bsi
    while dri+1 < len(hst) and hst[dri] < drm: dri = dri + 1
    [mhst,msrt] = [hst[bsi:dri],srt[bsi:dri]]
    if len(durs): return np.mean(durs),np.std(durs)
    return 0,0

