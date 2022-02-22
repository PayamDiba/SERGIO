import numpy as np
import pandas as pd

class grnParam(object):
    def __init__(self, k_act = [1,5], k_rep = [-5,-1], n = 2, decay = 0.8):
        try:
            len(n)
        except:
            n = [n] * 2

        try:
            len(decay)
            raise ValueError("decay rates should be the same for all genes")
        except:
            pass

        self.d_ = {'k+':k_act,
                   'k-':k_rep,
                   'n':n,
                   'decay':decay}

    def __call__(self):
        return self.d_

def getInterName(reg,tar):
    ret = ""
    for r in reg:
        ret += str(r.name_) +'-'
    return ret + str(tar.name_)

def parameterize_grn(df, param):
    param = param()

    ret = []
    names = ['index','reg','coop','tar','k','n','h','reg_decay','tar_decay']
    for _,r in df.iterrows():
        curr = [r.index, r.reg, r.coop, r.tar]
        kRange = param['k+'] if np.random.uniform() < 0.75 else param['k-']
        k = np.random.uniform(low = kRange[0], high = kRange[1])
        n = np.random.uniform(low = param['n'][0], high = param['n'][1])
        curr += [k,n,None,param['decay'],param['decay']]
        ret.append(curr)

    ret = pd.DataFrame(ret)
    ret.columns = names
    return ret
