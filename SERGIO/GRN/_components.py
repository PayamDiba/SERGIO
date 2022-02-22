from ._utils import getInterName
import numpy as np
from collections import defaultdict

class Gene(object):
    def __init__(self, name = None, decay = 0.8):
        self.name_ = name
        self.isMR_ = False
        self.level = None
        self.regs = []
        self.tars = []
        self.inInteractions = {} # 'reg names (seperated by comma if coop)' --> interactions
        self.prod_rates_ = None
        self.ss_conc_ = None
        self.sim_conc_ = defaultdict(list) # dict of cell_type_id to np.array of simulated conc
        self.decay_ = decay if decay else 0.8
        #self.outInteractions = {} # 'target names (seperated by comma if coop)' --> interactions
        #self.isTF_ = TF

    def _calc_prod(self, cTypes, regs_conc = 'ss'):
        if self.isMR_:
            return self.prod_rates_[cTypes] #np.arr

        ret = np.zeros(shape = (len(cTypes),))
        for i in self.inInteractions.keys():
            ret += self.inInteractions[i]._get_hill(cTypes, regs_conc = regs_conc)

        return ret

    def append_sim_conc(self, conc, cTypes):
        assert(len(conc) == len(cTypes))
        for cID, curr in zip(cTypes, conc):
            self.sim_conc_[cID] = np.append(self.sim_conc_[cID], max(0,curr)) #prevents negative expression

    def get_last_conc(self, cTypes):
        ret = [self.sim_conc_[ct][-1] for ct in cTypes]
        return np.array(ret)



#todo make a base interaction object
class SingleInteraction(object):
    def __init__(self, reg, tar, k = None, h = None, n = None):
        self.reg_ = reg #can be list if coop edge --> always be list
        self.tar_ = tar
        self.name_ = getInterName(reg,tar)
        self.k_ = k
        self.h_ = h
        self.n_ = n

    def _get_hill(self, cTypes, regs_conc = 'ss'):
        """
        cTypes : is list of cell type ids to compute hill for
        regs_conc (str): 'ss' (uses the steady state concentration of all regs),
        'sim' (uses the last simulated concentration of all regs)
        """
        if len(self.reg_) > 1:
            raise ValueError("Coop interactions are not yet implemented")

        if regs_conc == 'ss':
            x = self.reg_[0].ss_conc_[cTypes] #np.arr
        else:
            x = self.reg_[0].get_last_conc(cTypes) #np.arr

        ret = np.power(x, self.n_) / (np.power(x, self.n_) + np.power(self.h_, self.n_))
        if self.k_ > 0:
            return self.k_ * ret #np.arr of size No. provided cell types
        else:
            return np.abs(self.k_) * (1-ret) #np.arr of size No. provided cell types
