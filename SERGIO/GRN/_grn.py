"""
@author: Payam Dibaeinia
"""
from ._utils import getInterName
import numpy as np
import pandas as pd
from collections import defaultdict
from ._components import SingleInteraction

class GRN(object):
    def __init__(self, **kwargs):
        self.attr_ = self._build_attr(**kwargs)
        self.mrs_found_ = False

    def _build_attr(self,**kwargs):
        ret = {
            'genes' : {},
            'interactions' : {},
            'mrs' : set(),
            'level2gene' : defaultdict(list),
            'maxLevel' : -1,
        }
        for k in kwargs:
            ret[k] = {} if kwargs[k] == 'dict' else kwargs[k]

        return ret

    def setMRs(self):
        if not self.mrs_found_:
            mrs = set(self.attr_['genes'].keys()) - self.attr_['mrs']
            self.attr_['mrs'] = mrs
            self.mrs_found_ = True

    def get_mrs(self):
        return list(self.attr_['mrs'])

    def _add_gene(self, gene):
        if gene.name_ not in self.attr_['genes'].keys():
            self.attr_['genes'][gene.name_] = gene
            #self.attr_['mrs'].add(gene.name_)
            return gene
        else:
            #self.attr_['mrs'].add(gene.name_)
            return self.attr_['genes'][gene.name_]

    def add_interaction(self, interaction):
        self._check_inter_endPoints(interaction)
        if interaction.name_ not in self.attr_['interactions'].keys():
            self.attr_['interactions'][interaction.name_] = interaction


    def _check_inter_endPoints(self, interaction):
        interaction.tar_ = self._add_gene(interaction.tar_)
        #self.attr_['mrs'].remove(interaction.tar_.name_)
        self.attr_['mrs'].add(interaction.tar_.name_)
        interaction.reg_ = [self._add_gene(g) for g in interaction.reg_]

        interaction.tar_.regs += [interaction.reg_] #regs are added as list, regs of doubly interactions are added as lists of size > 1

        rList = interaction.reg_
        name = getInterName(rList,interaction.tar_)
        interaction.tar_.inInteractions[name] = interaction

        for g in interaction.reg_:
            g.tars += [interaction.tar_]


    def init(self, mr_profile, update_half_resp = True):
        for g in self.attr_['mrs']:
            self.attr_['genes'][g].isMR_ = True
        self._set_levels()
        self._set_MR_profile(mr_profile)
        self._estimate_params(half_resp = update_half_resp)

    def _set_levels(self):
        U = set()
        Z = set()
        V = set([g for g in self.attr_['genes'].keys()])

        currLayer = 0
        while U != V:
            currVerts = set(filter(lambda v: set([g.name_ for g in self.attr_['genes'][v].tars]).issubset(Z), V-U))

            for v in currVerts:
                self.attr_['genes'][v].level = currLayer
                self.attr_["level2gene"][currLayer] += [self.attr_['genes'][v]]
                U.add(v)

            currLayer += 1
            Z = Z.union(U)

        self.attr_['maxLevel'] = currLayer - 1

    def _set_MR_profile(self, mrProfile):
        """
        mrProfile: is an instance of MR object (# TODO: )
        """
        for mr_name in mrProfile.profile.keys():#Todo might be changed when MR object is implemented
            assert(mr_name in self.attr_['mrs'])
            self.attr_['genes'][mr_name].prod_rates_ = np.array(mrProfile.profile[mr_name]).flatten()

        self.nCellTypes_ = mrProfile.nTypes_

    def _estimate_steady_state(self, level):
        """
        get genes of this level
        for each gene:
            get prod rate
            divide by lambda and return
        prod rate -> gene method f(conc_all_its_regs = None) --> if None, uses regs.stead_state_estimate
        add steady_state_estimate to gene attribute
        Note that all conc should be accounting different cell types
        """
        genes = self.attr_["level2gene"][level]
        for g in genes:
            prod = g._calc_prod(cTypes = list(range(self.nCellTypes_)), regs_conc = 'ss') #np.arr of size No. cell types
            g.ss_conc_ = prod / g.decay_


    def _estimate_half_response(self, level):
        """
        get genes of this level
        for each gene:
            get ingoing interactions
            for each interaction:
                set h to the mean_types[reg.steady_state_estimate]

        need to do sth about doubly interactions
        """
        genes = self.attr_["level2gene"][level]
        for g in genes:
            for i in g.inInteractions.values():
                if isinstance(i, SingleInteraction):
                    i.h_ = np.mean(i.reg_[0].ss_conc_)
                else:
                    raise ValueError("not implemented")


    def _estimate_params(self, half_resp):
        self._estimate_steady_state(level = self.attr_['maxLevel'])
        for l in range(self.attr_['maxLevel'] - 1, -1, -1):
            if half_resp:
                self._estimate_half_response(level = l)

            try:
                self._estimate_steady_state(level = l)
            except:
                raise ValueError("Half response parameters are not defined. Consider setting update_half_resp = True in init.")

    def to_df(self, path):
        ret = []
        for i in self.attr_['interactions'].values():
            if len(i.reg_) > 1:
                raise ValueError("coop regulations is not implemented")
            curr = [j.name_ for j in i.reg_]
            curr += [i.tar_.name_, i.k_, i.n_, i.h_]
            ret.append(curr)

        ret = pd.DataFrame(ret)
        ret.columns = ['reg','tar','k','n','h']
        ret.to_csv(path, header = True, index = True)
