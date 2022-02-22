import numpy as np
import pandas as pd

class mrProfile(object):
    def __init__(self, MR_names, n_types):
        self.profile = dict.fromkeys(MR_names, None)
        for k in self.profile.keys():
            self.profile[k] = np.zeros(shape = n_types)
        self.mrNames_ = MR_names
        self.nTypes_ = n_types
        self.init_ = False

    def build_from_file(self, path, header = 0):
        if self.init_:
            print("profile is already initilized.")
            return

        prof = pd.read_csv(path, header = header, index_col = None)
        try:
            prof = prof.astype(float)
        except:
            pass

        prof.columns = ['index'] + list(prof.columns)[1:]
        prof.set_index('index', inplace=True)
        prof.index = prof.index.map(str)

        assert(set(list(prof.index)) == set(list(self.mrNames_)))
        assert(prof.shape[1] == self.nTypes_)
        for mr, p in prof.iterrows():
            self.profile[mr] = p.values

        self.init_ = True

    def build_rnd(self, range_dict = {'L':[1,2.5], 'H': [3.5,5]}):
        if self.init_:
            print("profile is already initilized.")
            return


        levels = list(range_dict.keys())
        for mr in self.mrNames_:
            for ct in range(self.nTypes_):
                currK = np.random.choice(levels,1)[0]
                currRange = range_dict[currK]
                self.profile[mr][ct] = np.random.uniform(low = currRange[0], high = currRange[1], size = 1)[0]

        self.init_ = True

    def build_rnd_complex(self, nMarker_per_type, nLow, nHigh, low = [1,2.5], high = [3.5,5]):
        if self.init_:
            print("profile is already initilized.")
            return


        total = nMarker_per_type * self.nTypes_ + nLow + nHigh
        assert(total <= len(self.mrNames_))
        if total == 0:
            raise ValueError("Specified parameters do not require complex randomization, use 'build_rnd' instead.")
        nMulti = len(self.mrNames_) - total
        conds = ["mu"] * nMulti + ["h"] * nHigh + ["l"] * nLow + ['mrk_{}'.format(i) for i in range(self.nTypes_)] * nMarker_per_type
        np.random.shuffle(conds)

        for mr,cond in zip(self.mrNames_,conds):
            if cond == 'mu':
                nh = np.random.choice(range(2,self.nTypes_-1), 1)[0]
                hInds = np.random.choice(range(self.nTypes_), nh, replace = False)
                for ct in range(self.nTypes_):
                    if ct in hInds:
                        self.profile[mr][ct] = np.random.uniform(low = high[0], high = high[1], size = 1)[0]
                    else:
                        self.profile[mr][ct] = np.random.uniform(low = low[0], high = low[1], size = 1)[0]

            elif cond == 'h':
                for ct in range(self.nTypes_):
                    self.profile[mr][ct] = np.random.uniform(low = high[0], high = high[1], size = 1)[0]

            elif cond == 'l':
                for ct in range(self.nTypes_):
                    self.profile[mr][ct] = np.random.uniform(low = low[0], high = low[1], size = 1)[0]

            else:
                mrkCT = int(float(cond.split('_')[1]))
                for ct in range(self.nTypes_):
                    if ct == mrkCT:
                        self.profile[mr][ct] = np.random.uniform(low = high[0], high = high[1], size = 1)[0]
                    else:
                        self.profile[mr][ct] = np.random.uniform(low = low[0], high = low[1], size = 1)[0]

        self.init_ = True





"""
things to implement
1- various pertubation methods: e.g. over-expr, down-reg, KO
"""
