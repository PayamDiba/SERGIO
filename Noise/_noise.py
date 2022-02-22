import pandas as pd
import numpy as np
from scipy.stats import wasserstein_distance
import cma
import matplotlib.pylab as plt
from sklearn.preprocessing import minmax_scale

class splatModel(object):
    def __init__(self):
        pass

    def outlier_effect(self, scData, outlier_prob, mean, scale):
        out_indicator = np.random.binomial(n = 1, p = outlier_prob, size = scData.shape[0])
        outlierGenesIndx = np.where(out_indicator == 1)[0]
        numOutliers = len(outlierGenesIndx)
        outFactors = np.random.lognormal(mean = mean, sigma = scale, size = numOutliers)
        for i, gIndx in enumerate(outlierGenesIndx):
            scData[gIndx,:] = scData[gIndx,:] * outFactors[i]

        return scData

    def lib_size_effect(self, scData, mean, scale):
        libFactors = np.random.lognormal(mean = mean, sigma = scale, size = scData.shape[1])
        normalizFactors = np.sum(scData, axis = 0 )
        binFactors = np.true_divide(libFactors, normalizFactors)
        binFactors = binFactors.reshape(1, scData.shape[1])
        binFactors = np.repeat(binFactors, scData.shape[0], axis = 0)
        return np.multiply(scData, binFactors)

    def dropout_indicator(self, scData, shape = 1, percentile = 65):
        scData_log = np.log(np.add(scData,1))
        log_mid_point = np.percentile(scData_log, percentile)
        prob_ber = np.true_divide (1, 1 + np.exp( -1*shape * (scData_log - log_mid_point) ))
        binary_ind = np.random.binomial(n = 1, p = prob_ber)
        return binary_ind

    def convert_to_UMIcounts (self, scData):
        return np.random.poisson (scData)

class scNoise(object):
    def __init__(self, model = 'splat'):
        if model == 'splat':
            self.sim_ = splatModel()
        else:
            raise ValueError("{} noise model is not implemented!".format(model))

        self.params = None

    def _add_outlier(self, expr, prob, mean, scale):
        return self.sim_.outlier_effect(expr, outlier_prob = prob, mean = mean, scale = scale)

    def _add_library_size(self, expr, mean, scale):
        return self.sim_.lib_size_effect(expr, mean = mean, scale = scale)

    def _add_dropout(self, expr, shape, percentile):
        ind = self.sim_.dropout_indicator(expr, shape = shape, percentile = percentile)
        return np.multiply(ind, expr)

    def _to_count(self, expr):
        return self.sim_.convert_to_UMIcounts(expr)

    def _add_noise(self,
                   expr,
                   prob_o = None,
                   mean_o = None,
                   scale_o = None,
                   mean_l = None,
                   scale_l = None,
                   shape_d = None,
                   perc_d = None,
                   to_count = True):

        ret = expr.copy()
        if prob_o and mean_o and scale_o:
            ret = self._add_outlier(ret, prob_o, mean_o, scale_o)

        if mean_l and scale_l:
            ret = self._add_library_size(ret, mean_l, scale_l)

        if shape_d and perc_d:
            ret = self._add_dropout(ret, shape_d, perc_d)

        if to_count:
            ret = self._to_count(ret)

        return ret

    def _stat_library_size(self, expr):
        return expr.sum(axis = 0)

    def _stat_zero_per_gene(self, expr):
        noneZeros = np.count_nonzero(expr, axis = 1)
        fracZero = (expr.shape[1] - noneZeros)/expr.shape[1]
        return fracZero

    def _stat_zero_per_cell(self, expr):
        noneZeros = np.count_nonzero(expr, axis = 0)
        fracZero = (expr.shape[0] - noneZeros)/expr.shape[0]
        return fracZero

    def _stat_mean_gene(self, expr):
        return np.mean(expr, axis = 1)

    def _stat_var_gene(self, expr):
        return np.var(expr, axis = 1)


    def _get_stat(self, expr, type):

        if type == 'library':
            return self._stat_library_size(expr)

        if type == 'zeroG':
            return self._stat_zero_per_gene(expr)

        if type == 'zeroC':
            return self._stat_zero_per_cell(expr)

        if type == 'meanG':
            return self._stat_mean_gene(expr)

        if type == 'varG':
            return self._stat_var_gene(expr)


    def _w1_dist(self, synth, real):
        return wasserstein_distance(synth, real)

    def _total_var(self, synth, real):
        """
        Calculates total variation
        """
        minVal = min(np.union1d(synth,real))
        maxVal = max(np.union1d(synth,real))

        pU,_ = np.histogram(synth, bins=50, range=(minVal,maxVal))
        pV,_ = np.histogram(real, bins=50, range=(minVal,maxVal))

        pU = np.true_divide(pU, np.sum(pU))
        pV = np.true_divide(pV, np.sum(pV))

        return 0.5 * np.sum(np.abs(pV-pU))

    def _quantile_dist(self, synth, real):
        synth_25 = np.percentile(synth,25)
        synth_75 = np.percentile(synth,75)
        real_25 = np.percentile(real,25)
        real_75 = np.percentile(real,75)
        synth_box = [i for i in synth if synth_25 <= i <= synth_75]
        real_box = [i for i in real if real_25 <= i <= real_75]

        synth_box = np.sort(synth_box)
        real_box = np.sort(real_box)

        concat = np.concatenate([synth_box, real_box]).flatten()
        concat = minmax_scale(concat,feature_range=(0,10))

        synth_box = concat[0:len(synth_box)].flatten()
        real_box = concat[len(synth_box):].flatten()

        synth_box = np.array_split(synth_box,10)
        real_box = np.array_split(real_box,10)
        synth_box = np.array([np.mean(i) for i in synth_box])
        real_box = np.array([np.mean(i) for i in real_box])

        return np.sum(np.abs(synth_box - real_box))

        #return np.abs(np.mean(synth_box) - np.mean(real_box))

    def _loss(self, synth_expr, real_expr, stats = ['library','zeroG','zeroC','meanG','varG']):
        loss = 0
        weights = [1,1,1,1,1]
        for s,w in zip(stats,weights):
            synth = self._get_stat(synth_expr, s)
            real = self._get_stat(real_expr, s)
            loss += w*self._quantile_dist(synth, real)
            #loss += w*self._w1_dist(synth, real) + w*self._total_var(synth, real) + w*self._median_dist(synth, real)
            #loss += w*self._total_var(synth, real)

        return loss

    def sample(self, expr, size):
        ind = np.random.choice(range(expr.shape[1]), size, replace = False)
        return expr[:,ind]

    def _eval_match(self, params, sample):
        curr_synth = self.sample(self.synth_expr_, sample) if sample else self.synth_expr_
        curr_real = self.sample(self.real_expr_, sample) if sample else self.real_expr_
        currNoisy = self._add_noise(expr = curr_synth,
                                    prob_o = params[0],
                                    mean_o = params[1],
                                    scale_o = params[2],
                                    mean_l = params[3],
                                    scale_l = params[4],
                                    shape_d = params[5],
                                    perc_d = params[6],
                                    to_count = True)

        loss = self._loss(currNoisy, curr_real)
        return loss

    def _optimize(self, tolX, tolIter, sample):
        #x0 = np.array([0.01, 0.8, 1, 5, 1, 5, 63])
        #low_bounds = [0.009, 0.75, 0.9, 1, 0.4, 2, 30]
        #up_bounds = [0.012, 5, 2.1, 10, 1.5, 8, 70]

        x0 = np.array([0.01, 0.8, 1, 5, 1, 4, 65])
        low_bounds = [0.009, 0.75, 0.01, 1, 0.3, 0.5, 30]
        up_bounds = [0.05, 10, 5, 10, 5, 10, 70]
        es = cma.CMAEvolutionStrategy(x0, 0.1, {'bounds':[low_bounds, up_bounds],'tolx':tolX, 'CMA_stds':[0.1,1,1,1,1,1,1]})
        while not es.stop() and es.result.iterations <= tolIter:
            solutions = es.ask()
            es.tell(solutions, [self._eval_match(s, sample) for s in solutions])
            es.disp()
        #es.optimize(self._eval_match)
        self.params = es.result.xbest

    def match(self, real, synth_clean, tol_param = 0.01, tol_iter = 200, sample = None):
        self.real_expr_ = real
        self.synth_expr_ = synth_clean
        self._optimize(tol_param, tol_iter, sample)
        self.noisy_expr = self._add_noise(expr = self.synth_expr_,
                                          prob_o = self.params[0],
                                          mean_o = self.params[1],
                                          scale_o = self.params[2],
                                          mean_l = self.params[3],
                                          scale_l = self.params[4],
                                          shape_d = self.params[5],
                                          perc_d = self.params[6],
                                          to_count = True)

    def get_noisy_expr(self):
        if self.params:
            return self.noisy_expr

    def plot(self):
        stats = {}
        for expr, key in zip([self.synth_expr_, self.real_expr_, self.noisy_expr], ['clean','real','noisy']):
            for s in ['library','zeroG','zeroC','meanG','varG']:
                stats[key,s] = self._get_stat(expr, s)

        fig,ax = plt.subplots(1,5,figsize=(16,5))
        ax[0].boxplot([stats[k,'library'] for k in ['real','noisy']])
        ax[1].boxplot([stats[k,'zeroG'] for k in ['real','noisy']])
        ax[2].boxplot([stats[k,'zeroC'] for k in ['real','noisy']])
        ax[3].boxplot([stats[k,'meanG'] for k in ['real','noisy']])
        ax[4].boxplot([stats[k,'varG'] for k in ['real','noisy']])

        ax[0].set_xticklabels(['real','noisy'], rotation = 45)
        ax[1].set_xticklabels(['real','noisy'], rotation = 45)
        ax[2].set_xticklabels(['real','noisy'], rotation = 45)
        ax[3].set_xticklabels(['real','noisy'], rotation = 45)
        ax[4].set_xticklabels(['real','noisy'], rotation = 45)

        ax[0].set_title('library size')
        ax[1].set_title('zeros per gene')
        ax[2].set_title('zeros per cell')
        ax[3].set_title('mean gene expression')
        ax[4].set_title('var gene expression')
        plt.tight_layout()
