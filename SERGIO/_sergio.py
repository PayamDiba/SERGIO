import numpy as np
import pandas as pd
from numba import jit
import logging
class sergio(object):
    def __init__(self, grn, diff_graph = None):
        """
        TODO: add functionality for dynamics simulations
        ,
        """

        self.grn_ = grn
        self.diff_graph_ = diff_graph
        self._init_conc()
        self.gNames_ = self._get_gene_names(grn)
        self.lambda_ = self._get_lambda(self.gNames_, grn)

    def _init_conc(self):
        cTypes = list(range(self.grn_.nCellTypes_))
        for g in self.grn_.attr_['genes'].values():
            g.append_sim_conc(conc = g.ss_conc_, cTypes = cTypes)

    def _get_gene_names(self, grn):
        return [k for k in grn.attr_['genes'].keys()]

    def _get_lambda(self, gnames, grn):
        ret = [grn.attr_['genes'][gn].decay_ for gn in gnames]
        return np.array(ret).reshape(-1,1)

    
    def _iter_ss(self, noise_ss, dt, cTypes):
        X = np.empty(shape = (len(self.gNames_),len(cTypes)))
        P = np.empty(shape = (len(self.gNames_),len(cTypes)))
        L = self.lambda_
        for ri,gn in enumerate(self.gNames_):
            gene = self.grn_.attr_['genes'][gn]
            X[ri] = gene.get_last_conc(cTypes)
            P[ri] = gene._calc_prod(cTypes, regs_conc = 'sim')

        P[P < 0] = 0 # numerical stability
        D = L*X
        rndP = np.random.normal(size = (len(self.gNames_),len(cTypes)))
        rndD = np.random.normal(size = (len(self.gNames_),len(cTypes)))

        newX = X + (P - D)*dt + (np.multiply(np.sqrt(P),rndP) + np.multiply(np.sqrt(D),rndD))*noise_ss*np.sqrt(dt)
        for gn,conc in zip(self.gNames_,newX):
            self.grn_.attr_['genes'][gn].append_sim_conc(conc.flatten(), cTypes)


    def _simulate_ss(self, nCells, noise_ss, dt = 0.01, safety_iter = 50, scale_iter = 10):
        """
        # TODO: make sure nCells is already np.array
        """
        # first do safety iterations
        cTypes = list(range(self.grn_.nCellTypes_))
        for _ in range(safety_iter):
            self._iter_ss(noise_ss, dt, cTypes)

        # next simulate required iterations
        req = nCells * scale_iter
        cTypes = list(range(self.grn_.nCellTypes_))
        nIter = 0
        while(cTypes):
            self._iter_ss(noise_ss, dt, cTypes)
            nIter += 1
            cTypes = np.where(req > nIter)[0].tolist()

    def simulate(self, nCells, noise_s, noise_u = None, safety_iter = 50, scale_iter = 10, dt = 0.01):
        if isinstance(nCells, int) or isinstance(nCells, float):
            self.nCells_ = np.array([nCells] * self.grn_.nCellTypes_)
        else:
            assert(len(nCells) == self.grn_.nCellTypes_)
            self.nCells_ = np.array(nCells)
        self.safety_iter_ = safety_iter
        self.scale_iter_ = scale_iter


        if not self.diff_graph_:
            self._simulate_ss(self.nCells_, noise_ss = noise_s, dt = dt, safety_iter = safety_iter, scale_iter = scale_iter)
        else:
            raise ValueError('not implement')

    def _get_rnd_scInd(self):
        ret = {}
        for ct in range(self.grn_.nCellTypes_):
            ret[ct] = self.safety_iter_ + np.random.randint(low = 0, high = self.nCells_[ct] * self.scale_iter_, size = self.nCells_[ct])
        return ret

    def getSimExpr(self):
        rndInd = self._get_rnd_scInd()
        expr = []
        cell_names = []
        gene_names = []
        for g in self.grn_.attr_['genes'].values():
            curr = []
            gene_names.append(g.name_)
            for ct in rndInd.keys():
                ind = rndInd[ct].tolist()
                curr += g.sim_conc_[ct][ind].tolist()
            expr.append(curr)


        cell_names += ['type_{}_cell_{}'.format(ct, i) for ct,n in enumerate(self.nCells_) for i in range(n)]
        expr = pd.DataFrame(expr)
        expr.columns = cell_names
        expr.index = gene_names

        return expr
    def _iter_ss_gpu(self, noise_ss, dt, cTypes):
        import cupy as cp
        X = cp.empty(shape = (len(self.gNames_),len(cTypes)))
        P = cp.empty(shape = (len(self.gNames_),len(cTypes)))
        L = cp.asarray(self.lambda_)
        for ri,gn in enumerate(self.gNames_):
            gene = self.grn_.attr_['genes'][gn]
            X[ri] = cp.asarray(gene.get_last_conc(cTypes))
            P[ri] = cp.asarray(gene._calc_prod(cTypes, regs_conc = 'sim'))

        P[P < 0] = 0 # numerical stability
        D = L*X
        rndP = cp.random.normal(size = (len(self.gNames_),len(cTypes)))
        rndD = cp.random.normal(size = (len(self.gNames_),len(cTypes)))

        newX = X + (P - D)*dt + (cp.multiply(cp.sqrt(P),rndP) + cp.multiply(cp.sqrt(D),rndD))*noise_ss*cp.sqrt(dt)
        for gn,conc in zip(self.gNames_,cp.asnumpy(newX)):
            self.grn_.attr_['genes'][gn].append_sim_conc(conc.flatten(), cTypes)


    def _simulate_ss_gpu(self, nCells, noise_ss, dt = 0.01, safety_iter = 50, scale_iter = 10):
        """
        # TODO: make sure nCells is already np.array
        """

        # first do safety iterations
        cTypes = list(range(self.grn_.nCellTypes_))
        for _ in range(safety_iter):
            self._iter_ss_gpu(noise_ss, dt, cTypes)

        # next simulate required iterations
        req = nCells * scale_iter
        cTypes = list(range(self.grn_.nCellTypes_))
        nIter = 0
        while(cTypes):
            self._iter_ss_gpu(noise_ss, dt, cTypes)
            nIter += 1
            cTypes = np.where(req > nIter)[0].tolist()

    def simulate_gpu(self, nCells, noise_s, noise_u = None, safety_iter = 50, scale_iter = 10, dt = 0.01):
        logging.warning("It is much slower than the CPU version, don't use it")
        if isinstance(nCells, int) or isinstance(nCells, float):
            self.nCells_ = np.array([nCells] * self.grn_.nCellTypes_)
        else:
            assert(len(nCells) == self.grn_.nCellTypes_)
            self.nCells_ = np.array(nCells)
        self.safety_iter_ = safety_iter
        self.scale_iter_ = scale_iter


        if not self.diff_graph_:
            self._simulate_ss_gpu(self.nCells_, noise_ss = noise_s, dt = dt, safety_iter = safety_iter, scale_iter = scale_iter)
        else:
            raise ValueError('not implement')

    """""""""""""""""""""""""""""""""""""""
    "" This part is to add technical noise
    """""""""""""""""""""""""""""""""""""""
    def outlier_effect(self, scData, outlier_prob, mean, scale):
        """
        scData: 2d np.array or pd.Dataframe contaiing genes as rows, cells as columns 
        """
        if isinstance(scData,pd.DataFrame):
            scData = scData.values
        out_indicator = np.random.binomial(n = 1, p = outlier_prob, size = len(self.gNames_))
        outlierGenesIndx = np.where(out_indicator == 1)[0]
        numOutliers = len(outlierGenesIndx)

        #### generate outlier factors ####
        outFactors = np.random.lognormal(mean = mean, sigma = scale, size = numOutliers)
        ##################################
        #if isinstance(scData, pd.DataFrame):
        #    scData = scData.values
        for i, gIndx in enumerate(outlierGenesIndx):
            scData[gIndx,:] = scData[gIndx,:] * outFactors[i]

        return scData


    def lib_size_effect(self, scData, mean, scale):
        """
        This functions adjusts the mRNA levels in each cell seperately to mimic
        the library size effect. To adjust mRNA levels, cell-specific factors are sampled
        from a log-normal distribution with given mean and scale.
        scData: the simulated data representing mRNA levels (concentrations);
        np.array (#bins * #genes * #cells)
        mean: mean for log-normal distribution
        var: var for log-normal distribution
        returns libFactors ( np.array(nBin, nCell) )
        returns modified single cell data ( np.array(nBin, nGene, nCell) )
        """

        ret_data = []

        libFactors = np.random.lognormal(mean = mean, sigma = scale, size = sum(self.nCells_))
        if len(libFactors)!=scData.shape[1]:
            raise ValueError('Number of cells in scData expected to be '+str(len(libFactors))+', having  '+str(scData.shape[1])+' instead')
        for cell_idx, cellFactors in enumerate(libFactors):
            cellExprMatrix = scData[:,cell_idx]#select a column of the dataframe
            normalizFactors = np.sum(cellExprMatrix, axis = 0 )
            cellFactors = np.true_divide(cellFactors, normalizFactors)
            cellFactors = np.repeat(cellFactors, len(self.gNames_), axis = 0)#create a matrix with repeated rows, where value is given by the technical noise term for each cell

            ret_data.append(np.multiply(cellExprMatrix, cellFactors))


        return libFactors, np.array(ret_data).T


    def dropout_indicator(self, scData, shape = 1, percentile = 65):
        """
        This is similar to Splat package
        Input:
        scData can be the output of simulator or any refined version of it
        (e.g. with technical noise)
        shape: the shape of the logistic function
        percentile: the mid-point of logistic functions is set to the given percentile
        of the input scData
        returns: np.array containing binary indactors showing dropouts
        """
        if isinstance(scData,pd.DataFrame):
            scData = scData.values
        scData = np.array(scData)
        scData_log = np.log(np.add(scData,1))
        log_mid_point = np.percentile(scData_log, percentile)
        prob_ber = np.true_divide (1, 1 + np.exp( -1*shape * (scData_log - log_mid_point) ))

        binary_ind = np.random.binomial( n = 1, p = prob_ber)

        return binary_ind

    def convert_to_UMIcounts (self, scData):
        """
        Input: scData can be the output of simulator or any refined version of it
        (e.g. with technical noise)
        """

        return np.random.poisson (scData)

    """""""""""""""""""""""""""""""""""""""""""""""""""""""""
    "" This part is to add technical noise to dynamics data
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""
    def outlier_effect_dynamics(self, U_scData, S_scData, outlier_prob, mean, scale):
        """
        This function
        """
        out_indicator = np.random.binomial(n = 1, p = outlier_prob, size = len(self.gNames_))
        outlierGenesIndx = np.where(out_indicator == 1)[0]
        numOutliers = len(outlierGenesIndx)

        #### generate outlier factors ####
        outFactors = np.random.lognormal(mean = mean, sigma = scale, size = numOutliers)
        ##################################

        U = np.concatenate(U_scData, axis = 1)
        S = np.concatenate(S_scData, axis = 1)
        for i, gIndx in enumerate(outlierGenesIndx):
            U[gIndx,:] = U[gIndx,:] * outFactors[i]
            S[gIndx,:] = S[gIndx,:] * outFactors[i]

        return np.split(U, self.nBins_, axis = 1), np.split(S, self.nBins_, axis = 1)


    def lib_size_effect_dynamics(self, U_scData, S_scData, mean, scale):
        """
        """

        #TODO make sure that having bins does not intefere with this implementation
        ret_data_U = []
        ret_data_S = []

        libFactors = np.random.lognormal(mean = mean, sigma = scale, size = (self.nBins_, self.nSC_))
        for binExprU, binExprS, binFactors in zip(U_scData, S_scData, libFactors):
            normalizFactors_U = np.sum(binExprU, axis = 0 )
            normalizFactors_S = np.sum(binExprS, axis = 0 )
            binFactors = np.true_divide(binFactors, normalizFactors_U + normalizFactors_S)
            binFactors = binFactors.reshape(1, self.nSC_)
            binFactors = np.repeat(binFactors, len(self.gNames_), axis = 0)

            ret_data_U.append(np.multiply(binExprU, binFactors))
            ret_data_S.append(np.multiply(binExprS, binFactors))


        return libFactors, np.array(ret_data_U), np.array(ret_data_S)


    def dropout_indicator_dynamics(self, U_scData, S_scData, shape = 1, percentile = 65):
        """
        """
        scData = np.array(U_scData) + np.array(S_scData)
        scData_log = np.log(np.add(scData,1))
        log_mid_point = np.percentile(scData_log, percentile)
        U_log = np.log(np.add(U_scData,1))
        S_log = np.log(np.add(S_scData,1))
        prob_ber_U = np.true_divide (1, 1 + np.exp( -1*shape * (U_log - log_mid_point) ))
        prob_ber_S = np.true_divide (1, 1 + np.exp( -1*shape * (S_log - log_mid_point) ))

        binary_ind_U = np.random.binomial( n = 1, p = prob_ber_U)
        binary_ind_S = np.random.binomial( n = 1, p = prob_ber_S)

        return binary_ind_U, binary_ind_S

    def convert_to_UMIcounts_dynamics (self, U_scData, S_scData):
        """
        Input: scData can be the output of simulator or any refined version of it
        (e.g. with technical noise)
        """

        return np.random.poisson (U_scData), np.random.poisson (S_scData)