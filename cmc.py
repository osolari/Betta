import numpy as np
from scipy.stats import pareto

class CMC:
    
    def __init__(self, shape = None, scale = None, c = 100, N = 10, n = 100, K = 50):
        
        if shape is None:
            self.shape = np.linspace( 1, 2, N, endpoint=False)
        else:
            self.shape = shape
            
        if scale is None:
            self.scale = np.repeat(1, len(shape))
        else:
            self.scale = scale
        
        assert isinstance(n, int), "n must be an integer."
        
        self.c = c
        self.N = len(self.shape)
        self.n = n
        self.K = K
        
    def take_sample(self):
        
        self.sample = np.asmatrix([np.random.pareto(self.shape[i], self.n)*self.scale[i] 
                                   for i in range(len(self.shape))])
        
    def cmc(self):
    
        z = 0
    
        def compZi(i, sample, c, shapei, scalei):
            
            ix = range(i) + range(i+1, sample.shape[0])
            return(pareto.sf(np.max(np.concatenate((c - np.sum(sample[ix,:], 0), 
                                                    np.max(sample[ix,:], 0))), 0), shapei, scale = scalei))
        
        self.Z = np.sum(np.concatenate([compZi(i, self.sample, self.c, self.shape[i], self.scale[i]) 
                         for i in range(self.sample.shape[0])]), 0)
        self.Zbar = np.mean(self.Z)
        
    def mc(self):
        
        self.muHat = (np.sum(self.sample, 0) > self.c).astype(int)
        self.muBar = np.mean(self.muHat)
        
    def compRelEff(self, mu, kappa = 0.01):
        
        def mccmc(obj):
            obj.take_sample()
            obj.cmc()
            obj.mc()
            
            return((obj.muBar, obj.Zbar))
        
        self.muHatZbars = np.asmatrix([mccmc(self) for i in range(self.K)])
        
        self.relEff = np.log(np.mean(np.abs(self.muHatZbars[:,1] - mu) > 
                              kappa * mu)) - np.log(np.mean(np.abs(self.muHatZbars[:,0] - mu) > 
                                                         kappa * mu))
        
