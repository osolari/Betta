from cmc import *
        
def compRelEffCurve(shape = np.linspace( 1, 2, 10, endpoint=False), 
                    n = np.linspace(2, 2000, 50).astype(np.int), kappa = .01, k = 20, c = 100, verbose = False):
    
    print(shape)
    
    def compMu(shape, c):
    
        tmp = CMC(shape=shape, n=np.int(1e7), c = c)
        tmp.take_sample()
        tmp.cmc()
        
        return(tmp.Zbar)
    
    print("Computing Mu,")
    MU = [compMu(shape, c) for i in range(1)]
    print(MU)
    mu = np.mean(MU)
    
    relEffs = np.zeros((k, len(n)))
    for i in range(len(n)):
        
        if verbose:
            print("Computing relative efficiency for n = {}".format(n[i]))
        
        tmp = CMC(shape = shape, n = n[i], c = c)
        
        for j in range(k):
            tmp.compRelEff(mu, kappa)
            relEffs[j, i] = tmp.relEff
    
    return(n, relEffs, shape, mu, c)


def compRelEffVarShape(kappa = 0.005):
    
    ALPHA = [np.linspace( 1, e, 10, endpoint=False) for e in range(2, 6)]
    
    return([compRelEffCurve(shape, kappa = kappa) for shape in ALPHA])


def compRelEffVarC(kappa = 0.005):
    
    C = np.linspace( 100, 1000, 10, endpoint=True)
    #C = np.linspace( 1, 5, 10)
    print "C:" + str(C)
    
    return([compRelEffCurve(shape = np.linspace( 1, 3, 10, endpoint=False), c = i, kappa = kappa) for i in C])


def compRelEffVarShapeVarMin(kappa = 0.005, c = 100, ALPHA = 
                             [(np.float(e)/2 + np.linspace( 0, 1, 10, endpoint=False)) for e in range(2, 6)]):
    
    return([compRelEffCurve(shape, kappa = kappa) for shape in ALPHA])


def processRatio(ratio):
    
    def rQuantile(x, q = [.25, .5, .75]):
        return(np.quantile(x[np.abs(x)< 1e5], q))
    
    def rMedian(x):
        try:
            return(np.median(x[np.abs(x) < 1e5]))
        except:
            return(np.nan)
    
    return(np.array([rMedian(ratio[:,i]) for i in range(ratio.shape[1])]))


def ratioCorrect(n, ratio):
    
    df = pd.DataFrame({'n':np.tile(n, ratio.shape[0]), 'ratio':ratio.flatten()})
    dfr = df[np.abs(df["ratio"]) < 3.5]
    md = smf.mixedlm("ratio ~ n",dfr, groups=dfr["n"])
    mdf = md.fit()
    
    return [dfr["n"], mdf.params[0] + mdf.params[1] * dfr["n"]]