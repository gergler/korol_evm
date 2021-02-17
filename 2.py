import warnings
warnings.filterwarnings("ignore")
from IPython.display import clear_output
clear_output(wait=True)

import sys, ROOT, math
from ROOT import gRandom, TCanvas, TFormula, TF1, TH1F, TStopwatch
from math import atan, tan, log

a = -1
b = 1

def f(x, p): return ((atan(x[0]))**2/(1+(x[0])**2) + 1/(1.0201-(x[0])**2))*p[0]

c = ROOT.TCanvas()
ff = TF1("myfunction", f, a, b, 1)
ff.SetNpx(10000)
ff.SetParameter(0, 1);
ff.Draw();
c.Draw();

def fonneuman(a, b):
    fmax = ff(1)
    while True:
        mu = gRandom.Uniform(0, fmax)
        r = gRandom.Uniform(a, b)
        if mu<=ff(r): return r

sw = ROOT.TStopwatch()
c1 = ROOT.TCanvas()
h1 = TH1F("h1","fonneuman method",50,-1,1)
sw.Start()
for i in range(0, 10000):
    h1.Fill(fonneuman(a, b))
sw.Stop
sw.Print()
h1.Draw();
h1.Fit(ff);
c1.Draw();

print('first moment = ', h1.GetMean())
print('mean square = ', h1.GetStdDev())

def g1(x): return ((atan(x))**2)/(1+(x)**2)
def g2(x): return 1/(1.0201-(x)**2)

с = 1.0201
rt = 1.01
d = 0.49505

b1 = (atan(b)**3)/3 - (atan(a)**3)/3
b2 = d*log(b + rt) - d*log(rt - b) - d*log(a + rt) + d*log(rt - a)

a1 = (b1)/(b1+b2)
a2 = (b2)/(b1+b2)

gf1 = ROOT.TF1("function", g1, a, b, 1)
gf1.SetParameter(0, 1/b1)
            
gf2 = ROOT.TF1("function", g2, a, b, 1)
gf2.SetParameter(0, 1/b2)

def F1(x, b1): return ((atan(x)**3)/3 - (atan(a)**3)/3)/b1
def F2(x, b2): return (d*log(x + rt) - d*log(rt - x) - d*log(a + rt) + d*log(rt - a))/b2

def inv1(x, b1): return tan(3*x*b1 + atan(a))
def inv2(x, b2): return rt*(-a*exp(x*b2/d) - a - rt*exp(x*b2/d) + rt)/(-a*exp(x*b2/d) + a - rt*exp(x*b2/d) -rt)

print('Inv1: ', inv1(0, b1))
print('Inv2: ', inv2(0, b2))

def fncomposition(a, b, a1, a2, b1, b2):
    while (True):
        nr = gRandom.Uniform(0,1)
        nk = gRandom.Uniform(0,1)
        if nk < a1:
            return inv1(nr, b1)
        else:
            return inv2(nr, b2)
        
c2 = ROOT.TCanvas("myCanvasName2","The Canvas Title", 1200, 500)    
h2 = TH1F("h2", "composition method", 50, -1, 1)
sw = TStopwatch()
sw.Start()
for i in range(0, 10000):
    r = fncomposition(a, b, a1, a2, b1, b2)
    h2.Fill(r)
sw.Stop()
sw.Print()
h2.Draw()
h2.Fit(ff)
c2.Draw()

print('First moment: ', h2.GetMean())
print('Mean square: ', h2.GetStdDev())
