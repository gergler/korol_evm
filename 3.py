import sys, ROOT, numpy
from ROOT import gRandom, TCanvas, TFormula, TF1,TF2, TH1F, TH2F, TStopwatch

ax = 0
bx = 0.5

ay = 1
by = 3

const = 0.5

def F(x, y): return (numpy.exp(x[0] + (x[1])**2))/(1 + (x[1])*numpy.sqrt(x[0])/2)

h2 = TH2F('h2', 'histogram tytle2', 100, 0, 0.5, 100, 1, 3);
for i in range(100):
    for j in range(100):
        for _ in range(int(F([0.001+0.5*i/100.0, 1.0+j/100.0],[]))):
            h2.Fill(0.001+0.5*i/100.0, 1+j/100.0);

c = ROOT.TCanvas("myCanvasName","The Canvas Title", 1200, 500)
h2.Draw('lego')
h2.GetXaxis().SetTitle('x')
h2.GetYaxis().SetTitle('y')
c.Draw();

# c2 = ROOT.TCanvas("myCanvasName2", "The Canvas Title", 1200, 500)  
# f1 = TF2("f1", F, 0.0, 0.5, 1.0, 3.0, 2)
# f1.SetParameter(0, 1)
# f1.Draw('surface')
# c2.Draw()

quantity = 10000
serial = 1000

arr = numpy.zeros(quantity)

integral = 0 
err = 0
cur = 0
cur_err = 0

for i in range (1, quantity + 1):
    if (i%serial == 0):
        err = numpy.sqrt((numpy.mean(arr[:i]**2)-numpy.mean(arr[:i])**2)/i) # sqrt from dispersion
        cur_err = numpy.sqrt((numpy.mean(arr[i-serial:i]**2)-numpy.mean(arr[i-serial:i])**2)/serial) # sqrt from dispersion for all
        print('Serial: %d, integral: %f, error: %f, cur_int/serial: %f, cur_error: %f' % (i, integral/i, err, cur/serial, cur_err))
        cur = 0
    x = gRandom.Uniform(ax, bx)
    y = gRandom.Uniform(ay, by)
    model = (bx-ax)*(by-ay)*F([x, y], 0)
    integral += model
    cur += model
    arr[i-1] = F([x, y], 0)*((bx-ax)*(by-ay)) # average method
    
def g1(x, y): return (numpy.exp(x[0] + (x[1])**2)/x[1] - (2 + x[1] * numpy.sqrt(x[0]))/x[1])
def g2(x, y): return 2*x[1]

def F2(x, y): return g1(x, y) + g2(x, y)

# def F2(x, y): return (numpy.exp(x[0] + (x[1])**2))/(1 + (x[1])*numpy.sqrt(x[0])/2) - 2 * x[1]

h3 = TH2F('h1', 'histogram tytle1', 100, 0, 0.5, 100, 1, 3);
for i in range(100):
    for j in range(100):
        for _ in range(int(F2([0.001+0.5*i/100.0, 1.0+j/100.0],[]))):
            h3.Fill(0.001+0.5*i/100.0, 1+j/100.0);

c4 = ROOT.TCanvas("myCanvasName","The Canvas Title", 1200, 500)
h3.Draw('lego')
h3.GetXaxis().SetTitle('x')
h3.GetYaxis().SetTitle('y')
c4.Draw();

# c3 = ROOT.TCanvas("myCanvasName3", "The Canvas Title", 1200, 500)  
# f2 = TF2("f2", F2, 0.0, 0.5, 1.0, 3.0, 1)
# f2.SetParameter(0, 1)
# f2.Draw('surface')
# c3.Draw()

analitic = (bx - ax)*(by - ay)*(by**2 - ay**2)
arr = numpy.zeros(quantity)
integral = 0 
err = 0
cur = 0
cur_err = 0

for i in range (1, quantity + 1):
    if (i%serial == 0):
        err = numpy.sqrt((numpy.mean(arr[:i]**2)-numpy.mean(arr[:i])**2)/i)
        cur_err = numpy.sqrt((numpy.mean(arr[i-serial:i]**2)-numpy.mean(arr[i-serial:i])**2)/serial)
        print('Serial: %d, integral: %f, error: %f, cur_int/serial: %f, cur_error: %f' % (i, integral/i + analitic, err, cur/serial + analitic, cur_err))
        cur = 0
    x = gRandom.Uniform(ax, bx)
    y = gRandom.Uniform(ay, by)
    model = (bx-ax)*(by-ay)*F2([x, y], 0)
    integral += model
    cur += model
    arr[i-1] = F2([x, y], 0)*((bx-ax)*(by-ay))
