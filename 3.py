import sys, ROOT, numpy
from ROOT import gRandom, TCanvas, TFormula, TF1,TF2, TH1F, TH2F, TStopwatch

ax = 0
bx = 0.5

ay = 1
by = 3

def F(x, y): return (numpy.exp(x[0] + (x[1])**2))/(1 + x[1]*numpy.sqrt(x[0])/2)

h1 = TH2F('h1', 'histogram', 100, 0, 0.5, 100, 1, 3);
for i in range(100):
    for j in range(100):
        for k in range(int(F([0.001+0.5*i/100.0, 1.0+j/100.0], []))):
            h1.Fill(0.001+0.5*i/100.0, 1.0+j/100.0);

# c1 = ROOT.TCanvas("myCanvasName1", "The Canvas Title", 1200, 500)  
# h1.Draw('lego')
# c1.Draw()

# c2 = ROOT.TCanvas("myCanvasName2", "The Canvas Title", 1200, 500)  
# f1 = TF2("f1", F, 0.0, 0.5, 1.0, 3.0, 2)
# f1.SetNpx(100)
# f1.Draw('colz')
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
        cur_err = numpy.sqrt((numpy.mean(arr[i-serial:i]**2)-numpy.mean(arr[i-serial:i])**2)/serial) # sqrt from dispersion normalized on serial
        print('Serial: %d, integral: %f, error: %f, cur_int/serial: %f, cur_error: %f' % (i, integral/i, err, cur/serial, cur_err))
        cur = 0
    x = gRandom.Uniform(ax, bx)
    y = gRandom.Uniform(ay, by)
    model = (bx-ax)*(by-ay)*F([x, y], 0)
    integral += model
    cur += model
    arr[i-1] = F([x, y], 0)*((bx-ax)*(by-ay))
 
Serial: 1000, integral: 611.258444, error: 37.139775, cur_int/serial: 611.258444, cur_error: 37.139775
Serial: 2000, integral: 575.719256, error: 24.880392, cur_int/serial: 540.180067, cur_error: 33.082532
Serial: 3000, integral: 577.240823, error: 20.353759, cur_int/serial: 580.283957, cur_error: 35.389154
Serial: 4000, integral: 582.451072, error: 17.922781, cur_int/serial: 598.081819, cur_error: 37.560612
Serial: 5000, integral: 578.137844, error: 15.865597, cur_int/serial: 560.884931, cur_error: 33.954857
Serial: 6000, integral: 592.946573, error: 14.802633, cur_int/serial: 666.990218, cur_error: 39.859935
Serial: 7000, integral: 594.763342, error: 13.812079, cur_int/serial: 605.663957, cur_error: 38.203830
Serial: 8000, integral: 600.576795, error: 13.030552, cur_int/serial: 641.270964, cur_error: 38.950444
Serial: 9000, integral: 593.339041, error: 12.169229, cur_int/serial: 535.437008, cur_error: 33.538270
Serial: 10000, integral: 596.776402, error: 11.566034, cur_int/serial: 627.712649, cur_error: 37.161468
