# James Osborn
# Plot field of view of telescope with satellite trails overlaid
# Observer and Field parameters defined in params.py

import datetime
from tqdm import tqdm
from astrosat import AstroSat
from astrosat import Plot

AS = AstroSat("params.yaml")

timeStep = 5  # 5 second time steps 

# get satellite positions
# use active satellites and visual objects (rocket bodies)
satTLEs = AS.get_TLEs('active')
satTLEs += AS.get_TLEs('visual')
sats = AS.get_satellites(satTLEs)

satDict = {}
for istep in tqdm(range(0, int(AS.parameters.duration), timeStep)):
    dateTemp = AS.parameters.date+datetime.timedelta(seconds=istep)
    for sat in sats:
        satDict = AS.process_satellite(sat, dateTemp, satDict)

# print satellite dictionary
sat_table = AS.print_satellite_dictionary(satDict)

# plot field
plot = Plot(AS)
plot.plot_satellites(satDict)
plot.plot_stars()
plot.plot_legend([AS.parameters.Mmin, int(round((AS.parameters.Mmin + AS.parameters.Mmax) / 3.)), 2 * int(round((AS.parameters.Mmin + AS.parameters.Mmax) / 3.))])
plot.make_plot()
plot.save_plot('skyView_%i_%i.png' % (AS.parameters.radius, AS.parameters.date.timestamp()))
input()
