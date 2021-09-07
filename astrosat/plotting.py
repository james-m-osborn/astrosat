
import numpy
import matplotlib.pylab as plt
import astropy.coordinates
import astropy.units as u
import scipy.interpolate
import os


class Plot:
    def __init__(self, AS):
        self.AS = AS

        # setup plots
        try:
            os.mkdir(AS.parameters.fdir)
        except FileExistsError:
            pass
        try:
            os.mkdir(AS.parameters.TLEdir)
        except FileExistsError:
            pass
        plt.style.use('dark_background')
        if AS.parameters.radius == 180:
            AS.markersize = 1.
        else:
            AS.markersize = 5

        AS.fig = plt.figure('Sky field', figsize=(7, 6))
        AS.axSky = AS.fig.add_subplot(111)  # , projection='polar')

    def plot_stars(self):
        # add stars to plot

        Mmin = self.AS.parameters.Mmin
        Mmax = self.AS.parameters.Mmax
        for i_d in range(len(self.AS.stars)):
            RAstar = self.AS.stars[i_d][1]
            DECstar = self.AS.stars[i_d][2]
            RAanglediff = (self.AS.parameters.RA - RAstar + 180 + 360) % 360 - 180
            DECanglediff = (self.AS.parameters.DEC - DECstar + 180 + 360) % 360 - 180

            # over sample as field is square
            if abs(RAanglediff) < self.AS.parameters.radius*1.5:
                if abs(DECanglediff) < self.AS.parameters.radius*1.5:
                    mag_star = self.AS.stars[i_d][3]
                    r = 1+1*(Mmin-mag_star)/float(Mmax-Mmin)
                    r = numpy.clip(r, 0, 1)

                    plt.plot(RAstar*12./180., DECstar, color=str(r), marker='o', markersize=r*5)

    def plot_satellites(self, sat_dict, label=0):
        Mmin = self.AS.parameters.Mmin
        Mmax = self.AS.parameters.Mmax
        # Add satellites to plot
        for i_sat in (sat_dict.keys()):
            mag_sat = sat_dict[i_sat]['MAG'][0]
            if mag_sat is not None:
                r = 1+1*(Mmin-mag_sat)/float(Mmax-Mmin)
                r = numpy.clip(r, 0, 1)

                if len(sat_dict[i_sat]['RA']) == 1:
                    plt.plot(sat_dict[i_sat]['RA'], sat_dict[i_sat]['DEC'], color=(r, 0., 0.), marker='o',
                                markersize=r * 5)

                else:
                    # find transit time through FoV
                    f = scipy.interpolate.interp1d(numpy.array(sat_dict[i_sat]['RA'])*180/12.,
                                                    sat_dict[i_sat]['DEC'], fill_value='extrapolate')
                    RAtemp = numpy.arange(self.AS.parameters.RA-self.AS.parameters.radius,
                                            self.AS.parameters.RA+self.AS.parameters.radius, 1./3600.)
                    DECextrap = f(RAtemp)
                    plt.plot(RAtemp*12./180., DECextrap, linestyle='--', color=(r, 0, 0), linewidth=r*5)

    def plot_legend(self, labels):
        Mmin = self.AS.parameters.Mmin
        Mmax = self.AS.parameters.Mmax
        # add legend to plot
        label_list = []
        for label in labels:
            r = 1+1*(Mmin-label)/float(Mmax-Mmin)
            r = numpy.clip(r, 0, 1.)
            p1, = plt.plot(-1, -400, 'o', color=str(r), markersize=r*5, label='mag %i' % label)

            label_list.append(p1)
        for label in labels:
            r = 1+1*(Mmin-label)/float(Mmax-Mmin)
            r = numpy.clip(r, 0, 1.)
            p2, = plt.plot(-1, -400, color=(r, 0., 0.), linestyle='--', linewidth=r*5, label='mag %i' % label)
            label_list.append(p2)

        plt.legend(handles=label_list, numpoints=1, bbox_to_anchor=(1.05, 1), loc='upper left')

    def make_plot(self):
        # make axis scale and labels
        self.AS.axSky.axis([self.AS.parameters.RA*12./180. - self.AS.parameters.radius*24./360.,
                            self.AS.parameters.RA*12./180. + self.AS.parameters.radius*24./360.,
                            self.AS.parameters.DEC-self.AS.parameters.radius, self.AS.parameters.DEC+self.AS.parameters.radius])
        labels = self.AS.axSky.xaxis.get_majorticklocs()
        xlabel_list = []
        for label in labels:
            angle = astropy.coordinates.Angle(float(label), u.hour)
            xlabel_list.append( angle.to_string(unit=u.hour, precision=0))
        labels = self.AS.axSky.yaxis.get_majorticklocs()
        ylabel_list = []
        for label in labels:
            angle = astropy.coordinates.Angle(float(label),u.degree)
            ylabel_list.append(angle.to_string(unit=u.degree, sep=':', precision=0))
        self.AS.axSky.set_xticklabels(xlabel_list, rotation=30)
        self.AS.axSky.set_yticklabels(ylabel_list)
        plt.subplots_adjust(left=0.15, right=0.75, top=0.9, bottom=0.2)
        self.AS.axSky.set_xlabel('RA')
        self.AS.axSky.set_ylabel('DEC')

    def save_plot(self, file_name):
        if not os.path.isdir(self.AS.parameters.plotDir):
            os.mkdir(self.AS.parameters.plotDir)

        self.AS.fig.savefig(self.AS.parameters.plotDir + file_name, bbox_inches='tight')
