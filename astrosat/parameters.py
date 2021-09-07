import numpy
import ephem
import datetime
import yaml


class Parameters:
    def __init__(self, parameter_file):

        self.parameter_file = parameter_file
        self.parameters_dict = self.load_parameters()

        self.satArea = float(self.parameters_dict['satArea'])
        self.albedo = float(self.parameters_dict['albedo'])
        # self.satType = self.parameters_dict['satType']
        self.radius = float(self.parameters_dict['radius'])
        self.duration = float(self.parameters_dict['duration'])
        self.Mmin = float(self.parameters_dict['Mmin'])
        self.Mmax = float(self.parameters_dict['Mmax'])
        self.plotDir = self.parameters_dict['plotDir']
        self.fdir = self.parameters_dict['fdir']
        self.TLEdir = self.parameters_dict['TLEdir']

        ra = self.parameters_dict['RA'].split(':')
        dec = self.parameters_dict['DEC'].split(':')
        self.RA = (numpy.array(ra).astype(float)*numpy.array([1, 1./60., 1./3600.])).sum()*180./12.
        self.DEC = (numpy.array(dec).astype(float)*numpy.array([1, 1./60., 1./3600.])).sum()

        self.date = datetime.datetime(year=self.parameters_dict['year'], month=self.parameters_dict['month'],
                                      day=self.parameters_dict['day'], hour=self.parameters_dict['hour'],
                                      minute=self.parameters_dict['minute'], second=self.parameters_dict['seconds'])

        self.obs = ephem.Observer()
        self.obs.lon = str(self.parameters_dict['lon'])
        self.obs.lat = str(self.parameters_dict['lat'])
        self.obs.elevation = self.parameters_dict['alt']

        self.verbose = self.parameters_dict['verbose']

    def load_parameters(self):
        with open(self.parameter_file, 'r') as stream:
            try:
                parameters_dict = yaml.safe_load(stream)
            except yaml.YAMLError as exc:
                print(exc)
        return parameters_dict
