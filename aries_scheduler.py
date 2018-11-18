from datetime import date
from datetime import timedelta
from datetime import datetime
from dateutil.relativedelta import relativedelta
import pandas as pd
import numpy as np
from bokeh.io import show, output_file
from bokeh.plotting import figure
from bokeh.models import ColumnDataSource, Range1d, HoverTool, Label, LabelSet, DatetimeTickFormatter, FixedTicker
from bokeh.layouts import column
from operator import itemgetter

verbosity = False
input_path = r'Capital_Inputs.xlsx'


rig_starts = {'HP 273':  datetime(2017, 12, 15,7,0,0), 'HP 329': datetime(2017, 9, 15,7,0,0), 
              'HP 330':  datetime(2017, 6, 20,7,0,0), 'RIG 4':  datetime(2019, 7,1,7,0,0), 
              'RIG 5':  datetime(2020, 1, 1,7,0,0),
              'RIG 6':  datetime(2020, 8, 15,7,0,0),
              'DUC':    None
             }
#'RIG 6':  datetime(2021, 1,1,7,0,0),

input_xl = pd.ExcelFile(input_path)
input_df = input_xl.parse('INPUTS')
scenario = input_df.loc[0, 'SCENARIO']
rig_df = input_xl.parse('SCHEDULE', parse_dates=[11, 12, 13, 14, 15, 16])
well_df = input_xl.parse('WELLS')
rig_dict = {}
pad_dict = {}
well_dict = {}




class Rig(object):

    def __init__(self, rig_name, start, verbose=False):
        self.rig_name = rig_name
        self.start = start
        self.num_pads = 0
        self.pad_list = []
        self.conductors = 14 #days
        self.mob_in = 7    #days
        self.mob_out = 7    #days
        self.verbose = verbose

        if verbose:
            print(self.rig_name, 'created. Project start:', self.start, '\nProject starts with setting conductors, then mob to location, then drill.')

    def info(self):
        """Print the name and start date of the rig."""

        print('Rig Name:\t', self.rig_name)
        print('Start Date:\t', self.start)

    def add_pad(self, Pad):
        """Assign a pad object to the rig. Increments pad count by 1."""

        self.pad_list.append(Pad)
        Pad.rig_name = self.rig_name
        self.num_pads += 1
        if self.verbose:
            print(Pad.pad_name, 'added to', self.rig_name, 'rig.')
    
    def drop_pad(self, Pad):
        """Drop a pad object from the rig. Decrements pad count by 1."""

        self.pad_list.remove(Pad)
        Pad.rig_name = None
        self.num_pads -= 1
        if self.verbose:
            print(Pad.pad_name, 'dropped from', self.rig_name, 'rig.')

    def pads(self):
        """Print the pads assigned to the rig."""

        if len(self.pad_list) == 0:
            print('No pads assigned to the', self.rig_name, 'rig.')
        else:
            for idx, pad in enumerate(self.pad_list):
                print(idx+1, '-', pad.pad_name)

    def pad_count(self):
        """Print the count of pad objects assigned to the rig."""
        
        if self.num_pads == 0:
            print('No pads assigned to the', self.rig_name, 'rig.')
        elif self.num_pads == 1:
            print('1 pad assigned to the', self.rig_name, 'rig.')
        else:
            print(self.num_pads, 'pads assigned to the', self.rig_name, 'rig.')

    def set_conductors(self, set_conductors_days):
        """Set the number of days to set conductors on the pad."""

        self.conductors = float(round(set_conductors_days, 1))
        for pad in self.pad_list:
            pad.conductors = self.conductors
        if self.verbose:
            print('Days to set conductors on the pad set to', str(self.conductors)+'.')

    def set_mob_in(self, mob_in_days):
        """Set the number of days to mob drilling rig to pad."""

        self.mob_in = float(round(mob_in_days, 1))
        for pad in self.pad_list:
            pad.conductors = self.conductors
        if self.verbose:
            print('Days to mob drilling rig to pad set to', str(self.mob_in)+'.')

    def set_mob_out(self, mob_out_days):
        """Set the number of days to mob drilling rig off pad."""

        self.mob_out = float(round(mob_out_days, 1))
        for pad in self.pad_list:
            pad.conductors = self.conductors
        if self.verbose:
            print('Days to mob drilling rig off pad set to', str(self.mob_out)+'.')

    def set_drill_time(self, drill_time=10.0):
        """Set the number of days to drill a well for all wells assigned to pad."""

        for pad in self.pad_list:
            for well in pad.well_list:
                well.set_drill_time(drill_time)

    def set_drill_time_by_depth(self, feet_per_day=1600):
        """Set the number of days to drill a well based on the depth of the well for all wells
           assigned to the pad. The number of feet drilled per day is multiplied by the well's
           depth to calculate the number of days drilled, rounded to the nearest tenth of a day."""

        for pad in self.pad_list:
            for well in pad.well_list:
                well.set_drill_time_by_depth(feet_per_day)

    def schedule(self):
        print('Rig Start:', self.start, '\n')
        for idp, pad in enumerate(self.pad_list):
            print('Mob to', pad.pad_name, 'pad', '('+str(idp+1)+'/'+str(self.num_pads)+')')
            for well in pad.well_list:
                print(well.well_name+':\tDrill: '+str(well.drill_date)+
                      '\tComplete: '+str(well.compl_date)+'\tStart: '+str(well.start_date))
            print('\n')

class Pad(object):

    def __init__(self, pad_name, Rig, verbose=False):
        self.pad_name = pad_name
        self.rig = Rig
        Rig.add_pad(self)
        self.num_wells= 0
        self.well_list = []
        self.build_pad = 30 #days
        self.conductors = self.rig.conductors
        self.mob_in = self.rig.mob_in
        self.mob_out = self.rig.mob_out
        self.log_pad = 0.5          #days/well
        self.build_facilities = 60  #days
        self.build_start = None
        self.build_finish = None
        self.drill_start = None
        self.drill_finish = None
        self.compl_start = None
        self.compl_finish = None
        self.prod_start = None
        self.prod_finish = None
        self.verbose = verbose
    
    def name(self):
        """Print the name of the pad."""

        print(self.pad_name)

    def add_well(self, Well):
        """Assign a well object to the pad. Increments well count by 1."""

        self.well_list.append(Well)
        Well.pad_name = self.pad_name
        self.num_wells += 1
        if self.verbose:
            print(Well.well_name, 'added to', self.pad_name, 'pad.')
    
    def drop_well(self, Well):
        """Drop a well object from the pad. Decrements well count by 1."""

        self.well_list.remove(Well)
        Well.pad_name = None
        self.num_wells -= 1
        if self.verbose:
            print(Well.well_name, 'dropped from', self.pad_name, 'pad.')

    def wells(self):
        """Print the wells assigned to the pad, sorted alphabetically."""

        if len(self.well_list) == 0:
            print('No wells assigned to the', self.pad_name, 'pad.')
        else:
            name_list = sorted([w.well_name for w in self.well_list])
            for idx, well in enumerate(name_list):
                print(idx+1, '-', well)

    def well_count(self):
        """Print the count of well objects assigned to the pad."""
        
        if self.num_wells == 0:
            print('No wells assigned to the', self.pad_name, 'pad.')
        elif self.num_wells == 1:
            print('1 well assigned to the', self.pad_name, 'pad.')
        else:
            print(self.num_wells, 'wells assigned to the', self.pad_name, 'pad.')

    def set_build_time(self, construct_days=30.0):
        """Set the number of days to construct the pad."""
        self.build_pad = float(round(construct_days, 1))
        if self.verbose:
            print('Days to build the', self.pad_name,
                'pad set to', str(self.build_pad)+'.')        
    
    def set_build_time_by_num_wells(self, construct_days_per_well=2.0):
        """Set the number of days to build the pad based on the well count.
           The number of wells assigned to the pad is multiplied by the number of days
           to construct for a single well to calculate the total number of days
           to build the pad, rounded to the nearest tenth of a day."""

        self.build_pad = float(round(construct_days_per_well * self.num_wells, 1))
        if self.verbose:
            print('Days to build the', self.pad_name,
                'pad set to', self.build_pad, 'based on well count.')

    def set_fac_time(self, fac_days=60.0):
        """Set the number of days to build facilities for the pad."""

        self.build_facilities = float(round(fac_days, 1))
        if self.verbose:
            print('Days to build facitilies for the', self.pad_name,
                'pad set to', str(self.build_facilities)+'.')
    
    def set_fac_time_by_num_wells(self, fac_days_per_well=2.0):
        """Set the number of days to build facilities based on the well count.
           The number of wells assigned to the pad is multiplied by the number of days
           to buld facilities for a single well to calculate the total number of days
           to build facilities, rounded to the nearest tenth of a day."""

        self.build_facilities = float(round(fac_days_per_well * self.num_wells, 1))
        if self.verbose:
            print('Days to build facitilies for the', self.pad_name,
                'pad set to', self.build_facilities, 'based on well count.')

    def set_drill_time(self, drill_time=10.0):
        """Set the number of days to drill a well for all wells assigned to pad."""

        for well in self.well_list:
            well.set_drill_time(drill_time)

    def set_drill_time_by_depth(self, feet_per_day=1600):
        """Set the number of days to drill a well based on the depth of the well for all wells
           assigned to the pad. The number of feet drilled per day is multiplied by the well's
           depth to calculate the number of days drilled, rounded to the nearest tenth of a day."""

        for well in self.well_list:
            well.set_drill_time_by_depth(feet_per_day)

    def set_logging_time(self, logging_days):
        """Set the number of days to log wells assigned to the pad."""

        self.log_pad = float(round(logging_days,1))
        if self.verbose:
            print('Days to log the', self.pad_name, 'pad set to', str(self.log_pad)+'.')

    def set_logging_time_per_well(self, logging_days_per_well):
        """Set the number of days to log wells assigned to the pad based on well count."""

        self.log_pad = float(round(logging_days_per_well*self.num_wells,1))
        if self.verbose:
            print('Days to log the', self.pad_name, 'pad set to', str(self.log_pad)+' based on well count.') 

    def set_compl_time(self, compl_days):
        """Set the number of days to complete wells assigned to the pad."""

        for well in self.well_list:
            well.set_compl_time(compl_days)

    def set_flowback_time(self, flowback_days):
        """Set the number of days post completion before first gas for wells assign to the pad."""

        for well in self.well_list:
            well.set_flowback_time(flowback_days)    

class Well(object):

    def __init__(self, propnum, well_name,
                 Pad, area, depth, verbose=False):
        self.propnum = propnum
        self.well_name = well_name
        self.pad = Pad
        Pad.add_well(self)
        self.area = area            #Ops Area
        self.depth = depth          #feet
        self.drill_time = 10.0
        self.compl_time = 3.         #days/well
        self.flowback_time = 30.     #days/well
        self.verbose = verbose
        self.drill_date = None
        self.compl_date = None
        self.start_date = None
        self.drill_cost = 1000.
        self.compl_cost = 2000.

    def info(self):
        """Print the name of the well."""

        print('Well Name:', self.well_name)
        print('Pad Name:', self.pad.pad_name)
        print('Rig Name:', self.pad.rig.rig_name)
    
    def set_drill_time(self, drill_time=10.0):
        """Set the number of days to drill a well."""

        self.drill_time = round(drill_time, 1)
        if self.verbose:
            print('Days to drill the', self.well_name, 'set to', str(self.drill_time)+'.')

    def set_drill_time_by_depth(self, feet_per_day=1600):
        """Set the number of days to drill a well based on the depth of the well.
           The number of feet drilled per day is multiplied by the well's depth to calculate
           the number of days drilled, rounded to the nearest tenth of a day."""

        self.drill_time = round(self.depth / feet_per_day, 1)
        if self.verbose:
            print('Days to drill the', self.well_name, 'set to', self.drill_time, 'based on depth.')

    def set_buffer(self, num_days=0):
        self.buffer = round(num_days, 1)

    def set_compl_time(self, compl_days):
        """Set the number of days to complete a well."""

        self.compl_time = float(round(compl_days, 1))
        if self.verbose:
            print('Days to complete the', self.well_name, 'set to', str(self.compl_time)+'.')

    def set_flowback_time(self, flowback_days):
        """Set the number of days post completion before first gas."""

        self.flowback_time = float(round(flowback_days, 1))
        if self.verbose:
            print('Days post completion before first gas for the', self.well_name, 
                  'set to', str(self.flowback_time)+'.')

def calc_drill_dates(Rig):
    for idp, pad in enumerate(Rig.pad_list):
        if pad.rig.rig_name == 'DUC':
            return None
        if pad.drill_finish == None:
            for idw, well in enumerate(pad.well_list):
                if idp == 0:
                    if idw == 0:
                        if pad.drill_start == None:
                            well.drill_date = Rig.start + timedelta(pad.conductors) + timedelta(pad.mob_in)
                            pad.drill_start = well.drill_date
                        else:
                            well.drill_date = pad.drill_start
                    elif idw == len(pad.well_list)-1:
                        prior_well = pad.well_list[idw-1]
                        well.drill_date = prior_well.drill_date + timedelta(prior_well.drill_time)
                        pad.drill_finish = well.drill_date + timedelta(well.drill_time) + timedelta(pad.mob_out)        
                    else:
                        prior_well = pad.well_list[idw-1]
                        well.drill_date = prior_well.drill_date + timedelta(prior_well.drill_time)
                else:
                    if idw == 0:
                        if pad.drill_start == None:
                            prior_pad_finish = Rig.pad_list[idp-1].drill_finish
                            well.drill_date = prior_pad_finish + timedelta(pad.mob_in)
                            pad.drill_start = well.drill_date
                        else:
                            well.drill_date = pad.drill_start
                    elif idw == len(pad.well_list)-1:
                        prior_well = pad.well_list[idw-1]
                        well.drill_date = prior_well.drill_date + timedelta(prior_well.drill_time)   
                        pad.drill_finish = well.drill_date + timedelta(well.drill_time) + timedelta(pad.mob_out)           
                    else:
                        prior_well = pad.well_list[idw-1]
                        well.drill_date = prior_well.drill_date + timedelta(prior_well.drill_time) 
        elif pad.drill_start != None and pad.drill_finish != None:
            drill_time = (pad.drill_finish - pad.drill_start).days / pad.num_wells
            for idw, well in enumerate(pad.well_list):
                if idw == 0:
                    well.drill_date = pad.drill_start
                elif idw == len(pad.well_list) - 1:
                    well.drill_date = pad.drill_finish - timedelta(drill_time)
                else:
                    well.drill_date = pad.well_list[idw-1].drill_date + timedelta(drill_time)
    #print('Drill dates calculated for', str(Rig.rig_name)+'.')

def calc_compl_dates(Rig):
    for pad in Rig.pad_list:
        for idw, well in enumerate(pad.well_list):
            if pad.compl_finish == None:
                if idw == 0:
                    if pad.compl_start == None:
                        well.compl_date = pad.drill_finish + timedelta(pad.log_pad) + timedelta(pad.build_facilities)
                        pad.compl_start = well.compl_date
                    else:
                        well.compl_date = pad.compl_start
                elif idw == len(pad.well_list)-1:
                    prior_well = pad.well_list[idw-1]
                    well.compl_date = prior_well.compl_date + timedelta(prior_well.compl_time)
                    pad.compl_finish = well.compl_date                
                else:
                    prior_well = pad.well_list[idw-1]
                    well.compl_date = prior_well.compl_date + timedelta(prior_well.compl_time)
            elif pad.compl_start != None and pad.compl_finish != None:
                compl_time = (pad.compl_finish - pad.compl_start).days / pad.num_wells
                for idw, well in enumerate(pad.well_list):
                    if idw == 0:
                        well.compl_date = pad.compl_start
                    elif idw == len(pad.well_list) - 1:
                        well.compl_date = pad.compl_finish - timedelta(compl_time) 
                    else:
                        well.compl_date = pad.well_list[idw-1].compl_date + timedelta(compl_time)                   
    #print('Completion dates calculated for', str(Rig.rig_name)+'.')

def calc_start_dates(Rig):
    for pad in Rig.pad_list:
        for idw, well in enumerate(pad.well_list):
            if pad.prod_finish == None:
                if idw == 0:
                    if pad.prod_start == None:
                        well.start_date = well.compl_date + timedelta(well.compl_time) + timedelta(well.flowback_time)
                        pad.prod_start = well.start_date
                    else:
                        well.start_date = pad.prod_start
                elif idw == len(pad.well_list)-1: 
                    well.start_date = well.compl_date + timedelta(well.compl_time) + timedelta(well.flowback_time)
                    pad.prod_finish = well.start_date
                else:
                    well.start_date = well.compl_date + timedelta(well.compl_time) + timedelta(well.flowback_time)
            elif pad.prod_start != None and pad.prod_finish != None:
                prod_time = (pad.prod_finish - pad.prod_start).days / pad.num_wells
                for idw, well in enumerate(pad.well_list):
                    if idw == 0:
                        well.start_date = pad.compl_start
                    elif idw == len(pad.well_list) - 1:
                        well.start_date = pad.start_finish
                    else:
                        well.start_date = pad.well_list[idw-1].start_date + timedelta(prod_time)                   
    #print('Production start dates calculated for', str(Rig.rig_name)+'.')

def calc_build_pad_dates(Rig):
    for pad in Rig.pad_list:
        if Rig.rig_name != 'DUC':
            if pad.build_start == None and pad.build_finish == None:
                pad.build_start = pad.drill_start - timedelta(pad.build_pad) - timedelta(1)
                pad.build_finish = pad.drill_start - timedelta(1)
            elif pad.build_start == None and pad.build_finish != None:
                pad.build_start = pad.build_finish - timedelta(pad.build_pad)
            elif pad.build_start != None and pad.build_finish == None:
                pad.build_finish = pad.build_start + timedelta(pad.build_pad)
            elif pad.build_start != None and pad.build_finish != None:
                continue

def calc_all_dates(rig_dict):
    for rig in rig_dict.values():
        if rig != 'DUC':
            calc_drill_dates(rig)
            calc_build_pad_dates(rig)
        calc_compl_dates(rig)
        calc_start_dates(rig)
        

def output_exact(rig_dict):
    with open(scenario+'_Exact_Dates.csv','w') as f:
        line = ','.join(['PROPNUM', 'DRILL', 'COMPL', 'START'])
        f.write(line+'\n')
        for rig in rig_dict.values():
            for _, pad in enumerate(rig.pad_list):
                for _, well in enumerate(pad.well_list):
                    if rig.rig_name != 'DUC':
                        line = ','.join([well.propnum, well.drill_date.strftime('%m/%d/%Y'),
                                            well.compl_date.strftime('%m/%d/%Y'),
                                            well.start_date.strftime('%m/%d/%Y')])
                    else:
                        line = ','.join([well.propnum, str(np.nan),
                                            well.compl_date.strftime('%m/%d/%Y'),
                                            well.start_date.strftime('%m/%d/%Y')])
                    f.write(line+'\n')

def output_duc_counts(rig_dict):
    output_list = [[], [], [], []]
    for rig in rig_dict.values():
        for _, pad in enumerate(rig.pad_list):
            for _, well in enumerate(pad.well_list):
                if rig.rig_name == 'DUC':
                    if 0 in output_list[0]:
                        output_list[3][output_list[0].index(0)] += 1
                    else:
                        output_list[0].append(0)
                        output_list[1].append(0)
                        output_list[2].append(0)
                        output_list[3].append(1)
                else:
                    drill_year = well.drill_date.year
                    if drill_year in output_list[0]:
                        output_list[1][output_list[0].index(drill_year)] += 1
                    else:
                        output_list[0].append(drill_year)
                        output_list[1].append(1)
                        output_list[2].append(0)
                        output_list[3].append(0)
                   
                    compl_year = well.compl_date.year
                    if compl_year in output_list[0]:
                        output_list[2][output_list[0].index(compl_year)] += 1
                    else:
                        output_list[0].append(compl_year)
                        output_list[2].append(1)
                        output_list[1].append(0)
                        output_list[3].append(0)
    write_list = []
    for i in range(len(output_list[0])):
        write_list.append([output_list[0][i], output_list[1][i], output_list[2][i], output_list[3][i]])
    output_list = sorted(write_list, key=itemgetter(0))
    for i in range(len(output_list)):
        if i == 0:
            continue
        else:
            output_list[i][3] = output_list[i-1][3]+output_list[i][1]-output_list[i][2]
    with open(scenario+'_DUC_Counts.csv','w') as f:
        line = ','.join(['YEAR', 'DRILL', 'COMPL', 'DUC'])
        f.write(line+'\n')
        for idx in range(len(output_list)):
            line = ','.join(map(str, [output_list[idx][0], output_list[idx][1], output_list[idx][2], output_list[idx][3]]))
            f.write(line+'\n')

def output_aries(rig_dict):
    with open(scenario+'_Aries_Dates.csv','w') as f:
        line = ','.join(['PROPNUM', 'DRILL', 'COMPL', 'START'])
        f.write(line+'\n')
        for rig in rig_dict.values():
            for _, pad in enumerate(rig.pad_list):
                for _, well in enumerate(pad.well_list):
                    if rig.rig_name != 'DUC':
                        if well.drill_date.day >= 15:
                            drill = well.drill_date + relativedelta(months=1)
                        else:
                            drill = well.drill_date
                    else:
                        drill = np.nan
                    if well.compl_date.day >= 15:
                        compl = well.compl_date + relativedelta(months=1)
                    else:
                        compl = well.compl_date
                    if well.start_date.day >= 15:
                        start = well.start_date + relativedelta(months=1)
                    else:
                        start = well.start_date
                    if rig.rig_name != 'DUC':
                        line = ','.join([well.propnum, str(well.drill_cost)+' 0 G '+drill.strftime('%m/%Y')+' AD PC 0',
                                            '0 '+str(well.compl_cost)+' G '+compl.strftime('%m/%Y')+' AD PC 0',
                                            '=\"'+start.strftime('%m/%Y')+'\"'])
                    else:
                        line = ','.join([well.propnum, str(drill),
                                            '0 '+str(well.compl_cost)+' G '+compl.strftime('%m/%Y')+' AD PC 0',
                                            '=\"'+start.strftime('%m/%Y')+'\"'])
                    f.write(line+'\n')

def create_df(rig_dict):
    df = pd.DataFrame()
    rigs = []
    pads = []
    start = []
    start_text = []
    end = []
    end_text = []
    colors = ['#d3d3d3', '#1b9e77','#d95f02','#7570b3']
    color = []
    drill_time = []
    depth = []
    compl_time = []
    well_count = []

    for rig in rig_dict.values():
        for pad in rig.pad_list:
            avg_drill = round(np.mean([well.drill_time for well in pad.well_list]), 1)
            avg_depth = round(np.mean([well.depth for well in pad.well_list]), 0)
            avg_compl = round(np.mean([well.compl_time for well in pad.well_list]), 1)
            if rig.rig_name != 'DUC':
                rigs.append(rig.rig_name+' - Pad')
                pads.append(pad.pad_name)
                start.append(pad.build_start)
                start_text.append(pad.build_start.strftime('%m/%d/%Y'))
                end.append(pad.build_finish)
                end_text.append(pad.build_finish.strftime('%m/%d/%Y'))
                color.append(colors[0])
                drill_time.append(avg_drill)
                depth.append(avg_depth)
                compl_time.append(avg_compl)
                well_count.append(pad.num_wells)

                rigs.append(rig.rig_name+' - Drill')
                pads.append(pad.pad_name)
                start.append(pad.drill_start)
                start_text.append(pad.drill_start.strftime('%m/%d/%Y'))
                end.append(pad.drill_finish)
                end_text.append(pad.drill_finish.strftime('%m/%d/%Y'))
                color.append(colors[1])
                drill_time.append(avg_drill)
                depth.append(avg_depth)
                compl_time.append(avg_compl)
                well_count.append(pad.num_wells)

            rigs.append(rig.rig_name+' - Complete')
            pads.append(pad.pad_name)
            start.append(pad.compl_start)
            start_text.append(pad.compl_start.strftime('%m/%d/%Y'))
            end.append(pad.compl_finish)
            end_text.append(pad.compl_finish.strftime('%m/%d/%Y'))
            color.append(colors[2])
            drill_time.append(avg_drill)
            depth.append(avg_depth)
            compl_time.append(avg_compl)
            well_count.append(pad.num_wells)

            rigs.append(rig.rig_name+' - Sales')
            pads.append(pad.pad_name)
            start.append(pad.prod_start)
            start_text.append(pad.prod_start.strftime('%m/%d/%Y'))
            end.append(pad.prod_finish)
            end_text.append(pad.prod_finish.strftime('%m/%d/%Y'))
            color.append(colors[3])
            drill_time.append(avg_drill)
            depth.append(avg_depth)
            compl_time.append(avg_compl)
            well_count.append(pad.num_wells)    

    df['rig'] = pd.Series(rigs)
    df['pad'] = pd.Series(pads)
    df['start'] = pd.Series(start)
    df['start_text'] = pd.Series(start_text)
    df['end'] = pd.Series(end)
    df['end_text'] = pd.Series(end_text)
    df['color'] = pd.Series(color)
    df['drill'] = pd.Series(drill_time)
    df['depth'] = pd.Series(depth)
    df['compl'] = pd.Series(compl_time)
    df['well_count'] = pd.Series(well_count)
    return df


def gantt_chart(df):
    output_file(scenario+'_Gantt_Chart.html', title=scenario+' Caerus D&C Schedule')
    height = 0.9
    plot_height = 300
    text_font_size = '7pt'
    text_font_style = 'bold'
    angle = 90
    angle_units = 'deg'
    x_offset = 10
    y_offset = 30
    fill_alpha = 0.75


    df1 = df[df.rig.str.contains('DUC')]
    source = ColumnDataSource(df1)
    rigs = list(df1.rig.unique())
    rigs.reverse()
    p1= figure(y_range=rigs, x_axis_type='datetime', x_range=Range1d(date(2017, 10, 1), date(2021, 12, 31)),
               plot_width=1500, plot_height=plot_height, toolbar_location='above', active_drag='pan',
               title='DUC Schedule')
    p1.hbar(y='rig', left='start', right='end', color='color', fill_alpha=fill_alpha, height=height, line_color='gray', source=source)
    hover = HoverTool(
        tooltips=[
                  ('Pad',               '@pad'          ),
                  ('Start',             '@start_text'   ), 
                  ('Finish',            '@end_text'     ),
                  ('Avg Compl Time',    '@compl'        ),
                  ('Well Count',        '@well_count')
        ]
    )
    labels = LabelSet(x='start', y='rig', text='pad', level='glyph', angle=angle, angle_units=angle_units,
              x_offset=x_offset, y_offset=-y_offset, source=source, render_mode='canvas', text_font_size=text_font_size, text_font_style=text_font_style, text_color='black')
    p1.add_tools(hover)
    p1.add_layout(labels)
    p1.ygrid.grid_line_color = None
    p1.xaxis.axis_label = "Date"
    p1.xaxis.formatter=DatetimeTickFormatter(
            months=["%b %Y"],
            days=["%b %e %Y"],
            years=["%b %Y"],
        )
    p1.outline_line_color = None

    df2 = df[df.rig.str.contains('330')]
    source = ColumnDataSource(df2)
    rigs = list(df2.rig.unique())
    rigs.reverse()
    p2 = figure(y_range=rigs, x_axis_type='datetime', x_range=Range1d(date(2017, 10, 1), date(2021, 12, 31)), 
               plot_width=1500, plot_height=plot_height, toolbar_location='above', active_drag='pan',
               title='HP 330 Rig Schedule')
    p2.hbar(y='rig', left='start', right='end', color='color', fill_alpha=fill_alpha, height=height, line_color='gray', source=source)
    hover = HoverTool(
        tooltips=[
                  ('Pad',               '@pad'          ),
                  ('Start',             '@start_text'   ), 
                  ('Finish',            '@end_text'     ),
                  ('Avg Drill Time',    '@drill'        ),
                  ('Avg Depth',         '@depth'        ),
                  ('Avg Compl Time',    '@compl'        ),
                  ('Well Count',        '@well_count')
        ]
    )
    labels = LabelSet(x='start', y='rig', text='pad', level='glyph', angle=angle, angle_units=angle_units,
              x_offset=x_offset, y_offset=-y_offset, source=source, render_mode='canvas', text_font_size=text_font_size, text_font_style=text_font_style, text_color='black')
    p2.add_tools(hover)
    p2.add_layout(labels)
    p2.ygrid.grid_line_color = None
    p2.xaxis.axis_label = "Date"
    p2.xaxis.formatter=DatetimeTickFormatter(
            months=["%b %Y"],
            days=["%b %e %Y"],
            years=["%b %Y"],
        )
    p2.outline_line_color = None

    df3 = df[df.rig.str.contains('329')]
    source = ColumnDataSource(df3)
    rigs = list(df3.rig.unique())
    rigs.reverse()
    p3 = figure(y_range=rigs, x_axis_type='datetime', x_range=Range1d(date(2017, 10, 1), date(2021, 12, 31)), 
               plot_width=1500, plot_height=plot_height, toolbar_location='above', active_drag='pan',
               title='HP 329 Rig Schedule')
    p3.hbar(y='rig', left='start', right='end', color='color', fill_alpha=fill_alpha, height=height, line_color='gray', source=source)
    hover = HoverTool(
        tooltips=[
                  ('Pad',               '@pad'          ),
                  ('Start',             '@start_text'   ), 
                  ('Finish',            '@end_text'     ),
                  ('Avg Drill Time',    '@drill'        ),
                  ('Avg Depth',         '@depth'        ),
                  ('Avg Compl Time',    '@compl'        ),
                  ('Well Count',        '@well_count')
        ]
    )
    labels = LabelSet(x='start', y='rig', text='pad', level='glyph', angle=angle, angle_units=angle_units,
              x_offset=x_offset, y_offset=-y_offset, source=source, render_mode='canvas', text_font_size=text_font_size, text_font_style=text_font_style, text_color='black')
    p3.add_tools(hover)
    p3.add_layout(labels)
    p3.ygrid.grid_line_color = None
    p3.xaxis.axis_label = "Date"
    p3.xaxis.formatter=DatetimeTickFormatter(
            months=["%b %Y"],
            days=["%b %e %Y"],
            years=["%b %Y"],
        )
    p3.outline_line_color = None

    df4 = df[df.rig.str.contains('273')]
    source = ColumnDataSource(df4)
    rigs = list(df4.rig.unique())
    rigs.reverse()
    p4 = figure(y_range=rigs, x_axis_type='datetime', x_range=Range1d(date(2017, 10, 1), date(2021, 12, 31)), 
               plot_width=1500, plot_height=plot_height, toolbar_location='above', active_drag='pan',
               title='HP 273 Rig Schedule')
    p4.hbar(y='rig', left='start', right='end', color='color', fill_alpha=fill_alpha, height=height, line_color='gray', source=source)
    hover = HoverTool(
        tooltips=[
                  ('Pad',               '@pad'          ),
                  ('Start',             '@start_text'   ), 
                  ('Finish',            '@end_text'     ),
                  ('Avg Drill Time',    '@drill'        ),
                  ('Avg Depth',         '@depth'        ),
                  ('Avg Compl Time',    '@compl'        ),
                  ('Well Count',        '@well_count')
        ]
    )
    labels = LabelSet(x='start', y='rig', text='pad', level='glyph', angle=angle, angle_units=angle_units,
              x_offset=x_offset, y_offset=-y_offset, source=source, render_mode='canvas', text_font_size=text_font_size, text_font_style=text_font_style, text_color='black')
    p4.add_tools(hover)
    p4.add_layout(labels)
    p4.ygrid.grid_line_color = None
    p4.xaxis.axis_label = "Date"
    p4.xaxis.formatter=DatetimeTickFormatter(
            months=["%b %Y"],
            days=["%b %e %Y"],
            years=["%b %Y"],
        )
    p4.outline_line_color = None

    df5 = df[df.rig.str.contains('4')]
    source = ColumnDataSource(df5)
    rigs = list(df5.rig.unique())
    rigs.reverse()
    p5 = figure(y_range=rigs, x_axis_type='datetime', x_range=Range1d(date(2017, 10, 1), date(2021, 12, 31)), 
               plot_width=1500, plot_height=plot_height, toolbar_location='above', active_drag='pan',
               title='Rig 4 Schedule')
    p5.hbar(y='rig', left='start', right='end', color='color', fill_alpha=fill_alpha, height=height, line_color='gray', source=source)
    hover = HoverTool(
        tooltips=[
                  ('Pad',               '@pad'          ),
                  ('Start',             '@start_text'   ), 
                  ('Finish',            '@end_text'     ),
                  ('Avg Drill Time',    '@drill'        ),
                  ('Avg Depth',         '@depth'        ),
                  ('Avg Compl Time',    '@compl'        ),
                  ('Well Count',        '@well_count')
        ]
    )
    labels = LabelSet(x='start', y='rig', text='pad', level='glyph', angle=angle, angle_units=angle_units,
              x_offset=x_offset, y_offset=-y_offset, source=source, render_mode='canvas', text_font_size=text_font_size, text_font_style=text_font_style, text_color='black')
    p5.add_tools(hover)
    p5.add_layout(labels)
    p5.ygrid.grid_line_color = None
    p5.xaxis.axis_label = "Date"
    p5.xaxis.formatter=DatetimeTickFormatter(
            months=["%b %Y"],
            days=["%b %e %Y"],
            years=["%b %Y"],
        )
    p5.outline_line_color = None

    df6 = df[df.rig.str.contains('5')]
    source = ColumnDataSource(df6)
    rigs = list(df6.rig.unique())
    rigs.reverse()
    p6 = figure(y_range=rigs, x_axis_type='datetime', x_range=Range1d(date(2017, 10, 1), date(2021, 12, 31)), 
               plot_width=1500, plot_height=plot_height, toolbar_location='above', active_drag='pan',
               title='Rig 5 Schedule')
    p6.hbar(y='rig', left='start', right='end', color='color', fill_alpha=fill_alpha, height=height, line_color='gray', source=source)
    hover = HoverTool(
        tooltips=[
                  ('Pad',               '@pad'          ),
                  ('Start',             '@start_text'   ), 
                  ('Finish',            '@end_text'     ),
                  ('Avg Drill Time',    '@drill'        ),
                  ('Avg Depth',         '@depth'        ),
                  ('Avg Compl Time',    '@compl'        ),
                  ('Well Count',        '@well_count')
        ]
    )
    labels = LabelSet(x='start', y='rig', text='pad', level='glyph', angle=angle, angle_units=angle_units,
              x_offset=x_offset, y_offset=-y_offset, source=source, render_mode='canvas', text_font_size=text_font_size, text_font_style=text_font_style, text_color='black')
    p6.add_tools(hover)
    p6.add_layout(labels)
    p6.ygrid.grid_line_color = None
    p6.xaxis.axis_label = "Date"
    p6.xaxis.formatter=DatetimeTickFormatter(
            months=["%b %Y"],
            days=["%b %e %Y"],
            years=["%b %Y"],
        )
    p6.outline_line_color = None

    df7 = df[df.rig.str.contains('6')]
    source = ColumnDataSource(df7)
    rigs = list(df7.rig.unique())
    rigs.reverse()
    p7 = figure(y_range=rigs, x_axis_type='datetime', x_range=Range1d(date(2017, 10, 1), date(2021, 12, 31)), 
               plot_width=1500, plot_height=plot_height, toolbar_location='above', active_drag='pan',
               title='Rig 6 Schedule')
    p7.hbar(y='rig', left='start', right='end', color='color', fill_alpha=fill_alpha, height=height, line_color='gray', source=source)
    hover = HoverTool(
        tooltips=[
                  ('Pad',               '@pad'          ),
                  ('Start',             '@start_text'   ), 
                  ('Finish',            '@end_text'     ),
                  ('Avg Drill Time',    '@drill'        ),
                  ('Avg Depth',         '@depth'        ),
                  ('Avg Compl Time',    '@compl'        ),
                  ('Well Count',        '@well_count')
        ]
    )
    labels = LabelSet(x='start', y='rig', text='pad', level='glyph', angle=angle, angle_units=angle_units,
              x_offset=x_offset, y_offset=-y_offset, source=source, render_mode='canvas', text_font_size=text_font_size, text_font_style=text_font_style, text_color='black')
    p7.add_tools(hover)
    p7.add_layout(labels)
    p7.ygrid.grid_line_color = None
    p7.xaxis.axis_label = "Date"
    p7.xaxis.formatter=DatetimeTickFormatter(
            months=["%b %Y"],
            days=["%b %e %Y"],
            years=["%b %Y"],
        )
    p7.outline_line_color = None
    
    p = column(p1, p2, p3, p4, p5,p6,p7)

    show(p)

for rigname in rig_df.RIG.unique():
    start = rig_starts[rigname]
    rig_dict[rigname] = Rig(rigname, start, verbosity)
for index, row in rig_df.iterrows():
    pad_dict[row['PAD']] = Pad(row['PAD'], rig_dict[row['RIG']], verbosity)
for index, row in well_df.iterrows():
    well_dict[row['PROPNUM']] = Well(row['PROPNUM'], row['WELLNAME'], pad_dict[row['PAD']], row['AREA'], row['DEPTH'], verbosity)
    well_dict[row['PROPNUM']].drill_cost = row['DRILL_COST']
    well_dict[row['PROPNUM']].compl_cost = row['COMPL_COST']

for rigname in rig_df.RIG.unique():
    df = rig_df[rig_df.RIG == rigname]
    if rigname != 'DUC':
        for index, row in df.iterrows():
            if input_df[(input_df['STYLE'] == row['STYLE']) & (input_df['UNIT'] == 'PER_WELL')]['CONDUCTORS'].values[0] == 0:
                rig_dict[row['RIG']].set_conductors(input_df[(input_df['STYLE'] == row['STYLE']) & (input_df['UNIT'] == 'FIXED')]['CONDUCTORS'].values[0])
            else:
                rig_dict[row['RIG']].set_conductors(input_df[(input_df['STYLE'] == row['STYLE']) & (input_df['UNIT'] == 'PER_WELL')]['CONDUCTORS'].values[0])
            if input_df[(input_df['STYLE'] == row['STYLE']) & (input_df['UNIT'] == 'PER_WELL')]['MOB_IN'].values[0] == 0:
                rig_dict[row['RIG']].set_mob_in(input_df[(input_df['STYLE'] == row['STYLE']) & (input_df['UNIT'] == 'FIXED')]['MOB_IN'].values[0])
            else:
                rig_dict[row['RIG']].set_mob_in(input_df[(input_df['STYLE'] == row['STYLE']) & (input_df['UNIT'] == 'PER_WELL')]['MOB_IN'].values[0])
            if input_df[(input_df['STYLE'] == row['STYLE']) & (input_df['UNIT'] == 'PER_WELL')]['DRILL'].values[0] == 0:
                pad_dict[row['PAD']].set_drill_time(input_df[(input_df['STYLE'] == row['STYLE']) & (input_df['UNIT'] == 'FIXED')]['DRILL'].values[0])
            else:
                pad_dict[row['PAD']].set_drill_time_by_depth(input_df[(input_df['STYLE'] == row['STYLE']) & (input_df['UNIT'] == 'PER_WELL')]['DRILL'].values[0])
            if input_df[(input_df['STYLE'] == row['STYLE']) & (input_df['UNIT'] == 'PER_WELL')]['MOB_OUT'].values[0] == 0:
                rig_dict[row['RIG']].set_mob_out(input_df[(input_df['STYLE'] == row['STYLE']) & (input_df['UNIT'] == 'FIXED')]['MOB_OUT'].values[0])
            else:
                rig_dict[row['RIG']].set_mob_out(input_df[(input_df['STYLE'] == row['STYLE']) & (input_df['UNIT'] == 'PER_WELL')]['MOB_OUT'].values[0])

        for index, row in df.iterrows():
            if not pd.isnull(row['CONDUCTORS']):
                pad_dict[row['PAD']].conductors = row['CONDUCTORS']
            if not pd.isnull(row['MOB_IN']):
                pad_dict[row['PAD']].mob_in = row['MOB_IN']
            if not pd.isnull(row['DRILL']):
                pad_dict[row['PAD']].set_drill_time(row['DRILL'])
            if not pd.isnull(row['MOB_OUT']):
                pad_dict[row['PAD']].mob_out = row['MOB_OUT']
            if not pd.isnull(row['DRILL_START']):
                pad_dict[row['PAD']].drill_start = row['DRILL_START']
            if not pd.isnull(row['DRILL_END']):
                pad_dict[row['PAD']].drill_finish = row['DRILL_END']

    for index, row in df.iterrows():
        if input_df[(input_df['STYLE'] == row['STYLE']) & (input_df['UNIT'] == 'PER_WELL')]['LOG'].values[0] == 0:
            pad_dict[row['PAD']].set_logging_time(input_df[(input_df['STYLE'] == row['STYLE']) & (input_df['UNIT'] == 'FIXED')]['LOG'].values[0])
        else:
            pad_dict[row['PAD']].set_logging_time_per_well(input_df[(input_df['STYLE'] == row['STYLE']) & (input_df['UNIT'] == 'PER_WELL')]['LOG'].values[0])

        if input_df[(input_df['STYLE'] == row['STYLE']) & (input_df['UNIT'] == 'PER_WELL')]['FAC'].values[0] == 0:
            pad_dict[row['PAD']].set_fac_time(input_df[(input_df['STYLE'] == row['STYLE']) & (input_df['UNIT'] == 'FIXED')]['FAC'].values[0])
        else:
            pad_dict[row['PAD']].set_fac_time_by_num_wells(input_df[(input_df['STYLE'] == row['STYLE']) & (input_df['UNIT'] == 'PER_WELL')]['FAC'].values[0])

        if input_df[(input_df['STYLE'] == row['STYLE']) & (input_df['UNIT'] == 'PER_WELL')]['FRAC'].values[0] == 0:
            pad_dict[row['PAD']].set_compl_time(input_df[(input_df['STYLE'] == row['STYLE']) & (input_df['UNIT'] == 'FIXED')]['FRAC'].values[0])
        else:
            pad_dict[row['PAD']].set_compl_time(input_df[(input_df['STYLE'] == row['STYLE']) & (input_df['UNIT'] == 'PER_WELL')]['FRAC'].values[0])
        
        if input_df[(input_df['STYLE'] == row['STYLE']) & (input_df['UNIT'] == 'PER_WELL')]['FLOWBACK'].values[0] == 0:
            pad_dict[row['PAD']].set_flowback_time(input_df[(input_df['STYLE'] == row['STYLE']) & (input_df['UNIT'] == 'FIXED')]['FLOWBACK'].values[0])
        else:
            pad_dict[row['PAD']].set_flowback_time(input_df[(input_df['STYLE'] == row['STYLE']) & (input_df['UNIT'] == 'PER_WELL')]['FLOWBACK'].values[0])
        if input_df[(input_df['STYLE'] == row['STYLE']) & (input_df['UNIT'] == 'PER_WELL')]['PAD'].values[0] == 0:
            pad_dict[row['PAD']].set_build_time(input_df[(input_df['STYLE'] == row['STYLE']) & (input_df['UNIT'] == 'FIXED')]['PAD'].values[0])
        else:
            pad_dict[row['PAD']].set_build_time_by_num_wells(input_df[(input_df['STYLE'] == row['STYLE']) & (input_df['UNIT'] == 'PER_WELL')]['PAD'].values[0])

    for index, row in df.iterrows():
        if not pd.isnull(row['LOG']):
            pad_dict[row['PAD']].set_logging_time(row['LOG'])
        if not pd.isnull(row['FAC']):
            pad_dict[row['PAD']].set_fac_time(row['FAC'])
        if not pd.isnull(row['FRAC']):
            pad_dict[row['PAD']].set_compl_time(row['FRAC'])
        if not pd.isnull(row['FLOWBACK']):
            pad_dict[row['PAD']].set_flowback_time(row['FLOWBACK'])
        if not pd.isnull(row['COMPL_START']):
            pad_dict[row['PAD']].compl_start = row['COMPL_START']
        if not pd.isnull(row['COMPL_END']):
            pad_dict[row['PAD']].compl_finish = row['COMPL_END']
        if not pd.isnull(row['PROD_START']):
            pad_dict[row['PAD']].prod_start = row['PROD_START']
        if not pd.isnull(row['PROD_END']):
            pad_dict[row['PAD']].prod_finish = row['PROD_END']
        if not pd.isnull(row['PAD_START']):
            pad_dict[row['PAD']].build_start = row['PAD_START']
        if not pd.isnull(row['PAD_END']):
            pad_dict[row['PAD']].build_finish = row['PAD_END']            

calc_all_dates(rig_dict)
output_aries(rig_dict)
output_exact(rig_dict)
output_duc_counts(rig_dict)
df = create_df(rig_dict)
gantt_chart(df)

df.to_csv('Gantt Chart data.csv')
