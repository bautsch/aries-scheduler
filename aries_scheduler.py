from datetime import date
from datetime import timedelta
from operator import itemgetter
from dateutil.relativedelta import relativedelta
import pandas as pd
import numpy as np
from bokeh.io import show
from bokeh.io import output_file
from bokeh.plotting import figure
from bokeh.models import ColumnDataSource
from bokeh.models import Range1d
from bokeh.models import HoverTool
from bokeh.models import LabelSet
from bokeh.models import DatetimeTickFormatter
from bokeh.layouts import column


class Schedule():
    def __init__(self, schedule_name, schedule_file_path,
                 output_file_path=None, gantt_start_date=date(2019, 1, 1),
                 verbose=False):
        self.schedule_name = schedule_name
        self.schedule_path = schedule_file_path
        self.output_path = output_file_path
        self.gantt_start_date = gantt_start_date
        self.gantt_end_date = gantt_start_date + relativedelta(years=+3)
        self.input_df = pd.read_excel(self.schedule_path,
                                      sheet_name='INPUTS')
        self.schedule_df = pd.read_excel(self.schedule_path,
                                         sheet_name='SCHEDULE',
                                         parse_dates=[11, 12, 13, 14, 15, 16])
        self.well_df = pd.read_excel(self.schedule_path,
                                     sheet_name='WELLS')
        self.verbosity = verbose
        self.rig_dict = {}
        self.pad_dict = {}
        self.well_dict = {}
        self.exact_dates = pd.DataFrame(columns=['WELL_ID', 'DRILL',
                                                 'COMPLETE', 'START'])
        self.gantt_df = None
        self.build_dictionaries()
        self.parse_inputs()
        self.parse_schedule_inputs()
        self.calc_all_dates()
        self.output_exact()
        self.output_aries()
        self.output_duc_counts()
        self.create_gantt_df()
        self.create_gantt_chart()
        print(self.schedule_name, 'schedule created.')

    def build_dictionaries(self):
        for rigname in self.schedule_df['RIG'].unique():
            if (rigname != 'DUC' and self.schedule_df.loc[
                    (self.schedule_df['RIG'] == rigname) &
                    (self.schedule_df['RIG_START'].notnull()),
                    'RIG_START'].size == 0):
                print('Enter a rig start date in the schedule.')
                return
            if (rigname != 'DUC' and self.schedule_df.loc[(self.schedule_df[
                    'RIG'] == rigname) &
                    (self.schedule_df['RIG_START'].notnull()),
                    'RIG_START'].size > 1):
                print('Too many rig start dates for', rigname, '.')
                return
            if rigname != 'DUC':
                rig_start = self.schedule_df.loc[(self.schedule_df[
                    'RIG'] == rigname) &
                    (self.schedule_df['RIG_START'].notnull()),
                    'RIG_START'].values[0]
                rig_start = pd.Timestamp(rig_start)
            else:
                rig_start = None
            self.rig_dict[rigname] = Rig(rigname,
                                         rig_start,
                                         self.verbosity)
        for _, row in self.schedule_df.iterrows():
            self.pad_dict[row['PAD']] = Pad(row['PAD'],
                                            self.rig_dict[row['RIG']],
                                            self.verbosity)
        for _, row in self.well_df.iterrows():
            self.well_dict[row['PROPNUM']] = Well(row['PROPNUM'],
                                                  row['WELLNAME'],
                                                  self.pad_dict[row['PAD']],
                                                  row['AREA'],
                                                  row['DEPTH'], self.verbosity)
            self.well_dict[row['PROPNUM']].drill_cost = row['DRILL_COST']
            self.well_dict[row['PROPNUM']].compl_cost = row['COMPL_COST']
        well_id = pd.Series([well.propnum for well in self.well_dict.values()])
        self.exact_dates['WELL_ID'] = well_id

    def parse_inputs(self):
        for rigname in self.schedule_df['RIG'].unique():
            schedule_df_subset = self.schedule_df[
                self.schedule_df.RIG == rigname]
            for _, row in schedule_df_subset.iterrows():
                if self.input_df[(self.input_df[
                        'STYLE'] == row['STYLE']) &
                                 (self.input_df['UNIT'] == 'PER_WELL')][
                                     'CONDUCTORS'].values[0] == 0:
                    self.rig_dict[
                        row['RIG']].set_conductors(
                            self.input_df[(self.input_df[
                                'STYLE'] == row[
                                    'STYLE']) & (self.input_df[
                                        'UNIT'] == 'FIXED')][
                                            'CONDUCTORS'].values[0])
                else:
                    self.rig_dict[row['RIG']].set_conductors(
                        self.input_df[(self.input_df['STYLE'] == row[
                            'STYLE']) & (self.input_df[
                                'UNIT'] == 'PER_WELL')][
                                    'CONDUCTORS'].values[0])
                if self.input_df[(self.input_df['STYLE'] == row[
                        'STYLE']) & (self.input_df['UNIT'] == 'PER_WELL')][
                            'MOB_IN'].values[0] == 0:
                    self.rig_dict[row['RIG']].set_mob_in(self.input_df[
                        (self.input_df['STYLE'] == row['STYLE']) &
                        (self.input_df['UNIT'] == 'FIXED')][
                            'MOB_IN'].values[0])
                else:
                    self.rig_dict[row['RIG']].set_mob_in(self.input_df[
                        (self.input_df['STYLE'] == row['STYLE']) &
                        (self.input_df['UNIT'] == 'PER_WELL')][
                            'MOB_IN'].values[0])
                if self.input_df[(self.input_df['STYLE'] == row['STYLE']) &
                                 (self.input_df['UNIT'] == 'PER_WELL')][
                                     'DRILL'].values[0] == 0:
                    self.pad_dict[row['PAD']].set_drill_time(self.input_df[
                        (self.input_df['STYLE'] == row['STYLE']) &
                        (self.input_df['UNIT'] == 'FIXED')][
                            'DRILL'].values[0])
                else:
                    self.pad_dict[row['PAD']].set_drill_time_by_depth(
                        self.input_df[(self.input_df['STYLE'] == row[
                            'STYLE']) & (self.input_df[
                                'UNIT'] == 'PER_WELL')][
                                    'DRILL'].values[0])
                if self.input_df[(self.input_df['STYLE'] == row['STYLE']) &
                                 (self.input_df['UNIT'] == 'PER_WELL')][
                                     'MOB_OUT'].values[0] == 0:
                    self.rig_dict[row['RIG']].set_mob_out(self.input_df[
                        (self.input_df['STYLE'] == row['STYLE']) &
                        (self.input_df['UNIT'] == 'FIXED')][
                            'MOB_OUT'].values[0])
                else:
                    self.rig_dict[row['RIG']].set_mob_out(self.input_df[
                        (self.input_df['STYLE'] == row['STYLE']) &
                        (self.input_df['UNIT'] == 'PER_WELL')][
                            'MOB_OUT'].values[0])

    def parse_schedule_inputs(self):
        for rigname in self.schedule_df.RIG.unique():
            schedule_df_subset = self.schedule_df[
                self.schedule_df.RIG == rigname]
            if rigname != 'DUC':
                for _, row in schedule_df_subset.iterrows():
                    if not pd.isnull(row['CONDUCTORS']):
                        self.pad_dict[row[
                            'PAD']].conductors = row['CONDUCTORS']
                    if not pd.isnull(row['MOB_IN']):
                        self.pad_dict[row['PAD']].mob_in = row['MOB_IN']
                    if not pd.isnull(row['DRILL']):
                        self.pad_dict[row['PAD']].set_drill_time(row['DRILL'])
                    if not pd.isnull(row['MOB_OUT']):
                        self.pad_dict[row['PAD']].mob_out = row['MOB_OUT']
                    if not pd.isnull(row['DRILL_START']):
                        self.pad_dict[row[
                            'PAD']].drill_start = row['DRILL_START']
                    if not pd.isnull(row['DRILL_END']):
                        self.pad_dict[row[
                            'PAD']].drill_finish = row['DRILL_END']
            else:
                for _, row in schedule_df_subset.iterrows():
                    self.pad_dict[row['PAD']].compl_start = row['COMPL_START']
            for _, row in schedule_df_subset.iterrows():
                if self.input_df[(self.input_df['STYLE'] == row['STYLE']) &
                                 (self.input_df['UNIT'] == 'PER_WELL')][
                                     'LOG'].values[0] == 0:
                    self.pad_dict[row['PAD']].set_logging_time(self.input_df[
                        (self.input_df['STYLE'] == row['STYLE']) &
                        (self.input_df['UNIT'] == 'FIXED')]['LOG'].values[0])
                else:
                    self.pad_dict[row['PAD']].set_logging_time_per_well(
                        self.input_df[(self.input_df['STYLE'] == row['STYLE'])
                                      & (self.input_df['UNIT'] == 'PER_WELL')][
                                          'LOG'].values[0])

                if self.input_df[(self.input_df['STYLE'] == row['STYLE'])
                                 & (self.input_df['UNIT'] == 'PER_WELL')][
                                     'FAC'].values[0] == 0:
                    self.pad_dict[row['PAD']].set_fac_time(self.input_df[
                        (self.input_df['STYLE'] == row['STYLE']) &
                        (self.input_df['UNIT'] == 'FIXED')][
                            'FAC'].values[0])
                else:
                    self.pad_dict[row['PAD']].set_fac_time_by_num_wells(
                        self.input_df[(self.input_df['STYLE'] == row['STYLE'])
                                      & (self.input_df[
                                          'UNIT'] == 'PER_WELL')][
                                              'FAC'].values[0])

                if self.input_df[(self.input_df['STYLE'] == row['STYLE'])
                                 & (self.input_df['UNIT'] == 'PER_WELL')][
                                     'FRAC'].values[0] == 0:
                    self.pad_dict[row['PAD']].set_compl_time(self.input_df[
                        (self.input_df['STYLE'] == row['STYLE']) &
                        (self.input_df['UNIT'] == 'FIXED')][
                            'FRAC'].values[0])
                else:
                    self.pad_dict[row['PAD']].set_compl_time(self.input_df[
                        (self.input_df['STYLE'] == row['STYLE']) &
                        (self.input_df['UNIT'] == 'PER_WELL')][
                            'FRAC'].values[0])

                if self.input_df[(self.input_df['STYLE'] == row['STYLE']) &
                                 (self.input_df['UNIT'] == 'PER_WELL')][
                                     'FLOWBACK'].values[0] == 0:
                    self.pad_dict[row['PAD']].set_flowback_time(self.input_df[
                        (self.input_df['STYLE'] == row['STYLE']) &
                        (self.input_df['UNIT'] == 'FIXED')][
                            'FLOWBACK'].values[0])
                else:
                    self.pad_dict[row['PAD']].set_flowback_time(self.input_df[
                        (self.input_df['STYLE'] == row['STYLE']) &
                        (self.input_df['UNIT'] == 'PER_WELL')][
                            'FLOWBACK'].values[0])
                if self.input_df[(self.input_df['STYLE'] == row['STYLE']) &
                                 (self.input_df['UNIT'] == 'PER_WELL')][
                                     'PAD'].values[0] == 0:
                    self.pad_dict[row['PAD']].set_build_time(self.input_df[(
                        (self.input_df['STYLE'] == row['STYLE']) &
                        (self.input_df['UNIT'] == 'FIXED'))]['PAD'].values[0])
                else:
                    self.pad_dict[row['PAD']].set_build_time_by_num_wells(
                        self.input_df[(self.input_df['STYLE'] == row['STYLE'])
                                      & (self.input_df['UNIT'] == 'PER_WELL')][
                                          'PAD'].values[0])
            for _, row in schedule_df_subset.iterrows():
                if not pd.isnull(row['LOG']):
                    self.pad_dict[row['PAD']].set_logging_time(row['LOG'])
                if not pd.isnull(row['FAC']):
                    self.pad_dict[row['PAD']].set_fac_time(row['FAC'])
                if not pd.isnull(row['FRAC']):
                    self.pad_dict[row['PAD']].set_compl_time(row['FRAC'])
                if not pd.isnull(row['FLOWBACK']):
                    self.pad_dict[row[
                        'PAD']].set_flowback_time(row['FLOWBACK'])
                self.pad_dict[row['PAD']].set_pod_size(row['POD_SIZE'])
                if not pd.isnull(row['COMPL_START']):
                    self.pad_dict[row['PAD']].compl_start = row['COMPL_START']
                if not pd.isnull(row['COMPL_END']):
                    self.pad_dict[row['PAD']].compl_finish = row['COMPL_END']
                if not pd.isnull(row['PROD_START']):
                    self.pad_dict[row['PAD']].prod_start = row['PROD_START']
                if not pd.isnull(row['PROD_END']):
                    self.pad_dict[row['PAD']].prod_finish = row['PROD_END']
                if not pd.isnull(row['PAD_START']):
                    self.pad_dict[row['PAD']].build_start = row['PAD_START']
                if not pd.isnull(row['PAD_END']):
                    self.pad_dict[row['PAD']].build_finish = row['PAD_END']

    def calc_drill_dates(self, rig):
        for idp, pad in enumerate(rig.pad_list):
            if pad.rig.rig_name == 'DUC':
                return None
            if pad.drill_finish is None:
                for idw, well in enumerate(pad.well_list):
                    if idp == 0:
                        if idw == 0:
                            if pad.drill_start is None:
                                well.drill_date = (rig.start
                                                   + timedelta(pad.conductors)
                                                   + timedelta(pad.mob_in))
                                pad.drill_start = well.drill_date
                                self.exact_dates.loc[self.exact_dates['WELL_ID'] == well.propnum,
                                                     'DRILL'] = well.drill_date
                            else:
                                well.drill_date = pad.drill_start
                                self.exact_dates.loc[self.exact_dates['WELL_ID'] == well.propnum,
                                                     'DRILL'] = well.drill_date
                        elif idw == len(pad.well_list)-1:
                            prior_well = pad.well_list[idw-1]
                            well.drill_date = (prior_well.drill_date
                                               + timedelta(
                                                   prior_well.drill_time))
                            pad.drill_finish = (well.drill_date
                                                + timedelta(well.drill_time)
                                                + timedelta(pad.mob_out))
                            self.exact_dates.loc[self.exact_dates['WELL_ID'] == well.propnum,
                                                 'DRILL'] = well.drill_date
                        else:
                            prior_well = pad.well_list[idw-1]
                            well.drill_date = (prior_well.drill_date
                                               + timedelta(
                                                   prior_well.drill_time))
                            self.exact_dates.loc[self.exact_dates['WELL_ID'] == well.propnum,
                                                 'DRILL'] = well.drill_date
                    else:
                        if idw == 0:
                            if pad.drill_start is None:
                                prior_pad_finish = rig.pad_list[
                                    idp-1].drill_finish
                                well.drill_date = (prior_pad_finish
                                                   + timedelta(pad.mob_in))
                                pad.drill_start = well.drill_date
                                self.exact_dates.loc[self.exact_dates['WELL_ID'] == well.propnum,
                                                     'DRILL'] = well.drill_date
                            else:
                                well.drill_date = pad.drill_start
                                self.exact_dates.loc[self.exact_dates['WELL_ID'] == well.propnum,
                                                     'DRILL'] = well.drill_date
                        elif idw == len(pad.well_list)-1:
                            prior_well = pad.well_list[idw-1]
                            well.drill_date = (prior_well.drill_date
                                               + timedelta(
                                                   prior_well.drill_time))
                            pad.drill_finish = (well.drill_date
                                                + timedelta(well.drill_time)
                                                + timedelta(pad.mob_out))
                            self.exact_dates.loc[self.exact_dates['WELL_ID'] == well.propnum,
                                                 'DRILL'] = well.drill_date
                        else:
                            prior_well = pad.well_list[idw-1]
                            well.drill_date = (prior_well.drill_date
                                               + timedelta(
                                                   prior_well.drill_time))
                            self.exact_dates.loc[self.exact_dates['WELL_ID'] == well.propnum,
                                                 'DRILL'] = well.drill_date
            elif pad.drill_start is not None and pad.drill_finish is not None:
                drill_time = ((pad.drill_finish - pad.drill_start).days
                              / pad.num_wells)
                for idw, well in enumerate(pad.well_list):
                    if idw == 0:
                        well.drill_date = pad.drill_start
                        self.exact_dates.loc[self.exact_dates['WELL_ID'] == well.propnum,
                                             'DRILL'] = well.drill_date
                    elif idw == len(pad.well_list) - 1:
                        well.drill_date = (pad.drill_finish
                                           - timedelta(drill_time))
                        self.exact_dates.loc[self.exact_dates['WELL_ID'] == well.propnum,
                                             'DRILL'] = well.drill_date
                    else:
                        well.drill_date = (pad.well_list[idw-1].drill_date
                                           + timedelta(drill_time))
                        self.exact_dates.loc[self.exact_dates['WELL_ID'] == well.propnum,
                                             'DRILL'] = well.drill_date

    def calc_compl_dates(self, rig):
        for pad in rig.pad_list:
            for idw, well in enumerate(pad.well_list):
                if rig.rig_name == 'DUC':
                    if pad.compl_finish is None:
                        if idw == 0:
                            well.compl_date = pad.compl_start
                            self.exact_dates.loc[self.exact_dates['WELL_ID'] == well.propnum,
                                                 'COMPLETE'] = well.compl_date

                        elif idw == len(pad.well_list)-1:
                            prior_well = pad.well_list[idw-1]
                            well.compl_date = (prior_well.compl_date
                                               + timedelta(
                                                   prior_well.compl_time))
                            pad.compl_finish = well.compl_date
                            self.exact_dates.loc[self.exact_dates['WELL_ID'] == well.propnum,
                                                 'COMPLETE'] = well.compl_date
                        else:
                            prior_well = pad.well_list[idw-1]
                            well.compl_date = (prior_well.compl_date
                                               + timedelta(
                                                   prior_well.compl_time))
                            self.exact_dates.loc[self.exact_dates['WELL_ID'] == well.propnum,
                                                 'COMPLETE'] = well.compl_date
                    else:
                        compl_time = ((pad.compl_finish - pad.compl_start).days
                                      / pad.num_wells)
                        for idw2, well2 in enumerate(pad.well_list):
                            if idw2 == 0:
                                well2.compl_date = pad.compl_start
                                self.exact_dates.loc[self.exact_dates['WELL_ID'] ==
                                                     well2.propnum,
                                                     'COMPLETE'] = well2.compl_date
                            elif idw2 == len(pad.well_list) - 1:
                                well2.compl_date = (pad.compl_finish
                                                    - timedelta(compl_time))
                                self.exact_dates.loc[self.exact_dates['WELL_ID'] ==
                                                     well2.propnum,
                                                     'COMPLETE'] = well2.compl_date
                            else:
                                well2.compl_date = (pad.well_list[
                                    idw2-1].compl_date +
                                                    timedelta(compl_time))
                                self.exact_dates.loc[self.exact_dates['WELL_ID'] ==
                                                     well2.propnum,
                                                 'COMPLETE'] = well2.compl_date

                elif pad.compl_finish is None:
                    if idw == 0:
                        if pad.compl_start is None:
                            well.compl_date = (pad.drill_finish
                                               + timedelta(pad.log_pad)
                                               + timedelta(
                                                   pad.build_facilities))
                            pad.compl_start = well.compl_date
                            self.exact_dates.loc[self.exact_dates['WELL_ID'] == well.propnum,
                                                 'COMPLETE'] = well.compl_date
                        else:
                            well.compl_date = pad.compl_start
                            self.exact_dates.loc[self.exact_dates['WELL_ID'] == well.propnum,
                                                 'COMPLETE'] = well.compl_date
                    elif idw == len(pad.well_list)-1:
                        prior_well = pad.well_list[idw-1]
                        well.compl_date = (prior_well.compl_date
                                           + timedelta(prior_well.compl_time))
                        pad.compl_finish = well.compl_date
                        self.exact_dates.loc[self.exact_dates['WELL_ID'] == well.propnum,
                                                 'COMPLETE'] = well.compl_date
                    else:
                        prior_well = pad.well_list[idw-1]
                        well.compl_date = (prior_well.compl_date
                                           + timedelta(prior_well.compl_time))
                        self.exact_dates.loc[self.exact_dates['WELL_ID'] == well.propnum,
                                             'COMPLETE'] = well.compl_date
                elif (pad.compl_start is not None and
                      pad.compl_finish is not None):
                    compl_time = ((pad.compl_finish - pad.compl_start).days
                                  / pad.num_wells)
                    for idw2, well2 in enumerate(pad.well_list):
                        if idw2 == 0:
                            well2.compl_date = pad.compl_start
                            self.exact_dates.loc[self.exact_dates['WELL_ID'] ==
                                                 well2.propnum,
                                                 'COMPLETE'] == well2.compl_date 
                        elif idw2 == len(pad.well_list) - 1:
                            well2.compl_date = (pad.compl_finish
                                                - timedelta(compl_time))
                            self.exact_dates.loc[self.exact_dates['WELL_ID'] ==
                                                 well2.propnum,
                                                 'COMPLETE'] == well2.compl_date
                        else:
                            well2.compl_date = (pad.well_list[
                                idw2-1].compl_date +
                                                timedelta(compl_time))
                            self.exact_dates.loc[self.exact_dates['WELL_ID'] ==
                                                 well2.propnum,
                                                 'COMPLETE'] == well2.compl_date

    def calc_start_dates(self, rig):
        for pad in rig.pad_list:
            well_ids = np.arange(1, len(pad.well_list)+1)
            if pad.pod_size > 0:
                pods = [w for w in well_ids
                        if (w % pad.pod_size == 0) and (w > 0)]
                if max(well_ids) > max(pods):
                    pods.append(max(well_ids))
            for idw, well in enumerate(pad.well_list):
                if pad.prod_finish is None:
                    if idw == 0:
                        if pad.prod_start is None:
                            if pad.pod_size == 0:
                                start = (well.compl_date
                                         + timedelta(well.compl_time)
                                         + timedelta(well.flowback_time))
                                well.start_date = start
                                pad.prod_start = well.start_date
                                self.exact_dates.loc[self.exact_dates['WELL_ID'] ==
                                                                      well.propnum,
                                                                      'START'] = well.start_date
                            else:
                                tmp_pods = [p - idw for p in pods]
                                idx_tmp = min(t for t in tmp_pods if t > 0)
                                idx_tmp = tmp_pods.index(idx_tmp)
                                idx_tmp = pods[idx_tmp] - 1
                                idx_tmp = pad.pod_size - 1
                                pod_compl = pad.well_list[idx_tmp].compl_date
                                pod_prod = (pod_compl
                                            + timedelta(well.compl_time)
                                            + timedelta(well.flowback_time))
                                well.start_date = pod_prod
                                pad.prod_start = well.start_date
                                self.exact_dates.loc[self.exact_dates['WELL_ID'] ==
                                                                      well.propnum,
                                                                      'START'] = well.start_date
                        else:
                            well.start_date = pad.prod_start
                            self.exact_dates.loc[self.exact_dates['WELL_ID'] ==
                                                                  well.propnum,
                                                                  'START'] = well.start_date
                    elif idw == len(pad.well_list)-1:
                        if pad.pod_size == 0:
                            well.start_date = (well.compl_date
                                               + timedelta(well.compl_time)
                                               + timedelta(well.flowback_time))
                            pad.prod_finish = well.start_date
                            self.exact_dates.loc[self.exact_dates['WELL_ID'] ==
                                                                  well.propnum,
                                                                  'START'] = well.start_date
                        else:
                            tmp_pods = [p - idw for p in pods]
                            idx_tmp = min(t for t in tmp_pods if t > 0)
                            idx_tmp = tmp_pods.index(idx_tmp)
                            idx_tmp = pods[idx_tmp] - 1
                            pod_compl = pad.well_list[idx_tmp].compl_date
                            pod_prod = (pod_compl
                                        + timedelta(well.compl_time)
                                        + timedelta(well.flowback_time))
                            well.start_date = pod_prod
                            pad.prod_finish = well.start_date
                            self.exact_dates.loc[self.exact_dates['WELL_ID'] ==
                                                                  well.propnum,
                                                                  'START'] = well.start_date
                    else:
                        if pad.pod_size == 0:
                            well.start_date = (well.compl_date
                                               + timedelta(well.compl_time)
                                               + timedelta(well.flowback_time))
                            self.exact_dates.loc[self.exact_dates['WELL_ID'] ==
                                                                  well.propnum,
                                                                  'START'] = well.start_date
                        else:
                            tmp_pods = [p - idw for p in pods]
                            idx_tmp = min(t for t in tmp_pods if t > 0)
                            idx_tmp = tmp_pods.index(idx_tmp)
                            idx_tmp = pods[idx_tmp] - 1
                            pod_compl = pad.well_list[idx_tmp].compl_date
                            pod_prod = (pod_compl
                                        + timedelta(well.compl_time)
                                        + timedelta(well.flowback_time))
                            well.start_date = pod_prod
                            self.exact_dates.loc[self.exact_dates['WELL_ID'] ==
                                                                  well.propnum,
                                                                  'START'] = well.start_date
                elif (pad.prod_start is not None and
                      pad.prod_finish is not None):
                    prod_time = ((pad.prod_finish - pad.prod_start).days
                                 / pad.num_wells)
                    for idw2, well2 in enumerate(pad.well_list):
                        if idw2 == 0:
                            well2.start_date = pad.compl_start
                            self.exact_dates.loc[self.exact_dates['WELL_ID'] ==
                                                                  well2.propnum,
                                                                  'START'] = well2.start_date
                        elif idw2 == len(pad.well_list) - 1:
                            well2.start_date = pad.start_finish
                            self.exact_dates.loc[self.exact_dates['WELL_ID'] ==
                                                                  well2.propnum,
                                                                  'START'] = well2.start_date
                        else:
                            well2.start_date = (pad.well_list[idw-1].start_date
                                                + timedelta(prod_time))
                            self.exact_dates.loc[self.exact_dates['WELL_ID'] ==
                                                                  well2.propnum,
                                                                  'START'] = well2.start_date
    def calc_build_pad_dates(self, rig):
        for pad in rig.pad_list:
            if rig.rig_name != 'DUC':
                if pad.build_start is None and pad.build_finish is None:
                    pad.build_start = (pad.drill_start -
                                       timedelta(1))
                    pad.build_finish = pad.drill_start - timedelta(1)
                elif (pad.build_start is None and
                      pad.build_finish is not None):
                    pad.build_start = (pad.build_finish -
                                       timedelta(pad.build_pad))
                elif (pad.build_start is not None and
                      pad.build_finish is None):
                    pad.build_finish = (pad.build_start +
                                        timedelta(pad.build_pad))
                elif (pad.build_start is not None and
                      pad.build_finish is not None):
                    continue

    def calc_all_dates(self):
        for rig in self.rig_dict.values():
            if rig != 'DUC':
                self.calc_drill_dates(rig)
                self.calc_build_pad_dates(rig)
            self.calc_compl_dates(rig)
            self.calc_start_dates(rig)
        print('Calculated all dates.')

    def output_exact(self):
        if self.output_path is None:
            file_name = self.schedule_name + ' Exact Dates.csv'
        else:
            file_name = (self.schedule_path +
                         self.schedule_name + ' Exact Dates.csv')
        with open(file_name, 'w') as write_file:
            line = ','.join(['PROPNUM', 'PAD', 'DRILL', 'COMPL', 'START'])
            write_file.write(line+'\n')
            for rig in self.rig_dict.values():
                for _, pad in enumerate(rig.pad_list):
                    for _, well in enumerate(pad.well_list):
                        if rig.rig_name != 'DUC':
                            line = ','.join([well.propnum,
                                             well.pad_name,
                                             well.drill_date.strftime(
                                                 '%m/%d/%Y'),
                                             well.compl_date.strftime(
                                                 '%m/%d/%Y'),
                                             well.start_date.strftime(
                                                 '%m/%d/%Y')])
                        else:
                            line = ','.join([well.propnum,
                                             well.pad_name,
                                             str(np.nan),
                                             well.compl_date.strftime(
                                                 '%m/%d/%Y'),
                                             well.start_date.strftime(
                                                 '%m/%d/%Y')])
                        write_file.write(line+'\n')
        print('Outputted exact dates.')

    def output_duc_counts(self):
        output_list = [[], [], [], []]
        for rig in self.rig_dict.values():
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
                            output_list[1][output_list[0].index(
                                drill_year)] += 1
                        else:
                            output_list[0].append(drill_year)
                            output_list[1].append(1)
                            output_list[2].append(0)
                            output_list[3].append(0)
                        compl_year = well.compl_date.year
                        if compl_year in output_list[0]:
                            output_list[2][output_list[0].index(
                                compl_year)] += 1
                        else:
                            output_list[0].append(compl_year)
                            output_list[2].append(1)
                            output_list[1].append(0)
                            output_list[3].append(0)
        write_list = []
        for i in range(len(output_list[0])):
            write_list.append([output_list[0][i], output_list[1][i],
                               output_list[2][i], output_list[3][i]])
        output_list = sorted(write_list, key=itemgetter(0))
        for i, _ in enumerate(output_list):
            if i == 0:
                continue
            else:
                output_list[i][3] = (output_list[i-1][3]
                                     + output_list[i][1] - output_list[i][2])
        if self.output_path is None:
            file_name = self.schedule_name+' DUC Counts.csv'
        else:
            file_name = self.schedule_path+self.schedule_name+' DUC Counts.csv'
        with open(file_name, 'w') as write_file:
            line = ','.join(['YEAR', 'DRILL', 'COMPL', 'DUC'])
            write_file.write(line+'\n')
            for idx, _ in enumerate(output_list):
                line = ','.join(map(str, [output_list[idx][0],
                                          output_list[idx][1],
                                          output_list[idx][2],
                                          output_list[idx][3]]))
                write_file.write(line+'\n')
        print('Outputted DUC counts.')

    def output_aries(self):
        if self.output_path is None:
            file_name = self.schedule_name + ' ARIES Dates.csv'
        else:
            file_name = (self.schedule_path +
                         self.schedule_name + ' ARIES Dates.csv')
        with open(file_name, 'w') as write_file:
            line = ','.join(['PROPNUM', 'DRILL', 'COMPL', 'START'])
            write_file.write(line+'\n')
            for rig in self.rig_dict.values():
                for _, pad in enumerate(rig.pad_list):
                    for _, well in enumerate(pad.well_list):
                        if rig.rig_name != 'DUC':
                            if well.drill_date.day >= 15:
                                drill = well.drill_date + relativedelta(
                                    months=1)
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
                            line = ','.join([well.propnum,
                                             str(well.drill_cost) +
                                             ' 0 G ' +
                                             drill.strftime('%m/%Y') +
                                             ' AD PC 0', '0 ' +
                                             str(well.compl_cost) +
                                             ' G ' + compl.strftime(
                                                 '%m/%Y') + ' AD PC 0',
                                             '=\"' + start.strftime(
                                                 '%m/%Y') + '\"'])
                        else:
                            line = ','.join([well.propnum, str(drill),
                                             '0 ' + str(well.compl_cost)
                                             + ' G ' + compl.strftime('%m/%Y')
                                             + ' AD PC 0', '=\"'
                                             + start.strftime('%m/%Y') + '\"'])
                        write_file.write(line+'\n')
        print('Outputted ARIES dates.')

    def create_gantt_df(self):
        temp_df = pd.DataFrame()
        rigs = []
        pads = []
        start = []
        start_text = []
        end = []
        end_text = []
        colors = ['#d3d3d3', '#1b9e77', '#d95f02', '#7570b3']
        color = []
        drill_time = []
        depth = []
        compl_time = []
        well_count = []
        for rig in self.rig_dict.values():
            for pad in rig.pad_list:
                avg_drill = round(np.mean(
                    [well.drill_time for well in pad.well_list]), 1)
                avg_depth = round(np.mean(
                    [well.depth for well in pad.well_list]), 0)
                avg_compl = round(np.mean(
                    [well.compl_time for well in pad.well_list]), 1)
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
        temp_df['rig'] = pd.Series(rigs)
        temp_df['pad'] = pd.Series(pads)
        temp_df['start'] = pd.Series(start)
        temp_df['start_text'] = pd.Series(start_text)
        temp_df['end'] = pd.Series(end)
        temp_df['end_text'] = pd.Series(end_text)
        temp_df['color'] = pd.Series(color)
        temp_df['drill'] = pd.Series(drill_time)
        temp_df['depth'] = pd.Series(depth)
        temp_df['compl'] = pd.Series(compl_time)
        temp_df['well_count'] = pd.Series(well_count)
        self.gantt_df = temp_df
        print('Created gantt chart dataframe.')

    def create_gantt_chart(self):
        if self.output_path is None:
            output_file(self.schedule_name+'.html')
        else:
            output_file(self.output_path+self.schedule_name+'.html',
                        title=self.schedule_name)
        height = 0.9
        plot_height = 300
        text_font_size = '7pt'
        text_font_style = 'bold'
        angle = 90
        angle_units = 'deg'
        x_offset = 10
        y_offset = 30
        fill_alpha = 0.75
        p_dict = {}
        for r in self.rig_dict.keys():
            dfr = self.gantt_df[self.gantt_df.rig.str.contains(r)]
            source = ColumnDataSource(dfr)
            rigs = list(dfr.rig.unique())
            rigs.reverse()
            p_dict[r] = figure(y_range=rigs, x_axis_type='datetime',
                               x_range=Range1d(self.gantt_start_date,
                                               self.gantt_end_date),
                               plot_width=1500, plot_height=plot_height,
                               toolbar_location='above', active_drag='pan',
                               title=r+' Schedule')
            p_dict[r].hbar(y='rig', left='start', right='end', color='color',
                           fill_alpha=fill_alpha, height=height,
                           line_color='gray', source=source)
            if r == 'DUC':
                hover = HoverTool(
                    tooltips=[
                        ('Pad', '@pad'),
                        ('Start', '@start_text'),
                        ('Finish', '@end_text'),
                        ('Avg Compl Time', '@compl'),
                        ('Well Count', '@well_count')
                    ]
                )
            else:
                hover = HoverTool(
                    tooltips=[
                        ('Pad', '@pad'),
                        ('Start', '@start_text'),
                        ('Finish', '@end_text'),
                        ('Avg Drill Time', '@drill'),
                        ('Avg Compl Time', '@compl'),
                        ('Well Count', '@well_count')
                    ]
                )
            labels = LabelSet(x='start', y='rig', text='pad', level='glyph',
                              angle=angle, angle_units=angle_units,
                              x_offset=x_offset, y_offset=-y_offset,
                              source=source, render_mode='canvas',
                              text_font_size=text_font_size,
                              text_font_style=text_font_style,
                              text_color='black')
            p_dict[r].add_tools(hover)
            p_dict[r].add_layout(labels)
            p_dict[r].ygrid.grid_line_color = None
            p_dict[r].xaxis.axis_label = "Date"
            p_dict[r].xaxis.formatter = DatetimeTickFormatter(
                months=["%b %Y"],
                days=["%b %e %Y"],
                years=["%b %Y"],
            )
            p_dict[r].outline_line_color = None
        p_columns = column(list(p_dict.values()))
        print('Created gantt chart.')
        show(p_columns)


class Rig():
    def __init__(self, rig_name, start, verbose=False):
        self.rig_name = rig_name
        self.start = start
        self.num_pads = 0
        self.pad_list = []
        self.conductors = 14
        self.mob_in = 7
        self.mob_out = 7
        self.verbose = verbose
        if verbose:
            print(self.rig_name, 'created. Project start:', self.start,
                  '\nProject starts with setting conductors, ',
                  'then mob to location, then drill.')

    def info(self):
        """Print the name and start date of the rig."""
        print('Rig Name:\t', self.rig_name)
        print('Start Date:\t', self.start)

    def add_pad(self, pad):
        """Assign a pad object to the rig. Increments pad count by 1."""
        self.pad_list.append(pad)
        pad.rig_name = self.rig_name
        self.num_pads += 1
        if self.verbose:
            print(pad.pad_name, 'added to', self.rig_name, 'rig.')

    def drop_pad(self, pad):
        """Drop a pad object from the rig. Decrements pad count by 1."""
        self.pad_list.remove(pad)
        pad.rig_name = None
        self.num_pads -= 1
        if self.verbose:
            print(pad.pad_name, 'dropped from', self.rig_name, 'rig.')

    def pads(self):
        """Print the pads assigned to the rig."""
        if not self.pad_list:
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
            print('Days to set conductors on the pad set to',
                  str(self.conductors)+'.')

    def set_mob_in(self, mob_in_days):
        """Set the number of days to mob drilling rig to pad."""
        self.mob_in = float(round(mob_in_days, 1))
        for pad in self.pad_list:
            pad.conductors = self.conductors
        if self.verbose:
            print('Days to mob drilling rig to pad set to',
                  str(self.mob_in)+'.')

    def set_mob_out(self, mob_out_days):
        """Set the number of days to mob drilling rig off pad."""
        self.mob_out = float(round(mob_out_days, 1))
        for pad in self.pad_list:
            pad.conductors = self.conductors
        if self.verbose:
            print('Days to mob drilling rig off pad set to',
                  str(self.mob_out)+'.')

    def set_drill_time(self, drill_time=10.0):
        """Set the number of days to drill a well for all
           wells assigned to pad."""
        for pad in self.pad_list:
            for well in pad.well_list:
                well.set_drill_time(drill_time)

    def set_drill_time_by_depth(self, feet_per_day=1600):
        """Set the number of days to drill a well based on the depth of the
           well for all wells assigned to the pad. The number of feet drilled
           per day is multiplied by the well's depth to calculate the number
           of days drilled, rounded to the nearest tenth of a day."""
        for pad in self.pad_list:
            for well in pad.well_list:
                well.set_drill_time_by_depth(feet_per_day)

    def schedule(self):
        print('Rig Start:', self.start, '\n')
        for idp, pad in enumerate(self.pad_list):
            print('Mob to', pad.pad_name, 'pad',
                  '(' + str(idp+1) + '/' + str(self.num_pads) + ')')
            for well in pad.well_list:
                print(well.well_name + ':\tDrill: ' + str(well.drill_date) +
                      '\tComplete: ' + str(well.compl_date) +
                      '\tStart: ' + str(well.start_date))
        print('\n')


class Pad():

    def __init__(self, pad_name, rig, verbose=False):
        self.pad_name = pad_name
        self.rig = rig
        self.rig.add_pad(self)
        self.num_wells = 0
        self.well_list = []
        self.build_pad = 30
        self.pod_size = 4
        self.conductors = self.rig.conductors
        self.mob_in = self.rig.mob_in
        self.mob_out = self.rig.mob_out
        self.log_pad = 0.5
        self.build_facilities = 60
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

    def add_well(self, well):
        """Assign a well object to the pad. Increments well count by 1."""
        self.well_list.append(well)
        well.pad_name = self.pad_name
        self.num_wells += 1
        if self.verbose:
            print(well.well_name, 'added to', self.pad_name, 'pad.')

    def drop_well(self, well):
        """Drop a well object from the pad. Decrements well count by 1."""
        self.well_list.remove(Well)
        well.pad_name = None
        self.num_wells -= 1
        if self.verbose:
            print(well.well_name, 'dropped from', self.pad_name, 'pad.')

    def wells(self):
        """Print the wells assigned to the pad, sorted alphabetically."""

        if not self.well_list:
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
            print(self.num_wells, 'wells assigned to the',
                  self.pad_name, 'pad.')

    def set_build_time(self, construct_days=30.0):
        """Set the number of days to construct the pad."""
        self.build_pad = float(round(construct_days, 1))
        if self.verbose:
            print('Days to build the', self.pad_name,
                  'pad set to', str(self.build_pad)+'.')

    def set_build_time_by_num_wells(self, construct_days_per_well=2.0):
        """Set the number of days to build the pad based on the well count.
           The number of wells assigned to the pad is multiplied by the number
           of days to construct for a single well to calculate the total
           number of days to build the pad, rounded to the nearest tenth of a
           day."""

        self.build_pad = float(round(construct_days_per_well *
                                     self.num_wells, 1))
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
           The number of wells assigned to the pad is multiplied by the number
           of days to buld facilities for a single well to calculate the total
           number of days to build facilities, rounded to the nearest tenth of
           a day."""
        self.build_facilities = float(round(fac_days_per_well *
                                            self.num_wells, 1))
        if self.verbose:
            print('Days to build facitilies for the', self.pad_name,
                  'pad set to', self.build_facilities, 'based on well count.')

    def set_drill_time(self, drill_time=10.0):
        """Set the number of days to drill a well for all wells assigned to
           pad."""
        for well in self.well_list:
            well.set_drill_time(drill_time)

    def set_drill_time_by_depth(self, feet_per_day=1600):
        """Set the number of days to drill a well based on the depth of the
           well for all wells assigned to the pad. The number of feet drilled
           per day is multiplied by the well's depth to calculate the number
           of days drilled, rounded to the nearest tenth of a day."""
        for well in self.well_list:
            well.set_drill_time_by_depth(feet_per_day)

    def set_logging_time(self, logging_days):
        """Set the number of days to log wells assigned to the pad."""
        self.log_pad = float(round(logging_days, 1))
        if self.verbose:
            print('Days to log the', self.pad_name,
                  'pad set to', str(self.log_pad)+'.')

    def set_logging_time_per_well(self, logging_days_per_well):
        """Set the number of days to log wells assigned to the pad based on
           well count."""
        self.log_pad = float(round(logging_days_per_well*self.num_wells, 1))
        if self.verbose:
            print('Days to log the', self.pad_name,
                  'pad set to', str(self.log_pad) + ' based on well count.')

    def set_compl_time(self, compl_days):
        """Set the number of days to complete wells assigned to the pad."""
        for well in self.well_list:
            well.set_compl_time(compl_days)

    def set_flowback_time(self, flowback_days):
        """Set the number of days post completion before first gas for wells
           assign to the pad."""

        for well in self.well_list:
            well.set_flowback_time(flowback_days)

    def set_pod_size(self, pod_size):
        """Set number of wells in each completion pod. Wells turn on in
           pod-sized batches at the same date post flowback."""

        self.pod_size = pod_size


class Well(object):

    def __init__(self, propnum, well_name,
                 pad, area, depth, verbose=False):
        self.propnum = propnum
        self.well_name = well_name
        self.pad = pad
        pad.add_well(self)
        self.area = area
        self.depth = depth
        self.drill_time = 10.0
        self.compl_time = 3.
        self.flowback_time = 30.
        self.buffer = 0.
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
            print('Days to drill the', self.well_name,
                  'set to', str(self.drill_time)+'.')

    def set_drill_time_by_depth(self, feet_per_day=1600):
        """Set the number of days to drill a well based on the depth ofd
           the well. The number of feet drilled per day is multiplied by the
           well's depth to calculate the number of days drilled, rounded to
           the nearest tenth of a day."""

        self.drill_time = round(self.depth / feet_per_day, 1)
        if self.verbose:
            print('Days to drill the', self.well_name,
                  'set to', self.drill_time, 'based on depth.')

    def set_buffer(self, num_days=0):
        self.buffer = round(num_days, 1)

    def set_compl_time(self, compl_days):
        """Set the number of days to complete a well."""
        self.compl_time = float(round(compl_days, 1))
        if self.verbose:
            print('Days to complete the', self.well_name,
                  'set to', str(self.compl_time)+'.')

    def set_flowback_time(self, flowback_days):
        """Set the number of days post completion before first gas."""
        self.flowback_time = float(round(flowback_days, 1))
        if self.verbose:
            print('Days post completion before first gas for the',
                  self.well_name,
                  'set to', str(self.flowback_time)+'.')
