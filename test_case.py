import aries_scheduler
from datetime import datetime

rigs = rig_starts = {'DUC':     None,
                     'RIG 01':  datetime(2019, 1, 1, 7, 0, 0), 
                     'RIG 02':  datetime(2019, 4, 1, 7, 0, 0),
                    }

name = 'Test Schedule'
path = 'test_case.xlsx'

test_schedule = aries_scheduler.Schedule(name, path, rigs)

