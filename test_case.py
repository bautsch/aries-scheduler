'''
Test scheduleer
'''

import aries_scheduler
# import pandas as pd

SCHEDULE_NAME = 'Test Schedule'
SCHEDULE_FILE_PATH = 'test_case.xlsx'

TEST_SCHEDULE = aries_scheduler.Schedule(SCHEDULE_NAME, SCHEDULE_FILE_PATH)

# schedule = pd.read_excel(SCHEDULE_FILE_PATH, sheet_name='SCHEDULE')
