'''
Test scheduler
'''

from datetime import date
import aries_scheduler
# import pandas as pd

SCHEDULE_NAME = 'Test Schedule'
SCHEDULE_FILE_PATH = 'test_case.xlsx'
GANTT_START_DATE = date(2019, 1, 1)

TEST_SCHEDULE = aries_scheduler.Schedule(schedule_name=SCHEDULE_NAME,
                                         schedule_file_path=SCHEDULE_FILE_PATH,
                                         gantt_start_date=GANTT_START_DATE)

# schedule = pd.read_excel(SCHEDULE_FILE_PATH, sheet_name='SCHEDULE')
