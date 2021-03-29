from definitions import  *
import os, sys
import datetime

def set_results_message(msg):
	global RESULT_MSG
	RESULT_MSG += msg


def get_results_message():
	global RESULT_MSG
	return RESULT_MSG
	
def write_daily_test(run_number, status):
	date = datetime.datetime.today().strftime('%d%m%Y')
	DAILY_TESTS_DIR = '/bioseq/bioSequence_scripts_and_constants/daily_tests'
	results_url = f'http://modelteller.tau.ac.il/results/{run_number}/output.html'
	with open(os.path.join(DAILY_TESTS_DIR, f'modelteller_{date}.txt'), "w") as f:
		f.write(f'{status},{results_url}')
	f.close()