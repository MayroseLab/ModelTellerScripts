from definitions import  *


def set_results_message(msg):
	global RESULT_MSG
	RESULT_MSG += msg


def get_results_message():
	global RESULT_MSG
	return RESULT_MSG