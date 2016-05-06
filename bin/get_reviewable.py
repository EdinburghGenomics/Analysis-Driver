
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from bin import automatic_review
from analysis_driver import rest_communication

def get_reviewable_runs():
    runs = rest_communication.get_documents('aggregate/all_runs', match={"proc_status":"finished"})