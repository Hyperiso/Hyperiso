import os

current_file = os.path.abspath(__file__)
current_dir = os.path.dirname(current_file)

db_path = os.path.abspath(os.path.join(current_dir, "..", "..", "../"))

