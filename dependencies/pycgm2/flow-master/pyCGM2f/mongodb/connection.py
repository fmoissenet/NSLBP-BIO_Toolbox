import pyCGM2f
from pymongo import MongoClient

class MongoConnection(object):
    def __init__(self):
        client = MongoClient(pyCGM2f.MONGODB_ADRESS, pyCGM2f.MONGODB_PORT)
        self.db = client['Gaitabase']

    def get_collection(self, name):
        return self.db[name]
