class Metabolite(object):
    def __init__(self, filename, source_id, mz, rt, intensity):
        self.filename = filename
        self.source_id = source_id
        self.mz = mz
        self.rt = rt
        self.intesity = intensity

class MetaboliteSet(object):
    def __init__(self, name):
        self.name = name
        self.files = []
        self.metabolites = []
        self.n_files = 0

    def add_metabolite(self, metabolite):
        self.metabolites.append(metabolite)

    def add_file(self, file):
        self.files.append(file)
        self.n_files += 1
