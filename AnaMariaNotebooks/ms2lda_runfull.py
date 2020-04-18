import numpy as np
import sys
import os
import csv

# Put this here as its now the only thing you need from the motifdb codebase
class FeatureMatcher(object):
    def __init__(self,db_features,other_features,bin_width=0.005):
        self.db_features = db_features
        self.other_features = other_features
        self.fmap = {}
        self.bin_width = bin_width
        self.augmented_features = {f:v for f,v in other_features.items()}

        self.match()
        self.match(ftype='loss')



    def match(self,ftype='fragment'):
        import bisect
        other_names = [f for f in self.other_features if f.startswith(ftype)]
        other_min_mz = [self.other_features[f][0] for f in self.other_features if f.startswith(ftype)]
        other_max_mz = [self.other_features[f][1] for f in self.other_features if f.startswith(ftype)]

        temp = list(zip(other_names,other_min_mz,other_max_mz))
        temp.sort(key = lambda x: x[1])
        other_names,other_min_mz,other_max_mz = zip(*temp)
        other_names = list(other_names)
        other_min_mz = list(other_min_mz)
        other_max_mz = list(other_max_mz)

        exact_match = 0
        new_ones = 0
        overlap_match = 0
        for f in [f for f in self.db_features if f.startswith(ftype)]:
            if f in other_names:
                self.fmap[f] = f;
                exact_match += 1
            else:
                fmz = float(f.split('_')[1])
                if fmz < other_min_mz[0] or fmz > other_max_mz[-1]:
                    self.fmap[f] = f
                    self.augmented_features[f] = (fmz-self.bin_width/2,fmz+self.bin_width/2)
                    new_ones += 1
                    continue
                fpos = bisect.bisect_right(other_min_mz,fmz)
                fpos -= 1
                if fmz <= other_max_mz[fpos]:
                    self.fmap[f] = other_names[fpos]
                    overlap_match += 1
                else:
                    self.fmap[f] = f
                    self.augmented_features[f] = (fmz-self.bin_width/2,fmz+self.bin_width/2)
                    new_ones += 1
        print("Finished matching ({}). {} exact matches, {} overlap matches, {} new features".format(ftype,exact_match,overlap_match,new_ones))

    def convert(self,dbspectra):
        for doc,spec in dbspectra.items():
            newspec = {}
            for f,i in spec.items():
                newspec[self.fmap[f]] = i
            dbspectra[doc] = newspec
        return dbspectra


"""Loading Motifs"""

# note to Ming: you'll need to set these correctly...

# motifdbcodepath = os.path.join(os.path.dirname(
#     os.path.realpath(__file__)), 'code/utilities')

#motifdbcodepath = '/Users/simon/git/motifdb/code/utilities'


# motifdbpath = os.path.join(os.path.dirname(
#     os.path.realpath(__file__)), 'motifs')

# #motifdbpath = '/Users/simon/git/motifdb/motifs'

# sys.path.append(motifdbcodepath)

# from motifdb_loader import load_db
# from motifdb_loader import MotifFilter

# NEW MOTIFDB CODE

import requests
server_url = 'http://ms2lda.org/motifdb/'
motifset_dict = requests.get(server_url+'list_motifsets/').json()
# db_list = ['gnps_binned_005']  # Can update this later with multiple motif sets
print(motifset_dict.keys())
db_list = [motifset_dict['Massbank library derived Mass2Motifs'], motifset_dict['GNPS library derived Mass2Motifs']]

# Obtain a token
client = requests.session()
token_output = client.get(server_url + 'initialise_api/').json()
token = token_output['token']
data = {'csrfmiddlewaretoken':token}
data['motifset_id_list'] = db_list
data['filter'] = 'True'

output = client.post(server_url + 'get_motifset/',data = data).json()
motifdb_spectra = output['motifs']
motifdb_metadata = output['metadata']
motifdb_features = set()
for m,spec in motifdb_spectra.items():
    for f in spec:
        motifdb_features.add(f)

# following should be mscluster or mzmine
input_format = 'mzmine'

input_csv_file = sys.argv[1]

# mgf file name
input_mgf_file = sys.argv[2]

# output prefix
output_prefix = sys.argv[3]



print(input_csv_file, input_mgf_file)


#ldacodepath = os.path.join(os.path.dirname(
#    os.path.realpath(__file__)), 'lda/code')

#ldacodepath = '/Users/anamaria/Documents/git/lda/code'
ldacodepath = 'D:/git/lda/code'

sys.path.append(ldacodepath)

from ms2lda_feature_extraction import LoadMGF, MakeBinnedFeatures

# Assuming only one mgf file...
name_field = "SCANS"

if input_format == 'mscluster':
    csv_id_field = None
    mgf_id_field = None
    input_csv_file = None
else:
    csv_id_col = 'row ID'
    mgf_id_field = 'SCANS'

name_field = 'SCANS'

if input_format == 'mscluster':
    l = LoadMGF(name_field=name_field)
else:  # mzmine pipeline, can do clever ID matching
    l = LoadMGF(name_field=name_field, peaklist=input_csv_file,
                csv_id_col=csv_id_col, id_field=mgf_id_field)

ms1, ms2, metadata = l.load_spectra([input_mgf_file])
print("Loaded {} spectra".format(len(ms1)))

m = MakeBinnedFeatures(bin_width=0.005)  #What value do you want here?? TODO: Parameterize
corpus, features = m.make_features(ms2)
corpus = corpus[list(corpus.keys())[0]]

# from motifdb_loader import FeatureMatcher
fm = FeatureMatcher(motifdb_features, features)
motifdb_spectra = fm.convert(motifdb_spectra)


# Add the motifdb features to avoid problems when loading the dict into vlda later
bin_width = m.bin_width
added = 0
for f in motifdb_features:
    if not f in features:
        word_mz = float(f.split('_')[1])
        word_mz_min = word_mz - bin_width / 2
        word_mz_max = word_mz + bin_width / 2
        features[f] = (word_mz_min, word_mz_max)
        added += 1

print("Added {} features".format(added))

from lda import VariationalLDA

# # Debug code!!! - remove once done!!

# molecules = corpus.keys()

# sub_corpus = {}

# for mol in molecules[:100]:

#     sub_corpus[mol] = corpus[mol]

# corpus = sub_corpus

# # End of bit to remove


K = 300  # number of *new* topics

vlda = VariationalLDA(corpus, K=K, normalise=1000.0,
                      fixed_topics=motifdb_spectra,
                      fixed_topics_metadata=motifdb_metadata)

# note that for real runs the number of iterations is recommended to be 1000 or higher
vlda.run_vb(initialise=True, n_its=10)

vd = vlda.make_dictionary(
    features=features, metadata=metadata, filename=output_prefix + '.dict')

from ms2lda_molnet_integration import write_output_files
write_output_files(vd, output_prefix, metadata,
                   overlap_thresh=0.3, p_thresh=0.1, X=5,motif_metadata = motifdb_metadata)



# Writing the report - ntoe that you might need to set the 'backend' argument
# for this method to work (see the method in lda.py) as it depends what on
# your system will render the pdf...
from lda import write_topic_report
write_topic_report(vd,output_prefix+'_topic_report.pdf')
sys.exit(0)


# TODO: what should this output format be???

with open(sys.argv[2], 'w') as tsvfile:

    fieldnames = ['filename', 'scan', 'motifs']

    writer = csv.DictWriter(tsvfile, fieldnames=fieldnames, delimiter='\t',)

    writer.writeheader()

    for output_object in output_list:

        writer.writerow(output_object)


#with open(sys.argv[3], 'w') as tsvfile:

#    fieldnames = ['filename', 'scan', 'motif_id',
#                  "motif_overlapscore", "motif_metadata"]

#    writer = csv.DictWriter(tsvfile, fieldnames=fieldnames, delimiter='\t',)

#    writer.writeheader()

#    for output_object in motif_separate_output_results:

#        writer.writerow(output_object)
