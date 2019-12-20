###methods useful for analysis of metabolomics files
import pandas as pd
import numpy as np
from scipy import stats
import csv
import matplotlib.pyplot as plt
import glob
import os
import re
import itertools
import GPy
import sys
sys.path.append("..")

from mzmine import pick_peaks, align, match_aligned_to_original, load_aligned_peaks

def create_std_dict(input_std_csv_file):

    """ Creates and returns a dictionary for a standards csv file obtained from ToxID."""

    standards = {}
    for _,row in input_std_csv_file.iterrows():
        compound_name = row['Compound Name']
        detected_mz = row['Detected m/z']
        expected_rt = row['Expected RT']
        actual_rt = row['Actual RT']
        polarity = row['Polarity']

        standards[compound_name] = (detected_mz, expected_rt, actual_rt, polarity)
    return standards

def select_std_from_mzmine(input_mzmine_file, input_std_dict, mz_range= 0.0003, rt_range= 0.5, polarity = "+"):

    """For 1 file: Matching the peaks picked with mzmine with the peaks from the standards csv files.
    Creating and returns a dictionary with the matched peaks: Key = Compound name, Value = (m/z (from mzmine), RT (from mzmine)) """

    mzmine_standards = {}
    name = input_mzmine_file.columns[3].split(".")[0]
    print(name)
    x=0
    p=0
    masses_detected_in_std_file = len(input_std_dict.values())
    for i in range(len(input_std_dict.values())):

        #print(i,':mz in std csv:',list(input_std_dict.values())[i][0], list(input_std_dict.keys())[i])
        if (list(input_std_dict.values())[i][3] == polarity):

            if (list(input_std_dict.values())[i][0] == 'N'):
                p+=1
            else:
                for _,row in input_mzmine_file.iterrows():
                    rowid = row['row ID']
                    mz = row['row m/z']
                    rt = row['row retention time']
                    n=0

                    if (mz <= float(list(input_std_dict.values())[i][0])+mz_range and mz >= float(list(input_std_dict.values())[i][0])-mz_range):
                        if (rt <= float(list(input_std_dict.values())[i][2])+rt_range and rt >= float(list(input_std_dict.values())[i][2])-rt_range):
                            #print("Found")
                            #print(mz)
                            #print(rt)
                            mzmine_standards[list(input_std_dict.keys())[i]]=(rowid,mz,rt)

                            #print(mz-list(standards.values())[i][0])
                            n+=1
                    if (n>=1):
                        x+=1
        #else:
            #print("wrong polarity")
    return (mzmine_standards, name)


def get_std_matches_for_files(original_files, output_dir, stds_csvs, original_csvs, rt_range = 0.5, mz_range = 0.0003):

    matches = {}
    for file_pos,o in enumerate(original_files):


        with open(stds_csvs[file_pos],'r') as s:
             for line in s:
                #to skip the comments read only the rows which start with a number
                if re.match(r"^\d+.*$",line):
                    name = line.split(',')[2]
                    polarity = line.split(',')[4]
                    mz = float(line.split(',')[6])
                    rt = float(line.split(',')[9])

                    with open(original_csvs[file_pos], 'r') as c:
                        reader = csv.reader(c)
                        header = next(reader)
                        for line in reader:

                            id_o = int(line[0])
                            mz_o = float(line[1])
                            rt_o = float(line[2])
                            int_o = float(line[3])

                            if not name in matches:
                                matches[name] = {}

                            if (mz_o <= mz+mz_range and mz_o >= mz-mz_range):
                                if (rt_o <= rt+rt_range and rt_o >= rt-rt_range):
                                    matches[name][o] = (id_o,mz_o,rt_o,int_o)



    output_file = os.path.join(output_dir,'stds_match_links.csv')
    with open(output_file,'w') as f:
        heads = ['stds_name'] + original_files
        writer = csv.writer(f)
        writer.writerow(heads)
        for matched_peak in matches:

            new_row = [matched_peak]

            for o in original_files:

                val = matches[matched_peak].get(o,None)

                if val:
                    new_row.append(val[0])
                else:
                    new_row.append('null')
            writer.writerow(new_row)

    return matches

def get_rt_diff_between_datasets(dataset1_mzmine_stds, dataset2_mzmine_stds, file_number = 3):

    dataset1_rt = []
    dataset2_rt = []
    stds = {}
    diff = []

    for i in range(file_number):
        dataset1 = set(dataset1_mzmine_stds[i])
        dataset2 = set(dataset2_mzmine_stds[i])

        for name in dataset1.intersection(dataset2):
            #print(name, (dataset1_mzmine_stds[i][name][1]-dataset2_mzmine_stds[i][name][1])*60)
            diff.append((dataset1_mzmine_stds[i][name][2]-dataset2_mzmine_stds[i][name][2])*60)
            dataset1_rt.append(dataset1_mzmine_stds[i][name][2])
            dataset2_rt.append(dataset2_mzmine_stds[i][name][2])

            stds[name] = (i+1, dataset1_mzmine_stds[i][name], dataset2_mzmine_stds[i][name])



    return (dataset1_rt, dataset2_rt, stds, diff)


def get_stats_on_diff(d1_stds_list,name1, d2_stds_list,name2, setbinwidth, zcut = 1):

    (d1_rt, d2_rt, stds_d1_d2, d1_d2_diff) = get_rt_diff_between_datasets(d1_stds_list,d2_stds_list)

    zscore = np.abs(stats.zscore(d1_d2_diff))
    d1_rt_mod = []
    d2_rt_mod = []

    for (d1, d2, z) in zip(d1_rt, d2_rt, zscore):
         if (z<zcut):
                d1_rt_mod.append(d1)
                d2_rt_mod.append(d2)

    d1_d2_diff_no_outliers = (np.subtract([x * 60 for x in d1_rt_mod],[x*60 for x in d2_rt_mod]))

    d1_d2_diff_std = np.std(d1_d2_diff_no_outliers)
    d1_d2_diff_mean = np.mean(d1_d2_diff_no_outliers)
    d1_d2_diff_max = np.max(abs(d1_d2_diff_no_outliers))

    print(d1_d2_diff_mean, d1_d2_diff_std,  d1_d2_diff_max)

    binwidth = setbinwidth
    plt.hist(d1_d2_diff_no_outliers, bins=np.arange(min(d1_d2_diff_no_outliers), max(d1_d2_diff_no_outliers) + binwidth, binwidth), alpha =0.5)
    plt.xlabel('RT '+ name1+ ' - RT ' +name2 + '(s)')
    plt.show()


def try_gp_regressions(t_dataset_mod, t_reference_mod):

    from sklearn import metrics

    k_lin = GPy.kern.Linear(1)
    k_exp = GPy.kern.Exponential(1)
    k_nn = GPy.kern.MLP(1)

    k_rbf = GPy.kern.RBF(input_dim=1, variance=2.25, lengthscale=1.5)
    k_mat32 = GPy.kern.Matern32(input_dim=1, variance=2., lengthscale=0.2)
    k_mat52 = GPy.kern.Matern52(1)

    k_per = GPy.kern.StdPeriodic(1, period=3.)
    k_cos = GPy.kern.Cosine(1)
    k_brwn = GPy.kern.Brownian(1)

    t_reference = np.array(t_reference_mod)
    t_dataset = np.array(t_dataset_mod)

    REF = t_reference[:,None]
    TOPRED = t_dataset[:,None]


    #Cross-validation
    training_data = TOPRED[:len(TOPRED)//2] #Need to put double / so it could work
    training_data_target = REF[:len(REF)//2]



    test_data = TOPRED[len(TOPRED)//2+1:]
    test_data_target = REF[len(REF)//2+1:]

    #write a vector with kernels which might work
    ks = [k_nn, k_nn+k_brwn, k_nn*k_brwn, k_nn+k_rbf, (k_cos+k_exp),(k_cos+k_exp)*k_rbf, (k_cos+k_lin)*(k_rbf), (k_cos+k_lin)+(k_per), (k_cos+k_lin)*(k_mat52)]
    ks_names = ['MLP','MLP+Brownian','MLP*Brownian','MLP+RBF','Cos+Exponential', '(Cos+Exponential)*RBF', '(Cos+Linear)*RBF', '(Cos+Linear)*Periodic', '(Cos+Linear)*Matern52']
    names = ['RBF','MLP','MLP+Brownian','MLP*Brownian','MLP+RBF','Cos+Exponential', '(Cos+Exponential)*RBF', '(Cos+Linear)*RBF', '(Cos+Linear)*Periodic', '(Cos+Linear)*Matern52']
    accuracy_list = []
    mae_list = []
    mse_list = []

    mtest = GPy.models.GPRegression(training_data,training_data_target,k_rbf)
    mtest.optimize()
    mtest.plot()
    plt.title('RBF')
    plt.show()

    gpr_predicted_data,_ = mtest.predict(test_data)
    plt.plot(test_data_target, gpr_predicted_data,  '.')
    plt.plot(test_data_target, test_data_target, '.r')
    plt.ylabel('Predicted data')
    plt.xlabel('True data')
    plt.show()

    accuracy0 = metrics.r2_score(test_data_target, gpr_predicted_data)
    mae0 = metrics.mean_absolute_error(test_data_target, gpr_predicted_data)
    mse0 = metrics.mean_squared_error(test_data_target, gpr_predicted_data)

    accuracy_list.append(accuracy0)
    mae_list.append(mae0)
    mse_list.append(mse0)

    print('Cross-Predicted Accuracy for','RBF',':', accuracy0)
    print('Mean absolute error for','RBF',':', mae0)
    print('Mean squared error for','RBF',':', mse0)

    mfinal = mtest
    kfinal = k_rbf
    namefinal = 'RBF'

    for k, name in zip(ks,ks_names):
        #regression and optimisation of model
        mtest = GPy.models.GPRegression(training_data,training_data_target,k)
        mtest.optimize()
        mtest.plot()
        plt.title(name)
        plt.show()

        gpr_predicted_data,_ = mtest.predict(test_data)
        plt.plot(test_data_target, gpr_predicted_data,  '.')
        plt.plot(test_data_target, test_data_target, '.r')
        plt.ylabel('Predicted data')
        plt.xlabel('True data')
        plt.show()

        accuracy = metrics.r2_score(test_data_target, gpr_predicted_data)
        mae = metrics.mean_absolute_error(test_data_target, gpr_predicted_data)
        mse = metrics.mean_squared_error(test_data_target, gpr_predicted_data)

        accuracy_list.append(accuracy)
        mae_list.append(mae)
        mse_list.append(mse)

        print('Cross-Predicted Accuracy for',name,':', accuracy)
        print('Mean absolute error for',name,':', mae)
        print('Mean squared error for',name,':', mse)

        if (accuracy >= accuracy0 and mae < mae0 and mse < mse0):
            accuracy0 = accuracy
            mae0 = mae
            mse0 = mse
            mfinal = mtest
            kfinal = k
            namefinal = name
    results_table = pd.DataFrame(data=[names, accuracy_list, mae_list, mse_list])
    print(namefinal)
    return mfinal, kfinal, results_table


def return_data_with_no_outliers(dataset1, dataset2, diff, zscore_cutoff):
    #Eliminating outliers using zscore
    zscore = np.abs(stats.zscore(diff))
    mod1 = []
    mod2 = []

    for (d1, d2, z) in zip(dataset1, dataset2, zscore):
         if (z<zscore_cutoff):
            mod1.append(d1)
            mod2.append(d2)

    return np.array(mod1), np.array(mod2)

def get_non_anomalies(dataset_1_rt, dataset_2_rt):
    from sklearn.ensemble import IsolationForest

    data = np.array([dataset_1_rt, dataset_2_rt]).transpose()
    df = pd.DataFrame(data)
    l = len(dataset_1_rt)

    clf = IsolationForest(behaviour = 'auto', max_samples=86, random_state = 1, contamination= 'auto')
    pred = clf.fit_predict(data)
    df['pred'] = pred
    df_no_anomalies = df[df['pred'] == 1]
    dataset_1_rt_mod = np.array(df_no_anomalies[0])
    dataset_2_rt_mod = np.array(df_no_anomalies[1])
    return dataset_1_rt_mod, dataset_2_rt_mod

def modify_rtdrift_in_mztab(output_dir, output_dir_new, model, name_tag, standards = False):

    if standards == False:
        filelen = len(glob.glob(os.path.join(output_dir,name_tag)))
    else:
        filelen = 6
    
   
    peak_file_list = glob.glob(os.path.join(output_dir, name_tag))
    peak_file_list.sort()
    peak_file_list_new = glob.glob(os.path.join(output_dir_new, name_tag))
    peak_file_list_new.sort()
    
    for i in range(filelen):
        
        with open(peak_file_list[i], 'r') as f:
            with open(peak_file_list_new[i],'w') as g:
                for line in f:
                    line = line.rstrip()
                    if line.startswith('SMH'):
                        tokens = line.split()
                        rt_pos = [tokens.index('retention_time')]
                        rt_pos.append(tokens.index('opt_assay[1]_peak_rt'))
                    if not line.startswith('SML'):
                        g.write(line + '\n')
                    else:
                        tokens = line.split()
                        old_rt = float(tokens[rt_pos[0]])

                        new_rt,_ = model.predict(np.array([[old_rt]]))
                        print(new_rt[0][0])
                        new_rt = new_rt[0][0]
                        print(old_rt)
                        new_rt += old_rt
                        print(new_rt)

                        tokens[rt_pos[0]] = (new_rt)
                        tokens[rt_pos[1]] = (new_rt)
                        new_line = '\t'.join([str(t) for t in tokens])
                        g.write(new_line + '\n')

def modify_rt_in_mztab(output_dir, output_dir_new, model, name_tag = '*.mzTab', standards = False):

    if standards == False:
        filelen = len(glob.glob(os.path.join(output_dir,name_tag)))
    else:
        filelen = 6

    for i in range(filelen):
        with open(glob.glob(os.path.join(output_dir,name_tag))[i], 'r') as f:
            with open(glob.glob(os.path.join(output_dir_new,name_tag))[i],'w') as g:
                for line in f:
                    line = line.rstrip()
                    if line.startswith('SMH'):
                        tokens = line.split()
                        rt_pos = [tokens.index('retention_time')]
                        rt_pos.append(tokens.index('opt_assay[1]_peak_rt'))
                    if not line.startswith('SML'):
                        g.write(line + '\n')
                    else:
                        tokens = line.split()
                        old_rt = float(tokens[rt_pos[0]])
                        print(old_rt)
                        #new_rt = brr.predict(np.array(old_rt).reshape(-1,1))[0]
                        new_rt, = model.predict(np.array(old_rt).reshape(-1,1))[0]
                        print(new_rt[0])

                        tokens[rt_pos[0]] = new_rt[0]
                        tokens[rt_pos[1]] = new_rt[0]
                        new_line = '\t'.join([str(t) for t in tokens])
                        g.write(new_line + '\n')

def get_good_peaks(file, matches, stds_matches):
    good_peaks = 0
    for match in matches:
        for name in stds_matches:
            if check_peak_is_good(name, matches[match], file, stds_matches):
                good_peaks+=1
    return good_peaks

def check_peak_is_good(name, peak, files, stds_matches):

    for file in files:
        if file in peak:
            peak_id = peak[file][0]
            if file in stds_matches[name]:
                peak_id2 = stds_matches[name][file][0]
                if (peak_id != peak_id2):
                    return False
            else:
                return False
        else:
            return False

    return True

def get_total_stds(file, stds_matches):
    total = 0
    for name in stds_matches:
        if file in stds_matches[name]:
            total+=1

    return total

def get_peaks_for_files(output_dir, files, rt_window, xml_template, stds_matches):
    MZMINE_COMMAND = '/Users/anamaria/Documents/git/MZmine-2.40.1/startMZmine_MacOSX.command'
    
    print(files)
    dict_peaks = {}
    dict_peaks[str(files)] = []

    for rt in rt_window:
        print(rt)
        align(output_dir, xml_template = xml_template, mzmine_command = MZMINE_COMMAND, rt_window = rt)
        aligned_peaks,original_file,f_idx_dict = load_aligned_peaks(output_dir)
        #print(f_idx_dict)
        matches = match_aligned_to_original(aligned_peaks,files,output_dir,f_idx_dict,write_file = True)
        good_peaks = get_good_peaks(files, matches, stds_matches)

        dict_peaks[str(files)].append((len(aligned_peaks),good_peaks))


    return dict_peaks
